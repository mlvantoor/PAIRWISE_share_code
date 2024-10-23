###################################################################################################
##
## fetch "Beta-Lactamase DataBase" (BLDB)
## match BLDB database with sequencing results from Resfinder
##
## version: 2024-10-23
##
## Note: this code is not entirely stable as it depends on the resfinder file having the same structure throughout all sheets. In the files that I received, that has not necessarily been the case. 
##
####################################################################################################

if(!'rvest' %in% row.names(installed.packages())){
  install.packages('rvest', dependencies=T)
}
if(!'plyr' %in% row.names(installed.packages())){
  install.packages('plyr', dependencies=T)
}
if(!'tidyr' %in% row.names(installed.packages())){
  install.packages('tidyr', dependencies=T)
}
if(!'readxl' %in% row.names(installed.packages())){
  install.packages('readxl', dependencies=T)
}

# load required packages
library(rvest) # web-scraping
library(plyr) # data management
library(tidyr) # I'm generally a base-R person, but this is great for converting tables from wide to long & vice-versa
library(readxl) # for reading .xls/xlsx-files into R

###################################################################################################
## fetch the BLDB data base

# links to the individual tables on BLDB (there is one table per class)
bldb.links <- data.frame(
  link=c('http://bldb.eu/BLDB.php?prot=A', 'http://bldb.eu/BLDB.php?prot=B1', 
         'http://bldb.eu/BLDB.php?prot=B2', 'http://bldb.eu/BLDB.php?prot=C', 
         'http://bldb.eu/BLDB.php?prot=D'),
  class=c('A', 'B1', 'B2', 'C', 'D')
)

# scrape tables from BLDB
bldb <- ddply(bldb.links, 'link', function(x){
  tab <- x$link %>%
    read_html() %>% 
    html_table()
  tab <- tab[[2]]
  tab <- tab[which(tab[,1]==x$class),]
  return(tab)
})

# remove links from the information on protein function
bldb$Functionalinformation <- gsub(' view', '', bldb$Functionalinformation)
bldb$Functionalinformation <- gsub('view', '', bldb$Functionalinformation)
bldb$Functionalinformation <- gsub('Narrow', 'narrow spectrum', bldb$Functionalinformation)
bldb$Functionalinformation <- gsub('NOTE: .*', '', bldb$Functionalinformation)

# simplify protein naming (convert to all lower cases and remove hyphens)
# this is to simplify matching between Resfinder data and data base
bldb$protein.smpl <- tolower(bldb$Proteinname)
bldb$protein.smpl <- gsub('-', '', bldb$protein.smpl)

###################################################################################################
## import & shape raw Resfinder data

# choose your corresponding resfinder file
resf <- file.choose()

# sequencing runs are stored in separate sheets
runs <- excel_sheets(resf)

# read sheets
# NOTE: I'm keeping sample ID and 
seq_results <- rbind.fill(lapply(runs[1:10], function(run, resfinder.file){
  print(run) # print sheet
  
  # import data & header separately & stick things together
  tmp <- read_excel(resfinder.file, sheet=run, skip=2, col_names=F, col_types='text')
  header1 <- as.vector(t(read_excel(resfinder.file, sheet=run, col_names=F, col_types='text')[1,]))
  header1[1] <- 'isolate'
  header2 <- as.vector(t(read_excel(resfinder.file, sheet=run, col_names=F, col_types='text')[2,]))
  header <- c(header1[1:4], header2[5:length(header2)])
  names(tmp) <- header
  
  tmp$`Z-score` <- as.numeric(tmp$`Z-score`)
  
  # there are multiple rows for each isolate, but the information about isolate id, species, and z-score are not complete
  # this bit of code fills those gaps 
  while(any(is.na(tmp$isolate))){
    start <- which(is.na(tmp$isolate))[1]
    tmp$isolate[start] <- tmp$isolate[max(which(!is.na(tmp$isolate[1:start])))]
    tmp$`Species (MALDI-TOF)`[start] <- tmp$`Species (MALDI-TOF)`[max(which(!is.na(tmp$`Species (MALDI-TOF)`[1:start])))]
    tmp$`Confirmed species`[start] <- tmp$`Confirmed species`[max(which(!is.na(tmp$`Confirmed species`[1:start])))]
    tmp$`Z-score`[start] <- tmp$`Z-score`[max(which(!is.na(tmp$`Z-score`[1:start])))]
  }
  
  tmp.long <- ddply(tmp, 'isolate', function(y){
    id <- unique(y$isolate); print(id)
    
    # this bit is for my own use and splits the isolate ID into separate sample and isolate IDs.
    # I'm commenting out because it might not work for all sample naming schemes
    # split <- strsplit(id, '')[[1]]
    # y$sample_id <- paste(split[1:4], collapse='')
    # if(!y$sample_id[1] %in% se.samples$sample_id){return(NULL)}
    
    # y$sub_id <- paste(split[5:length(split)], collapse='')
    y$species <- ifelse(any(!is.na(y$`Confirmed species`)), unique(y$`Confirmed species`[!is.na(y$`Confirmed species`)]), NA)
    y$z_score <- ifelse(any(!is.na(y$`Z-score`)), unique(y$`Z-score`[!is.na(y$`Z-score`)]), NA)
    
    # which groups are present in the data?
    keep <- c('isolate', 'species', header2[!is.na(header2)])
    keep <- keep[keep %in% names(y)]
    
    # convert to long table, with one row per gene and isolate
    y <- tidyr:::pivot_longer(y[,keep], cols=keep[-c(1:2)], names_to='group', values_to='gene')
    return(y[!is.na(y$gene),])
  })
  
  return(tmp.long)
}, resfinder.file=resf))


# from now on, we will only look at beta-lactams
seq <- seq_results[!is.na(seq_results$gene) & seq_results$group=='β-lactams',]
seq$gene.smpl <- ifelse(grepl('bla-', seq$gene), gsub('bla-', '', seq$gene),
                        ifelse(grepl('bla', seq$gene), gsub('bla', '', seq$gene), seq$gene))
seq$gene.smpl <- tolower(gsub(' \\(.*', '', seq$gene.smpl))
seq$gene.smpl <- gsub('-', '', seq$gene.smpl)


# if we can find information, include, else suggest that a sample needs to be checked manually
seq$functional <- unlist(lapply(seq$gene.smpl, function(x){
  if(x %in% bldb$protein.smpl){
    return(bldb$Functionalinformation[bldb$protein.smpl==x])
  }else{
    return('manual check')
  }
}))

# check whether genes are natural (N) or have been acquired (A)
seq$natural.aquired <- unlist(lapply(seq$gene.smpl, function(x){
  if(x %in% bldb$protein.smpl){
    return(bldb$`Natural (N) orAcquired (A)`[bldb$protein.smpl==x])
  }else{
    return('manual check')
  }
}))

table(seq$functional=='manual check')

table(seq$functional[seq$natural.aquired=='A'])


