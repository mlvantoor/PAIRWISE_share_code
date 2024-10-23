###################################################################################################
##
## import sample information from cultivation & AST result files
##
## version: 2024-10-23
##
## Note: this code is not entirely stable as it depends on the cultivation results files having the same, consistent structure - which is not always a given. 
##
####################################################################################################

if(!'plyr' %in% row.names(installed.packages())){
  install.packages('plyr', dependencies=T)
}
if(!'readxl' %in% row.names(installed.packages())){
  install.packages('readxl', dependencies=T)
}

# load required packages
library(plyr) # data management
library(readxl) # for reading .xls/xlsx-files into R

# import function for sample information
import.smpl <- function(file.cult, sheet, keep.col=c('Serial-number', 'date collected', 'detailed location', 'species', 'sample type', 'ESBL-A')){
  print(sheet)
  headers <- unlist(as.vector(read_excel(file.cult, sheet=sheet, n_max=1, col_names=F, col_types='text')[1,]))
  samples <- read_excel(file.cult, sheet=sheet, skip=5, col_names=F, col_types='text', trim_ws=F)
  samples <- samples[,which(headers %in% keep.col)]
  names(samples) <- c('sample_id', 'sample_date', 'sample_site', 'sample_origin', 'sample_location', 'ESBL-A')
  samples <- samples[!is.na(samples$sample_id),]
  samples$sample_opportunity <- sheet
  return(samples)
}

# vector of cultivation result files
res.files <- list.files('', full.names=T, pattern='.xlsx')

# import all sheets within the vector of files, excluding only the ference sheet that does not contain data
all.data <- rbind.fill(lapply(res.files, function(f){
  print(f)
  sheets <- excel_sheets(f)
  sheets <- sheets[-which(sheets=='Ab Reference sheet')]
  d <- rbind.fill(lapply(sheets, import.smpl, file.cult=f))
  return(d)
}))


# if you want to format dates as Date-objects, note that the origin for dates in excel appears to be 1899-12-30 - see below
# whether the following line of code works will depend on the format of the date in the result files
all.data$date <- as.Date(as.numeric(all.data$sample_date), origin='1899-12-30')
