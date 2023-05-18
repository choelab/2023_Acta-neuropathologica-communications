DelayedArray:::set_verbose_block_processing(TRUE)

# Passing a higher value will make some computations faster but use more memory. Adjust with caution!
options(DelayedArray.block.size=1000e6)

base::gc()
base::rm(list = ls()) # Clear the environment
options(warn=-1)

load.lib <- c("DEP","ggplot2","ggrepel","ggpubr","scales","dplyr","tidyverse","devtools","data.table","httr","stringr")

if (!requireNamespace(load.lib, quietly = TRUE))
    install.packages(load.lib)

sapply(load.lib,require,character.only = TRUE)

require.lib <- c("pbapply","reshape2","stringr","ggvenn","stringi","BiocManager","RColorBrewer","UniprotR","enrichR","xlsx")

if (!requireNamespace(require.lib, quietly = TRUE))
    install.packages(require.lib)

bioconductor.lib <- c("EnhancedVolcano","STRINGdb")

if (!requireNamespace(bioconductor.lib, quietly = TRUE))
    BiocManager::install(bioconductor.lib)

if (!requireNamespace("nichenetr", quietly = TRUE))
    devtools::install_github("saeyslab/nichenetr")

httr::set_config(config(ssl_verifypeer = 0L))

#################################################################
uniprot_mapping <- function(ids) {
  uri <- 'https://rest.uniprot.org/uniprotkb/search?size=1&&query='
  #uri <- 'https://www.uniprot.org/uniprotkb/search?query='
  #idStr <- paste(ids, collapse="+or+")
  format <- '&fields=accession%2Cgene_names&format=tsv'
  #fullUri <- paste0(uri,idStr,format)
  fullUri <- paste0(uri,ids,format)
  #dat <- read.delim(fullUri)
  dat <- GET(fullUri,
             accept_json(),
             add_headers(Accept = 'application/json'))
  return(dat)
}

require(pbapply)
pboptions(type="txt", char=":")

require(reshape2)

get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

#################################################################


path <- file.path(get_wd())
rawPD<-list.files(paste0(path,"/rawdata"), pattern = "_Proteins.txt$")


# Load data exported from the PD software
datafromPD<-lapply(rawPD, function(file){
  return(fread(paste0(path,"/",file)))
})

require(stringr)
for(i in seq_len(length(datafromPD))){
  colnames(datafromPD[[i]])<-str_replace_all(colnames(datafromPD[[i]]),"[ ]",".")
  colnames(datafromPD[[i]])<-str_replace_all(colnames(datafromPD[[i]]),"[-]","_")
}

counts<-lapply(datafromPD, function(data){
  cts<-data %>% base::subset(select = c("Accession",grep("^Abundances.", colnames(data), value=TRUE))) 
  return(cts)
})
