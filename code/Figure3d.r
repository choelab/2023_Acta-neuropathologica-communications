source(paste0(getwd(),"/code/proteomic_function.R"))

path <- file.path(get_wd())
rawPD<-list.files(paste0(path,"/rawdata_from_PD"), pattern = "_Proteins.txt$")

