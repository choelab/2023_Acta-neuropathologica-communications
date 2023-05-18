source("~/R/proteomic_function.R")
path <- file.path("/","mnt","hadoop_tmp","working","220107_APPc_abcam_BSB_SepHippoSyn")
rawPD<-list.files(path, pattern = "_Proteins.txt$")


datafromPD<-lapply(rawPD, function(file){
  return(fread(paste0(path,"/",file)))
})


require(stringr)
for(i in seq_len(length(datafromPD))){
  colnames(datafromPD[[i]])<-str_replace_all(colnames(datafromPD[[i]]),"[ ]",".")
  colnames(datafromPD[[i]])<-str_replace_all(colnames(datafromPD[[i]]),"[-]","_")
}
