###############################
# title : Figure 2
# author : Jaemyung, Jang (piloter2@kbir.re.kr)
# Date : March 18, 2023
###############################

source(paste0(getwd(),"/proteomic_function.R"))

path <- file.path("/","mnt","hadoop_tmp","working","190722_MS hipp Synaptosome")
rawPD<-list.files(path, pattern = "_Proteins.txt$")


# Load data
datafromPD<-lapply(rawPD, function(file){
  return(fread(paste0(path,"/",file)))
})

require(stringr)
for(i in seq_len(length(datafromPD))){
  colnames(datafromPD[[i]])<-str_replace_all(colnames(datafromPD[[i]]),"[ ]",".")
  colnames(datafromPD[[i]])<-str_replace_all(colnames(datafromPD[[i]]),"[-]","_")
}

# Venn diagram
counts<-lapply(datafromPD, function(data){
  cts<-data %>% base::subset(select = c("Accession",grep("^Abundances.", colnames(data), value=TRUE))) 
  return(cts)
})

k=1
Syn.count<-list("WT"=subset(counts[[k]], select = c("Accession",grep("WT$", colnames(counts[[k]]), value=TRUE))),
                "Tg6799"=subset(counts[[k]], select = c("Accession",grep("Tg6799$", colnames(counts[[k]]), value=TRUE))))
                

