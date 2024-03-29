###############################
# title : Figure 3g in Abnormal accumulation of extracellular vesicles in hippocampal dystrophic axons and regulation by the primary cilia gene intraflagellar transport homolog 88 in Alzheimer’s disease
# author : Jaemyung, Jang (piloter2@kbri.re.kr)
# kenel : R 4.2.0
# Date : May 18, 2023
###############################

source(paste0(getwd(),"/code/proteomic_function.R"))

k = 2 # select number - 220107_APPc_abcam_SepHippoSyn_Proteins.txt ; from Dr. Baek

# Figure 3g - Venn diagram
counts<-lapply(datafromPD, function(data){
  cts<-data %>% base::subset(select = c("Accession",grep("^Abundances.", colnames(data), value=TRUE))) 
  return(cts)
})

Syn.count<-list("C57BL/6"=subset(counts[[k]], select = c("Accession",grep(".Control.WT$", colnames(counts[[k]]), value=TRUE))),
                "5xFAD"=subset(counts[[k]], select = c("Accession",grep(".Sample.Tg6799$", colnames(counts[[k]]), value=TRUE))))
                
Syn.types <- list( "WT" = Syn.count[['C57BL/6']]$Accession[unique(c(unlist(apply(Syn.count[['C57BL/6']], 2, function(x) which(!is.na(x)))[-1])))],
                   "5xFAD"=Syn.count[['5xFAD']]$Accession[unique(c(unlist(apply(Syn.count[['5xFAD']], 2, function(x) which(!is.na(x)))[-1])))])
                
require(ggvenn)

p <- ggvenn(Syn.types,
       fill_alpha = 0.1,
       fill_color =c( '#fde725ff', '#21908dff') #"#440154ff"
      )+ ggtitle("APP(abcam)-IP of septum-hippocampus synaptosome") 

ggsave(paste0(path, "/results/Figure3g_venn.pdf"), p, width = 5, height =5 , units = "in", device = "pdf")

# Matching gene symbols in the PD software and the unitProt
cts <- lapply(datafromPD, function(dt){
  dts <- dt %>% subset( select = c("Accession", grep("^Abundance.Ratio.log2.|^Abundance.Ratio.Adj.P_Value.|^Abundance.Ratio.Weight.", colnames(dt), value=TRUE)))
  return(dts)
})

resultUNIPROT <- pbapply::pblapply(datafromPD[[k]]$Accession,function(ids){ #[-grep("ProteinCenter:sp_incl_isoforms",datafromPD$Accession)]
    data <- uniprot_mapping(ids)
    content(data, as= "text", encoding = 'UTF-8')
    res<-unlist(str_split(unlist(str_split(unlist(str_split(content(data, as= "text", encoding = 'UTF-8'),"\\t")),"\\n"))[4]," "))[1]
    return(res)
  })
  
  results_from_uniprot <- data.frame('Accession' = datafromPD[[k]]$Accession,  #[-grep("ProteinCenter:sp_incl_isoforms",datafromPD$Accession)]
                                     'Gene.names' = unlist(resultUNIPROT))
  require(stringi)
  result<-lapply(datafromPD[[k]]$Description, function(desc){ #[-grep("ProteinCenter:sp_incl_isoforms",datafromPD$Accession)]
    str_sub(desc, 
            str_locate(desc,pattern = "GN=")[2]+1,
            str_locate(desc,pattern = "PE=")[1]-2)
  })
  
  result<-make.unique(unlist(result))
  
  results_from_datasheet <- data.frame('Protein.IDs' = datafromPD[[k]]$"Accession",  'Gene.names' = result) #[-grep("ProteinCenter:sp_incl_isoforms",datafromPD$Accession)]
  results_from_uniprot <- na.omit(results_from_uniprot)
  results_from_uniprot$Gene.names[which(results_from_uniprot$Gene.names=="")] <- result[which(datafromPD[[k]]$"Accession" %in% results_from_uniprot$Accession[which(results_from_uniprot$Gene.names=="")])]
 
data.merge<-merge(cts[[k]], results_from_uniprot, by = "Accession")
data_unique <- make_unique(data.merge, "Gene.names", "Accession", delim = ";")
colnames(data_unique)

data.merged <- subset(data_unique, select = c('Gene.names','Accession','Abundance.Ratio.log2.Sample..Control','Abundance.Ratio.Adj.P_Value.Sample..Control','Abundance.Ratio.Weight.Sample..Control'))
colnames(data.merged) <- c("Gene.names","Accession", "Abundance.Ratio.log2","Abundance.Ratio.Adj.P_Value","Abundance.Ratio.Weight")

# Figure 3g - Volcano plot
selAccession <- data.merged %>% 
           dplyr::filter(Abundance.Ratio.log2 < -0.25 & Abundance.Ratio.Adj.P_Value < 10e-1) %>% dplyr::pull(Accession)


selAccession.1 <- data.merged %>% 
           dplyr::filter(Abundance.Ratio.log2 > 0.25 & Abundance.Ratio.Adj.P_Value < 10e-1) %>% dplyr::pull(Accession)


# options(repr.plot.width = 16, repr.plot.height = 16)
keyvals.colour <- ifelse(data.merged$Accession %in% selAccession, 'royalblue', 'lightgrey')
keyvals.colour[data.merged$Accession %in% selAccession.1] <- "darkred"
keyvals.colour[is.na(keyvals.colour)] <- 'lightgrey'

names(keyvals.colour)[keyvals.colour == 'lightgrey'] <- 'ns'
names(keyvals.colour)[keyvals.colour == 'royalblue'] <- 'down-regulated'
names(keyvals.colour)[keyvals.colour == 'darkred'] <- 'up-regulated'

data.merged$group <- names(keyvals.colour)

require(EnhancedVolcano)

pv<-EnhancedVolcano(data.merged,
    lab = data.merged$Gene.names,
    x = 'Abundance.Ratio.log2',
    y = 'Abundance.Ratio.Adj.P_Value',
    title = "Proteomic analysis for APP(abcam)-IP of septum-hippocampus synaptosome",
    subtitle = "Differential expression between 5xFAD and C56BL/6",
    xlab = bquote(~Log[2]~ 'fold change'),
    ylab = bquote(~Log[10]~ 'adjust P-value'),
    pCutoff = 0.1,
    FCcutoff = 0.25,
    pointSize = 0.1,
    labSize = 4.0,
    colAlpha = 1,
    border = 'full',
    borderWidth = 1.0,
    legendPosition = 'none',
    legendLabSize = 4,
    legendIconSize = 2.4,
    drawConnectors = TRUE,
    widthConnectors = 0.001)

ggsave(paste0(path, "/results/Figure3g_vc.pdf"), pv, width = 7.5, height =7.5 , units = "in", device = "pdf")


# Extended Data Supplementary Figure 3 - Gene Set Enrichment Analysis
require(enrichR)

  listEnrichrSites()
  websiteLive <- TRUE
  dbs <- listEnrichrDbs()

  dbs <- "SynGO_2022"
  if (websiteLive) {
      enriched <- enrichr(data.merged$Gene.names[data.merged$Accession %in% selAccession], dbs)
      enriched.1 <- enrichr(data.merged$Gene.names[data.merged$Accession %in% selAccession.1], dbs)
  }

enriches.syngo<-enriched[[2]][grep("Vesicle|Presynap",enriched[[2]]$Term),]
enriches.syngo1<-enriched.1[[2]][grep("Vesicle|Presynap",enriched.1[[2]]$Term),]

pdf(paste0(path,"/results/FigureS3_SynGO-upregulated.pdf"), width=10, height=5)
    plotEnrich(enriches.syngo, showTerms = 25, numChar = 85, y = "Count", orderBy = "P.value")
dev.off()

pdf(paste0(path,"/results/FigureS3_SynGO-downregulated.pdf"), width=10, height=6)
 plotEnrich(enriches.syngo1, showTerms = 25, numChar = 85, y = "Count", orderBy = "P.value")
 dev.off()
