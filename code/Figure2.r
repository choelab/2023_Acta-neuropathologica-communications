###############################
# title : Figure 2
# author : Jaemyung, Jang (piloter2@kbri.re.kr)
# kenel : R 4.3.0
# Date : March 18, 2023
###############################

source(paste0(getwd(),"/code/proteomic_function.R"))

path <- file.path(get_wd(),"/rawdata_from_PD")
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
                
Syn.types <- list( "C57BL/6" = Syn.count[['WT']]$Accession[unique(c(unlist(apply(Syn.count[['WT']], 2, function(x) which(!is.na(x)))[-1])))],
                   "5xFAD"=Syn.count[['Tg6799']]$Accession[unique(c(unlist(apply(Syn.count[['Tg6799']], 2, function(x) which(!is.na(x)))[-1])))])
                
require(ggvenn)

p <- ggvenn(Syn.types,
       fill_alpha = 0.1,
       fill_color =c( '#fde725ff', '#21908dff') ) + #"#440154ff"
        ggtitle("Hippocampal synaptosome between C57BL/6 and 5xFAD")

ggsave(paste0(path, "/results/Figure_venn.pdf"), p, width = 5, height =5 , units = "in", device = "pdf")

cts <- lapply(datafromPD, function(dt){
  dts <- dt %>% subset( select = c("Accession", grep("^Abundance.Ratio\\.\\(log2\\)|^Abundance.Ratio.Adj..P_Value|^Abundance.Ratio.Weight", colnames(dt), value=TRUE)))
  #dts <- dt %>% subset( select = c("Accession", grep("^Abundance.Ratio\\.\\(log2\\)|^Abundance.Ratio.P_Value|^Abundance.Ratio.Weight", colnames(dt), value=TRUE)))
  return(dts)
})

# Gene symbol

  k = 1 # select number number from "print(rawPD)"
  
  resultUNIPROT <- pbapply::pblapply(datafromPD[[k]]$Accession,function(ids){ #[-grep("ProteinCenter:sp_incl_isoforms",datafromPD$Accession)]
    data <- uniprot_mapping(ids)
    content(data, as= "text", encoding = 'UTF-8')
    res<-unlist(str_split(unlist(str_split(unlist(str_split(content(data, as= "text", encoding = 'UTF-8'),"\\t")),"\\n"))[4]," "))[1]
    return(res)
  })
  
  results_from_uniprot <- data.frame('Accession' = datafromPD[[k]]$Accession,  #[-grep("ProteinCenter:sp_incl_isoforms",datafromPD$Accession)]
                                     'Gene.names' = unlist(resultUNIPROT))
  require(stringi)
  # str_sub(datafromPD[[k]]$Description[1], #[-grep("ProteinCenter:sp_incl_isoforms",datafromPD$Accession)]
  #         str_locate(datafromPD[[k]]$Description[1],pattern = "GN=")[2]+1, #[-grep("ProteinCenter:sp_incl_isoforms",datafromPD$Accession)]
  #         str_locate(datafromPD[[k]]$Description[1],pattern = "PE=")[1]-2) #[-grep("ProteinCenter:sp_incl_isoforms",datafromPD$Accession)]
  # seq_len(length(dataPD$Description))
  result<-lapply(datafromPD[[k]]$Description, function(desc){ #[-grep("ProteinCenter:sp_incl_isoforms",datafromPD$Accession)]
    # str_sub(grep("GN=",unlist(str_split(desc, " ")), value = TRUE),start = 4L)
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

data.merged <- subset(data_unique, select = c('Gene.names','Accession','Abundance.Ratio.(log2):.(Sample)./.(Control)','Abundance.Ratio.Adj..P_Value:.(Sample)./.(Control)','Abundance.Ratio.Weight:.(Sample)./.(Control)'))
colnames(data.merged) <- c("Gene.names","Accession", "Abundance.Ratio.log2","Abundance.Ratio.Adj.P_Value","Abundance.Ratio.Weight")

# Scatter Plot
require(RColorBrewer)
require(ggrepel)

    selAccession <- data.merged %>% 
            dplyr::filter(Abundance.Ratio.log2 < -0.25 & Abundance.Ratio.Adj.P_Value < .1) %>% dplyr::pull(Accession)



    selAccession.1 <- data.merged %>% 
            dplyr::filter(Abundance.Ratio.log2 > 0.25 & Abundance.Ratio.Adj.P_Value < .1) %>% dplyr::pull(Accession)


    # options(repr.plot.width = 16, repr.plot.height = 16)
    keyvals.colour <- ifelse(data.merged$Accession %in% selAccession, 'royalblue', 'lightgrey')
    keyvals.colour[data.merged$Accession %in% selAccession.1] <- "darkred"
    keyvals.colour[is.na(keyvals.colour)] <- 'lightgrey'

    names(keyvals.colour)[keyvals.colour == 'lightgrey'] <- 'ns'
    names(keyvals.colour)[keyvals.colour == 'royalblue'] <- 'down-regulated'
    names(keyvals.colour)[keyvals.colour == 'darkred'] <- 'up-regulated'

    data.merged$group <- names(keyvals.colour)

pscatter <- ggplot(data.merged, aes(x=Abundance.Ratio.log2, y=Abundance.Ratio.Weight, color =group)) + 
  geom_point() +
  geom_text_repel(data = data.merged[ which(data.merged$Gene.names %in% selgene),],
                    #   data.merged [(abs(data.merged$Abundance.Ratio.log2) > 0.25) 
                    #                   &  ( data.merged$Abundance.Ratio.Adj.P_Value < 0.1 )],
                  aes(label=Gene.names), 
                  arrow = arrow(length = unit(0.000001, "npc")),
                  colour = "black",
                  max.overlaps = 25, min.segment.length=0.1)+
  #xlab("Log2 Fold change between single-cell and bulk-cell analysis")+
  #ylab("Log2 Fold change between 100-cell and bulk-cell analysis")+
  scale_colour_manual(values = c("royalblue","lightgrey","darkred"))+ #brewer.pal(3,"Paired")
  
  theme_bw()+
  theme(legend.position=c(0.9, 0.2),
        axis.text.x = element_text(size = 12, hjust = 0, colour = "black"),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_blank(),
        legend.background = element_blank()) + 
   theme(legend.position = "none")

ggsave(paste0(path, "/results/Figure2_scatter_renew.pdf"), pscatter, width = 5, height = 5 , units = "in", device = "pdf")

# Pathway analysis

require(enrichR)

  listEnrichrSites()
  websiteLive <- TRUE
  dbs <- listEnrichrDbs()

  dbs <- "GO_Biological_Process_2021"
  if (websiteLive) {
      enriched <- enrichr(data.merged$Gene.names[data.merged$Accession %in% selAccession], dbs)
      enriched.1 <- enrichr(data.merged$Gene.names[data.merged$Accession %in% selAccession.1], dbs)
  }

  enriches<-enriched[grep("synaptic vesicle|synaptic transmission|signal release|neuron projection morphogenesis",enriched$Term),]  %>% dplyr::filter(P.value < 0.1) 
  enriches.1<-enriched.1[grep("synaptic vesicle|synaptic transmission|signal release|neuron projection morphogenesis",enriched.1$Term),]  %>% dplyr::filter(P.value < 0.1)

  enrichDF<-rbind(cbind(enriches, group = "negative"),cbind(enriches.1, group = "positive"))
  enrichDF$NumGenes<-as.numeric(str_split(enrichDF$Overlap, pattern = "/", simplify = TRUE)[,1])
  enrichDF <- enrichDF %>% dplyr::filter(NumGenes > 1)
  enrichDF <- enrichDF[-grep("neuromuscular ", enrichDF$Term),]

selgene<-nichenetr::convert_human_to_mouse_symbols(unique(unlist(str_split(enrichDF$Genes,";")))) %>% na.omit()

p<-cowplot::plot_grid(
    ggplot(enrichDF %>% dplyr::filter(group == "negative"), aes(x=Term, y=NumGenes, fill = Adjusted.P.value)) + 
    geom_bar(stat = "identity", position = position_dodge(width = 0.1), width=0.5)+
    #scale_fill_gradient(low = "yellow", high = "red", na.value = NA)+
    scale_fill_distiller(palette = "Blues")+
    theme_bw() + 
#    facet_grid( group ~ .) + 
    #  ylim(0, 3) + 
    coord_flip() +
     theme(axis.text.x = element_text(face="bold", #color="#993333", , angle=45
                           size=12),
          axis.text.y = element_text(face="bold", #color="#993333", , angle=45
                            size=10),
                        #    size=13.5),
                           legend.background=element_rect(fill = alpha("white", 0)),
          legend.key=element_rect(fill = alpha("white", .5)),
                           legend.position="right")        ,

    ggplot(enrichDF %>% dplyr::filter(group == "positive"), aes(x=Term, y=NumGenes, fill = Adjusted.P.value)) + 
    geom_bar(stat = "identity", position = position_dodge(width = 0.1), width=0.5)+
#   scale_fill_gradient2(low = "grey", high = "royalblue", na.value = NA)+
    scale_fill_distiller(palette = "OrRd")+
    theme_bw() + 
#    facet_grid( group ~ .) + 
    #  ylim(0, 11) + 
    coord_flip() +   
     theme(axis.text.x = element_text(face="bold", #color="#993333", , angle=45
                           size=12),
          axis.text.y = element_text(face="bold", #color="#993333", , angle=45
                        #    size=9),
                           size=10),
          legend.background=element_rect(fill = alpha("white", 0)),
          legend.key=element_rect(fill = alpha("white", .5)),
                           legend.position="right"),
        rel_heights = c(dim(enrichDF %>% dplyr::filter(group == "negative"))[1],dim(enrichDF %>% dplyr::filter(group == "positive"))[1]),
        ncol=1)

print(p)

ggsave(paste0(path, "/results/Figure2_pathway_summary.pdf"), p, width = 10, height = 10 , units = "in", device = "pdf")
