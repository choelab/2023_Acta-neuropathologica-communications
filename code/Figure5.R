###############################
# title : Figure 5 in Abnormal accumulation of extracellular vesicles in hippocampal dystrophic axons and regulation by the primary cilia gene intraflagellar transport homolog 88 in Alzheimer’s disease
# author : Jaemyung, Jang (piloter2@kbri.re.kr)
# kenel : R 4.3.0
# Date : May 18, 2023
###############################

source(paste0(getwd(),"/code/proteomic_function.R"))

k = 4 # select number - 230316_220315_N_Lysate_siIft_Proteins.txt ; from Dr. Yeo

# Figure 5 - Venn Diagram

new.count<-list("Expr"=subset(counts[[k]], select = c("Accession",grep("Sample.siIft88_Abeta$", colnames(counts[[k]]), value=TRUE))),
                "Ctrl"=subset(counts[[k]], select = c("Accession",grep("Sample.Con_Abeta$", colnames(counts[[k]]), value=TRUE))))
                
new.types <- list( "Expr" = new.count[['Expr']]$Accession[unique(c(unlist(apply(new.count[['Expr']], 2, function(x) which(!is.na(x)))[-1])))],
                   "Ctrl"=new.count[['Ctrl']]$Accession[unique(c(unlist(apply(new.count[['Ctrl']], 2, function(x) which(!is.na(x)))[-1])))])               

p <- ggvenn(new.types,
       fill_alpha = 0.1,
       fill_color =c( '#fde725ff', '#21908dff') ) + #"#440154ff"
        ggtitle("Exosomes from normal neuron treated with siIft88 and veh")

ggsave(paste0(path, "/results/Figure5_venn.pdf"), p, width = 5, height =5 , units = "in", device = "pdf")

# Matching gene symbols in the PD software and the unitProt

cts <- lapply(datafromPD, function(dt){
  dts <- dt %>% subset( select = c("Accession", grep("^Abundance.Ratio\\.log2|^Abundance.Ratio.Adj.P_Value|^Abundance.Ratio.Weight", colnames(dt), value=TRUE)))
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

data.merged <- subset(data_unique, select = c('Gene.names','Accession','Abundance.Ratio.log2.siIft88_Abeta..siIft88','Abundance.Ratio.Adj.P_Value.siIft88_Abeta..siIft88','Abundance.Ratio.Weight.siIft88_Abeta..siIft88'))
colnames(data.merged) <- c("Gene.names","Accession", "Abundance.Ratio.log2","Abundance.Ratio.Adj.P_Value","Abundance.Ratio.Weight")

# Figure 5 - Scatter Plot

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
  geom_text_repel(data = data.merged[(abs(data.merged$Abundance.Ratio.log2) > 0.25) 
                                      &  ( data.merged$Abundance.Ratio.Adj.P_Value < 0.5 ) ,],
                  aes(label=Gene.names), 
                  arrow = arrow(length = unit(0.0001, "npc")),
                  colour = "black",max.overlaps = 25, min.segment.length=0.1)+
  scale_colour_manual(values = c("royalblue","lightgrey","darkred"))+ #brewer.pal(3,"Paired")
  
  theme_bw()+
  theme(legend.position=c(0.9, 0.2),
        axis.text.x = element_text(size = 12, hjust = 0, colour = "black"),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_blank(),
        legend.background = element_blank()) + 
   theme(legend.position = "none")

ggsave(paste0(path, "/results/Figure5_scatter.pdf"), pscatter, width = 5, height = 5 , units = "in", device = "pdf")

# Figure 5 - Gene Set Enrichment Analysis

require(enrichR)

  listEnrichrSites()
  websiteLive <- TRUE
  dbs <- listEnrichrDbs()

  dbs <- "GO_Biological_Process_2021"
  
  if (websiteLive) {
      enriched <- enrichr(data.merged$Gene.names[data.merged$Accession %in% selAccession], dbs)
      enriched.1 <- enrichr(data.merged$Gene.names[data.merged$Accession %in% selAccession.1], dbs)
  }

  enriches<-enriched[grep("cili|endocy|phago|exocy",enriched$Term),]  %>% dplyr::filter(P.value < 0.1) 
  enriches.1<-enriched.1[grep("cili|endocy|phago|exocy",enriched.1$Term),]  %>% dplyr::filter(P.value < 0.1)

  enrichDF<-rbind(cbind(enriches, group = "negative"),cbind(enriches.1, group = "positive"))
  enrichDF$NumGenes<-as.numeric(str_split(enrichDF$Overlap, pattern = "/", simplify = TRUE)[,1])
  enrichDF <- enrichDF %>% dplyr::filter(NumGenes > 1)
  enrichDF <- enrichDF[-grep("neuromuscular ", enrichDF$Term),]

selgene<-nichenetr::convert_human_to_mouse_symbols(unique(unlist(str_split(enrichDF$Genes,";")))) %>% na.omit()

pPathway<-cowplot::plot_grid(
    ggplot(enrichDF %>% dplyr::filter(group == "negative"), aes(x=Term, y=NumGenes, fill = Adjusted.P.value)) + 
      geom_bar(stat = "identity", position = position_dodge(width = 0.1), width=0.5)+
      #scale_fill_gradient(low = "yellow", high = "red", na.value = NA)+
      scale_fill_distiller(palette = "Blues")+
      theme_bw() + 
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
      scale_fill_distiller(palette = "OrRd")+
      theme_bw() + 
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

ggsave(paste0(path, "/results/Figure5_pathway_summary.pdf"), pPathway, width = 10, height = 10 , units = "in", device = "pdf")
