###############################
# title : Figure 6 in Abnormal accumulation of extracellular vesicles in hippocampal dystrophic axons and regulation by the primary cilia gene intraflagellar transport homolog 88 in Alzheimer’s disease
# author : Jaemyung, Jang (piloter2@kbri.re.kr)
# kenel : R 4.3.0
# Date : March 18, 2023
###############################

DelayedArray:::set_verbose_block_processing(TRUE)

# Passing a higher value will make some computations faster but use more memory. Adjust with caution!
options(DelayedArray.block.size=1000e6)
options(stringsAsFactors = FALSE)

base::gc()
base::rm(list = ls()) # Clear the environment
options(warn=-1)

path_base <- getwd()
source(paste0(path_base,"/code/single_cell_function.R"))

Idents(object = hipp_05m) <- sprintf("%02d",hipp_05m@meta.data$`SCT_snn_res.0.3`)
hipp_05m@meta.data$cell_type_age[which(Idents(hipp_05m)=="01")] <- "OLG"
hipp_05m@meta.data$cell_type_age[which(Idents(hipp_05m)=="02")] <- "ASC"
hipp_05m@meta.data$cell_type_age[which(Idents(hipp_05m)=="03")] <- "microglia"
hipp_05m@meta.data$cell_type_age[which(Idents(hipp_05m)=="04")] <- "Neuron A"
hipp_05m@meta.data$cell_type_age[which(Idents(hipp_05m)=="05")] <- "Neuron B"
hipp_05m@meta.data$cell_type_age[which(Idents(hipp_05m)=="06")] <- "Endo"
hipp_05m@meta.data$cell_type_age[which(Idents(hipp_05m)=="07")] <- "OPC"
hipp_05m@meta.data$cell_type_age[which(Idents(hipp_05m)=="08")] <- "Fib"
hipp_05m@meta.data$cell_type_age[which(Idents(hipp_05m)=="09")] <- "NEUT"
hipp_05m@meta.data$cell_type_age[which(Idents(hipp_05m)=="10")] <- "OLG"
hipp_05m@meta.data$cell_type_age[which(Idents(hipp_05m)=="11")] <- "EPC"
#hipp_05m@meta.data$cell_type_age[which(Idents(hipp_05m)=="12")] <- "VSMC" # Since the total number was less than 100, we excluded them from further analysis

#Extended Data Suppmentary Figure 7a
pS1a<-DimPlot(subset(hipp_05m, subset = genetype == "C57BL/6", ident =names(which(table(Idents(hipp_05m))>100))), group.by = "cell_type_age", label =T, label.size = 2.5, pt.size = 0.01,raster=FALSE,
             cols =  c(brewer.pal(5,"Paired"),brewer.pal(7,"Pastel1"))) + 
  theme(legend.position = "none",  axis.title = element_blank(), 
        plot.title = element_text(size = 12, face = "bold",hjust = 0), axis.line = element_blank(), 
        axis.text = element_blank(),  axis.ticks = element_blank()) + 
  ggtitle(paste0("WT hippocampus : ",comma(length(colnames(subset(hipp_05m, subset = genetype == "C57BL/6", ident =names(which(table(Idents(hipp_05m))>100))))),format = "d")," cells")) +
  NoLegend() + NoAxes()

pS1b<-DimPlot(subset(hipp_05m, subset = genetype == "5xFAD (C57BL/6)", ident = names(which(table(Idents(hipp_05m))>100))), group.by = "cell_type_age", label =T, label.size = 2.5, pt.size = 0.01,raster=FALSE,
             cols =  c(brewer.pal(5,"Paired"),brewer.pal(7,"Pastel1"))) + 
  theme(legend.position = "none",  axis.title = element_blank(), 
        plot.title = element_text(size = 12, face = "bold",hjust = 0), axis.line = element_blank(), 
        axis.text = element_blank(),  axis.ticks = element_blank()) + 
  ggtitle(paste0("Tg6799 hippocampus : ",comma(length(colnames(subset(hipp_05m, subset = genetype == "5xFAD (C57BL/6)", ident =names(which(table(Idents(hipp_05m))>100))))),format = "d")," cells")) +
  NoLegend() + NoAxes()

Idents(hipp_05m) <- "cell_type_age"
hipp_05m.markers <- FindAllMarkers(hipp_05m,
                                    only.pos = TRUE, logfc.threshold = 0.25, 
                                    assay = "RNA")

hipp_05m.markers.top20 <- hipp_05m.markers %>% group_by(cluster) #%>%  filter (p_val_adj < 0.1) %>% top_n(n = 5, wt = avg_log2FC) 
hipp_05m.markers.arrange <- rbind(hipp_05m.markers.top20 %>% dplyr::filter(cluster == "ASC"),
      hipp_05m.markers.top20 %>% dplyr::filter(cluster == "Endo"),
      hipp_05m.markers.top20 %>% dplyr::filter(cluster == "EPC"),
      hipp_05m.markers.top20 %>% dplyr::filter(cluster == "Fib"),
      hipp_05m.markers.top20 %>% dplyr::filter(cluster == "MG"),
      hipp_05m.markers.top20 %>% dplyr::filter(cluster == "Neuron A"),
      hipp_05m.markers.top20 %>% dplyr::filter(cluster == "Neuron B"),
      hipp_05m.markers.top20 %>% dplyr::filter(cluster == "NEUT"),
      hipp_05m.markers.top20 %>% dplyr::filter(cluster == "OLG"),
      hipp_05m.markers.top20 %>% dplyr::filter(cluster == "OPC"))
hipp_05m.markers.filter <- setdiff(hipp_05m.markers.arrange$gene, grep("^Rps|^Rpl|^mt.|^Hba\\.|^Hbb\\.|$Rik|^Gm",hipp_05m.markers.top20$gene, value = TRUE))

#Extended Data Suppmentary Figure 7b
pS2<-DoHeatmap(subset(hipp_05m,downsample = 1000), features = hipp_05m.markers.filter,
          group.by = "cell_type_age",
          # group.colors = c(brewer.pal(5,"Paired"),brewer.pal(7,"Pastel1")),
          size = 3, assay = "SCT",angle = 45)+
  scale_fill_gradientn(colors = c("royalblue","white", "darkred"))

#Extended Data Suppmentary Figure 7c
require(Nebulosa)
pS3<-plot_density(hipp_05m, features = c("Meg3","Atp1b1"), reduction = "umap" ,
               pal = "inferno", 
               method = c("wkde"), size = 0.025, adjust=0.25)

# subsetting neuron subtype
neuron <- subset(hipp_05m, ident = c("04","05"))
   
neuron <- neuron %>%
     RunICA(verbose = FALSE, features = unique(c("Prox1","Dkk3","Calb2","Ociad2","Cacng5","Fibcd1",
                                                 "Pdzd2","Tox3","Iyd","Coch","Wfs1","Dcn",
                                                 "Fgf2","Maob","Plcb4","Scgn","Vcan","Ntsr2", "Srgap2", #CA2
                                                 "Chga","Cntn6","Drd2","Ephb6","Gal","Glp1r","Grm8","Hpcal1","Nmb","Prrx1","Rac3","Vat1",
                                                 "Lct","Trhr","Gsg1l","Spata13","Stra6","Cpne7","Grp","Nr2f2",
                                                 "Csf2rb2","Nmb","Thbs2","Tm4sf1",                            #Mossy cells
                                                 "Cadm2","Chgb","Dagla","Epha5","Epha7","Inhba","Ing2","Ppp1r1a","Rimbp2","Mgll","Sertad4","Syt7",
                                                 "Atp2b4","Cadm1","Crym","Cpne7","Crtac1","Efnb2","Etv1","Fam70a","Gfra1","Gpr123","Hasp1","Kcnd3","Magee2","Necab1","Negr1","Nr2f2","Nrbp2","Odz4","Rab26","Timp2","Wdr6","Resp18","Ypel2",
                                                 "Gad1","Gad2","Slc17a7","Slc32a1","Lamp5","Ndnf","Sncg","Vip","Sst","Chodl","Pvalb","Cux2","Rorb","Fezf2","Sulf1","Foxp2","Nxph4",
                                                 "Wfs1","Dsp","Gad2","Pou3f1","Thsd4","Dcn","Ly6g6e","C1ql2","Ddit4l","Ndnf",
                                                 "Atpfou2b1","Hpca","Spink8","Arpc2","Ptn",
                                                 "Ly6e","Atp1b1","Slc24a2","Prkca","Cpe",
                                                 "Dcn","Nnat","Timp2","Cpne7","Ndufc2",
                                                 "Serpini1","Ly6g6e","Cbin2","Fezf2","Nxph3",
                                                 "Cplx2","Ncdn","Olfm1","Synpr","Sem5a",
                                                 "Gad2","Gad1","Slc32a1","Nap1l5","Slc6a1",
                                                 "Tshz2","Arpp21","Ddit4l","Nrep","Meis2",
                                                 "Ndnf","Diablo","Cacna2d2","Marcks")), nics = 10) %>% #, approx=FALSE
     RunUMAP(reduction = "ica", umap.method = "umap-learn", dims=1:10) %>%
     FindNeighbors(reduction = "ica", dims=1:10) %>%
     FindClusters(resolution =seq(0.1,1,0.1), algorithm = 2)

Idents(object = neuron) <- neuron@meta.data$'SCT_snn_res.0.2' 
neuron@meta.data$cell_type_age[which(Idents(neuron)=='0')] <- "Neuron:sub-u"
neuron@meta.data$cell_type_age[which(Idents(neuron)=='1')] <- "Neuron:sub-1"
neuron@meta.data$cell_type_age[which(Idents(neuron)=='2')] <- "Neuron:sub-2"
neuron@meta.data$cell_type_age[which(Idents(neuron)=='3')] <- "Neuron:sub-3"
neuron@meta.data$cell_type_age[which(Idents(neuron)=='4')] <- "Neuron:sub-4"
neuron@meta.data$cell_type_age[which(Idents(neuron)=='5')] <- "Neuron:sub-5"
neuron@meta.data$cell_type_age[which(Idents(neuron)=='6')] <- "Neuron:sub-6"
neuron@meta.data$cell_type_age[which(Idents(neuron)=='7')] <- "Neuron:sub-7"
neuron@meta.data$cell_type_age[which(Idents(neuron)=='8')] <- "Neuron:sub-8"

neuron.markers <- FindAllMarkers(neuron,
                                 only.pos = TRUE, logfc.threshold = 0.25,
                                 assay = "RNA")

neuron.markers.top20 <- neuron.markers %>% group_by(cluster) %>%  filter (p_val_adj <1) %>% top_n(n = 20, wt = avg_log2FC)

# Figure 6a
p6a<-DimPlot(neuron, 
             group.by = "cell_type_age",
             reduction = "umap",
             # group.by = "ident",
             label =T, label.size = 2.5, pt.size = 0.001,raster=FALSE,
             cols =  c(brewer.pal(11,"Paired")))+
  theme(axis.title = element_blank(), 
        plot.title = element_text(size = 12, face = "bold",hjust = 0), axis.line = element_blank(), 
        axis.text = element_blank(),  axis.ticks = element_blank()) + 
  theme(legend.position = "none") +
  ggtitle(paste0("neuron subtype : ", comma(length (colnames(neuron)),format="d")," cells"))
  
# Figure 6b
p6b<-DotPlot(neuron, 
             features = c("Thy1",neuron.markers.top20$gene),
             assay = "SCT", 
             group.by = "cell_type_age",
             cols = "RdBu", 
             dot.scale	= 10)+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.9, hjust=1))+ RotatedAxis() + coord_flip()+
  theme(legend.title = element_text(size = 12, hjust = 0),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10))

# Figure 6c
require("scProportionTest")

hipp_05m@meta.data$cell_type_age[which(WhichCells(hipp_05m) %in% WhichCells(neuron, idents = '0'))] <- "Neuron:Undef."
hipp_05m@meta.data$cell_type_age[which(WhichCells(hipp_05m) %in% WhichCells(neuron, idents = '1'))] <- "Neuron:Gad1"
hipp_05m@meta.data$cell_type_age[which(WhichCells(hipp_05m) %in% WhichCells(neuron, idents = '2'))] <- "Neuron:Olfm"
hipp_05m@meta.data$cell_type_age[which(WhichCells(hipp_05m) %in% WhichCells(neuron, idents = '3'))] <- "Neuron:Hpca"
hipp_05m@meta.data$cell_type_age[which(WhichCells(hipp_05m) %in% WhichCells(neuron, idents = '4'))] <- "Neuron:Chgb"
hipp_05m@meta.data$cell_type_age[which(WhichCells(hipp_05m) %in% WhichCells(neuron, idents = '5'))] <- "Neuron:Gad2"
hipp_05m@meta.data$cell_type_age[which(WhichCells(hipp_05m) %in% WhichCells(neuron, idents = '6'))] <- "Neuron:Calb2"
hipp_05m@meta.data$cell_type_age[which(WhichCells(hipp_05m) %in% WhichCells(neuron, idents = '7'))] <- "Neuron:Sst"
hipp_05m@meta.data$cell_type_age[which(WhichCells(hipp_05m) %in% WhichCells(neuron, idents = '8'))] <- "Neuron:Slc32a1"

  prop_test <- sc_utils(hipp_05m)
  prop_test <- permutation_test(
    hipp_05m, cluster_identity = "cell_type_age",
    sample_1 = names(table(hipp_05m$genetype))[2], sample_2 = names(table(hipp_05m$genetype))[1],
    sample_identity = "genetype"
  )

prop_test@results$permutation <- prop_test@results$permutation[-18,]
p6c <- permutation_plot(prop_test,log2FD_threshold = 0.3,FDR_threshold = 0.1)

ggsave(paste0(path_base,"/results/Figure6a_DimPlot.pdf"), p6a, units = "in", width =5, height =5, device = "pdf")
ggsave(paste0(path_base,"/results/Figure6b_DotPlot.pdf"), p6b, units = "in", width =5, height =10, device = "pdf")
ggsave(paste0(path_base,"/results/Figure6c_permutation.pdf"), p6c, units = "in", width =10, height =5, device = "pdf")

neuron.DE <- lapply(levels(Idents(neuron)), function(ID){#lapply(names(which(table(Idents(hipp_05m))>50)), function(ID){
  obj<-subset(neuron, idents = ID)
  Idents(obj) <- "genetype"
  if (length(names(table(obj$genetype))) < 2){
    marker <- c()
  }
  else if(table(obj$genetype)[1] < 3 | table(obj$genetype)[2] < 3) {
    
    marker <- c()
    
  } else {
    
    marker <- FindMarkers(obj, ident.1 = names(table(neuron$genetype))[1], ident.2 = names(table(neuron$genetype))[2],assay = "RNA" ,test.use = "MAST", logfc.threshold = 0.0001)
    
  }
  return(marker)
})

require(DOSE)
require(AnnotationDbi)
require(org.Mm.eg.db)
require(clusterProfiler)

Mm <- org.Mm.eg.db

neuron.GO.mf <- list()
genenlist <- list()

for(id in c("4","7")) {

    de<-c()
    selectedGenes <- c()
    gene_list <-c()
    gene_lists <-c()
    
    selectedGene <- neuron.DE[[id]] %>% dplyr::filter(p_val < 0.01  & abs(avg_log2FC) > 0.25)
    selectedGenes <- selectedGene[!grepl("^mt\\.", rownames(selectedGene)),]    
    gene_lists<- selectedGenes %>% dplyr::select(avg_log2FC)

    de<-AnnotationDbi::select(Mm,
                              keys = rownames(gene_lists),
                              columns = c("ENTREZID", "SYMBOL"),
                              keytype = "SYMBOL")
    
    geneLST<-data.frame("SYMBOL"=rownames(gene_lists), "avg_log2FC" = gene_lists$avg_log2FC)
    gene_df <- de %>% mutate(rank = rank(de$ENTREZID,  ties.method = "random")) # %>% 
    gene_dfs <-merge(gene_df, geneLST, by = "SYMBOL", all=TRUE) %>%   arrange(desc(rank))
    gene_dfs<-gene_dfs[!is.na(gene_dfs$ENTREZID),]
    
    if(length(gene_dfs$ENTREZID)>0) {      
      gene.list<-gene_dfs$'avg_log2FC'
      names(gene.list)<-gene_dfs$ENTREZID
      
      
      gene.list = sort(gene.list, decreasing = TRUE)
      gene.list <- na.omit(gene.list)
      genenlist[[id]] <- gene.list

       gse.bp <- gseGO(geneList=gene.list, 
                      ont ="BP",
                      keyType = "ENTREZID",
                      #             nPerm = 10000,
                      minGSSize = 3,
                      maxGSSize = 1000,
                      pvalueCutoff = 0.1,
                      pAdjustMethod = "none",
                      verbose = TRUE,
                      OrgDb = Mm)
       neuron.GO.bp[[id]] <- setReadable(gse.bp, 'org.Mm.eg.db', 'ENTREZID')
    } 
    
}

neuron.DE.A <- merge(neuron.DE[["4"]] %>% dplyr::filter(p_val < 0.01),neuron.DE[["7"]] %>% dplyr::filter(p_val < 0.01), by=0,all = TRUE) 

  common<-intersect(neuron.GO.bp[["4"]]@result$Description,neuron.GO.bp[["7"]]@result$Description)
  commons<-setdiff(grep("synaptic|vesicl|phago|secret|exocy|cili|endo", common, value=TRUE),
                  grep("insulin|hormone|glial|peptide|endothe|tubule", common, value=TRUE))
                  
  id = "7" # or "4"
  pathGenes <- unique(unlist(str_split(neuron.GO.bp[[id]]@result$core_enrichment[neuron.GO.bp[[id]]@result$Description %in% commons],"/")))

    selAccession.n <- neuron.DE.A %>% 
            dplyr::filter(avg_log2FC.x < -0.25 & avg_log2FC.y < -0.25) %>% dplyr::pull(Row.names)
    selAccession.p <- neuron.DE.A %>% 
            dplyr::filter(avg_log2FC.x > 0.25 & avg_log2FC.y > 0.25) %>% dplyr::pull(Row.names)

    # options(repr.plot.width = 16, repr.plot.height = 16)
    keyvals.colour <- ifelse(neuron.DE.A$Row.names %in% selAccession.n, 'royalblue', 'lightgrey')
    keyvals.colour[neuron.DE.A$Row.names %in% selAccession.p] <- "darkred"
    keyvals.colour[is.na(keyvals.colour)] <- 'lightgrey'

    names(keyvals.colour)[keyvals.colour == 'lightgrey'] <- 'ns'
    names(keyvals.colour)[keyvals.colour == 'royalblue'] <- 'down-regulated'
    names(keyvals.colour)[keyvals.colour == 'darkred'] <- 'up-regulated'

    neuron.DE.A$group <- names(keyvals.colour)
    

# Figure 6d
require(ggrepel)
p6d<-ggplot(neuron.DE.A, aes(x=avg_log2FC.x, y=avg_log2FC.y, color =group)) +
      geom_point()+
      geom_label_repel(aes(label = ifelse(Row.names %in% pathGenes, Row.names, ""),
                          label.size = 20), #ifelse(abs(avg_log2FC.x) > 0.25  &  abs(avg_log2FC.y) > 0.25, Row.names, "")),
                        arrow = arrow(length = unit(0.00001, "npc")),
                        colour = "black",max.overlaps = 100, min.segment.length=0.1)+
                      #  max.overlaps = 100, min.segment.length = 0.05) +
      xlab ("5xFAD vs C57BL/6 in Sst-expressing neurons (Log2 Fold Change, p<0.01)") +
      ylab ("5xFAD vs C57BL/6 in Thy1-expressing neurons (Log2 Fold Change, p<0.01)") +
      #ggtitle (paste0("CRPC patients vs healthy control in ",names(table(subset(tumor,idents=id)[['Population']]))," cluster (adjust P-value < 1)")) +
      scale_colour_manual(values = c("royalblue","lightgrey","darkred"))+
      theme_bw()+
        theme(#legend.position=c(0.9, 0.2),
                axis.text.x = element_text(size = 12, hjust = 0, colour = "black"),
                axis.text.y = element_text(size = 12),
                axis.title = element_text(size = 16),
                legend.title = element_blank(),
                legend.background = element_blank()) + 
      theme(legend.position = "none")
      
dtf1<- neuron.GO.bp[["4"]]@result[neuron.GO.bp[["4"]]@result$Description %in% commons, ]
dtf3<- neuron.GO.bp[["7"]]@result[neuron.GO.bp[["7"]]@result$Description %in% commons, ]
dtf<- rbind(cbind(dtf1, "group"="Thy1-expressing neurons"),
            cbind(dtf3, "group"="Sst-expressing neurons"))

# Figure 6e
p6e<-ggplot(dtf, aes(x = Description, y = NES, fill = group)) +
        geom_col(position = "dodge", colour = "black") +
        scale_fill_brewer(palette="Pastel1") +
        xlab("") + 
        ylab("Normalized Enrichment Score (NES)") + 
        theme_bw() + 
        theme(legend.position=c(0.8, 0.9),
                axis.text.x = element_text(size = 16, hjust = 0, colour = "black"),
                axis.text.y = element_text(size = 16),
                axis.title = element_text(size = 16),
                legend.title = element_blank(),
                legend.background = element_blank()) + 
        coord_flip()
        
ggsave(paste0(path_base,"/results/Figure1g3_4ways.pdf"), p6d, units = "in", width =7.5, height =7.5, device = "pdf")
ggsave(paste0(path_base,"/results/Figure1g3_pathway.pdf"), p6e, units = "in", width =10, height =10, device = "pdf")
