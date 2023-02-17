# Author : Jaemyung
# Contact : piloter2@kbri.re.kr
# Date : Feb. 17, 2023

# Pass TRUE if you want to see progress output on some of Monocle 3's operations
DelayedArray:::set_verbose_block_processing(TRUE)

# Passing a higher value will make some computations faster but use more memory. Adjust with caution!
options(DelayedArray.block.size=1000e6)

base::gc()
base::rm(list = ls()) # Clear the environment
options(warn=-1)

if(!"RCy3" %in% installed.packages()){
    install.packages("BiocManager")
    BiocManager::install("RCy3")
}

load.lib <- c("Seurat","dplyr","tidyverse","data.table","scales","RColorBrewer","ggplot2")
load.libs <- c("DOSE","GO.db","GSEABase","org.Mm.eg.db","clusterProfiler","AnnotationDbi","SingleCellExperiment","RCy3")
#,"dplyr","tidyr",  "ggplot2","stringr",  "RColorBrewer","rWikiPathways","RCy3"
# install.packages(load.lib)
# BiocManager::install(load.libs)

#p_load(load.libs, update = TRUE, character.only = TRUE)
sapply(load.lib,require,character.only = TRUE)
sapply(load.libs,require,character.only = TRUE)


category.bp<-list()
for(id in c("4","6","7")) {
  print(id)
  category.bp[[id]] = grep("cili|exocy|secretion",neuron.GO.bp[[id]]@result$Description, value = TRUE) #unique(org.GO.bp@result$Description[grep(paste0(pathGenes,collapse = "$|"),org.GO.bp@result$core_enrichment)])#unique(hipp.GO@result$Description[grep(paste0(pathGenes[[id]],collapse = "|"),hipp.GO@result$core_enrichment)])
  print(category.bp[[id]]) 
}

id <- "7"
egobp <- neuron.GO.bp[[id]]
## extract a dataframe with results from object of type enrichResult
egobp.results.df <- egobp@result[which(egobp@result$Description %in% category.bp[[id]]),]
egobp.results.dfs <- egobp.results.df[-grep("hormone|insulin",egobp.results.df$Description),]
## create a new column for term size from BgRatio
#egobp.results.df$term.size <- gsub("/(\\d+)", "", egobp.results.df$BgRatio)
#egobp.results.dfs$geneID <- gsub("/", ",", egobp.results.dfs$core_enrichment)
# head(neuron.DE[[id]])
de.genes <- neuron.DE[[id]][,c(1,2,5)]
colnames(de.genes) <- c("P.Value","log2FC","adj.P.Value")
de.genes$Gene <- rownames(de.genes)
#loadTableData(de.genes)
#loadTableData(egobp.results.dfs)

egobp.results.filtered <- egobp.results.dfs[,c(2:3,5)]
egobp.results.filtered$geneID <- str_split(egobp.results.dfs$core_enrichment,"/")

view(egobp.results.filtered)

cytoscapePing ()
cytoscapeVersionInfo ()

#gD <- igraph::simplify(igraph::graph.data.frame(egobp.results.filtered, directed=FALSE))
relations <- data.frame(source = egobp.results.filtered$geneID[[1]],
                        target=rep(egobp.results.filtered$Description[1],length(egobp.results.filtered$geneID[[1]])))

for(i in 2:length(egobp.results.filtered$Description)) {
    relations <- rbind(relations,
                        data.frame(source = egobp.results.filtered$geneID[[i]],
                                   target = rep(egobp.results.filtered$Description[i],length(egobp.results.filtered$geneID[[i]]))))
}

colnames(relations)

actors <- egobp.results.filtered[,c(1:3)]

de.genes.filtered <- de.genes[which(de.genes$Gene %in% relations$source),]
colnames(de.genes.filtered) <- c("P.Value","log2FC","adj.P.Value","source")

head(relations)
head(de.genes.filtered)

relation <- merge(relations, de.genes.filtered, by = "source", all = TRUE )

# colnames(de.genes.filtered) <- c("P.Value","log2FC","adj.P.Value","name")
ig <- igraph::graph_from_data_frame(relation, directed=FALSE)
# print(ig, e = TRUE, v = TRUE)

createNetworkFromIgraph(ig,"Figure5_GRN")

style.name = "Figure5"
createVisualStyle(style.name)
setVisualStyle(style.name)

setNodeShapeDefault("ellipse", style.name) #remember to specify your style.name!
setNodeSizeDefault(60, style.name)
setNodeColorDefault("#AAAAAA", style.name)
setEdgeLineWidthDefault(2, style.name)
setNodeLabelMapping('name', style.name)

expr.network <- getTableColumns(table = "edge", columns = 'log2FC') 

loadTableData(
  de.genes.filtered[,c(1:3)],
  data.key.column = "row.names",
  table = "node",
  table.key.column = "name")

loadTableData(
  egobp.results.dfs[,c(2,5)],
  data.key.column = "Description",
  table = "node",
  table.key.column = "name")

setNodeShapeBypass(node.names =unique(egobp.results.dfs$Description) , "ROUND_RECTANGLE")
setNodeColorBypass(node.names =unique(egobp.results.dfs$Description), '#FF55AA')

min.logFC = min(expr.network[,1],na.rm=TRUE)
max.logFC = max(expr.network[,1],na.rm=TRUE)
data.values = c(-abs(min.logFC),0,abs(min.logFC))
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
# display.brewer.all(length(data.values), colorblindFriendly=TRUE, type="div") 
selectNodes(nodes=rownames(de.genes.filtered), by.col = "name")
setNodeColorMapping(table.column= 'log2FC', data.values, node.colors, style.name=style.name)

selectNodes(nodes=unique(egobp.results.dfs$Description), by.col = "name")

data.values2 = c(-abs(max(egobp.results.dfs$NES)),0,abs(max(egobp.results.dfs$NES)))
node.colors2 <- c(rev(brewer.pal(length(data.values2), "RdBu")))
setNodeColorMapping(table.column= 'NES', data.values2, node.colors2, style.name=style.name)
setNodeSizeMapping(table.column= 'NES', style.name=style.name)
