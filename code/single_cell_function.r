DelayedArray:::set_verbose_block_processing(TRUE)

# Passing a higher value will make some computations faster but use more memory. Adjust with caution!
options(DelayedArray.block.size=1000e6)

base::gc()
base::rm(list = ls()) # Clear the environment
options(warn=-1)
options(stringsAsFactors = FALSE)

load.lib <- c("Seurat","ggplot2","ggrepel","ggpubr","scales","dplyr","tidyverse","devtools","data.table","SingleCellExperiment")

if (!requireNamespace(load.lib, quietly = TRUE))
    install.packages(load.lib)

sapply(load.lib,require,character.only = TRUE)

require.lib <- c("cowplot","BiocManager","RColorBrewer","ggpubr","enrichR","ggrepel")

if (!requireNamespace(require.lib, quietly = TRUE))
    install.packages(require.lib)

bioconductor.lib <- c("EnhancedVolcano","Nebulosa","DOSE","AnnotationDbi","org.Mm.eg.db","clusterProfiler")

if (!requireNamespace(bioconductor.lib, quietly = TRUE))
    BiocManager::install(bioconductor.lib)

if (!requireNamespace("nichenetr", quietly = TRUE))
    devtools::install_github("saeyslab/nichenetr")

if (!requireNamespace("scProportionTest", quietly = TRUE))
    devtools::install_github("rpolicastro/scProportionTest")
