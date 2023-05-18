DelayedArray:::set_verbose_block_processing(TRUE)

# Passing a higher value will make some computations faster but use more memory. Adjust with caution!
options(DelayedArray.block.size=1000e6)

base::gc()
base::rm(list = ls()) # Clear the environment
options(warn=-1)
options(stringsAsFactors = FALSE)

load.lib <- c("Seurat","ggplot2","ggrepel","ggpubr","scales","dplyr","tidyverse","devtools","data.table","SingleCellExperiment","RColorBrewer")

if (!requireNamespace(load.lib, quietly = TRUE))
    install.packages(load.lib)

sapply(load.lib,require,character.only = TRUE)
