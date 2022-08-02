library(Seurat)
library(cowplot)
library(harmony)
library(tidyverse)
library(glue)
library(Seurat.utils)
require(MarkdownHelpers)
require(MarkdownReports)
require(Stringendo)
require(ggExpress)

in_data_dir <- '/Volumes/groups/knoblich/Organoid_Research/SCN1A/Experiments/PRT_43/scSeq/R13889/cellranger/199216_cr700_20220725'

OutDir = '/Users/ramsey.najm/GitHub/Projects/RN_SCN1A_PRT43/'

wA4=12
hA4=21

samples <- list.files(glue("{in_data_dir}/outs/per_sample_outs/"))


# create seurat object:
seurat_list <- lapply(samples, function(s){
  
  cur_data <- Read10X(glue("{in_data_dir}/outs/per_sample_outs/{s}/count/sample_filtered_feature_bc_matrix"))
  
  cur_seurat <- CreateSeuratObject(
    counts = cur_data[["Gene Expression"]],
    project = s
  )
  
  cur_seurat <- RenameCells(cur_seurat, add.cell.id = s)
  return(cur_seurat)
})

# merge seurat object
seurat_obj <- merge(x=seurat_list[[1]], y=seurat_list[2:length(seurat_list)])

head(seurat_obj@meta.data)
table(seurat_obj@meta.data$orig.ident)