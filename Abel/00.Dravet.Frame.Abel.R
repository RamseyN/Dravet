######################################################################
# 00.Dravet.Frame.Abel.R
######################################################################
# source('~/GitHub/Projects/Dravet/Abel/00.Dravet.Frame.Abel.R')
stop()
rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)


# Functions ------------------------
source('~/GitHub/Packages/Seurat.pipeline/Load.packages.local.R')

# Parameters ------------------------
source('~/GitHub/Projects/Dravet/Abel/Parameters.Dravet.R')
source('~/GitHub/Projects/Dravet/Abel/GeneLists.Dravet.R')

# Setup ------------------------
OutDir <- OutDirOrig <- PasteOutdirFromFlags(  "~/Dropbox (VBC)/Abel.IMBA/AnalysisD/Dravet"
                                               # , flag.nameiftrue(p$'premRNA')
                                               # , flag.nameiftrue(p$"dSample.Organoids")
                                               , flag.names_list(p$'variables.2.regress')
)

setup_MarkdownReports(OutDir = OutDir, scriptname = '00.Dravet.Frame.Abel.R')
oo()
md.LogSettingsFromList(p)



# Convert to RDS or Read in ------------------------
Continue = F
if (Continue) {
  # combined.obj <- read_rds('~/Dropbox (VBC)/Abel.IMBA/AnalysisD/CON/......Rds.gz'); say()
  combined.obj <- read_rds('~/Dropbox (VBC)/Group Folder Knoblich/Users/Abel.ds/CON.ds/sc4.Crimson.Dlx.2021.10.30/00.CON.Frame.sc4.BC.Dlx.Crimson.R_2021_11_24-07h.only.AAV/combined.obj__2021.11.24_08.17.Rds.gz'); say()
  {
    recall.all.genes()
    all.genes <- list.fromNames(rownames(combined.obj))
    recall.parameters(overwrite = T)
    recall.meta.tags.n.datasets()
    set.mm()
    # qMarkerCheck.BrainOrg()
    GetClusteringRuns()
    clUMAP( "RNA_snn_res.0.`3")
    getClusterNames()
  }
  
} else {
  "Connect to: smb://storage.imp.ac.at/groups/knoblich"
  # InputDir = "/Volumes/demultiplexed/173770/173770_Crimson_Dlx_Fus/"
  InputDir = "/Volumes/knoblich/Organoid_Research/SCN1A/Experiments/PRT_43/scSeq/R13889/cellranger/199216_cr700_20220725/outs/per_sample_outs"
  
  MSeq.obj.fname = list.files(InputDir, pattern = "*.Rds")
  MseqTabelExists = length(MSeq.obj.fname)
  
  if (!MseqTabelExists) Convert10Xfolders(InputDir, min.cells = 3, depth = 3)
  ls.Seurat <- LoadAllSeurats(InputDir, file.pattern =  "*.Rds", string.remove1 = F)
  (samples.short <- samples <- names(ls.Seurat))
  ls.input <- ls.VarGenes <- list.fromNames(samples.short)
  (n.datasets = length(ls.Seurat))
}

# create seurat object:

{
  InputDir = "/Volumes/knoblich/Organoid_Research/SCN1A/Experiments/PRT_43/scSeq/R13889/cellranger/199216_cr700_20220725/"
  
  (samples <- list.files(glue::glue("{InputDir}/outs/per_sample_outs/")))
  ls.Seurat <- lapply(samples, function(s){
    
    cur_data <- Read10X(glue::glue("{InputDir}/outs/per_sample_outs/{s}/count/sample_filtered_feature_bc_matrix"))
    
    cur_seurat <- CreateSeuratObject(
      counts = cur_data[["Gene Expression"]],
      project = s
    )
    
    cur_seurat <- RenameCells(cur_seurat, add.cell.id = s)
    return(cur_seurat)
  })
  
}
ls.Seurat.bac <- ls.Seurat

## meta.tags ------------------------
# samples.short = "d80.AUT"
meta.tags <- list(
  'library' = samples,
  'genotype' = substr(samples,1,2),
  'batch' = substr(samples,nchar(samples),nchar(samples))
  
  
)
  n.datasets=l(samples)
##  ------------------------

# Metadata ---------------------------
i = 1
for(i in 1:n.datasets ) {
  ls.Seurat[[i]] <- add.meta.fraction(obj = ls.Seurat[[i]], col.name = "percent.mito", gene.symbol.pattern = "^MT\\.|^MT-")
  ls.Seurat[[i]] <- add.meta.fraction(obj = ls.Seurat[[i]], col.name = "percent.ribo", gene.symbol.pattern = "^RPL|^RPS")
  ls.Seurat[[i]] <- add.meta.fraction(obj = ls.Seurat[[i]], col.name = "percent.AC.GenBank", gene.symbol.pattern = "^AC[0-9]{6}\\.")
  ls.Seurat[[i]] <- add.meta.fraction(obj = ls.Seurat[[i]], col.name = "percent.AL.EMBL", gene.symbol.pattern = "^AL[0-9]{6}\\.")
  ls.Seurat[[i]] <- add.meta.fraction(obj = ls.Seurat[[i]], col.name = "percent.LINC", gene.symbol.pattern = "^LINC0")
  ls.Seurat[[i]] <- add.meta.fraction(obj = ls.Seurat[[i]], col.name = "percent.MALAT1", gene.symbol.pattern = "^MALAT1")
  ls.Seurat[[i]] <- add.meta.fraction(obj = ls.Seurat[[i]], col.name = "percent.antisense", gene.symbol.pattern = "-AS[0-9]")
  ls.Seurat[[i]] <- add.meta.fraction(obj = ls.Seurat[[i]], col.name = "percent.RabV", gene.symbol.pattern = "RabV.N2c$")
  ls.Seurat[[i]] <- add.meta.fraction(obj = ls.Seurat[[i]], col.name = "percent.AAV", gene.symbol.pattern = "^AAV-")
  # More in source('~/GitHub/Projects/SEO/elements/Add.Metadata.R')
}





for(i in 1:n.datasets ) {
  total_expr = Matrix::colSums(GetAssayData(object = ls.Seurat[[i]]))
  META = ls.Seurat[[i]]@meta.data
  
  HGA_Markers <- intersect(genes.ls$'HGA_MarkerGenes', rownames(ls.Seurat[[i]]))
  HGA_Markers = GetAssayData(object = ls.Seurat[[i]])[HGA_Markers, ]
  META$'log10.HGA_Markers'  <- log10(Matrix::colSums(HGA_Markers) + 1) / (total_expr * (1 - META$'percent.ribo') * (1 - META$'percent.mito'))
  
  RabV_Markers <- intersect(genes.ls$'RabV', rownames(ls.Seurat[[i]]))
  RabV_Markers = GetAssayData(object = ls.Seurat[[i]])[RabV_Markers, ]
  META$'log10.RabV_Markers'  <- log10(Matrix::colSums(RabV_Markers) + 1) / (total_expr * (1 - META$'percent.ribo') * (1 - META$'percent.mito'))
  
  META$'log10.nCount_RNA' <- log10(as.named.vector(META[, 'nCount_RNA', drop = F]))
  META$'log10.nFeature_RNAs' <- log10(as.named.vector(META[, 'nFeature_RNA', drop = F]))
  nCells = nrow(META)
  
  for (variable in names(meta.tags)) {
    META[[variable]] <- rep(meta.tags[[variable]][i], nCells)
  }
  ls.Seurat[[i]]@meta.data <- META
}

  


# Filtering / Prepare for CCA ------------------------------------------------------------------------------------------
names(ls.Seurat) <- samples
if (TRUE) PlotFilters(ls.obj = ls.Seurat)
if (TRUE) source("~/GitHub/Packages/Seurat.pipeline/elements/Filtering.plots.3D.multiplex.R"); create_set_Original_OutDir()
Nr.Cells.Before.Filtering = unlapply(ls.Seurat, ncol); # names(Nr.Cells.Before.Filtering) = suffices

# isave.RDS(ls.Seurat) # if (!exists("ls.Seurat.bac")) ls.Seurat.bac <- ls.Seurat # ls.Seurat <- ls.Seurat.bac
if (TRUE) source("~/GitHub/Packages/Seurat.pipeline/elements/Filter.Dataset.R"); create_set_Original_OutDir()



# Downsample per dataset  ------------------------------------------------------------------------------------------
# if (p$"dSample.Organoids") source("~/GitHub/Projects/TSC2/s3/s3.elements/dSample.Organoids.R")
unlapply(ls.Seurat, ncol)

# ls.Seurat <- ls.Seurat.bac


# parallelize.new = T
# if (parallelize.new) parallel.computing.by.future(maxMemSize = 12000 * 1024^2, workers_ = 1)


tic(); ls.Seurat <- future.apply::future_lapply(X = ls.Seurat, FUN = NormalizeData, normalization.method = "LogNormalize", scale.factor = 10000)
# foreach(i=1:n.datasets) %dopar% { NormalizeData(object = ls.Seurat[[i]], normalization.method = "LogNormalize", scale.factor = 10000) };
toc();


tic(); ls.Seurat <- future.apply::future_lapply(ls.Seurat, FUN = FindVariableFeatures, mean.function = 'FastExpMean', dispersion.function = 'FastLogVMR', nfeatures =10000)
# foreach(i=1:n.datasets) %dopar% { FindVariableFeatures(object = ls.Seurat[[i]], mean.function = 'FastExpMean', dispersion.function = 'FastLogVMR', nfeatures =10000) };
toc();




# VariableFeaturePlot ---------------
PlotVariableFeatures = T
if (PlotVariableFeatures) source('~/GitHub/Packages/Seurat.pipeline/elements/Plot.variable.Genes.R'); create_set_Original_OutDir()




# merge seurat object
combined.obj <- merge(x=ls.Seurat[[1]], y=ls.Seurat[2:length(ls.Seurat)])

# rm(list = c( 'obj', 'ls.obj',  'ls.Seurat'))


{
  TotalReadCounts <- sparseMatrixStats::rowSums2(combined.obj@assays$RNA@counts)
  names(TotalReadCounts) <- rownames(combined.obj@assays$RNA@counts)
  combined.obj@misc$'TotalReadCounts' <- TotalReadCounts
  
  TotalExpressingCells <- sparseMatrixStats::rowSums2(combined.obj@assays$RNA@counts > 0)
  names(TotalExpressingCells) <- rownames(combined.obj@assays$RNA@counts)
  combined.obj@misc$'TotalExpressingCells' <- TotalExpressingCells
  
}

# Some Metadata calculations ---------------------------------------------------------------------------
# if (F) write_clip (rownames(combined.obj))

combined.obj <- calc.q99.Expression.and.set.all.genes(obj = combined.obj, quantileX = .99)

# all.genes <- combined.obj@misc$all.genes
combined.obj@misc$'meta.tags' <- meta.tags

## Perform integrated analysis --------------
# if (parallelize.new) parallel.computing.by.future(maxMemSize = 12000 * 1024^2, workers_ = 1)

for (assayX in c("RNA")) {
  tic(); combined.obj <- ScaleData(combined.obj, assay = assayX, verbose = T, vars.to.regress = p$'variables.2.regress'); toc()
}; say()

if (parallelize.new) parallel.computing.by.future(maxMemSize = 12000 * 1024^2, workers_ = 4)

# GetAssayData(object = combined.obj, assay = assay.def, slot = "data")
# t-SNE and Clustering ---------------------------------------------------------------------------
combined.obj <- FindVariableFeatures(object = combined.obj, mean.function = 'FastExpMean', dispersion.function = 'FastLogVMR', nfeatures =10000)

tic(); combined.obj <- RunPCA(combined.obj, npcs = p$'n.PC', verbose = T); toc()
# tic(); combined.obj <- RunUMAP(combined.obj, reduction = "pca", dims = 1:p$'n.PC'); toc()
tic(); combined.obj <- SetupReductionsNtoKdimensions(obj = combined.obj, nPCs = p$'n.PC', dimensions = 3:2, reduction = "umap"); toc()

tic(); combined.obj <- schex::make_hexbin(combined.obj, nbins = 15, dimension_reduction = "UMAP"); toc()

tic(); combined.obj <- FindNeighbors(combined.obj, reduction = "pca", dims = 1:p$'n.PC'); toc()
tic(); combined.obj <- FindClusters(combined.obj, resolution = p$'snn_res'); toc()
# tic(); combined.obj <- RunTSNE(combined.obj, reduction = "pca", dims = 1:p$'n.PC'); toc()
isave.RDS(combined.obj, inOutDir = T)

say()
# Visualization ---------------------------------------------------------------------------

# mm = combined.obj@meta.data
# combined.obj@meta.data$'v.project' = ppp(mm$experiment,mm$RNA.model)
# combined.obj@meta.data$'RNA.model'   = translate(mm$RNA.model, oldvalues = "WT_mRNA", newvalues = "mRNA")
# if (F) {source("~/GitHub/Projects/TSC2/s3/s3.elements/zz.plots.4.ORC.R")}

# Basic Stats ------------------------------------------------------------------------
DefaultAssay(combined.obj) <- "RNA"
if (TRUE) source('~/GitHub/Packages/Seurat.pipeline/elements/Plots.stats.R'); create_set_Original_OutDir()

p$'also.tSNE' = F
if (TRUE) source('~/GitHub/Packages/Seurat.pipeline/elements/Gene.expression.gene.lists.R'); create_set_Original_OutDir()

# Second part of heavy analysis ------------------------------------------------------------------------

"Do THIS"
"Do THIS"
"Do THIS"
"Do THIS"
"Do THIS"
"Do THIS"
if (FALSE) source('~/GitHub/Projects/CON/sc4.Crimson/sc4.GEX/Remove.NonNeural.Clusters.R')
"Do THIS"
"Do THIS"
"Do THIS"





Part.2 = F
if (Part.2) {
  p$"Reorder.Dim" = -2
  n.datasets=1
  if (TRUE) source("~/GitHub/Packages/Seurat.pipeline/elements/Renumber.Clusters.R"); create_set_Original_OutDir()
  
  
  if (TRUE) source('~/GitHub/Packages/Seurat.pipeline/elements/Differential.gene.expression.Loop.R'); create_set_Original_OutDir()
  # if (F) source("~/GitHub/Projects/TSC2/s3/s3.elements.mseq/Differential.gene.expression.Batch.comparison.R"); create_set_Original_OutDir()
  
  
  
  "weird: batch was not set"
  combined.obj$batch <- substr(combined.obj$orig.ident,nchar(samples),nchar(samples))
  
  "weird: batch was not set"
  
  plotCellFractionsBarplots <- T
  if (TRUE) source('~/GitHub/Packages/Seurat.pipeline/elements/Cell.cycle.scoring.R'); create_set_Original_OutDir()
  
  if (TRUE) source('~/GitHub/Packages/Seurat.pipeline/elements/PCA.R'); create_set_Original_OutDir()
  
  if (TRUE) source("~/GitHub/Packages/Seurat.pipeline/elements/Plot.3D.umaps.R"); create_set_Original_OutDir()
  combined.obj <- RecallReduction(obj = combined.obj, dim = 2, reduction = "umap")
  
}



if (TRUE) { combined.obj@misc$p <- p; isave.RDS(obj = combined.obj, inOutDir = T) } # Add a parameter list to a Seurat object's misc slot # seuSaveRds(object = combined.obj, use_Original_OutDir = T)


memory.biggest.objects()
say()


