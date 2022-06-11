### Genome Research
### Michael Vermeulen
### Mosaic loss of chromosome Y in aged human microglia
### sn/scRNA-seq processing 


### This is the template used to process all single-nuclei and single-cell data. 
### Here, preprocessing of GSE174332 primary motor cortex data is provided. 
### All other preprocessing scripts are available upon request. 

### once each cell-type is processed individually and inspected, the cell-type files are merged. 


## load packages 

library(Biobase)
library(GEOquery)
library(Seurat)
library(ggpubr)
library(tidyverse)
library(scibetR)
library(harmony)
library(viridis)
library(tibble)
library(reshape2)
library(scDblFinder)
library(SingleR)
library(celldex)
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
edb <- EnsDb.Hsapiens.v86

setwd("H:/LOY/scLOY/scLOY/")
source("E:/LOY/Example/PACKAGE/call_LOY.R")
source("E:/LOY/scLOY/scripts/loy_functions.R")
source('H:/LOY/scLOY/scLOY/R/load.R')
source('H:/LOY/scLOY/scLOY/data_raw/build_data_files.R')
source("E:/LOY/scLOY/scripts/annotation_functions.R")


## Download GEO metadata if available.
gsm <-getGEO("GSE174332", destdir = "H:/LOY/scLOY/data/GEO")
gsm[[1]] -> gsm
pheno<-pData(phenoData(gsm))


# downloaded snRNAseq matricies from: 
# https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE174332&format=file


# load raw data from dir
files <- list.dirs(path = "E:/LOY/scLOY/data/BRAIN_GSE174332", full.names = T); files <- files[-1]
lapply(X = files, FUN = function(x){
  basename(x) %>% gsub(pattern = "_filtered_feature_bc_matrix", replacement = "") -> lib
  message(lib)
  Matrix::readMM(file = paste0(x,"/counts_fil.mtx")) -> mat
  data.table::fread(file = paste0(x,"/col_metadata.tsv"), data.table = F) -> col
  data.table::fread(file = paste0(x,"/row_metadata.tsv"), data.table = F) -> row
  dimnames(mat) <- list(row$Gene,col$Barcode)
  rownames(col) <- col$Barcode
  CreateSeuratObject(counts = mat, min.cells = 0, min.features = 200, meta.data = col) -> o
  subset(o, CellClass %in% c("Glia")) -> o
  o@meta.data$library_id <- lib
  return(o)
}) -> l

l[[1]] -> o
for(i in 2:length(l)){
  message(i)
  merge(o, l[[i]]) -> o
}

# Get expressed PAR, chrY, ribosomal and mito genes
rownames(o)[grep(x = rownames(o), pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")] -> ribo
ribo[grep(ribo,pattern = "K",invert = T)] -> ribo
rownames(o)[grep(x = rownames(o), pattern = "^MT-")] -> mito
rownames(o)[(rownames(o) %in% Y_scLOY_genes$gene_name)] -> Y
Y_scLOY_genes$gene_name[Y_scLOY_genes$gene_name %in% rownames(o)] -> male_genes

# Get % PercentageFeatureSet of mito and ribo genes
o[["percent.mt"]] <- Seurat::PercentageFeatureSet(o, pattern = "^MT-")
o[["percent.rb"]] <- Seurat::PercentageFeatureSet(o, features = ribo)

# print out percentile distribution of QC metrics 
quantile(o@meta.data$nCount_RNA, probs = seq(0,1,0.01))
quantile(o@meta.data$nFeature_RNA, probs = seq(0,1,0.01))
quantile(o@meta.data$percent.mt, probs = seq(0,1,0.01))
quantile(o@meta.data$percent.rb, probs = seq(0,1,0.01))

par(mfrow=c(2,2))
hist(o@meta.data$nCount_RNA,breaks = 50)
hist(o@meta.data$nFeature_RNA,breaks = 50)
hist(o@meta.data$percent.mt,breaks = 50)
hist(o@meta.data$percent.rb,breaks = 50)

# filter seurat object for QC metrics
subset(o, nCount_RNA > 1000  & nCount_RNA < 7500 &
         nFeature_RNA > 800 &  nFeature_RNA < 3500 & 
         percent.mt < 5 & 
         percent.rb < 3 ) -> o

## Run doubletfinder
## find doublet and singlet nuclei/cells
as.SingleCellExperiment(o) -> sce
sce <- scDblFinder(sce, samples = "orig.ident")
sce$scDblFinder.class -> o@meta.data$scDblFinder.class
rm(sce)
table(o@meta.data$scDblFinder.class)

# Seurat normalization pipeline
o <- Seurat::NormalizeData(o, scale.factor = 10000)
o <- Seurat::FindVariableFeatures(o, selection.method = "vst", nfeatures = 3000)

# dissociation gene module
# Genes commonly observed in disturbed microglia 
diss <- c("RGS1","DUSP1","AC022217.3","DUSP1","HSPA1A","HSPA1B","DDIT4","SLC2A3","FOS","JUN",
          "AC131944.1","HIF1A-AS3","HIST1H2AC","HIST1H2BD","HIST1H2BG","HIST2H2BE","RGS2","RHOB",
          "HIST1H2BJ","TEX14","UBC","HSP90AA1","AL118516.1","JUNB","HIST2H2AA4","HIST1H1C","BTG1","BTG2","DUSP5")
list(diss[diss %in% rownames(o)]) -> diss
Seurat::AddModuleScore(object = o, features = diss, name = "diss") -> o

# PAR and Y modules
## add PAR and MSY scores for each nuclei/cell
list(PAR_scLOY_genes$gene_name[PAR_scLOY_genes$gene_name %in% rownames(o)]) -> P
Seurat::AddModuleScore(object = o, features = P, name = "PAR_score") -> o
list(Y_scLOY_genes$gene_name[Y_scLOY_genes$gene_name %in% rownames(o)]) -> Y
Seurat::AddModuleScore(object = o, features = Y, name = "Y_score") -> o

# Scale data and run PCA
o <- Seurat::ScaleData(o, vars.to.regress = c("nCount_RNA","percent.mt","percent.rb"))
v <- Seurat::VariableFeatures(o)
v <- v[!(v %in% c(ribo,mito,male_genes,"XIST","JPX"))]  # remove ribo/mito and sex specific genes from clustering
o <- Seurat::RunPCA(o, npcs = 50, features = v, assay = "RNA")
ElbowPlot(o, ndims = 50)

# Dimension reduction
o <- Seurat::RunUMAP(o,reduction = "pca", dims = 1:15)
o <- Seurat::FindNeighbors(o,reduction = "pca", dims = 1:15)
o <- Seurat::FindClusters(o,resolution = 0.1)
DimPlot(o, group.by = "RNA_snn_res.0.1", label = T)
DimPlot(o, group.by = "scDblFinder.class", label = T)
Idents(o) <- "RNA_snn_res.0.1"

Idents(o) <- "RNA_snn_res.0.1"

### Run harmony batch correction
o <- harmony::RunHarmony(o,
                         group.by.vars = c("Batch"),
                         assay.use = "RNA",
                         plot_convergence = T )
o <- Seurat::RunUMAP(o,reduction = "harmony", dims = 1:20)
o <- Seurat::FindNeighbors(o,reduction = "harmony", dims = 1:20)
o <- Seurat::FindClusters(o,resolution = seq(0.1,2,0.1))

# visualize UMAP
DimPlot(o, group.by = "RNA_snn_res.0.5", label = T, raster = T) -> dim
DimPlot(o, group.by = "Sample_ID", label = T, raster = T) -> dim


## annotate clusters
Idents(o) <- "RNA_snn_res.0.5"

FindMarkers(o, ident.1 = 3, ident.2 = NULL, group.by = "RNA_snn_res.0.1", only.pos = T, logfc.threshold = 0.75)

plot_brain(s_obj = o, seurat_res = 0.5, output = "")

o@meta.data <- dplyr::left_join(x = o@meta.data %>% tibble::rownames_to_column("X"),
                                y = as.data.frame(main) %>% tibble::rownames_to_column("X") %>% dplyr::select(X,labels_main),
                                by = "X") %>% tibble::column_to_rownames("X")

o@meta.data <- dplyr::left_join(x = o@meta.data %>% tibble::rownames_to_column("X"),
                                y = as.data.frame(fine) %>% tibble::rownames_to_column("X") %>% dplyr::select(X,labels),
                                by = "X") %>% tibble::column_to_rownames("X")


FeaturePlot(o, features = "CD74", order = T) + scale_color_viridis_c()
FeaturePlot(o, features = "C3", order = T) + scale_color_viridis_c()
FeaturePlot(o, features = "SLC17A7", order = T) + scale_color_viridis_c()
FeaturePlot(o, features = "GAD1", order = T) + scale_color_viridis_c()
FeaturePlot(o, features = "VWF", order = T) + scale_color_viridis_c()
FeaturePlot(o, features = "AQP4", order = T) + scale_color_viridis_c()
FeaturePlot(o, features = "VCAN", order = T) + scale_color_viridis_c()


o@meta.data$cell_ident_MV <- o@meta.data$RNA_snn_res.0.1
levels(o@meta.data$cell_ident_MV)[levels(o@meta.data$cell_ident_MV)%in%c(4)] <- "Microglia"
levels(o@meta.data$cell_ident_MV)[levels(o@meta.data$cell_ident_MV)%in%c(0,1)] <- "Oligo"
levels(o@meta.data$cell_ident_MV)[levels(o@meta.data$cell_ident_MV)%in%c(2)] <- "Astro"
levels(o@meta.data$cell_ident_MV)[levels(o@meta.data$cell_ident_MV)%in%c(3)] <- "OPC"
DimPlot(o, group.by = "cell_ident_MV", label = T)

# call loss of Y 
o@meta.data$donor_organism.biomaterial_core.biomaterial_id <- o@meta.data$library_id
call_sample_sex(o) -> o
LOY_genes_to_exclude_expr(o) -> remove
call_LOY(o,genes_to_exclude = c(remove) %>% unique()) -> o

# add cell_broad_MV
NA ->  o@meta.data$cell_broad_MV

o@meta.data$cohort <- "Brain ALS/FTLD (GSE174332)"

# sex + cluster specific PAR and Y modules 
module_score_by_cluster(o, ident = "cell_ident_MV") -> o

list(PAR_scLOY_genes$gene_name[PAR_scLOY_genes$gene_name %in% rownames(o)]) -> P
Seurat::AddModuleScore(object = o, features = P, name = "PAR_score") -> o
list(Y_scLOY_genes$gene_name[Y_scLOY_genes$gene_name %in% rownames(o)]) -> Y
Seurat::AddModuleScore(object = o, features = Y, name = "Y_score") -> o

# metadata cleaning 
o@meta.data$Condition -> o@meta.data$diagnosis
ifelse(o@meta.data$Condition %in% c("ALS","FTLD"), F, T) -> o@meta.data$neuro_degen
o@meta.data$Condition -> o@meta.data$neuro_degen_type
o@meta.data[o@meta.data$neuro_degen_type=="PN",]$neuro_degen_type <- "Control"

o@meta.data$smoking <- NA
o@meta.data$donor_organism.age <- NA
o@meta.data$tissue <- "Brain"
o@meta.data$specific_tissue <- "Motor cortex"
o@meta.data$GEO <- "GSE174332"
o@meta.data$sequencing_method2 <- "single-nuclei"
o@meta.data$sequencing_method <- "10X 3' v3"
o@meta.data$introns <- "included"
o@meta.data$cancer_diagnosis <- NA

subset(o, nCount_RNA > 2500) -> o.
LOY_dataset_summary(o.@meta.data) %>% View()

check_object_meta(o@meta.data)

saveRDS(o, file = "E:/LOY/scLOY/processed_seurat/BRAIN/GSE174332_BRAIN_GLIA.RDS")

### CREATE SUBSETS FOR EACH MAJOR CELL TYPE 

############# Neuron ################
DimPlot(o, group.by = "cell_ident_MV", label = T)
subset(o, cell_ident_MV == "Neuron") -> o.

# Scale data and run PCA
o. <- Seurat::NormalizeData(o., scale.factor = 10000)
o. <- Seurat::FindVariableFeatures(o., selection.method = "vst", nfeatures = 2000)
o. <- Seurat::ScaleData(o., vars.to.regress = c("nCount_RNA","percent.mt","percent.rb"))
v <- Seurat::VariableFeatures(o.)
v <- v[!(v %in% c(ribo,mito,male_genes,"XIST"))]
o. <- Seurat::RunPCA(o., npcs = 50, features = v, assay = "RNA")
ElbowPlot(o., ndims = 50)

o. <- Seurat::RunUMAP(o.,reduction = "pca", dims = 1:15)
o. <- Seurat::FindNeighbors(o.,reduction = "pca", dims = 1:15)
o. <- Seurat::FindClusters(o.,resolution = 0.1)
DimPlot(o., group.by = "donor_organism.biomaterial_core.biomaterial_id", label = T)
Idents(o.) <- "RNA_snn_res.0.1"

### 
o. <- harmony::RunHarmony(o. , theta = c(2), lambda  = c(1),
                          group.by.vars = c("donor_organism.biomaterial_core.biomaterial_id"),
                          assay.use = "RNA",
                          plot_convergence = T )
o. <- Seurat::RunUMAP(o.,reduction = "harmony", dims = 1:30)
o. <- Seurat::FindNeighbors(o.,reduction = "harmony", dims = 1:30)
o. <- Seurat::FindClusters(o.,resolution = seq(0.1,2,0.1))
DimPlot(o., group.by = "donor_organism.biomaterial_core.biomaterial_id", label = T)
DimPlot(o., group.by = "RNA_snn_res.0.1", label = T)

FeaturePlot(o., features = "GAD1", order = T)
FeaturePlot(o., features = "SLC17A7", order = T)
FeaturePlot(o., features = "SST", order = T)
FeaturePlot(o., features = "VIP", order = T)
FeaturePlot(o., features = "PVALB", order = T)

FindMarkers(o., ident.1 = 0, ident.2 = NULL, group.by = "RNA_snn_res.0.1", only.pos = T)

o.@meta.data$cell_subtype_MV <- o.@meta.data$RNA_snn_res.0.1
levels(o.@meta.data$cell_subtype_MV)[levels(o.@meta.data$cell_subtype_MV)%in%c(1,11,4,8,10)] <- "Inhib neuron"
levels(o.@meta.data$cell_subtype_MV)[levels(o.@meta.data$cell_subtype_MV)%in%c(13,2,7,0,5,6,3,9,12)] <- "Ex neuron"
DimPlot(o., group.by = "cell_subtype_MV", label = T)

saveRDS(o.,"E:/LOY/scLOY/processed_seurat/CELL_TYPE_SPECIFIC/NEURON/")

############# Astrocyte ################
DimPlot(o, group.by = "cell_ident_MV", label = T)
subset(o, cell_ident_MV == "Astro") -> o.

# Scale data and run PCA
o. <- Seurat::NormalizeData(o., scale.factor = 10000)
o. <- Seurat::FindVariableFeatures(o., selection.method = "vst", nfeatures = 2000)
o. <- Seurat::ScaleData(o., vars.to.regress = c("nCount_RNA","percent.mt","percent.rb"))
v <- Seurat::VariableFeatures(o.)
v <- v[!(v %in% c(ribo,mito,male_genes,"XIST"))]
o. <- Seurat::RunPCA(o., npcs = 50, features = v, assay = "RNA")
ElbowPlot(o., ndims = 50)

o. <- Seurat::RunUMAP(o.,reduction = "pca", dims = 1:15)
o. <- Seurat::FindNeighbors(o.,reduction = "pca", dims = 1:15)
o. <- Seurat::FindClusters(o.,resolution = 0.1)
DimPlot(o., group.by = "donor_organism.biomaterial_core.biomaterial_id", label = T) + NoLegend()
Idents(o.) <- "RNA_snn_res.0.1"

### 
o. <- harmony::RunHarmony(o. , theta = c(2), lambda  = c(1),
                          group.by.vars = c("donor_organism.biomaterial_core.biomaterial_id"),
                          assay.use = "RNA",
                          plot_convergence = T )
o. <- Seurat::RunUMAP(o.,reduction = "harmony", dims = 1:15)
o. <- Seurat::FindNeighbors(o.,reduction = "harmony", dims = 1:15)
o. <- Seurat::FindClusters(o.,resolution = seq(0.1,2,0.1))
DimPlot(o., group.by = "donor_organism.biomaterial_core.biomaterial_id", label = T)
DimPlot(o., group.by = "RNA_snn_res.0.1", label = T)

FeaturePlot(o., features = "GFAP", order = T)
FeaturePlot(o., features = "", order = T)
FeaturePlot(o., features = "SST", order = T)
FeaturePlot(o., features = "VIP", order = T)
FeaturePlot(o., features = "PVALB", order = T)

FindMarkers(o., ident.1 = 0, ident.2 = NULL, group.by = "RNA_snn_res.0.1", only.pos = T)

o.@meta.data$cell_subtype_MV <- "Astrocyte"

saveRDS(o.,"E:/LOY/scLOY/processed_seurat/CELL_TYPE_SPECIFIC/ASTROCYTE/GSE174332_BRAIN_MOTOR_CORTEX_ASTROCYTE.RDS")

############# Microglia ################
DimPlot(o, group.by = "cell_ident_MV", label = T)
subset(o, cell_ident_MV == "Microglia") -> o.

# Scale data and run PCA
o. <- Seurat::NormalizeData(o., scale.factor = 10000)
o. <- Seurat::FindVariableFeatures(o., selection.method = "vst", nfeatures = 2000)
o. <- Seurat::ScaleData(o., vars.to.regress = c("nCount_RNA","percent.mt","percent.rb"))
v <- Seurat::VariableFeatures(o.)
v <- v[!(v %in% c(ribo,mito,male_genes,"XIST"))]
o. <- Seurat::RunPCA(o., npcs = 50, features = v, assay = "RNA")
ElbowPlot(o., ndims = 50)

o. <- Seurat::RunUMAP(o.,reduction = "pca", dims = 1:20)
o. <- Seurat::FindNeighbors(o.,reduction = "pca", dims = 1:20)
o. <- Seurat::FindClusters(o.,resolution = 0.1)
DimPlot(o., group.by = "donor_organism.biomaterial_core.biomaterial_id", label = T)
Idents(o.) <- "RNA_snn_res.0.1"

### 
o. <- harmony::RunHarmony(o. , theta = c(2), lambda  = c(1),
                          group.by.vars = c("donor_organism.biomaterial_core.biomaterial_id"),
                          assay.use = "RNA",
                          plot_convergence = T )
o. <- Seurat::RunUMAP(o.,reduction = "harmony", dims = 1:30)
o. <- Seurat::FindNeighbors(o.,reduction = "harmony", dims = 1:30)
o. <- Seurat::FindClusters(o.,resolution = seq(0.1,2,0.1))
DimPlot(o., group.by = "donor_organism.biomaterial_core.biomaterial_id", label = T)
DimPlot(o., group.by = "RNA_snn_res.0.1", label = T)

FeaturePlot(o., features = "CD74", order = T)
FeaturePlot(o., features = "CD163", order = T)
FeaturePlot(o., features = "C3", order = T)
FeaturePlot(o., features = "CX3CR1", order = T)
FeaturePlot(o., features = "FCN1", order = T)

FindMarkers(o., ident.1 = 0, ident.2 = NULL, group.by = "RNA_snn_res.0.1", only.pos = T)

o.@meta.data$cell_subtype_MV <- o.@meta.data$RNA_snn_res.0.1
levels(o.@meta.data$cell_subtype_MV)[levels(o.@meta.data$cell_subtype_MV)%in%c(0,2)] <- "Microglia"
levels(o.@meta.data$cell_subtype_MV)[levels(o.@meta.data$cell_subtype_MV)%in%c(3)] <- "CAM"
levels(o.@meta.data$cell_subtype_MV)[levels(o.@meta.data$cell_subtype_MV)%in%c(1)] <- "Cycling"
DimPlot(o., group.by = "cell_subtype_MV", label = T)

saveRDS(o.,"E:/LOY/scLOY/processed_seurat/CELL_TYPE_SPECIFIC/MICROGLIA/GSE174332_MOTOR_CORTEX_MICROGLIA.RDS")

############# Oligo ################
DimPlot(o, group.by = "cell_ident_MV", label = T)
subset(o, cell_ident_MV == "Oligo") -> o.

# Scale data and run PCA
o. <- Seurat::NormalizeData(o., scale.factor = 10000)
o. <- Seurat::FindVariableFeatures(o., selection.method = "vst", nfeatures = 2000)
o. <- Seurat::ScaleData(o., vars.to.regress = c("nCount_RNA","percent.mt","percent.rb"))
v <- Seurat::VariableFeatures(o.)
v <- v[!(v %in% c(ribo,mito,male_genes,"XIST"))]
o. <- Seurat::RunPCA(o., npcs = 50, features = v, assay = "RNA")
ElbowPlot(o., ndims = 50)

o. <- Seurat::RunUMAP(o.,reduction = "pca", dims = 1:20)
o. <- Seurat::FindNeighbors(o.,reduction = "pca", dims = 1:20)
o. <- Seurat::FindClusters(o.,resolution = 0.1)
DimPlot(o., group.by = "donor_organism.biomaterial_core.biomaterial_id", label = T)
Idents(o.) <- "RNA_snn_res.0.1"

### 
o. <- harmony::RunHarmony(o. , theta = c(2), lambda  = c(1),
                          group.by.vars = c("donor_organism.biomaterial_core.biomaterial_id"),
                          assay.use = "RNA",
                          plot_convergence = T )
o. <- Seurat::RunUMAP(o.,reduction = "harmony", dims = 1:20)
o. <- Seurat::FindNeighbors(o.,reduction = "harmony", dims = 1:20)
o. <- Seurat::FindClusters(o.,resolution = seq(0.1,2,0.1))
DimPlot(o., group.by = "donor_organism.biomaterial_core.biomaterial_id", label = T)
DimPlot(o., group.by = "RNA_snn_res.0.1", label = T)

FeaturePlot(o., features = "MBP", order = T)
FeaturePlot(o., features = "ACTA2", order = T)


FindMarkers(o., ident.1 = 0, ident.2 = NULL, group.by = "RNA_snn_res.0.1", only.pos = T)

o.@meta.data$cell_subtype_MV <- o.@meta.data$RNA_snn_res.0.1
levels(o.@meta.data$cell_subtype_MV)[levels(o.@meta.data$cell_subtype_MV)%in%c()] <- "Oligo"
DimPlot(o., group.by = "cell_subtype_MV", label = T)

saveRDS(o.,"E:/LOY/scLOY/processed_seurat/CELL_TYPE_SPECIFIC/OLIGO/GSE174332_MOTOR_CORTEX_OLIGO.RDS")

############# OPC ################
o<-readRDS("E:/LOY/scLOY/processed_seurat/CELL_TYPE_SPECIFIC/")
subset(o, cell_ident_MV == "OPC") -> o.

# Scale data and run PCA
o. <- Seurat::NormalizeData(o., scale.factor = 10000)
o. <- Seurat::FindVariableFeatures(o., selection.method = "vst", nfeatures = 2000)
o. <- Seurat::ScaleData(o., vars.to.regress = c("nCount_RNA","percent.mt","percent.rb"))
v <- Seurat::VariableFeatures(o.)
v <- v[!(v %in% c(ribo,mito,male_genes,"XIST"))]
o. <- Seurat::RunPCA(o., npcs = 50, features = v, assay = "RNA")
ElbowPlot(o., ndims = 50)

o. <- Seurat::RunUMAP(o.,reduction = "pca", dims = 1:20)
o. <- Seurat::FindNeighbors(o.,reduction = "pca", dims = 1:20)
o. <- Seurat::FindClusters(o.,resolution = 0.1)
DimPlot(o., group.by = "donor_organism.biomaterial_core.biomaterial_id", label = T)
Idents(o.) <- "RNA_snn_res.0.1"

### 
o. <- harmony::RunHarmony(o. , theta = c(2), lambda  = c(1),
                          group.by.vars = c("donor_organism.biomaterial_core.biomaterial_id"),
                          assay.use = "RNA",
                          plot_convergence = T )
o. <- Seurat::RunUMAP(o.,reduction = "harmony", dims = 1:20)
o. <- Seurat::FindNeighbors(o.,reduction = "harmony", dims = 1:20)
o. <- Seurat::FindClusters(o.,resolution = seq(0.1,2,0.1))
DimPlot(o., group.by = "donor_organism.biomaterial_core.biomaterial_id", label = T)
DimPlot(o., group.by = "RNA_snn_res.0.1", label = T)

FeaturePlot(o., features = "VCAN", order = T)
FeaturePlot(o., features = "MBP", order = T)


FindMarkers(o., ident.1 = 0, ident.2 = NULL, group.by = "RNA_snn_res.0.1", only.pos = T)

o.@meta.data$cell_subtype_MV <- o.@meta.data$RNA_snn_res.0.1
levels(o.@meta.data$cell_subtype_MV)[levels(o.@meta.data$cell_subtype_MV)%in%c(0,1)] <- "OPC"
levels(o.@meta.data$cell_subtype_MV)[levels(o.@meta.data$cell_subtype_MV)%in%c()] <- "Doublet"
DimPlot(o., group.by = "cell_subtype_MV", label = T)

saveRDS(o.,"E:/LOY/scLOY/processed_seurat/CELL_TYPE_SPECIFIC/OPC/GSE174332_MOTOR_CORTEX_OPC.RDS")



