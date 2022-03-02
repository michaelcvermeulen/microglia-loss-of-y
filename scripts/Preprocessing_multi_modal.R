### Preprocessing multimodal data for 

library(Seurat)
library(magrittr)
library(data.table)
library(Signac)
library(Biobase)
library(GEOquery)
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


### add link 
Read10X_h5("e:/LOY/scLOY/data/10X_MULTIOME_LYMPH/lymph_node_lymphoma_14k_filtered_feature_bc_matrix.h5") -> h5

# fragment file
frag <- "e:/LOY/scLOY/data/10X_MULTIOME_LYMPH/lymph_node_lymphoma_14k_atac_fragments.tsv.gz"


# take Gene expression from the multi modal object and create a Seurat object
h5$`Gene Expression` -> expr
CreateSeuratObject(counts = expr, min.cells = 10, min.features = 200 ) -> o
o[["percent.mt"]] <- PercentageFeatureSet(o, pattern = "^MT-")
NormalizeData(o) -> o

# This sample is a male, all cells are male
o@meta.data$donor_organism.sex <- "male"


# find suitable chrY genes and run LOY function
LOY_genes_to_exclude_expr(o) -> remove
call_LOY(s_obj = o, genes_to_exclude = remove) -> o

# remove low quality cells
subset(o, nCount_RNA > 1000 & nFeature_RNA > 800) -> o

## Preprocess chromatin peaks next 
h5$Peaks -> peaks

grange.counts <- StringToGRanges(rownames(peaks), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
peaks <- peaks[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = edb)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"

peaks[grepl(x = rownames(peaks), pattern = "chrY"),] %>% colSums() -> chrY

dplyr::left_join( x = tibble::rownames_to_column(o@meta.data,"CB"),
                  y = tibble::rownames_to_column(as.data.frame(chrY),"CB") ,
                  by = "CB") %>% tibble::column_to_rownames("CB") -> o@meta.data

o@meta.data %>% group_by(LOY) %>% summarise(median = median(chrY))

chrom_assay <- Signac::CreateChromatinAssay(
  counts = peaks,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag,
  annotation = annotations
)

#o$blacklist_fraction <- FractionCountsInRegion(
#  object = o,
#  assay = 'peaks',
#  regions = blacklist_hg38_unified
#)

colnames(chrom_assay)[colnames(chrom_assay) %in% colnames(o)] -> keep

subset(chrom_assay, cells = keep ) -> chrom_assay
subset(o, cells = keep) -> o

o[["ATAC"]] <- chrom_assay

dplyr::mutate(o@meta.data, chrY_ratio_ATAC = chrY/nCount_ATAC) -> o@meta.data

ifelse(o@meta.data$LOY == "LOY" & o@meta.data$chrY == 0, "LOY", "NORMAL") -> o@meta.data$multimodal_LOY

FindMarkers(subset(o, nCount_RNA > 1500), ident.1 = "LOY", ident.2 = "NORMAL", group.by = "multimodal_LOY")


DefaultAssay(o) <- "ATAC"
o <- NucleosomeSignal(o)
o <- TSSEnrichment(o)


VlnPlot(
  object = o,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

o <- subset(
  x = o,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 20000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 1.5 &
    TSS.enrichment > 1
)

### install MACS3


# RNA analysis
DefaultAssay(o) <- "RNA"

# ribo, mito
rownames(o)[grep(x = rownames(o), pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")] -> ribo
ribo[grep(ribo,pattern = "K",invert = T)] -> ribo
rownames(o)[grep(x = rownames(o), pattern = "^MT-")] -> mito
rownames(o)[(rownames(o) %in% Y_scLOY_genes$gene_name)] -> Y
Y_scLOY_genes$gene_name[Y_scLOY_genes$gene_name %in% rownames(o)] -> male_genes

o[["percent.mt"]] <- Seurat::PercentageFeatureSet(o, pattern = "^MT-")
o[["percent.rb"]] <- Seurat::PercentageFeatureSet(o, features = ribo)

o <- Seurat::NormalizeData(o, scale.factor = 10000)
o <- Seurat::FindVariableFeatures(o, selection.method = "vst", nfeatures = 2500)
# dissociation gene module
diss <- c("RGS1","DUSP1","AC022217.3","DUSP1","HSPA1A","HSPA1B","DDIT4","SLC2A3","FOS","JUN",
          "AC131944.1","HIF1A-AS3","HIST1H2AC","HIST1H2BD","HIST1H2BG","HIST2H2BE","RGS2","RHOB",
          "HIST1H2BJ","TEX14","UBC","HSP90AA1","AL118516.1","JUNB","HIST2H2AA4","HIST1H1C","BTG1","BTG2","DUSP5")
Seurat::AddModuleScore(object = o, features = diss, name = "diss") -> o

# PAR and Y modules
list(PAR_scLOY_genes$gene_name[PAR_scLOY_genes$gene_name %in% rownames(o)]) -> P
Seurat::AddModuleScore(object = o, features = P, name = "PAR_score") -> o
list(Y_scLOY_genes$gene_name[Y_scLOY_genes$gene_name %in% rownames(o)]) -> Y
Seurat::AddModuleScore(object = o, features = Y, name = "Y_score") -> o

# Scale data and run PCA
o <- Seurat::ScaleData(o, vars.to.regress = c("nCount_RNA","percent.mt","percent.rb"))
v <- Seurat::VariableFeatures(o)
v <- v[!(v %in% c(ribo,mito,male_genes,"XIST","JPX"))]
o <- Seurat::RunPCA(o, npcs = 50, features = v, assay = "RNA")
ElbowPlot(o, ndims = 50)

o <- Seurat::RunUMAP(o,reduction = "pca", dims = 1:20)
o <- Seurat::FindNeighbors(o,reduction = "pca", dims = 1:20)
o <- Seurat::FindClusters(o,resolution = 0.1)
DimPlot(o, group.by = "RNA_snn_res.0.1", label = T)
Idents(o) <- "RNA_snn_res.0.1"

FindMarkers(o, ident.1 = 1, ident.2 = NULL, group.by = "RNA_snn_res.0.1")

DimPlot(subset(o, nCount_RNA > 1500), group.by = "LOY")

FindMarkers(subset(o, nCount_RNA > 1500 & nFeature_RNA >= 1000 & RNA_snn_res.0.1 == 1), group.by = "LOY",
            ident.1 = "LOY", ident.2 = "NORMAL", test.use = "MAST", latent.vars = c("nCount_RNA","nFeature_RNA","percent.rb"))

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(o) <- "ATAC"
o <- RunTFIDF(o)
o <- FindTopFeatures(o, min.cutoff = "q25")
o <- RunSVD(o)
o <- RunUMAP(o, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

library(BSgenome.Hsapiens.UCSC.hg38)
o <- RegionStats(o, genome = BSgenome.Hsapiens.UCSC.hg38)

gene.activities <- GeneActivity(o)

o[['ATAC_gene_activity']] <- CreateAssayObject(counts = gene.activities)
o <- NormalizeData(
  object = o,
  assay = 'ATAC_gene_activity',
  normalization.method = 'LogNormalize',
  scale.factor = median(o$nCount_RNA)
)

o <- FindMultiModalNeighbors(o, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
o <- RunUMAP(o, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
o <- FindClusters(o, graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution = 0.1)

DefaultAssay(o) <- "ATAC"
fc <- FoldChange(subset(o, RNA_snn_res.0.1 == 1 & nCount_RNA > 1500), ident.1 = "4", ident.2 = "0", group.by = "LOY_levels")
fc <- dplyr::mutate(fc, mean_0_4 = ( mean_4 / mean_0 ))
fc[grepl( x = rownames(fc), pattern = "chrX" ),] -> dat
subset(o, RNA_snn_res.0.1 == 1 & nCount_RNA > 1500) -> o.; o.@meta.data$LOY %>% table()

t <- rownames(dat[dat$log2_fold_change > 0.5 & dat$mean_NORMAL > 0.2, ])

closest_genes_t <- ClosestFeature(o, regions = t)
ClosestFeature(o, regions = "chr15-33158297-33158868")


closest_genes_cd14mono <- ClosestFeature(pbmc, regions = open_cd14mono)


p3 <- DimPlot(o, reduction = "wnn.umap",label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

o <- LinkPeaks(
  object = o,
  peak.assay = "ATAC",
  expression.assay = "SCT",
  genes.use = c("LYZ", "MS4A1")
)


p1 <- CoveragePlot(
  object = o,
  region = "MS4A1",
  features = "MS4A1",
  expression.assay = "SCT",group.by = "LOY",
  extend.upstream = 500,
  extend.downstream = 10000
)




DefaultAssay(o) <- 'ATAC'
Idents(o) <- "LOY"
da_peaks <- FindMarkers(
  object = o,
  ident.1 = "LOY",
  ident.2 = "NORMAL",
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)


head(da_peaks)


DefaultAssay(o) <- 'RNA'
NormalizeData(o) -> o
Idents(o) <- "LOY"
da_peaks <- FindMarkers(
  object = o,
  ident.1 = "LOY",
  ident.2 = "NORMAL",
  min.pct = 0.05
)


gene_plot <- AnnotationPlot(
  object = o,
  region = "UTY"
)

tile_plot <- TilePlot(
  object = o,
  region = "chr2-87011729-87035519",
  idents = c("0", "1")
)

Idents(o) <- "LOY"
CoveragePlot(
  object = o,
  region = "SLC25A6",
  expression.assay = "RNA", features = "SLC25A6",
  extend.upstream = 40000,
  extend.downstream = 20000
)

Idents(o) <- "LOY"
CoveragePlot(
  object = o,
  region = "CD99",
  expression.assay = "RNA", features = "CD99",
  extend.upstream = 40000,
  extend.downstream = 20000
)

Idents(o) <- "LOY"
CoveragePlot(
  object = subset(o, nCount_RNA > 2500),
  region = "UTY",
  expression.assay = "RNA", features = "UTY",
  extend.upstream = 20000,
  extend.downstream = 20000
)

TilePlot(object = o, region = "UTY", group.by = "LOY", tile.size = 5000 )
TilePlot(object = subset(o, nCount_RNA > 2500), region = "UTY", group.by = "LOY", tile.size = 2500 ) + scale_fill_continuous()

DefaultAssay(o) <- "ATAC"
Idents(o) <- "LOY"
CoveragePlot(
  object = subset(o, nCount_RNA > 2500),
  region = "USP9Y", window = 500,
  expression.assay = "RNA", features = "USP9Y",
  extend.upstream = 100000,
  extend.downstream = 100000
)

DefaultAssay(o) <- "ATAC"
Idents(o) <- "LOY"
CoveragePlot(
  object = subset(o, nCount_RNA > 3000),
  region = "P2RY8",
  expression.assay = "RNA", features = "P2RY8",
  extend.upstream = 2000,
  extend.downstream = 2000
)

DefaultAssay(o) <- "ATAC"
Idents(o) <- "RNA_snn_res.0.1"
CoveragePlot(
  object = o,
  region = "MKI67",
  expression.assay = "RNA", features = "MKI67",
  extend.upstream = 10000,
  extend.downstream = 20000
)

library(BSgenome.Hsapiens.UCSC.hg38)
o <- RegionStats(o, genome = BSgenome.Hsapiens.UCSC.hg38)
# link peaks to genes
o <- LinkPeaks(
  object = o,
  peak.assay = "ATAC",
  expression.assay = "RNA",
  genes.use = c("UTY")
)

Idents(o) <- "LOY"
p1 <- CoveragePlot(
  object = subset(o, nCount_RNA > 1500),
  region = "UTY",
  features = "UTY",
  expression.assay = "RNA",
  extend.upstream = 10000,
  extend.downstream = 10000
)

p2 <- CoveragePlot(
  object = subset(o, nCount_RNA > 1500),
  region = "USP9Y",
  features = "USP9Y",
  expression.assay = "RNA",
  extend.upstream = 10000,
  extend.downstream = 10000
)

p3 <- CoveragePlot(
  object = subset(o, nCount_RNA > 1500),
  region = "LINC00278",
  features = "LINC00278",
  expression.assay = "RNA",
  extend.upstream = 10000,
  extend.downstream = 10000
)


#### NORM ACCESS
o[["ATAC_gene_activity"]]@counts -> counts
o[["ATAC_gene_activity"]]@data -> data

counts %>% as.matrix() -> counts
counts[rownames(counts) %in% expressed_Y_genes(o),] -> Y_counts
Y_counts %>% colSums() %>% as.data.frame() -> Y
left_join(x = tibble::rownames_to_column(o@meta.data,"CB"), y = tibble::rownames_to_column(Y,"CB"), by = "CB" ) %>%
  tibble::column_to_rownames("CB") -> o@meta.data

ggboxplot(o@meta.data[o@meta.data$nCount_RNA > 2500,], x = "LOY", y = ".")
o@meta.data[o@meta.data$nCount_RNA > 2500,] %>% dplyr::group_by(LOY) %>% summarise(med = mean(.))

o[["ATAC_gene_activity"]]@data -> data
data %>% as.matrix() -> data
data[rownames(data) %in% expressed_Y_genes(o),] -> Y_counts
Y_counts %>% colMeans() %>% as.data.frame() -> Y
names(Y) <- "Norm_Access"
left_join(x = tibble::rownames_to_column(o@meta.data,"CB"), y = tibble::rownames_to_column(Y,"CB"), by = "CB" ) %>%
  tibble::column_to_rownames("CB") -> o@meta.data