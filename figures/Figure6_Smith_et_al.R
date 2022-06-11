### Genome Research
### Michael Vermeulen
### Mosaic loss of chromosome Y in aged human microglia

### This is the template used to process all microglia data
### Here, preprocessing and microglia annotation of GSE160936 SSC and EC brain tissue data is provided. 
### All other preprocessing + annotation scripts are available upon request. 

library(celldex)
library(monocle3)
library(monocle)
library(scDblFinder)
library(SingleR)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)
library(ggpubr)
library(nichenetr)
library(dittoSeq)
library(escape)
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
edb <- EnsDb.Hsapiens.v86

setwd("H:/LOY/scLOY/scLOY/")
source("E:/LOY/Example/PACKAGE/call_LOY.R")
source("E:/LOY/scLOY/scripts/loy_functions.R")
source('H:/LOY/scLOY/scLOY/R/load.R')
source('H:/LOY/scLOY/scLOY/data_raw/build_data_files.R')
###


# color palette
vega <- c('#1f77b4','#aec7e8','#ff7f0e',
          '#ffbb78','#2ca02c','#98df8a',
          '#d62728','#ff9896','#9467bd',
          '#c5b0d5','#8c564b','#c49c94',
          '#e377c2','#f7b6d2','#7f7f7f',
          '#c7c7c7','#bcbd22','#dbdb8d',
          '#17becf','#9edae5')

cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

# microglia dissociation associated genes
diss <- c("RGS1","DUSP1","AC022217.3","DUSP1","HSPA1A","HSPA1B","DDIT4","SLC2A3","FOS","JUN",
          "AC131944.1","HIF1A-AS3","HIST1H2AC","HIST1H2BD","HIST1H2BG","HIST2H2BE","RGS2","RHOB",
          "HIST1H2BJ","TEX14","UBC","HSP90AA1","AL118516.1","JUNB","HIST2H2AA4","HIST1H1C","BTG1","BTG2","DUSP5")

# microglia gene sets from previous publications
# https://github.com/michaelcvermeulen/microglia-loss-of-y/blob/main/data/microglia_gene_sets.csv
read.csv("h:/LOY/scLOY/data/microglia_gene_sets.csv") -> gene_sets


## read in preprocessed GSE160936 seurat object.
## processed using Preprocessing_scrnaseq.R script template. 
## https://github.com/michaelcvermeulen/microglia-loss-of-y/blob/main/scripts/Preprocessing_scrnaseq.R
readRDS(file = "E:/LOY/scLOY/processed_seurat/BRAIN/GSE160936_AD_GLIAL_introns_included.RDS") -> o

### ALL #############
subset(o, 
       nCount_RNA >= 1000 & 
         nFeature_RNA >= 800 & 
         cell_ident_MV == "Microglia" ) -> o

### run seurat integration
ifnb.list <- SplitObject(o, split.by = "title")

ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2500)
  return(x)
})

rownames(o)[grep(x = rownames(o), pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")] -> ribo
ribo[grep(ribo,pattern = "K",invert = T)] -> ribo
rownames(o)[grep(x = rownames(o), pattern = "^MT-")] -> mito
rownames(o)[(rownames(o) %in% Y_scLOY_genes$gene_name)] -> Y
Y_scLOY_genes$gene_name[Y_scLOY_genes$gene_name %in% rownames(o)] -> male_genes

features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 2500)
features <- features[!(features %in% c("XIST",ribo,mito,male_genes,"JPX"))]
anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
combined <- IntegrateData(anchorset = anchors)
DefaultAssay(combined) <- "integrated"
combined -> o
rm(combined)

rownames(o)[grep(x = rownames(o), pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")] -> ribo
ribo[grep(ribo,pattern = "K",invert = T)] -> ribo

o <- Seurat::ScaleData(o, vars.to.regress = c("nCount_RNA","percent.mt","percent.rb"))
v <- Seurat::VariableFeatures(o)
v <- v[!(v %in% c(ribo,mito,male_genes,"XIST","JPX"))]
o <- Seurat::RunPCA(o, npcs = 50, features = v, assay = "integrated")
ElbowPlot(o,ndims = 50 )

o <- Seurat::RunUMAP(o,reduction = "pca", dims = 1:20)
o <- Seurat::FindNeighbors(o,reduction = "pca", dims = 1:20)
o <- Seurat::FindClusters(o,resolution = seq(0.1,2,0.1))
DimPlot(o, group.by = "integrated_snn_res.0.5", label = T, cols = dittoColors())
DimPlot(o, group.by = "title")

FindMarkers(o, ident.1 = 10, ident.2 = NULL, group.by = "integrated_snn_res.0.5",
            only.pos = T, logfc.threshold = 0.75)

Idents(o) <- "integrated_snn_res.0.5"
RenameIdents(o, '0' = '0', '7' = '1', '4' = '2', 
             '6' = '3', '1' = '4', '5' = '5',
             '3' = '6', '2' = '7', '9' = '8', '8' = '9', '10' = '10',
             '11' = '11') -> o
Idents(o) -> o@meta.data$integrated_snn_res.0.5_rearranged
DimPlot(o, group.by = "integrated_snn_res.0.5_rearranged", label = T)
Idents(o) <- "integrated_snn_res.0.3"

saveRDS(o, "e:/LOY/scLOY/processed_integrated/MICROGLIA/IMPERIAL_CCA_JULY6.RDS")
### HOMEOSTATIC GENES ##########
output <- "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/FINAL/IMPERIAL/FULL"
dir.create(paste0(output,"/SUBCLUSTERS"), recursive = T)
dir.create(paste0(output,"/UMAP"), recursive = T)
dir.create(paste0(output,"/BARPLOTS"), recursive = T)

DefaultAssay(o) <- "RNA"
p <- plot_density(o, c("CX3CR1","P2RY12","CSF1R"), reduction = "umap",
                  method = "wkde", adjust = 3,
                  pal = "inferno", joint = T,
                  combine = F)[[4]];
ggsave(plot = p, dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/HOMEO_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(),
       dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/HOMEO.pdf"))

p <- plot_density(o, c("CX3CR1","P2RY12","CSF1R"), reduction = "umap",
                  method = "wkde", adjust = 2,
                  pal = "inferno", joint = T,
                  combine = F)[[1]];
ggsave(plot = p, dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/CX3CR1_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(), dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/CX3CR1.pdf"))

p <- plot_density(o, c("CX3CR1","P2RY12","CSF1R"), reduction = "umap",
                  method = "wkde", adjust = 2,
                  pal = "inferno", joint = T,
                  combine = F)[[2]];
ggsave(plot = p, dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/P2RY12_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(),
       dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/P2RY12.pdf"))

p <- plot_density(o, c("CX3CR1","P2RY12","CSF1R"), reduction = "umap",
                  method = "wkde", adjust = 2,
                  pal = "inferno", joint = T,
                  combine = F)[[3]];
ggsave(plot = p,  dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/CSF1R_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(),
       dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/CSF1R.pdf"))


### GENE SETS #############
### find homeostatic and AD associated clusters
### pull these clusters

a <- gene_sets$SalaFrigerio_ARM
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "ARM", assay = "RNA") -> o
plot_density(o, "ARM1", reduction = "umap", pal = "inferno", adjust = 1) -> p
FeaturePlot(o, features = "ARM1", order = T, pt.size = 0.9) + scale_color_viridis_c(option = "B")
ggsave(plot = p + NoAxes(), dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/ARM_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(), dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/ARM.pdf"))

a <- gene_sets$Cytokine_Response
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "Cytokine_Response", assay = "RNA") -> o
plot_density(o, "Cytokine_Response1", reduction = "umap", pal = "inferno") -> p
ggsave(plot = p + NoAxes(), dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/Cytokine_Response_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(), dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/Cytokine_Response.pdf"))

plot_density(o, features = c("CCL3","CCL4","IL1B"), method = "wkde", adjust = 2,
             reduction = "umap",
             pal = "inferno", joint = T, combine = F) -> p
p[[4]] + NoAxes() -> p
ggsave(plot = p + NoAxes(), dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/CCL3_CCL4_IL1B_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(), dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/CCL3_CCL4_IL1B.pdf"))


a <- gene_sets$Prinz_CAM
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "CAM", assay = "RNA") -> o
plot_density(o, "CAM1", reduction = "umap", pal = "inferno") -> p
ggsave(plot = p, dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/CAM_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(),
       dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/CAM.pdf"))

a <- gene_sets$DAM_UP_SOBUE
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "DAM_UP_SOBUE", assay = "RNA") -> o
plot_density(o, "DAM_UP_SOBUE1", reduction = "umap", pal = "inferno") -> p
ggsave(plot = p + NoAxes(), dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/DAM_UP_SOBUE_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(),
       dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/DAM_UP_SOBUE.pdf"))

a <- gene_sets$DAM_UP_KRASEMANN_mouse %>% nichenetr::convert_mouse_to_human_symbols()
a[!is.na(a)] -> a
names(a) <- NULL
a %>% as.vector() -> a
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "DAM_UP_KRASEMAN", assay = "RNA") -> o
plot_density(o, "DAM_UP_KRASEMAN1", reduction = "umap", pal = "inferno") -> p
ggsave(plot = p+ NoAxes(), dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/DAM_UP_KRASEMAN_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(),
       dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/DAM_UP_KRASEMAN.pdf"))

FeaturePlot(o, features = "AGE_UP_MICROGLIA1", pt.size = 0.8, order = T) + scale_color_viridis_c(option = "B") -> p
ggsave(plot = p+ NoAxes(), dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/DAM_UP_KRASEMAN_LEGEND_FEATUREPLOT.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(), dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/DAM_UP_KRASEMAN_LEGEND_FEATUREPLOT.pdf"))


a <- gene_sets$AGE_UP_MICROGLIA
a[!is.na(a)] -> a
names(a) <- NULL
a %>% as.vector() -> a
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "AGE_UP_MICROGLIA", assay = "RNA") -> o
plot_density(o, "AGE_UP_MICROGLIA1", reduction = "umap", pal = "inferno") -> p
ggsave(plot = p+ NoAxes(), dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/AGE_UP_MICROGLIA_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(), dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/AGE_UP_MICROGLIA.pdf"))

FeaturePlot(o, features = "AGE_UP_MICROGLIA1", pt.size = 0.8, order = T) + scale_color_viridis_c(option = "B") -> p
ggsave(plot = p+ NoAxes(), dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/AGE_UP_MICROGLIA_LEGEND_FEATUREPLOT.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(), dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/AGE_UP_MICROGLIA_FEATUREPLOT.pdf"))

a <- gene_sets$HOMEOSTATIC_KRASEMANN_mouse %>% nichenetr::convert_mouse_to_human_symbols()
a[!is.na(a)] -> a
names(a) <- NULL
a %>% as.vector() -> a
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "HOMEOSTATIC_KRASEMANN", assay = "RNA") -> o
plot_density(o, "HOMEOSTATIC_KRASEMANN1", reduction = "umap", pal = "inferno") -> p
ggsave(dpi = "retina", units = "in", width = 4.5, height = 3.5,
       plot = p+ NoAxes(),
       filename = paste0(output,"/SUBCLUSTERS/HOMEOSTATIC_KRASEMANN_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(),
       dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/HOMEOSTATIC_KRASEMANN.pdf"))

a <- gene_sets$AD_traj_late
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "AD_TRAJ_LATE", assay = "RNA") -> o
plot_density(o, "AD_TRAJ_LATE1", reduction = "umap", pal = "inferno") -> p
ggsave(dpi = "retina", units = "in", width = 4.5, height = 3.5,
       plot = p+ NoAxes(),
       filename = paste0(output,"/SUBCLUSTERS/AD_TRAJ_LATE_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(),
       dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/AD_TRAJ_LATE.pdf"))

a <- gene_sets$AD_traj_conversion
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "AD_TRAJ_CONVERSION", assay = "RNA") -> o
plot_density(o, "AD_TRAJ_CONVERSION1", reduction = "umap", pal = "inferno") -> p
ggsave(dpi = "retina", units = "in", width = 4.5, height = 3.5,
       plot = p+ NoAxes(),
       filename = paste0(output,"/SUBCLUSTERS/AD_TRAJ_CONVERSION_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(),
       dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/AD_TRAJ_CONVERSION.pdf"))

a <- gene_sets$AD_traj_early_end
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "AD_TRAJ_EARLY_END", assay = "RNA") -> o
plot_density(o, "AD_TRAJ_EARLY_END1", reduction = "umap", pal = "inferno") -> p
ggsave(dpi = "retina", units = "in", width = 4.5, height = 3.5,
       plot = p+ NoAxes(),
       filename = paste0(output,"/SUBCLUSTERS/AD_TRAJ_EARLY_END_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(),
       dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/AD_TRAJ_EARLY_END.pdf"))

a <- gene_sets$AD_traj_early
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "AD_TRAJ_EARLY", assay = "RNA") -> o
plot_density(o, "AD_TRAJ_EARLY1", reduction = "umap", pal = "inferno") -> p
ggsave(dpi = "retina", units = "in", width = 4.5, height = 3.5,
       plot = p+ NoAxes(),
       filename = paste0(output,"/SUBCLUSTERS/AD_TRAJ_EARLY_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(),
       dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/AD_TRAJ_EARLY.pdf"))

a <- gene_sets$SalaFrigerio_HM
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "SalaFrigerio_HM", assay = "RNA") -> o
plot_density(o, "HM1", reduction = "umap", pal = "inferno") -> p
ggsave(plot = p+ NoAxes(), dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/HM_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(), dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/HM.pdf"))

plot_density(o, features = "CDKN2A", reduction = "umap", pal = "inferno", method = "wkde", adjust = 2) -> p
ggsave(plot = p+ NoAxes(), dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/SENES_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(), dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/SENES.pdf"))

a <- gene_sets$HAM_UP_Srinivasan_2021
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "HAM_UP_Srinivasan_2021", assay = "RNA") -> o
plot_density(subset(o, diagnosis == "AD"), "HAM_UP_Srinivasan_20211", reduction = "umap", pal = "inferno") -> p
ggsave(plot = p+ NoAxes(), dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/HM_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(), dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/HM.pdf"))





###### ANNOTATE #########


FindMarkers(o, ident.1 = 10, ident.2 = NULL, group.by = "integrated_snn_res.0.5", only.pos = T, logfc.threshold = 0.5)

o@meta.data$MG_subtype <- o@meta.data$integrated_snn_res.0.5
levels(o@meta.data$MG_subtype)[levels(o@meta.data$MG_subtype)%in%c(8)] <- "CAM"
levels(o@meta.data$MG_subtype)[levels(o@meta.data$MG_subtype)%in%c(0,1,4)] <- "Homeostatic"
levels(o@meta.data$MG_subtype)[levels(o@meta.data$MG_subtype)%in%c(6,2)] <- "Pro-inflammatory"
levels(o@meta.data$MG_subtype)[levels(o@meta.data$MG_subtype)%in%c(5,3)] <- "DAM"
levels(o@meta.data$MG_subtype)[levels(o@meta.data$MG_subtype)%in%c(9)] <- "CRM"
levels(o@meta.data$MG_subtype)[levels(o@meta.data$MG_subtype)%in%c(7)] <- "Cycling"
levels(o@meta.data$MG_subtype)[levels(o@meta.data$MG_subtype)%in%c(10)] <- "Monocyte"
levels(o@meta.data$MG_subtype)[levels(o@meta.data$MG_subtype)%in%c()] <- "Plasma"
levels(o@meta.data$MG_subtype)[levels(o@meta.data$MG_subtype)%in%c()] <- "Mast"
levels(o@meta.data$MG_subtype)[levels(o@meta.data$MG_subtype)%in%c()] <- "Doublet"
levels(o@meta.data$MG_subtype)[levels(o@meta.data$MG_subtype)%in%c()] <- "B"
DimPlot(o, group.by = "MG_subtype", label = T)

o@meta.data$MG_subtype_C <- o@meta.data$integrated_snn_res.0.5
levels(o@meta.data$MG_subtype_C)[levels(o@meta.data$MG_subtype_C)%in%c(8)] <- "CAM"
levels(o@meta.data$MG_subtype_C)[levels(o@meta.data$MG_subtype_C)%in%c(0)] <- "H1"
levels(o@meta.data$MG_subtype_C)[levels(o@meta.data$MG_subtype_C)%in%c(1)] <- "H2"
levels(o@meta.data$MG_subtype_C)[levels(o@meta.data$MG_subtype_C)%in%c(4)] <- "H3"
levels(o@meta.data$MG_subtype_C)[levels(o@meta.data$MG_subtype_C)%in%c(2)] <- "I1"
levels(o@meta.data$MG_subtype_C)[levels(o@meta.data$MG_subtype_C)%in%c(6)] <- "I2"
levels(o@meta.data$MG_subtype_C)[levels(o@meta.data$MG_subtype_C)%in%c(10)] <- "Monocyte"
levels(o@meta.data$MG_subtype_C)[levels(o@meta.data$MG_subtype_C)%in%c(7)] <- "Cycling"
levels(o@meta.data$MG_subtype_C)[levels(o@meta.data$MG_subtype_C)%in%c(5)] <- "DAM1"
levels(o@meta.data$MG_subtype_C)[levels(o@meta.data$MG_subtype_C)%in%c(9)] <- "CRM"
levels(o@meta.data$MG_subtype_C)[levels(o@meta.data$MG_subtype_C)%in%c(3)] <- "DAM2"
DimPlot(o, group.by = "MG_subtype_C", label = T)

Idents(o) <- "MG_subtype_C"
FindAllMarkers(o, logfc.threshold = 0.3, min.pct = 0.1, test.use = "MAST",
               latent.vars = c("nCount_RNA","nFeature_RNA","percent.rb")) -> DE
data.table::fwrite(x = DE, file = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/FINAL/IMPERIAL/FULL/DE_final.txt",
                   sep = "\t")

factor(o@meta.data$MG_subtype_C, levels = c("H1","H2","H3","I1","I2",
                                            "DAM1","DAM2","CRM","Cycling","CAM","Monocyte")) -> o@meta.data$ord

colors <- c(dittoColors()[21],dittoColors()[5],"#89c2d9",
            dittoColors()[1], dittoColors()[17],
            dittoColors()[7], dittoColors()[23],
            "red",
            dittoColors()[11],
            "#7209B7",
            dittoColors()[16])

dittoDimPlot(o, var = "ord", color.panel = colors,  
             do.label = T, labels.highlight = T, labels.size = 4) -> p
p %>% ggpar(title = "") -> p
ggsave(plot = p + NoAxes() + NoLegend(), dpi = "retina", units = "in", width = 6, height = 5,
       filename = paste0(output,"/UMAP/NO_AXES_NO_LEG_UMAP.pdf"))

dittoDimPlot(o, var = "ord", color.panel = colors,  
             do.label = F, labels.highlight = T, labels.size = 4) -> p
p %>% ggpar(title = "") -> p
ggsave(plot = p + NoAxes(), dpi = "retina", units = "in", width = 6, height = 5,
       filename = paste0(output,"/UMAP/NO_AXES_UMAP.pdf"))

dittoDimPlot(o, var = "ord", color.panel = colors,split.by = "diagnosis", size = 0.8,
             do.label = F, labels.highlight = T, labels.size = 4, do.contour = T) -> p
p %>% ggpar(title = "") + NoAxes() + NoLegend() -> p
ggsave(plot = p + NoAxes(), dpi = "retina", units = "in", width = 5, height = 3,
       filename = paste0(output,"/UMAP/AD_SPLIT.pdf"))

dittoDimPlot(o, var = "ord", color.panel = colors,split.by = "specific_tissue", size = 0.8,
             do.label = F, labels.highlight = T, labels.size = 4, do.contour = T) -> p
p %>% ggpar(title = "") + NoAxes() + NoLegend() -> p
ggsave(plot = p + NoAxes(), dpi = "retina", units = "in", width = 5, height = 3,
       filename = paste0(output,"/UMAP/TISSUE_SPLIT.pdf"))

dittoDimPlot(o, var = "title", size = 0.8,
             do.label = F) -> p
p %>% ggpar(title = "") + NoAxes() -> p
ggsave(plot = p + NoAxes(), dpi = "retina", units = "in", width = 6, height = 4,
       filename = paste0(output,"/UMAP/SAMPLE_SPLIT.pdf"))

################# BAR PLOT  # 

DimPlot(o, group.by = "ord", label = T, cols = dittoColors())
DimPlot(o, split.by = "diagnosis", label = T, group.by = "integrated_snn_res.0.5", cols = dittoColors())

dittoBarPlot(o, group.by = "title", color.panel = colors, var.labels.reorder = c(6,7,8,9,10,4,5,2,3,1,11),
             var = "ord") %>% ggpar(xlab = "", title = "") -> a
ggsave(plot = a + NoLegend(), dpi = "retina", units = "in", width = 5, height = 3,
       filename = paste0(output,"/BARPLOTS/SAMPLE_PROP_BAR.pdf"))


dittoBarPlot(o, group.by = "title", var = "LOY", var.labels.reorder = c(2,1),
             color.panel = c("grey","red")) %>% ggpar(xlab = "") -> b 
ggarrange(plotlist = list(a,b), ncol = 1, nrow = 2, align = "hv")

dittoBarPlot(o, group.by = "integrated_snn_res.0.5", var = "LOY", var.labels.reorder = c(2,1),
             color.panel = c("grey","red")) %>% ggpar(xlab = "") -> c

dittoBarPlot(o, group.by = "diagnosis", var = "LOY", var.labels.reorder = c(2,1),
             color.panel = c("grey","red")) %>% ggpar(xlab = "") -> d

dittoBarPlot(o, group.by = "diagnosis", var.labels.reorder = c(6,7,8,9,10,4,5,2,3,1,11), color.panel = colors,
             var = "ord") %>% ggpar(xlab = "", title = "") -> f
ggsave(plot = e + NoLegend(), dpi = "retina", units = "in", width = 1.5, height = 4,
       filename = paste0(output,"/BARPLOTS/DIAGNOSIS_PROP_BAR.pdf"))

dittoBarPlot(o, group.by = "donor_organism.sex", var.labels.reorder = c(6,7,8,9,10,4,5,2,3,1,11), color.panel = colors,
             var = "ord") %>% ggpar(xlab = "", title = "", ylab = "") -> e
ggsave(plot = e + NoLegend(), dpi = "retina", units = "in", width = 1.7, height = 4,
       filename = paste0(output,"/BARPLOTS/SEX_PROP_BAR.pdf"))

ggarrange(plotlist = list(f + NoLegend(),e + NoLegend()), ncol = 2, nrow = 1, align = "hv") -> p
ggsave(plot = p, dpi = "retina", units = "in", width = 3, height = 4,
       filename = paste0(output,"/BARPLOTS/SEX_DIAGNOSIS_PROP_BAR.pdf"))

paste0(o@meta.data$donor_organism.sex, "_", o@meta.data$diagnosis) -> o@meta.data$diagnosis_sex
dittoBarPlot(o, group.by = "donor_organism.sex", var.labels.reorder = c(1,2,4,5,6,7,8,9,10,11,3),
             var = "integrated_snn_res.0.5") %>% ggpar(xlab = "") -> e

o@meta.data$diagnosis_sex <- factor(o@meta.data$diagnosis_sex, levels = c("male_AD","female_AD","male_Non-disease control","female_Non-disease control"))
dittoBarPlot(o, group.by = "diagnosis_sex", var.labels.reorder = c(1,2,4,5,6,7,8,9,10,11,3),
             var = "integrated_snn_res.0.5") %>% ggpar(xlab = "") -> e



### GENE PANELS

plot_density(o, "CX3CR1", reduction = "umap", pal = "inferno", adjust = 1) + NoAxes() -> p1
plot_density(o, "CD163", reduction = "umap", pal = "inferno", adjust = 1) + NoAxes() -> p2
plot_density(o, "CD83", reduction = "umap", pal = "inferno", adjust = 1) + NoAxes() -> p3
plot_density(o, "TGFBI", reduction = "umap", pal = "inferno", adjust = 1) + NoAxes() -> p9
plot_density(o, "P2RY12", reduction = "umap", pal = "inferno", adjust = 1) + NoAxes() -> p4
plot_density(o, "CTNNA2", reduction = "umap", pal = "inferno", adjust = 1) + NoAxes() -> p5
plot_density(o, "CCL3", reduction = "umap", pal = "inferno", adjust = 1) + NoAxes() -> p6
plot_density(o, "VCAN", reduction = "umap", pal = "inferno", adjust = 1) + NoAxes() -> p7
plot_density(o, "BRIP1", reduction = "umap", pal = "inferno", adjust = 1) + NoAxes() -> p8
ggarrange(plotlist = list(p1,p4,p2,p3,p9,p6,p7,p8,p5), ncol = 3, nrow = 3) -> p

ggsave(plot = p, dpi = "retina", units = "in", width = 8, height = 7,
       filename = paste0(output,"/UMAP/MARKERS_PANEL.pdf"))

DotPlot(o, group.by = "ord", 
        features = c("C3","CD74","CX3CR1","P2RY12",
                     "NAV3","CTNNA2",
                     "CD163","TGFBI",
                     "CD83","LPL","APOE","TREM2",
                     "CCL3","IL1B",
                     "BRIP1","CENPP","BRCA2",
                     "F13A1","LYVE1",
                     "VCAN")) %>% ggpar(x.text.angle = 45, rotate = T)

dittoDotPlot(o, vars = c("C3","CD74","CX3CR1","P2RY12",
                         "NAV3","CTNNA2",
                         "CD163","TGFBI","DPYD",
                         "CD83","LPL","APOE","TREM2","MYO1E",
                         "CCL3","CCL4","IL1B",
                         "BRIP1","CENPP","BRCA2",
                         "F13A1","LYVE1",
                         "VCAN"), group.by = "ord", 
             min.color = "white",
             max.color = "red") %>% ggpar(x.text.angle = 45, rotate = T) -> p
ggsave(plot = p, dpi = "retina", units = "in", width = 5, height =5,
       filename = paste0(output,"/UMAP/DOTPLOT_PANEL.pdf"))


a <- gene_sets$AGE_DOWN_MICROGLIA
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "AGE_DOWN_MICROGLIA", assay = "RNA") -> o

a <- gene_sets$HAM_DOWN_Srinivasan_2021
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "HAM_DOWN_Srinivasan_2021", assay = "RNA") -> o

a <- gene_sets$AD_traj_early
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "AD_traj_early", assay = "RNA") -> o

a <- gene_sets$KEGG_CELL_CYCLE
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "KEGG_cell_cycle", assay = "RNA") -> o

o@meta.data %>% dplyr::select(c("HOMEOSTATIC_KRASEMANN1",
                                "SalaFrigerio_HM1",
                                "ARM1",
                                "HAM_UP_Srinivasan_20211",
                                "HAM_DOWN_Srinivasan_20211",
                                "DAM_UP_KRASEMAN1",
                                "AGE_UP_MICROGLIA1",
                                "AGE_DOWN_MICROGLIA1",
                                "Cytokine_Response1",
                                "AD_TRAJ_LATE1",
                                "AD_traj_early1",
                                "CAM1",
                                "diss1","ord"
)) -> df



dittoDotPlot(o, vars = c("HOMEOSTATIC_KRASEMANN1",
                         "SalaFrigerio_HM1",
                         "AD_traj_early1",
                         "AGE_DOWN_MICROGLIA1",
                         "HAM_DOWN_Srinivasan_20211",
                         "ARM1",
                         "HAM_UP_Srinivasan_20211",
                         
                         "DAM_UP_KRASEMAN1",
                         "AGE_UP_MICROGLIA1",
                         
                         "Cytokine_Response1",
                         "AD_TRAJ_LATE1",
                         
                         "CAM1",
                         "diss1", "KEGG_cell_cycle1"
), group.by = "ord",  
min.color = "white", size = 3,
max.color = "red") %>% ggpar(x.text.angle = 45, rotate = T) + 
  scale_color_gradient2(limits = c(-1,2), oob = scales::squish, low = "blue", high = "red", mid = "grey98")-> p
ggsave(plot = p + NoLegend(), dpi = "retina", units = "in", width = 5, height =4,
       filename = paste0(output,"/UMAP/DOTPLOT_MODULE_PANEL.pdf"))
ggsave(plot = p, dpi = "retina", units = "in", width = 5, height =5,
       filename = paste0(output,"/UMAP/DOTPLOT_MODULE_PANEL_LEGEND.pdf"))

reshape2::dcast(p$data, formula = var ~ grouping, value.var = "color") -> df
rownames(df) <- df$var
df[,-1] -> df
pheatmap::pheatmap(df, cluster_cols = T, 
                   breaks = seq(-1, 2, length.out = 11), color = RColorBrewer::brewer.pal(n = 11, name = "Reds")) 



y.labels = c("Stressed | Thrupp",
             "CAM | Prinz",
             "MG AD TRAJ Late | Gerrits",
             "Cytokine response MG | ",
             "Age associated UP | Galatro",
             "DAM UP | Krasemann",
             "HAM UP | Srinivasan",
             "Activated MG | SalaFrigerio",
             "Homeostatic MG | SalaFrigerio",
             "Homeostatic MG | Krasemann")

ggsave(plot = p, dpi = "retina", units = "in", width = 5, height =5,
       filename = paste0(output,"/UMAP/DOTPLOT_PANEL.pdf"))







### Astrocyte LOY in Smith AM

library(celldex)
library(monocle3)
library(monocle)
library(slingshot)
library(scDblFinder)
library(SingleR)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)
library(ggpubr)
library(nichenetr)
library(dittoSeq)
library(escape)
library(Nebulosa)
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
edb <- EnsDb.Hsapiens.v86

setwd("H:/LOY/scLOY/scLOY/")
source("E:/LOY/Example/PACKAGE/call_LOY.R")
source("E:/LOY/scLOY/scripts/loy_functions.R")
source('H:/LOY/scLOY/scLOY/R/load.R')
source('H:/LOY/scLOY/scLOY/data_raw/build_data_files.R')
###

vega <- c('#1f77b4','#aec7e8','#ff7f0e',
          '#ffbb78','#2ca02c','#98df8a',
          '#d62728','#ff9896','#9467bd',
          '#c5b0d5','#8c564b','#c49c94',
          '#e377c2','#f7b6d2','#7f7f7f',
          '#c7c7c7','#bcbd22','#dbdb8d',
          '#17becf','#9edae5')

cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

diss <- c("RGS1","DUSP1","AC022217.3","DUSP1","HSPA1A","HSPA1B","DDIT4","SLC2A3","FOS","JUN",
          "AC131944.1","HIF1A-AS3","HIST1H2AC","HIST1H2BD","HIST1H2BG","HIST2H2BE","RGS2","RHOB",
          "HIST1H2BJ","TEX14","UBC","HSP90AA1","AL118516.1","JUNB","HIST2H2AA4","HIST1H1C","BTG1","BTG2","DUSP5")

read.csv("h:/LOY/scLOY/data/microglia_gene_sets.csv") -> gene_sets

readRDS(file = "e:/LOY/scLOY/processed_seurat/BRAIN/GSE160936_AD_GLIAL_introns_included.RDS") -> tmp
#DimPlot(o, group.by = "RNA_snn_res.1.1", reduction = "umap")

### EC #############
subset(tmp, 
       nCount_RNA >= 1000 & 
         nFeature_RNA >= 800 & 
         specific_tissue == "EC" & donor_organism.sex == "male" & cell_ident_MV == "Astro") -> o

### run seurat integration
ifnb.list <- SplitObject(o, split.by = "title")

ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  return(x)
})

rownames(o)[grep(x = rownames(o), pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")] -> ribo
ribo[grep(ribo,pattern = "K",invert = T)] -> ribo
rownames(o)[grep(x = rownames(o), pattern = "^MT-")] -> mito
rownames(o)[(rownames(o) %in% Y_scLOY_genes$gene_name)] -> Y
Y_scLOY_genes$gene_name[Y_scLOY_genes$gene_name %in% rownames(o)] -> male_genes

features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 2000)
features <- features[!(features %in% c("XIST",ribo,mito,"JPX"))]
anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
combined <- IntegrateData(anchorset = anchors)
DefaultAssay(combined) <- "integrated"
combined -> o
rm(combined)

rownames(o)[grep(x = rownames(o), pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")] -> ribo
ribo[grep(ribo,pattern = "K",invert = T)] -> ribo

o <- Seurat::ScaleData(o, vars.to.regress = c("nCount_RNA","percent.mt","percent.rb"))
v <- Seurat::VariableFeatures(o)
v <- v[!(v %in% c(ribo,mito,"XIST","JPX"))]
o <- Seurat::RunPCA(o, npcs = 50, features = v, assay = "integrated")
ElbowPlot(o,ndims = 50 )

o <- Seurat::RunUMAP(o,reduction = "pca", dims = 1:12)
o <- Seurat::FindNeighbors(o,reduction = "pca", dims = 1:12)
o <- Seurat::FindClusters(o,resolution = seq(0.1,2,0.1))
DimPlot(o, group.by = "integrated_snn_res.0.5", label = T, cols = dittoColors())
DimPlot(o, group.by = "integrated_snn_res.0.3", label = T, cols = dittoColors())
DimPlot(o, group.by = "title")

FindMarkers(o, ident.1 = 8, ident.2 = NULL, group.by = "integrated_snn_res.0.5",
            only.pos = T, logfc.threshold = 0.75)


Idents(o) <- "integrated_snn_res.0.5"

dir.create("e:/LOY/scLOY/processed_integrated/ASTROCYTE/")
saveRDS(o, "e:/LOY/scLOY/processed_integrated/ASTROCYTE/IMPERIAL_CCA_EC_MALE_JULY12.RDS")

### HOMEOSTATIC GENES ##########
output <- "EC_ASTRO"
dir.create(paste0(output,"/SUBCLUSTERS"), recursive = T)
dir.create(paste0(output,"/UMAP"), recursive = T)
dir.create(paste0(output,"/BARPLOTS"), recursive = T)

DefaultAssay(o) <- "RNA"
p <- plot_density(o, c("GFAP","DPP10","SLC1A2"), reduction = "umap",
                  method = "wkde", adjust = 3,
                  pal = "inferno", joint = T,
                  combine = F)[[4]];
ggsave(plot = p, dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/HOMEO_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(),
       dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/HOMEO.pdf"))

p <- plot_density(o, c("GFAP","DPP10","SLC1A2"), reduction = "umap",
                  method = "wkde", adjust = 2,
                  pal = "inferno", joint = T,
                  combine = F)[[1]];
ggsave(plot = p, dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/GFAP_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(), dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/GFAP.pdf"))

p <- plot_density(o, c("GFAP","DPP10","SLC1A2"), reduction = "umap",
                  method = "wkde", adjust = 2,
                  pal = "inferno", joint = T,
                  combine = F)[[2]];
ggsave(plot = p, dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/DPP10_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(),
       dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/DPP10.pdf"))

p <- plot_density(o, c("GFAP","DPP10","SLC1A2"), reduction = "umap",
                  method = "wkde", adjust = 2,
                  pal = "inferno", joint = T,
                  combine = F)[[3]];
ggsave(plot = p,  dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/SLC1A2_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(),
       dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/SLC1A2.pdf"))


### GENE SETS #############
### find homeostatic and AD associated clusters
### pull these clusters

a <- gene_sets$SalaFrigerio_ARM
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "ARM", assay = "RNA") -> o
plot_density(o, "ARM1", reduction = "umap", pal = "inferno", adjust = 1) -> p
FeaturePlot(o, features = "ARM1", order = T, pt.size = 0.9) + scale_color_viridis_c(option = "B")
ggsave(plot = p + NoAxes(), dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/ARM_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(), dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/ARM.pdf"))

a <- gene_sets$Cytokine_Response
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "Cytokine_Response", assay = "RNA") -> o
plot_density(o, "Cytokine_Response1", reduction = "umap", pal = "inferno") -> p
ggsave(plot = p + NoAxes(), dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/Cytokine_Response_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(), dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/Cytokine_Response.pdf"))

plot_density(o, features = c("CCL3","CCL4","IL1B"), method = "wkde", adjust = 2,
             reduction = "umap",
             pal = "inferno", joint = T, combine = F) -> p
p[[4]] + NoAxes() -> p
ggsave(plot = p + NoAxes(), dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/CCL3_CCL4_IL1B_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(), dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/CCL3_CCL4_IL1B.pdf"))


a <- gene_sets$Prinz_CAM
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "CAM", assay = "RNA") -> o
plot_density(o, "CAM1", reduction = "umap", pal = "inferno") -> p
ggsave(plot = p, dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/CAM_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(),
       dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/CAM.pdf"))

a <- gene_sets$DAM_UP_SOBUE
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "DAM_UP_SOBUE", assay = "RNA") -> o
plot_density(o, "DAM_UP_SOBUE1", reduction = "umap", pal = "inferno") -> p
ggsave(plot = p + NoAxes(), dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/DAM_UP_SOBUE_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(),
       dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/DAM_UP_SOBUE.pdf"))

a <- gene_sets$DAM_UP_KRASEMANN_mouse %>% nichenetr::convert_mouse_to_human_symbols()
a[!is.na(a)] -> a
names(a) <- NULL
a %>% as.vector() -> a
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "DAM_UP_KRASEMAN", assay = "RNA") -> o
plot_density(o, "DAM_UP_KRASEMAN1", reduction = "umap", pal = "inferno") -> p
ggsave(plot = p+ NoAxes(), dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/DAM_UP_KRASEMAN_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(),
       dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/DAM_UP_KRASEMAN.pdf"))

FeaturePlot(o, features = "AGE_UP_MICROGLIA1", pt.size = 0.8, order = T) + scale_color_viridis_c(option = "B") -> p
ggsave(plot = p+ NoAxes(), dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/DAM_UP_KRASEMAN_LEGEND_FEATUREPLOT.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(), dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/DAM_UP_KRASEMAN_LEGEND_FEATUREPLOT.pdf"))


a <- gene_sets$AGE_UP_MICROGLIA
a[!is.na(a)] -> a
names(a) <- NULL
a %>% as.vector() -> a
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "AGE_UP_MICROGLIA", assay = "RNA") -> o
plot_density(o, "AGE_UP_MICROGLIA1", reduction = "umap", pal = "inferno") -> p
ggsave(plot = p+ NoAxes(), dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/AGE_UP_MICROGLIA_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(), dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/AGE_UP_MICROGLIA.pdf"))

FeaturePlot(o, features = "AGE_UP_MICROGLIA1", pt.size = 0.8, order = T) + scale_color_viridis_c(option = "B") -> p
ggsave(plot = p+ NoAxes(), dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/AGE_UP_MICROGLIA_LEGEND_FEATUREPLOT.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(), dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/AGE_UP_MICROGLIA_FEATUREPLOT.pdf"))

a <- gene_sets$HOMEOSTATIC_KRASEMANN_mouse %>% nichenetr::convert_mouse_to_human_symbols()
a[!is.na(a)] -> a
names(a) <- NULL
a %>% as.vector() -> a
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "HOMEOSTATIC_KRASEMANN", assay = "RNA") -> o
plot_density(o, "HOMEOSTATIC_KRASEMANN1", reduction = "umap", pal = "inferno") -> p
ggsave(dpi = "retina", units = "in", width = 4.5, height = 3.5,
       plot = p+ NoAxes(),
       filename = paste0(output,"/SUBCLUSTERS/HOMEOSTATIC_KRASEMANN_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(),
       dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/HOMEOSTATIC_KRASEMANN.pdf"))

a <- gene_sets$AD_traj_late
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "AD_TRAJ_LATE", assay = "RNA") -> o
plot_density(o, "AD_TRAJ_LATE1", reduction = "umap", pal = "inferno") -> p
ggsave(dpi = "retina", units = "in", width = 4.5, height = 3.5,
       plot = p+ NoAxes(),
       filename = paste0(output,"/SUBCLUSTERS/AD_TRAJ_LATE_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(),
       dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/AD_TRAJ_LATE.pdf"))

a <- gene_sets$AD_traj_conversion
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "AD_TRAJ_CONVERSION", assay = "RNA") -> o
plot_density(o, "AD_TRAJ_CONVERSION1", reduction = "umap", pal = "inferno") -> p
ggsave(dpi = "retina", units = "in", width = 4.5, height = 3.5,
       plot = p+ NoAxes(),
       filename = paste0(output,"/SUBCLUSTERS/AD_TRAJ_CONVERSION_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(),
       dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/AD_TRAJ_CONVERSION.pdf"))

a <- gene_sets$AD_traj_early_end
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "AD_TRAJ_EARLY_END", assay = "RNA") -> o
plot_density(o, "AD_TRAJ_EARLY_END1", reduction = "umap", pal = "inferno") -> p
ggsave(dpi = "retina", units = "in", width = 4.5, height = 3.5,
       plot = p+ NoAxes(),
       filename = paste0(output,"/SUBCLUSTERS/AD_TRAJ_EARLY_END_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(),
       dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/AD_TRAJ_EARLY_END.pdf"))

a <- gene_sets$AD_traj_early
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "AD_TRAJ_EARLY", assay = "RNA") -> o
plot_density(o, "AD_TRAJ_EARLY1", reduction = "umap", pal = "inferno") -> p
ggsave(dpi = "retina", units = "in", width = 4.5, height = 3.5,
       plot = p+ NoAxes(),
       filename = paste0(output,"/SUBCLUSTERS/AD_TRAJ_EARLY_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(),
       dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/AD_TRAJ_EARLY.pdf"))

a <- gene_sets$SalaFrigerio_HM
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "SalaFrigerio_HM", assay = "RNA") -> o
plot_density(o, "HM1", reduction = "umap", pal = "inferno") -> p
ggsave(plot = p+ NoAxes(), dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/HM_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(), dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/HM.pdf"))

plot_density(o, features = "CDKN2A", reduction = "umap", pal = "inferno", method = "wkde", adjust = 2) -> p
ggsave(plot = p+ NoAxes(), dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/SENES_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(), dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/SENES.pdf"))

a <- gene_sets$HAM_UP_Srinivasan_2021
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "HAM_UP_Srinivasan_2021", assay = "RNA") -> o
plot_density(subset(o, diagnosis == "AD"), "HAM_UP_Srinivasan_20211", reduction = "umap", pal = "inferno") -> p
ggsave(plot = p+ NoAxes(), dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/HM_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(), dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/HM.pdf"))

a <- gene_sets$AGE_DOWN_MICROGLIA
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "AGE_DOWN_MICROGLIA", assay = "RNA") -> o

a <- gene_sets$HAM_DOWN_Srinivasan_2021
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "HAM_DOWN_Srinivasan_2021", assay = "RNA") -> o

a <- gene_sets$AD_traj_early
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "AD_traj_early", assay = "RNA") -> o

a <- gene_sets$KEGG_CELL_CYCLE
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "KEGG_cell_cycle", assay = "RNA") -> o

dittoDotPlot(o, vars = c("HOMEOSTATIC_KRASEMANN1",
                         "SalaFrigerio_HM1",
                         "AD_traj_early1",
                         "AGE_DOWN_MICROGLIA1",
                         "HAM_DOWN_Srinivasan_20211",
                         "ARM1",
                         "HAM_UP_Srinivasan_20211",
                         
                         "DAM_UP_KRASEMAN1",
                         "AGE_UP_MICROGLIA1",
                         
                         "Cytokine_Response1",
                         "AD_TRAJ_LATE1",
                         
                         "CAM1",
                         "diss1", "KEGG_cell_cycle1"
), group.by = "integrated_snn_res.0.5",  
min.color = "white", size = 4,
max.color = "red") %>% ggpar(x.text.angle = 45, rotate = T) + 
  scale_color_gradient2(limits = c(-1,2), oob = scales::squish, low = "blue", high = "red", mid = "grey98")-> p



###### ANNOTATE #########

subset(o, integrated_snn_res.0.4 == "4", invert = T) -> o
dittoDimPlot(o, var = "integrated_snn_res.0.4", 
             do.label = T, labels.highlight = T, labels.size = 4) -> p
FindMarkers(o, ident.1 = 4, ident.2 = NULL, group.by = "integrated_snn_res.0.4", only.pos = T, logfc.threshold = 0.5)

o@meta.data$MG_subtype2 <- o@meta.data$integrated_snn_res.0.4
levels(o@meta.data$MG_subtype2)[levels(o@meta.data$MG_subtype2)%in%c(2,3)] <- "GFAPhi"
levels(o@meta.data$MG_subtype2)[levels(o@meta.data$MG_subtype2)%in%c(1,0)] <- "GFAPlow"
DimPlot(o, group.by = "MG_subtype2", label = T)

o@meta.data$MG_subtype_C2 <- o@meta.data$integrated_snn_res.0.4
levels(o@meta.data$MG_subtype_C2)[levels(o@meta.data$MG_subtype_C2)%in%c(2)] <- "GFAPhi1"
levels(o@meta.data$MG_subtype_C2)[levels(o@meta.data$MG_subtype_C2)%in%c(3)] <- "GFAPhi2"
levels(o@meta.data$MG_subtype_C2)[levels(o@meta.data$MG_subtype_C2)%in%c(1)] <- "GFAPlow1"
levels(o@meta.data$MG_subtype_C2)[levels(o@meta.data$MG_subtype_C2)%in%c(0)] <- "GFAPlow2"
DimPlot(o, group.by = "MG_subtype_C2", label = T)

Idents(o) <- "MG_subtype_C2"
FindAllMarkers(o, logfc.threshold = 0.2, min.pct = 0.1, test.use = "MAST",
               latent.vars = c("nCount_RNA","nFeature_RNA","percent.rb")) -> DE
data.table::fwrite(x = DE, file = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/FINAL/IMPERIAL/EC_ASTRO/DE_final.txt",
                   sep = "\t")

factor(o@meta.data$MG_subtype_C2, levels = c("GFAPhi1","GFAPhi2","GFAPlow1","GFAPlow2")) -> o@meta.data$ord2

colors <- c(dittoColors()[21],dittoColors()[5],
            dittoColors()[1], dittoColors()[17])
           

dittoDimPlot(o, var = "ord2", color.panel = colors,  
             do.label = T, labels.highlight = T, labels.size = 4) -> p
p %>% ggpar(title = "") -> p
ggsave(plot = p + NoAxes() + NoLegend(), dpi = "retina", units = "in", width = 6, height = 5,
       filename = paste0(output,"/UMAP/NO_AXES_NO_LEG_UMAP.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(), dpi = "retina", units = "in", width = 3, height = 2.5,
       filename = paste0(output,"/UMAP/NO_AXES_NO_LEG_UMAP_SMALL.pdf"))

dittoDimPlot(o, var = "ord2", color.panel = colors,  
             do.label = T, labels.highlight = T, labels.size = 4) -> p
p %>% ggpar(title = "") -> p
ggsave(plot = p + NoAxes() + NoLegend(), dpi = "retina", units = "in", width = 6, height = 5,
       filename = paste0(output,"/UMAP/NO_AXES_NO_LEG_UMAP.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(), dpi = "retina", units = "in", width = 3, height = 2.5,
       filename = paste0(output,"/UMAP/NO_AXES_NO_LEG_UMAP_SMALL.pdf"))

dittoDimPlot(o, var = "ord2", color.panel = colors,  
             do.label = F, labels.highlight = T, labels.size = 4) -> p
p %>% ggpar(title = "") -> p
ggsave(plot = p + NoAxes(), dpi = "retina", units = "in", width = 6, height = 5,
       filename = paste0(output,"/UMAP/NO_AXES_UMAP.pdf"))
ggsave(plot = p + NoAxes(), dpi = "retina", units = "in", width = 3, height = 2.5,
       filename = paste0(output,"/UMAP/NO_AXES_UMAP_SMALL.pdf"))

dittoDimPlot(o, var = "ord2", color.panel = colors,split.by = "diagnosis", size = 0.8,
             do.label = F, labels.highlight = T, labels.size = 4, do.contour = T) -> p
p %>% ggpar(title = "") + NoAxes() + NoLegend() -> p
ggsave(plot = p + NoAxes(), dpi = "retina", units = "in", width = 5, height = 3,
       filename = paste0(output,"/UMAP/AD_SPLIT.pdf"))

#dittoDimPlot(o, var = "ord", color.panel = colors,split.by = "specific_tissue", size = 0.8,
#             do.label = F, labels.highlight = T, labels.size = 4, do.contour = T) -> p
#p %>% ggpar(title = "") + NoAxes() + NoLegend() -> p
#ggsave(plot = p + NoAxes(), dpi = "retina", units = "in", width = 5, height = 3,
#       filename = paste0(output,"/UMAP/TISSUE_SPLIT.pdf"))

dittoDimPlot(o, var = "title", size = 0.8,
             do.label = F) -> p
p %>% ggpar(title = "") + NoAxes() -> p
ggsave(plot = p + NoAxes(), dpi = "retina", units = "in", width = 6, height = 4,
       filename = paste0(output,"/UMAP/SAMPLE_SPLIT.pdf"))

dittoDimPlot(subset(o, nCount_RNA > 2500), var = "LOY", size = 1, split.by = "diagnosis", color.panel = c("red","grey"),
             do.label = F, order = "decreasing") -> p
p %>% ggpar(title = "") + NoAxes() -> p
ggsave(plot = p + NoAxes(), dpi = "retina", units = "in", width = 7.5, height = 3.5,
       filename = paste0(output,"/UMAP/DIAGNOSIS_LOY_2500.pdf"))


dittoDimPlot(o, var = "LOY", size = 1, split.by = "title", color.panel = c("red","grey"),
             do.label = F, order = "decreasing") -> p
p %>% ggpar(title = "") + NoAxes() -> p
ggsave(plot = p + NoAxes(), dpi = "retina", units = "in", width = 8, height = 8,
       filename = paste0(output,"/UMAP/TITLE_LOY.pdf"))

dittoDimPlot(o, var = "ord2", size = 1, split.by = "LOY", color.panel = colors,
             do.label = F, do.contour = T) -> p
p %>% ggpar(title = "") + NoAxes() -> p
ggsave(plot = p + NoAxes(), dpi = "retina", units = "in", width = 7.5, height = 3.5,
       filename = paste0(output,"/UMAP/LOY_CLUSTERS_CONTOUR.pdf"))




################# BAR PLOT  # 

DimPlot(o, group.by = "ord2", label = T, cols = dittoColors())
DimPlot(o, split.by = "diagnosis", label = T, group.by = "integrated_snn_res.0.5", cols = dittoColors())

order <- NULL

dittoBarPlot(o, group.by = "title", color.panel = colors, var.labels.reorder = NULL,
             var = "ord2") %>% ggpar(xlab = "", title = "") -> a
ggsave(plot = a + NoLegend(), dpi = "retina", units = "in", width = 5, height = 3,
       filename = paste0(output,"/BARPLOTS/SAMPLE_PROP_BAR.pdf"))

o@meta.data[o@meta.data$diagnosis=="AD" & o@meta.data$nCount_RNA > 2500,] %>% rownames() -> cells
dittoBarPlot(o, group.by = "title", var = "LOY", var.labels.reorder = c(2,1),
             cells.use = cells,
             color.panel = c("grey","red")) %>% ggpar(xlab = "", legend = "none", title = "") -> b 
o@meta.data[o@meta.data$diagnosis=="Non-disease control" & o@meta.data$nCount_RNA > 2500,] %>% rownames() -> cells
dittoBarPlot(o, group.by = "title", var = "LOY", var.labels.reorder = c(2,1),
             cells.use = cells,
             color.panel = c("grey","red")) %>% ggpar(xlab = "", legend = "none", title = "") -> c 
ggarrange(plotlist = list(b ,c + rremove("ylab") + rremove("y.text") ), ncol = 2, nrow = 1, align = "hv") -> a 
ggsave(plot = a + NoLegend(), dpi = "retina", units = "in", width = 5, height = 3,
       filename = paste0(output,"/BARPLOTS/SAMPLE_PROP_BAR_BY_DIAGNOSIS_2500.pdf"))

#ggarrange(plotlist = list(a,b), ncol = 1, nrow = 2, align = "hv")
ggsave(plot = b + NoLegend(), dpi = "retina", units = "in", width = 3, height = 3,
       filename = paste0(output,"/BARPLOTS/SAMPLE_PROP_LOY_BAR.pdf"))

dittoBarPlot(o, group.by = "ord2", var = "LOY", var.labels.reorder = c(2,1), x.reorder = order,
             color.panel = c("grey","red")) %>% ggpar(xlab = "") -> c
ggsave(plot = c + NoLegend(), dpi = "retina", units = "in", width = 3, height = 3,
       filename = paste0(output,"/BARPLOTS/CLUSTER_PROP_LOY_BAR.pdf"))

c$data %>% dplyr::filter(label == "LOY") -> dat
factor(dat$grouping, levels = dat[order(-dat$percent),]$grouping) -> dat$grouping
k <- colors
dat %>% ggbarplot(x = "grouping", y = "percent", fill = "grouping", palette = k) %>% 
  ggpar(legend = "none", xlab = "", ylab = "Percent LOY", x.text.angle = 45, ylim = c(0,0.25)) -> p
ggsave(plot = p + NoLegend(), dpi = "retina", units = "in", width = 4, height = 3,
       filename = paste0(output,"/BARPLOTS/CLUSTER_PROP_LOY_BAR_ranked.pdf"))

dittoBarPlot(o, group.by = "diagnosis", var = "LOY", var.labels.reorder = c(2,1),
             color.panel = c("grey","red")) %>% ggpar(xlab = "") -> d
ggsave(plot = d + NoLegend(), dpi = "retina", units = "in", width = 1, height = 3,
       filename = paste0(output,"/BARPLOTS/DIAGNOSIS_PROP_LOY_BAR.pdf"))

dittoBarPlot(subset(o, nCount_RNA > 2500 & nFeature_RNA > 1000), group.by = "title", var = "LOY", var.labels.reorder = c(2,1),
             color.panel = c("grey","red")) %>% ggpar(xlab = "") -> b 
#ggarrange(plotlist = list(a,b), ncol = 1, nrow = 2, align = "hv")
ggsave(plot = b + NoLegend(), dpi = "retina", units = "in", width = 3, height = 3,
       filename = paste0(output,"/BARPLOTS/SAMPLE_PROP_LOY_BAR_2500.pdf"))

dittoBarPlot(subset(o, nCount_RNA > 2500 & nFeature_RNA > 1000), group.by = "ord", var = "LOY", var.labels.reorder = c(2,1), x.reorder = order,
             color.panel = c("grey","red")) %>% ggpar(xlab = "") -> c
ggsave(plot = c + NoLegend(), dpi = "retina", units = "in", width = 3, height = 3,
       filename = paste0(output,"/BARPLOTS/CLUSTER_PROP_LOY_BAR_2500.pdf"))

dittoBarPlot(subset(o, nCount_RNA > 2500 & nFeature_RNA > 1000), group.by = "diagnosis", var = "LOY", var.labels.reorder = c(2,1),
             color.panel = c("grey","red")) %>% ggpar(xlab = "") -> d
ggsave(plot = d + NoLegend(), dpi = "retina", units = "in", width = 1, height = 3,
       filename = paste0(output,"/BARPLOTS/DIAGNOSIS_PROP_LOY_BAR_2500.pdf"))


dittoBarPlot(o, group.by = "diagnosis", var.labels.reorder = order, color.panel = colors,
             var = "ord2") %>% ggpar(xlab = "", title = "") -> f
ggsave(plot = e + NoLegend(), dpi = "retina", units = "in", width = 1.5, height = 4,
       filename = paste0(output,"/BARPLOTS/DIAGNOSIS_PROP_BAR.pdf"))

dittoBarPlot(o, group.by = "donor_organism.sex", var.labels.reorder = order, color.panel = colors,
             var = "ord2") %>% ggpar(xlab = "", title = "", ylab = "") -> e
ggsave(plot = e + NoLegend(), dpi = "retina", units = "in", width = 1.7, height = 4,
       filename = paste0(output,"/BARPLOTS/SEX_PROP_BAR.pdf"))

ggarrange(plotlist = list(f + NoLegend(),e + NoLegend()), ncol = 2, nrow = 1, align = "hv") -> p
ggsave(plot = p, dpi = "retina", units = "in", width = 3, height = 4,
       filename = paste0(output,"/BARPLOTS/SEX_DIAGNOSIS_PROP_BAR.pdf"))

paste0(o@meta.data$donor_organism.sex, "_", o@meta.data$diagnosis) -> o@meta.data$diagnosis_sex
dittoBarPlot(o, group.by = "donor_organism.sex", var.labels.reorder = order,
             var = "integrated_snn_res.0.5") %>% ggpar(xlab = "") -> e

o@meta.data$diagnosis_sex <- factor(o@meta.data$diagnosis_sex, levels = c("male_AD","female_AD","male_Non-disease control","female_Non-disease control"))
dittoBarPlot(o, group.by = "diagnosis_sex", var.labels.reorder = order,
             var = "integrated_snn_res.0.5") %>% ggpar(xlab = "") -> e



### GENE PANELS

plot_density(o, "CX3CR1", reduction = "umap", pal = "inferno", adjust = 1) + NoAxes() -> p1
plot_density(o, "CD163", reduction = "umap", pal = "inferno", adjust = 1) + NoAxes() -> p2
plot_density(o, "CD83", reduction = "umap", pal = "inferno", adjust = 1) + NoAxes() -> p3
plot_density(o, "TGFBI", reduction = "umap", pal = "inferno", adjust = 1) + NoAxes() -> p9
plot_density(o, "P2RY12", reduction = "umap", pal = "inferno", adjust = 1) + NoAxes() -> p4
plot_density(o, "CTNNA2", reduction = "umap", pal = "inferno", adjust = 1) + NoAxes() -> p5
plot_density(o, "CCL3", reduction = "umap", pal = "inferno", adjust = 1) + NoAxes() -> p6
plot_density(o, "VCAN", reduction = "umap", pal = "inferno", adjust = 1) + NoAxes() -> p7
plot_density(o, "BRIP1", reduction = "umap", pal = "inferno", adjust = 1) + NoAxes() -> p8
ggarrange(plotlist = list(p1,p4,p2,p3,p9,p6,p7,p8,p5), ncol = 3, nrow = 3) -> p

ggsave(plot = p, dpi = "retina", units = "in", width = 8, height = 7,
       filename = paste0(output,"/UMAP/MARKERS_PANEL.pdf"))

DotPlot(o, group.by = "ord", 
        features = c("C3","CD74","CX3CR1","P2RY12",
                     "NAV3","CTNNA2",
                     "CD163","TGFBI","DPYD","F13A1",
                     "CD83","LPL","APOE","TREM2","CCL3","IL1B",
                     "GPNMB",
                     "BRIP1","CENPP","BRCA2",
                     "F13A1","LYVE1",
                     "VCAN")) %>% ggpar(x.text.angle = 45, rotate = T)

dittoDotPlot(o, vars = c("C3","CD74","CX3CR1","P2RY12",
                         "GPM6A","CTNNA2",
                         "CD163","TGFBI","DPYD",
                         "CD83","LPL","APOE","TREM2","MYO1E","CCL3","CCL4","IL1B",
                         "BRIP1","CENPP","BRCA2",
                         "F13A1","LYVE1",
                         "VCAN"), group.by = "ord2", 
             min.color = "white",
             max.color = "red") %>% ggpar(x.text.angle = 45, rotate = T) -> p
ggsave(plot = p, dpi = "retina", units = "in", width = 5, height =5,
       filename = paste0(output,"/UMAP/DOTPLOT_PANEL.pdf"))




o@meta.data %>% dplyr::select(c("HOMEOSTATIC_KRASEMANN1",
                                "SalaFrigerio_HM1",
                                "ARM1",
                                "HAM_UP_Srinivasan_20211",
                                "HAM_DOWN_Srinivasan_20211",
                                "DAM_UP_KRASEMAN1",
                                "AGE_UP_MICROGLIA1",
                                "AGE_DOWN_MICROGLIA1",
                                "Cytokine_Response1",
                                "AD_TRAJ_LATE1",
                                "AD_traj_early1",
                                "CAM1",
                                "diss1","ord"
)) -> df



dittoDotPlot(o, vars = c("HOMEOSTATIC_KRASEMANN1",
                         "SalaFrigerio_HM1",
                         "AD_traj_early1",
                         "AGE_DOWN_MICROGLIA1",
                         "HAM_DOWN_Srinivasan_20211",
                         "ARM1",
                         "HAM_UP_Srinivasan_20211",
                         
                         "DAM_UP_KRASEMAN1",
                         "AGE_UP_MICROGLIA1",
                         
                         "Cytokine_Response1",
                         "AD_TRAJ_LATE1",
                         
                         "CAM1",
                         "diss1", "KEGG_cell_cycle1"
), group.by = "ord2",  
min.color = "white", size = 3,
max.color = "red") %>% ggpar(x.text.angle = 45, rotate = T) + 
  scale_color_gradient2(limits = c(-1,2), oob = scales::squish, low = "blue", high = "red", mid = "grey98")-> p
ggsave(plot = p + NoLegend(), dpi = "retina", units = "in", width = 5, height =4,
       filename = paste0(output,"/UMAP/DOTPLOT_MODULE_PANEL.pdf"))
ggsave(plot = p, dpi = "retina", units = "in", width = 5, height =5,
       filename = paste0(output,"/UMAP/DOTPLOT_MODULE_PANEL_LEGEND.pdf"))

dittoDotPlot(o, vars = c("HOMEOSTATIC_KRASEMANN1",
                         "SalaFrigerio_HM1",
                         "AD_traj_early1",
                         "AGE_DOWN_MICROGLIA1",
                         "HAM_DOWN_Srinivasan_20211",
                         "ARM1",
                         "HAM_UP_Srinivasan_20211",
                         
                         "DAM_UP_KRASEMAN1",
                         "AGE_UP_MICROGLIA1",
                         
                         "Cytokine_Response1",
                         "AD_TRAJ_LATE1",
                         
                         "CAM1",
                         "diss1", "KEGG_cell_cycle1"
), group.by = "LOY",  
min.color = "white", size = 3,
max.color = "red") %>% ggpar(x.text.angle = 45, rotate = T) + 
  scale_color_gradient2(limits = c(-1,2), oob = scales::squish, low = "blue", high = "red", mid = "grey98")-> p
ggsave(plot = p + NoLegend(), dpi = "retina", units = "in", width = 3, height =4,
       filename = paste0(output,"/UMAP/DOTPLOT_MODULE_PANEL_LOY.pdf"))
ggsave(plot = p, dpi = "retina", units = "in", width = 3, height =5,
       filename = paste0(output,"/UMAP/DOTPLOT_MODULE_PANEL_LEGEND_LOY.pdf"))


reshape2::dcast(p$data, formula = var ~ grouping, value.var = "color") -> df
rownames(df) <- df$var
df[,-1] -> df
pheatmap::pheatmap(df, cluster_cols = T, 
                   breaks = seq(-1, 2, length.out = 11), color = RColorBrewer::brewer.pal(n = 11, name = "Reds")) 



y.labels = c("Stressed | Thrupp",
             "CAM | Prinz",
             "MG AD TRAJ Late | Gerrits",
             "Cytokine response MG | ",
             "Age associated UP | Galatro",
             "DAM UP | Krasemann",
             "HAM UP | Srinivasan",
             "Activated MG | SalaFrigerio",
             "Homeostatic MG | SalaFrigerio",
             "Homeostatic MG | Krasemann")

ggsave(plot = p, dpi = "retina", units = "in", width = 5, height =5,
       filename = paste0(output,"/UMAP/DOTPLOT_PANEL.pdf"))



readRDS(file = "e:/LOY/scLOY/processed_seurat/BRAIN/GSE160936_AD_GLIAL_introns_included.RDS") -> tmp
dittoDimPlot(object = tmp, var = "cell_ident_MV", do.label = F) %>% ggpar(title = "") -> p
p + NoAxes() -> p
ggsave(plot = p, dpi = "retina", units = "in", width = 7, height =5,
       filename = paste0(output,"/UMAP/FULL_COHORT_UMAP.pdf"))




#### TRAJECTORY



readRDS(file = "e:/LOY/scLOY/processed_seurat/BRAIN/GSE160936_AD_GLIAL_introns_included.RDS") -> tmp
#DimPlot(o, group.by = "RNA_snn_res.1.1", reduction = "umap")

### SSC #############
subset(tmp, 
       nCount_RNA >= 1000 & 
         nFeature_RNA >= 800 & 
         specific_tissue == "SSC" & donor_organism.sex == "male" & cell_ident_MV == "Astro") -> o

### run seurat integration
ifnb.list <- SplitObject(o, split.by = "title")

ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  return(x)
})

rownames(o)[grep(x = rownames(o), pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")] -> ribo
ribo[grep(ribo,pattern = "K",invert = T)] -> ribo
rownames(o)[grep(x = rownames(o), pattern = "^MT-")] -> mito
rownames(o)[(rownames(o) %in% Y_scLOY_genes$gene_name)] -> Y
Y_scLOY_genes$gene_name[Y_scLOY_genes$gene_name %in% rownames(o)] -> male_genes

features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 2000)
features <- features[!(features %in% c("XIST",ribo,mito,"JPX"))]
anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
combined <- IntegrateData(anchorset = anchors)
DefaultAssay(combined) <- "integrated"
combined -> o
rm(combined)

rownames(o)[grep(x = rownames(o), pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")] -> ribo
ribo[grep(ribo,pattern = "K",invert = T)] -> ribo

o <- Seurat::ScaleData(o, vars.to.regress = c("nCount_RNA","percent.mt","percent.rb"))
v <- Seurat::VariableFeatures(o)
v <- v[!(v %in% c(ribo,mito,"XIST","JPX"))]
o <- Seurat::RunPCA(o, npcs = 50, features = v, assay = "integrated")
ElbowPlot(o,ndims = 50 )

o <- Seurat::RunUMAP(o,reduction = "pca", dims = 1:12)
o <- Seurat::FindNeighbors(o,reduction = "pca", dims = 1:12)
o <- Seurat::FindClusters(o,resolution = seq(0.1,2,0.1))
DimPlot(o, group.by = "integrated_snn_res.0.5", label = T, cols = dittoColors())
DimPlot(o, group.by = "integrated_snn_res.0.3", label = T, cols = dittoColors())
DimPlot(o, group.by = "title")

FindMarkers(o, ident.1 = 8, ident.2 = NULL, group.by = "integrated_snn_res.0.5",
            only.pos = T, logfc.threshold = 0.75)


Idents(o) <- "integrated_snn_res.0.5"

saveRDS(o, "e:/LOY/scLOY/processed_integrated/ASTROCYTE/IMPERIAL_CCA_SSC_MALE_JULY12.RDS")

### HOMEOSTATIC GENES ##########
output <- "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/FINAL/IMPERIAL/SSC_ASTRO"
dir.create(paste0(output,"/SUBCLUSTERS"), recursive = T)
dir.create(paste0(output,"/UMAP"), recursive = T)
dir.create(paste0(output,"/BARPLOTS"), recursive = T)

DefaultAssay(o) <- "RNA"
p <- plot_density(o, c("GFAP","DPP10","SLC1A2"), reduction = "umap",
                  method = "wkde", adjust = 3,
                  pal = "inferno", joint = T,
                  combine = F)[[4]];
ggsave(plot = p, dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/HOMEO_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(),
       dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/HOMEO.pdf"))

p <- plot_density(o, c("GFAP","DPP10","SLC1A2"), reduction = "umap",
                  method = "wkde", adjust = 2,
                  pal = "inferno", joint = T,
                  combine = F)[[1]];
ggsave(plot = p, dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/GFAP_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(), dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/GFAP.pdf"))

p <- plot_density(o, c("GFAP","DPP10","SLC1A2"), reduction = "umap",
                  method = "wkde", adjust = 2,
                  pal = "inferno", joint = T,
                  combine = F)[[2]];
ggsave(plot = p, dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/DPP10_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(),
       dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/DPP10.pdf"))

p <- plot_density(o, c("GFAP","DPP10","SLC1A2"), reduction = "umap",
                  method = "wkde", adjust = 2,
                  pal = "inferno", joint = T,
                  combine = F)[[3]];
ggsave(plot = p,  dpi = "retina", units = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/SLC1A2_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(),
       dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/SLC1A2.pdf"))


### GENE SETS #############
### find homeostatic and AD associated clusters
### pull these clusters

a <- gene_sets$SalaFrigerio_ARM
list(a[a %in% rownames(o)]) -> a
Seurat::AddModuleScore(object = o, features = a, name = "ARM", assay = "RNA") -> o
plot_density(o, "ARM1", reduction = "umap", pal = "inferno", adjust = 1) -> p
FeaturePlot(o, features = "ARM1", order = T, pt.size = 0.9) + scale_color_viridis_c(option = "B")
ggsave(plot = p + NoAxes(), dpi = "retina", unitsR = "in", width = 4.5, height = 3.5,
       filename = paste0(output,"/SUBCLUSTERS/ARM_LEGEND.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(), dpi = "retina", units = "in", width = 4, height = 4,
       filename = paste0(output,"/SUBCLUSTERS/ARM.pdf"))

###### ANNOTATE #########

subset(o, integrated_snn_res.0.4 == "7", invert = T) -> o
dittoDimPlot(o, var = "integrated_snn_res.0.4", 
             do.label = T, labels.highlight = T, labels.size = 4) -> p
FindMarkers(o, ident.1 = 7, ident.2 = NULL, group.by = "integrated_snn_res.0.4", only.pos = T, logfc.threshold = 0.5)

o@meta.data$MG_subtype2 <- o@meta.data$integrated_snn_res.0.4
levels(o@meta.data$MG_subtype2)[levels(o@meta.data$MG_subtype2)%in%c(5,6,4,2,3)] <- "GFAPhi"
levels(o@meta.data$MG_subtype2)[levels(o@meta.data$MG_subtype2)%in%c(0,1)] <- "GFAPlow"
DimPlot(o, group.by = "MG_subtype2", label = T)

o@meta.data$MG_subtype_C2 <- o@meta.data$integrated_snn_res.0.4
levels(o@meta.data$MG_subtype_C2)[levels(o@meta.data$MG_subtype_C2)%in%c(2,3)] <- "GFAPhi1"
levels(o@meta.data$MG_subtype_C2)[levels(o@meta.data$MG_subtype_C2)%in%c(4,6)] <- "GFAPhi2"
levels(o@meta.data$MG_subtype_C2)[levels(o@meta.data$MG_subtype_C2)%in%c(5)] <- "GFAPhi3"
levels(o@meta.data$MG_subtype_C2)[levels(o@meta.data$MG_subtype_C2)%in%c(1)] <- "GFAPlow1"
levels(o@meta.data$MG_subtype_C2)[levels(o@meta.data$MG_subtype_C2)%in%c(0)] <- "GFAPlow2"
DimPlot(o, group.by = "MG_subtype_C2", label = T)

Idents(o) <- "MG_subtype_C2"
FindAllMarkers(o, logfc.threshold = 0.2, min.pct = 0.1, test.use = "MAST",
               latent.vars = c("nCount_RNA","nFeature_RNA","percent.rb")) -> DE
data.table::fwrite(x = DE, file = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/FINAL/IMPERIAL/SSC_ASTRO/DE_final.txt",
                   sep = "\t")

factor(o@meta.data$MG_subtype_C2, levels = c("GFAPhi1","GFAPhi2","GFAPhi3","GFAPlow1","GFAPlow2")) -> o@meta.data$ord2

colors <- c(dittoColors()[21],dittoColors()[5],dittoColors()[13],
            dittoColors()[1], dittoColors()[17])


dittoDimPlot(o, var = "ord2", color.panel = colors,  
             do.label = T, labels.highlight = T, labels.size = 4) -> p
p %>% ggpar(title = "") -> p
ggsave(plot = p + NoAxes() + NoLegend(), dpi = "retina", units = "in", width = 6, height = 5,
       filename = paste0(output,"/UMAP/NO_AXES_NO_LEG_UMAP.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(), dpi = "retina", units = "in", width = 3, height = 2.5,
       filename = paste0(output,"/UMAP/NO_AXES_NO_LEG_UMAP_SMALL.pdf"))

dittoDimPlot(o, var = "ord2", color.panel = colors,  
             do.label = T, labels.highlight = T, labels.size = 4) -> p
p %>% ggpar(title = "") -> p
ggsave(plot = p + NoAxes() + NoLegend(), dpi = "retina", units = "in", width = 6, height = 5,
       filename = paste0(output,"/UMAP/NO_AXES_NO_LEG_UMAP.pdf"))
ggsave(plot = p + NoAxes() + NoLegend(), dpi = "retina", units = "in", width = 3, height = 2.5,
       filename = paste0(output,"/UMAP/NO_AXES_NO_LEG_UMAP_SMALL.pdf"))

dittoDimPlot(o, var = "ord2", color.panel = colors,  
             do.label = F, labels.highlight = T, labels.size = 4) -> p
p %>% ggpar(title = "") -> p
ggsave(plot = p + NoAxes(), dpi = "retina", units = "in", width = 6, height = 5,
       filename = paste0(output,"/UMAP/NO_AXES_UMAP.pdf"))
ggsave(plot = p + NoAxes(), dpi = "retina", units = "in", width = 3, height = 2.5,
       filename = paste0(output,"/UMAP/NO_AXES_UMAP_SMALL.pdf"))

dittoDimPlot(o, var = "ord2", color.panel = colors,split.by = "diagnosis", size = 0.8,
             do.label = F, labels.highlight = T, labels.size = 4, do.contour = T) -> p
p %>% ggpar(title = "") + NoAxes() + NoLegend() -> p
ggsave(plot = p + NoAxes(), dpi = "retina", units = "in", width = 5, height = 3,
       filename = paste0(output,"/UMAP/AD_SPLIT.pdf"))

#dittoDimPlot(o, var = "ord", color.panel = colors,split.by = "specific_tissue", size = 0.8,
#             do.label = F, labels.highlight = T, labels.size = 4, do.contour = T) -> p
#p %>% ggpar(title = "") + NoAxes() + NoLegend() -> p
#ggsave(plot = p + NoAxes(), dpi = "retina", units = "in", width = 5, height = 3,
#       filename = paste0(output,"/UMAP/TISSUE_SPLIT.pdf"))

dittoDimPlot(o, var = "title", size = 0.8,
             do.label = F) -> p
p %>% ggpar(title = "") + NoAxes() -> p
ggsave(plot = p + NoAxes(), dpi = "retina", units = "in", width = 6, height = 4,
       filename = paste0(output,"/UMAP/SAMPLE_SPLIT.pdf"))

dittoDimPlot(subset(o, nCount_RNA > 2500), var = "LOY", size = 1, split.by = "diagnosis", color.panel = c("red","grey"),
             do.label = F, order = "decreasing") -> p
p %>% ggpar(title = "") + NoAxes() -> p
ggsave(plot = p + NoAxes(), dpi = "retina", units = "in", width = 7.5, height = 3.5,
       filename = paste0(output,"/UMAP/DIAGNOSIS_LOY_2500.pdf"))


dittoDimPlot(o, var = "LOY", size = 1, split.by = "title", color.panel = c("red","grey"),
             do.label = F, order = "decreasing") -> p
p %>% ggpar(title = "") + NoAxes() -> p
ggsave(plot = p + NoAxes(), dpi = "retina", units = "in", width = 8, height = 8,
       filename = paste0(output,"/UMAP/TITLE_LOY.pdf"))

dittoDimPlot(o, var = "ord2", size = 1, split.by = "LOY", color.panel = colors,
             do.label = F, do.contour = T) -> p
p %>% ggpar(title = "") + NoAxes() -> p
ggsave(plot = p + NoAxes(), dpi = "retina", units = "in", width = 7.5, height = 3.5,
       filename = paste0(output,"/UMAP/LOY_CLUSTERS_CONTOUR.pdf"))




################# BAR PLOT  # 

DimPlot(o, group.by = "ord2", label = T, cols = dittoColors())
DimPlot(o, split.by = "diagnosis", label = T, group.by = "integrated_snn_res.0.5", cols = dittoColors())

order <- NULL

dittoBarPlot(o, group.by = "title", color.panel = colors, var.labels.reorder = NULL,
             var = "ord2") %>% ggpar(xlab = "", title = "") -> a
ggsave(plot = a + NoLegend(), dpi = "retina", units = "in", width = 5, height = 3,
       filename = paste0(output,"/BARPLOTS/SAMPLE_PROP_BAR.pdf"))

o@meta.data[o@meta.data$diagnosis=="AD" & o@meta.data$nCount_RNA > 2500,] %>% rownames() -> cells
dittoBarPlot(o, group.by = "title", var = "LOY", var.labels.reorder = c(2,1),
             cells.use = cells,
             color.panel = c("grey","red")) %>% ggpar(xlab = "", legend = "none", title = "") -> b 
o@meta.data[o@meta.data$diagnosis=="Non-disease control" & o@meta.data$nCount_RNA > 2500,] %>% rownames() -> cells
dittoBarPlot(o, group.by = "title", var = "LOY", var.labels.reorder = c(2,1),
             cells.use = cells,
             color.panel = c("grey","red")) %>% ggpar(xlab = "", legend = "none", title = "") -> c 
ggarrange(plotlist = list(b ,c + rremove("ylab") + rremove("y.text") ), ncol = 2, nrow = 1, align = "hv") -> a 
ggsave(plot = a + NoLegend(), dpi = "retina", units = "in", width = 5, height = 3,
       filename = paste0(output,"/BARPLOTS/SAMPLE_PROP_BAR_BY_DIAGNOSIS_2500.pdf"))

#ggarrange(plotlist = list(a,b), ncol = 1, nrow = 2, align = "hv")
ggsave(plot = b + NoLegend(), dpi = "retina", units = "in", width = 3, height = 3,
       filename = paste0(output,"/BARPLOTS/SAMPLE_PROP_LOY_BAR.pdf"))

dittoBarPlot(subset(o, nCount_RNA > 2500), group.by = "ord2", var = "LOY", var.labels.reorder = c(2,1), x.reorder = order,
             color.panel = c("grey","red")) %>% ggpar(xlab = "") -> c
ggsave(plot = c + NoLegend(), dpi = "retina", units = "in", width = 3, height = 3,
       filename = paste0(output,"/BARPLOTS/CLUSTER_PROP_LOY_BAR.pdf"))

c$data %>% dplyr::filter(label == "LOY") -> dat
factor(dat$grouping, levels = dat[order(-dat$percent),]$grouping) -> dat$grouping
k <- colors[c(3,4,5,1,2)]
dat %>% ggbarplot(x = "grouping", y = "percent", fill = "grouping", palette = k) %>% 
  ggpar(legend = "none", xlab = "", ylab = "Percent LOY", x.text.angle = 45, ylim = c(0,0.25)) -> p
ggsave(plot = p + NoLegend(), dpi = "retina", units = "in", width = 4, height = 3,
       filename = paste0(output,"/BARPLOTS/CLUSTER_PROP_LOY_BAR_ranked.pdf"))

dittoBarPlot(o, group.by = "diagnosis", var = "LOY", var.labels.reorder = c(2,1),
             color.panel = c("grey","red")) %>% ggpar(xlab = "") -> d
ggsave(plot = d + NoLegend(), dpi = "retina", units = "in", width = 1, height = 3,
       filename = paste0(output,"/BARPLOTS/DIAGNOSIS_PROP_LOY_BAR.pdf"))

dittoBarPlot(subset(o, nCount_RNA > 2500 & nFeature_RNA > 1000), group.by = "title", var = "LOY", var.labels.reorder = c(2,1),
             color.panel = c("grey","red")) %>% ggpar(xlab = "") -> b 
#ggarrange(plotlist = list(a,b), ncol = 1, nrow = 2, align = "hv")
ggsave(plot = b + NoLegend(), dpi = "retina", units = "in", width = 3, height = 3,
       filename = paste0(output,"/BARPLOTS/SAMPLE_PROP_LOY_BAR_2500.pdf"))

dittoBarPlot(subset(o, nCount_RNA > 2500 & nFeature_RNA > 1000), group.by = "ord", var = "LOY", var.labels.reorder = c(2,1), x.reorder = order,
             color.panel = c("grey","red")) %>% ggpar(xlab = "") -> c
ggsave(plot = c + NoLegend(), dpi = "retina", units = "in", width = 3, height = 3,
       filename = paste0(output,"/BARPLOTS/CLUSTER_PROP_LOY_BAR_2500.pdf"))

dittoBarPlot(subset(o, nCount_RNA > 2500 & nFeature_RNA > 1000), group.by = "diagnosis", var = "LOY", var.labels.reorder = c(2,1),
             color.panel = c("grey","red")) %>% ggpar(xlab = "") -> d
ggsave(plot = d + NoLegend(), dpi = "retina", units = "in", width = 1, height = 3,
       filename = paste0(output,"/BARPLOTS/DIAGNOSIS_PROP_LOY_BAR_2500.pdf"))


dittoBarPlot(o, group.by = "diagnosis", var.labels.reorder = order, color.panel = colors,
             var = "ord2") %>% ggpar(xlab = "", title = "") -> f
ggsave(plot = e + NoLegend(), dpi = "retina", units = "in", width = 1.5, height = 4,
       filename = paste0(output,"/BARPLOTS/DIAGNOSIS_PROP_BAR.pdf"))

dittoBarPlot(o, group.by = "donor_organism.sex", var.labels.reorder = order, color.panel = colors,
             var = "ord2") %>% ggpar(xlab = "", title = "", ylab = "") -> e
ggsave(plot = e + NoLegend(), dpi = "retina", units = "in", width = 1.7, height = 4,
       filename = paste0(output,"/BARPLOTS/SEX_PROP_BAR.pdf"))

ggarrange(plotlist = list(f + NoLegend(),e + NoLegend()), ncol = 2, nrow = 1, align = "hv") -> p
ggsave(plot = p, dpi = "retina", units = "in", width = 3, height = 4,
       filename = paste0(output,"/BARPLOTS/SEX_DIAGNOSIS_PROP_BAR.pdf"))

paste0(o@meta.data$donor_organism.sex, "_", o@meta.data$diagnosis) -> o@meta.data$diagnosis_sex
dittoBarPlot(o, group.by = "donor_organism.sex", var.labels.reorder = order,
             var = "integrated_snn_res.0.5") %>% ggpar(xlab = "") -> e

o@meta.data$diagnosis_sex <- factor(o@meta.data$diagnosis_sex, levels = c("male_AD","female_AD","male_Non-disease control","female_Non-disease control"))
dittoBarPlot(o, group.by = "diagnosis_sex", var.labels.reorder = order,
             var = "integrated_snn_res.0.5") %>% ggpar(xlab = "") -> e



### GENE PANELS

plot_density(o, "CX3CR1", reduction = "umap", pal = "inferno", adjust = 1) + NoAxes() -> p1
plot_density(o, "CD163", reduction = "umap", pal = "inferno", adjust = 1) + NoAxes() -> p2
plot_density(o, "CD83", reduction = "umap", pal = "inferno", adjust = 1) + NoAxes() -> p3
plot_density(o, "TGFBI", reduction = "umap", pal = "inferno", adjust = 1) + NoAxes() -> p9
plot_density(o, "P2RY12", reduction = "umap", pal = "inferno", adjust = 1) + NoAxes() -> p4
plot_density(o, "CTNNA2", reduction = "umap", pal = "inferno", adjust = 1) + NoAxes() -> p5
plot_density(o, "CCL3", reduction = "umap", pal = "inferno", adjust = 1) + NoAxes() -> p6
plot_density(o, "VCAN", reduction = "umap", pal = "inferno", adjust = 1) + NoAxes() -> p7
plot_density(o, "BRIP1", reduction = "umap", pal = "inferno", adjust = 1) + NoAxes() -> p8
ggarrange(plotlist = list(p1,p4,p2,p3,p9,p6,p7,p8,p5), ncol = 3, nrow = 3) -> p

ggsave(plot = p, dpi = "retina", units = "in", width = 8, height = 7,
       filename = paste0(output,"/UMAP/MARKERS_PANEL.pdf"))

DotPlot(o, group.by = "ord", 
        features = c("C3","CD74","CX3CR1","P2RY12",
                     "NAV3","CTNNA2",
                     "CD163","TGFBI","DPYD","F13A1",
                     "CD83","LPL","APOE","TREM2","CCL3","IL1B",
                     "GPNMB",
                     "BRIP1","CENPP","BRCA2",
                     "F13A1","LYVE1",
                     "VCAN")) %>% ggpar(x.text.angle = 45, rotate = T)

dittoDotPlot(o, vars = c("C3","CD74","CX3CR1","P2RY12",
                         "GPM6A","CTNNA2",
                         "CD163","TGFBI","DPYD",
                         "CD83","LPL","APOE","TREM2","MYO1E","CCL3","CCL4","IL1B",
                         "BRIP1","CENPP","BRCA2",
                         "F13A1","LYVE1",
                         "VCAN"), group.by = "ord2", 
             min.color = "white",
             max.color = "red") %>% ggpar(x.text.angle = 45, rotate = T) -> p
ggsave(plot = p, dpi = "retina", units = "in", width = 5, height =5,
       filename = paste0(output,"/UMAP/DOTPLOT_PANEL.pdf"))




o@meta.data %>% dplyr::select(c("HOMEOSTATIC_KRASEMANN1",
                                "SalaFrigerio_HM1",
                                "ARM1",
                                "HAM_UP_Srinivasan_20211",
                                "HAM_DOWN_Srinivasan_20211",
                                "DAM_UP_KRASEMAN1",
                                "AGE_UP_MICROGLIA1",
                                "AGE_DOWN_MICROGLIA1",
                                "Cytokine_Response1",
                                "AD_TRAJ_LATE1",
                                "AD_traj_early1",
                                "CAM1",
                                "diss1","ord"
)) -> df



dittoDotPlot(o, vars = c("HOMEOSTATIC_KRASEMANN1",
                         "SalaFrigerio_HM1",
                         "AD_traj_early1",
                         "AGE_DOWN_MICROGLIA1",
                         "HAM_DOWN_Srinivasan_20211",
                         "ARM1",
                         "HAM_UP_Srinivasan_20211",
                         
                         "DAM_UP_KRASEMAN1",
                         "AGE_UP_MICROGLIA1",
                         
                         "Cytokine_Response1",
                         "AD_TRAJ_LATE1",
                         
                         "CAM1",
                         "diss1", "KEGG_cell_cycle1"
), group.by = "ord2",  
min.color = "white", size = 3,
max.color = "red") %>% ggpar(x.text.angle = 45, rotate = T) + 
  scale_color_gradient2(limits = c(-1,2), oob = scales::squish, low = "blue", high = "red", mid = "grey98")-> p
ggsave(plot = p + NoLegend(), dpi = "retina", units = "in", width = 5, height =4,
       filename = paste0(output,"/UMAP/DOTPLOT_MODULE_PANEL.pdf"))
ggsave(plot = p, dpi = "retina", units = "in", width = 5, height =5,
       filename = paste0(output,"/UMAP/DOTPLOT_MODULE_PANEL_LEGEND.pdf"))

dittoDotPlot(o, vars = c("HOMEOSTATIC_KRASEMANN1",
                         "SalaFrigerio_HM1",
                         "AD_traj_early1",
                         "AGE_DOWN_MICROGLIA1",
                         "HAM_DOWN_Srinivasan_20211",
                         "ARM1",
                         "HAM_UP_Srinivasan_20211",
                         
                         "DAM_UP_KRASEMAN1",
                         "AGE_UP_MICROGLIA1",
                         
                         "Cytokine_Response1",
                         "AD_TRAJ_LATE1",
                         
                         "CAM1",
                         "diss1", "KEGG_cell_cycle1"
), group.by = "LOY",  
min.color = "white", size = 3,
max.color = "red") %>% ggpar(x.text.angle = 45, rotate = T) + 
  scale_color_gradient2(limits = c(-1,2), oob = scales::squish, low = "blue", high = "red", mid = "grey98")-> p
ggsave(plot = p + NoLegend(), dpi = "retina", units = "in", width = 3, height =4,
       filename = paste0(output,"/UMAP/DOTPLOT_MODULE_PANEL_LOY.pdf"))
ggsave(plot = p, dpi = "retina", units = "in", width = 3, height =5,
       filename = paste0(output,"/UMAP/DOTPLOT_MODULE_PANEL_LEGEND_LOY.pdf"))


reshape2::dcast(p$data, formula = var ~ grouping, value.var = "color") -> df
rownames(df) <- df$var
df[,-1] -> df
pheatmap::pheatmap(df, cluster_cols = T, 
                   breaks = seq(-1, 2, length.out = 11), color = RColorBrewer::brewer.pal(n = 11, name = "Reds")) 



y.labels = c("Stressed | Thrupp",
             "CAM | Prinz",
             "MG AD TRAJ Late | Gerrits",
             "Cytokine response MG | ",
             "Age associated UP | Galatro",
             "DAM UP | Krasemann",
             "HAM UP | Srinivasan",
             "Activated MG | SalaFrigerio",
             "Homeostatic MG | SalaFrigerio",
             "Homeostatic MG | Krasemann")

ggsave(plot = p, dpi = "retina", units = "in", width = 5, height =5,
       filename = paste0(output,"/UMAP/DOTPLOT_PANEL.pdf"))



readRDS(file = "e:/LOY/scLOY/processed_seurat/BRAIN/GSE160936_AD_GLIAL_introns_included.RDS") -> tmp
dittoDimPlot(object = tmp, var = "cell_ident_MV", do.label = F) %>% ggpar(title = "") -> p
p + NoAxes() -> p
ggsave(plot = p, dpi = "retina", units = "in", width = 7, height =5,
       filename = paste0(output,"/UMAP/FULL_COHORT_UMAP.pdf"))




##### SUPP FIG LOY PROPS

readRDS(file = "e:/LOY/scLOY/processed_seurat/BRAIN/GSE160936_AD_GLIAL_introns_included.RDS") -> tmp
#DimPlot(o, group.by = "RNA_snn_res.1.1", reduction = "umap")


output <- ""
subset(tmp, 
       nCount_RNA >= 2500 & 
         nFeature_RNA >= 1000 & 
         donor_organism.sex == "male") -> o

subset(o, cell_ident_MV == "Neuron") -> o.
dittoPlot(o., var = "nCount_RNA", group.by = "LOY", plots = c("vlnplot","boxplot"), color.panel = c("red","grey"),
          color.by = "LOY") %>% ggpar(xlab = "", title = "", ylab = "", legend = "")
ggsave()




