### combine all included nuclei using harmony
### plot using TSNE 
library(Biobase)
library(GEOquery)
library(Seurat)
library(ggpubr)
library(tidyverse)
library(scibetR)
library(viridis)
library(reshape2)
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
edb <- EnsDb.Hsapiens.v86


devtools::source_url("https://raw.githubusercontent.com/michaelcvermeulen/microglia-loss-of-y/main/scripts/LOYtools/functions.R")
list.files(path = "E:/LOY/scLOY/processed_seurat/MERGED_COHORT/", full.names = T) -> files


## list of included datasets
inc <- c("GSE160936","GSE178265","GSE174367",
"GSE148822","GSE174332","GSE137444","GSE157783",
"GSE152058","syn12514624","syn21670836","GSE160189","GSE126836","GSE163577","GSE119212","GSE157760",
"PRJNA530977","GSE157827","GSE153807","syn18485175","GSE144136","GSE164485")

  

files[grepl( x = files, pattern = paste(paste(inc, collapse="|")))] -> files


lapply(X = 1:length(files), FUN = function(z){
  
  readRDS(file = files[z]) -> o
}) -> l
l[[1]] -> o
for(z in c(2:6)){   message(z); merge(o,l[[z]]) -> o }
saveRDS(object = o, file = "e:/LOY/scLOY/processed_seurat/MERGED_COHORT_ALIGNED/dataset1-6.RDS")

l[[7]] -> o
for(z in c(8:11)){   message(z); merge(o,l[[z]]) -> o }
saveRDS(object = o, file = "e:/LOY/scLOY/processed_seurat/MERGED_COHORT_ALIGNED/dataset7-12.RDS")

l[[13]] -> o
for(z in c(14:18)){   message(z); merge(o,l[[z]]) -> o }
saveRDS(object = o, file = "e:/LOY/scLOY/processed_seurat/MERGED_COHORT_ALIGNED/dataset8-18.RDS")

l[[19]] -> o
for(z in c(20:21)){   message(z); merge(o,l[[z]]) -> o }
saveRDS(object = o, file = "e:/LOY/scLOY/processed_seurat/MERGED_COHORT_ALIGNED/dataset19-21.RDS")

o <- subset(o, nCount_RNA >= 1500 & nFeature_RNA >= 1000)

rownames(o)[grep(x = rownames(o), pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")] -> ribo
ribo[grep(ribo,pattern = "K",invert = T)] -> ribo
rownames(o)[grep(x = rownames(o), pattern = "^MT-")] -> mito
rownames(o)[(rownames(o) %in% Y_scLOY_genes$gene_name)] -> Y
Y_scLOY_genes$gene_name[Y_scLOY_genes$gene_name %in% rownames(o)] -> male_genes
PAR_scLOY_genes$gene_name[PAR_scLOY_genes$gene_name %in% rownames(o)] %>% unique() -> PAR_genes

o <- Seurat::NormalizeData(o, scale.factor = 10000)
o <- Seurat::FindVariableFeatures(o, selection.method = "vst", nfeatures = 1500)
diss <- c("RGS1","DUSP1","AC022217.3","DUSP1","HSPA1A","HSPA1B","DDIT4","SLC2A3","FOS","JUN",
          "AC131944.1","HIF1A-AS3","HIST1H2AC","HIST1H2BD","HIST1H2BG","HIST2H2BE","RGS2","RHOB",
          "HIST1H2BJ","TEX14","UBC","HSP90AA1","AL118516.1","JUNB","HIST2H2AA4","HIST1H1C","BTG1","BTG2","DUSP5")
list(diss[diss %in% rownames(o)]) -> diss
Seurat::AddModuleScore(object = o, features = diss, name = "diss") -> o

o <- Seurat::ScaleData(o, vars.to.regress = c("nCount_RNA","percent.mt","percent.rb"))
v <- Seurat::VariableFeatures(o)
v <- v[!(v %in% c(ribo,mito,"XIST","JPX"))]

o <- Seurat::RunPCA(o, npcs = 50, features = v, assay = "RNA")
Seurat::ElbowPlot(o, ndims = 50)
o <- Seurat::RunUMAP(o,reduction = "pca", dims = 1:20)
o <- Seurat::FindNeighbors(o,reduction = "pca", dims = 1:20)
o <- Seurat::FindClusters(o,resolution = 0.1)
DimPlot(o, group.by = "orig.ident")
DimPlot(o, group.by = "RNA_snn_res.0.1",label=T,raster=T)
DimPlot(o, group.by = "donor",label=T,raster=T)
DimPlot(o, group.by = "GEO",label=T,raster=T)

o <- harmony::RunHarmony(o,
                         group.by.vars = c("GEO","specific_tissue","donor"),
                         assay.use = "RNA",
                         plot_convergence = T )
o <- Seurat::RunUMAP(o,reduction = "harmony", dims = 1:20)
o <- Seurat::FindNeighbors(o,reduction = "harmony", dims = 1:20)
o <- Seurat::FindClusters(o,resolution = seq(0.1,2,0.1))


