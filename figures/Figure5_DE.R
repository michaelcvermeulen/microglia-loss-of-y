library(ggrepel)
library(Seurat)
library(ggpubr)
library(tidyverse)
library(viridis)
library(tibble)
library(reshape2)
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
edb <- EnsDb.Hsapiens.v86

setwd("H:/LOY/scLOY/scLOY/")
source("E:/LOY/Example/PACKAGE/call_LOY.R")
source("E:/LOY/scLOY/scripts/loy_functions.R")
source("E:/LOY/scLOY/scripts/loy_summary_functions.R")
source('H:/LOY/scLOY/scLOY/R/load.R')
source('H:/LOY/scLOY/scLOY/data_raw/build_data_files.R')

devtools::source_url("https://raw.githubusercontent.com/michaelcvermeulen/microglia-loss-of-y/main/scripts/LOYtools/functions.R")

list() -> l
#### remake heatmap for figure 3


### GSE160936
readRDS(file = "e:/LOY/scLOY/processed_seurat/MERGED_COHORT/GSE160936_cleaned_filtered.RDS") -> o
o@meta.data$`brain region (somatosensory or entorhinal cortex):ch1` -> o@meta.data$brain_region
o@meta.data$donor_organism.biomaterial_core.biomaterial_id <- o@meta.data$title
gsub(x = o@meta.data$donor_organism.biomaterial_core.biomaterial_id , pattern = " SSC", replacement = "") -> o@meta.data$donor_organism.biomaterial_core.biomaterial_id
gsub(x = o@meta.data$donor_organism.biomaterial_core.biomaterial_id , pattern = " EC", replacement = "") -> o@meta.data$donor_organism.biomaterial_core.biomaterial_id

list( ) -> l
# oo@meta.data[oo@meta.data$cell_ident_MV%in%c("CAM","Microglia"),]$cell_ident_MV <- "Microglia/CAM"
DE1 <- run_DE_multi(s_obj = subset(o, brain_region =="SSC"),
                    sample = o@meta.data$donor_organism.biomaterial_core.biomaterial_id %>% unique(),
                    latent = c("nCount_RNA","nFeature_RNA","percent.mt","donor_organism.biomaterial_core.biomaterial_id"),
                    cell_type = c("Microglia"), nUMI = 2000, nFeature = 1000, logfc.threshold = 0, min.pct = 0.05)
DE1$ID <- "GSE160936 SSC"
rlist::list.append(l,DE1) -> l


DE2 <- run_DE_multi(s_obj = subset(o, brain_region =="EC"),
                    sample = o@meta.data$donor_organism.biomaterial_core.biomaterial_id %>% unique(),
                    latent = c("nCount_RNA","nFeature_RNA","percent.mt","donor_organism.biomaterial_core.biomaterial_id"),
                    cell_type = c("Microglia"), nUMI = 2000, nFeature = 1000, logfc.threshold = 0, min.pct = 0.05)
DE2$ID <- "GSE106936 EC"
rlist::list.append(l,DE2) -> l

DE3 <- run_DE_multi(s_obj = subset(o),
                    sample = o@meta.data$donor_organism.biomaterial_core.biomaterial_id %>% unique(),
                    latent = c("nCount_RNA","nFeature_RNA","percent.mt",
                               "donor_organism.biomaterial_core.biomaterial_id","brain_region"),
                    cell_type = c("Microglia"), nUMI = 2000, nFeature = 1000, logfc.threshold = 0, min.pct = 0.05)
DE3$ID <- "GSE106936 FULL"
rlist::list.append(l,DE3) -> l

saveRDS(l, file = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/FINAL/DE/SEPT_21/GSE106936_DE.RDS")


### syn12514624
readRDS(file = "e:/LOY/scLOY/processed_seurat/MERGED_COHORT/syn12514624_cleaned_filtered.RDS") -> o

list( ) -> l
# oo@meta.data[oo@meta.data$cell_ident_MV%in%c("CAM","Microglia"),]$cell_ident_MV <- "Microglia/CAM"
DE1 <- run_DE_multi(s_obj = subset(o),
                    sample = o@meta.data$donor_organism.biomaterial_core.biomaterial_id %>% unique(),
                    latent = c("nCount_RNA","nFeature_RNA","percent.mt","donor_organism.biomaterial_core.biomaterial_id"),
                    cell_type = c("Microglia"), nUMI = 2000, nFeature = 1000, logfc.threshold = 0, min.pct = 0.05)
DE1$ID <- "syn12514624 DLFPC"
rlist::list.append(l,DE1) -> l
saveRDS(l, file = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/FINAL/DE/SEPT_21/syn12514624_DE.RDS")

### GSE178265
readRDS(file = "e:/LOY/scLOY/processed_seurat/MERGED_COHORT/GSE178265_cleaned_filtered.RDS") -> o

list( ) -> l
# oo@meta.data[oo@meta.data$cell_ident_MV%in%c("CAM","Microglia"),]$cell_ident_MV <- "Microglia/CAM"
DE1 <- run_DE_multi(s_obj = subset(o, specific_tissue == "Substantia nigra pars compacta"),
                    sample = o@meta.data$donor_organism.biomaterial_core.biomaterial_id %>% unique(),
                    latent = c("nCount_RNA","nFeature_RNA","percent.mt","donor_organism.biomaterial_core.biomaterial_id"),
                    cell_type = c("Microglia"), nUMI = 2000, nFeature = 1000, logfc.threshold = 0, min.pct = 0.05)
DE1$ID <- "GSE178265 SN"
rlist::list.append(l,DE1) -> l

# oo@meta.data[oo@meta.data$cell_ident_MV%in%c("CAM","Microglia"),]$cell_ident_MV <- "Microglia/CAM"
DE1 <- run_DE_multi(s_obj = subset(o, specific_tissue == "Caudate nucleus"),
                    sample = o@meta.data$donor_organism.biomaterial_core.biomaterial_id %>% unique(),
                    latent = c("nCount_RNA","nFeature_RNA","percent.mt","donor_organism.biomaterial_core.biomaterial_id"),
                    cell_type = c("Microglia"), nUMI = 2000, nFeature = 1000, logfc.threshold = 0, min.pct = 0.05)
DE1$ID <- "GSE178265 CN"
rlist::list.append(l,DE1) -> l
saveRDS(l, file = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/FINAL/DE/SEPT_21/GSE178265_DE.RDS")


### GSE174367
readRDS(file = "e:/LOY/scLOY/processed_seurat/MERGED_COHORT/GSE174367_cleaned_filtered.RDS") -> o

list( ) -> l
# oo@meta.data[oo@meta.data$cell_ident_MV%in%c("CAM","Microglia"),]$cell_ident_MV <- "Microglia/CAM"
DE1 <- run_DE_multi(s_obj = subset(o),
                    sample = o@meta.data$donor_organism.biomaterial_core.biomaterial_id %>% unique(),
                    latent = c("nCount_RNA","nFeature_RNA","percent.mt","orig.ident"),
                    cell_type = c("Microglia"), nUMI = 2000, nFeature = 1000, logfc.threshold = 0, min.pct = 0.05)
DE1$ID <- "GSE174367 PFC"
rlist::list.append(l,DE1) -> l
saveRDS(l, file = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/FINAL/DE/SEPT_21/GSE174367_DE_OLD.RDS")


# GSE148822

readRDS(file = "e:/LOY/scLOY/processed_seurat/MERGED_COHORT/GSE148822_cleaned_filtered.RDS") -> o
list( ) -> l
# oo@meta.data[oo@meta.data$cell_ident_MV%in%c("CAM","Microglia"),]$cell_ident_MV <- "Microglia/CAM"

o@meta.data$donor -> o@meta.data$donor_organism.biomaterial_core.biomaterial_id
DE1 <- run_DE_multi(s_obj = subset(o, specific_tissue == "occipitotemporal cortex"),
                    sample = o@meta.data$donor_organism.biomaterial_core.biomaterial_id %>% unique(),
                    latent = c("nCount_RNA","nFeature_RNA","percent.mt","donor_organism.biomaterial_core.biomaterial_id"),
                    cell_type = c("Microglia"), nUMI = 2000, nFeature = 1000, logfc.threshold = 0, min.pct = 0.05)
DE1$ID <- "GSE148822 OTC"
rlist::list.append(l,DE1) -> l

o@meta.data$donor -> o@meta.data$donor_organism.biomaterial_core.biomaterial_id
DE1 <- run_DE_multi(s_obj = subset(o, specific_tissue == "occipital cortex"),
                    sample = o@meta.data$donor_organism.biomaterial_core.biomaterial_id %>% unique(),
                    latent = c("nCount_RNA","nFeature_RNA","percent.mt","donor_organism.biomaterial_core.biomaterial_id"),
                    cell_type = c("Microglia"), nUMI = 2000, nFeature = 1000, logfc.threshold = 0, min.pct = 0.05)
DE1$ID <- "GSE148822 OC"
rlist::list.append(l,DE1) -> l


saveRDS(l, file = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/FINAL/DE/SEPT_21/GSE148822_DE.RDS")


# GSE137444
readRDS(file = "e:/LOY/scLOY/processed_seurat/MERGED_COHORT/GSE137444_cleaned_filtered.RDS") -> o
list( ) -> l

DE1 <- run_DE_multi(s_obj = subset(o),
                    sample = o@meta.data$donor_organism.biomaterial_core.biomaterial_id %>% unique(),
                    latent = c("nCount_RNA","nFeature_RNA","percent.mt","donor_organism.biomaterial_core.biomaterial_id"),
                    cell_type = c("Microglia"), nUMI = 2000, nFeature = 1000, logfc.threshold = 0, min.pct = 0.05)
DE1$ID <- "GSE137444 TC"
rlist::list.append(l,DE1) -> l
# oo@meta.data[oo@meta.data$cell_ident_MV%in%c("CAM","Microglia"),]$cell_ident_MV <- "Microglia/CAM"
saveRDS(l, file = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/FINAL/DE/SEPT_21/GSE1374442_DE.RDS")


# GSE174332
readRDS(file = "e:/LOY/scLOY/processed_seurat/MERGED_COHORT/GSE174332_cleaned_filtered.RDS") -> o
list( ) -> l

o@meta.data[o@meta.data$cell_ident_MV=="Cycling microglia",]$cell_ident_MV <- "Microglia"
DE1 <- run_DE_multi(s_obj = subset(o),
                    sample = o@meta.data$donor_organism.biomaterial_core.biomaterial_id %>% unique(),
                    latent = c("nCount_RNA","nFeature_RNA","percent.mt","donor_organism.biomaterial_core.biomaterial_id"),
                    cell_type = c("Microglia"), nUMI = 2000, nFeature = 1000, logfc.threshold = 0, min.pct = 0.05)
DE1$ID <- "GSE174332 MC"
rlist::list.append(l,DE1) -> l
# oo@meta.data[oo@meta.data$cell_ident_MV%in%c("CAM","Microglia"),]$cell_ident_MV <- "Microglia/CAM"
saveRDS(l, file = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/FINAL/DE/SEPT_21/GSE174332_DE.RDS")

# GSE157783
readRDS(file = "e:/LOY/scLOY/processed_seurat/MERGED_COHORT/GSE157783_cleaned_filtered.RDS") -> o
list( ) -> l

DE1 <- run_DE_multi(s_obj = subset(o),
                    sample = o@meta.data$donor_organism.biomaterial_core.biomaterial_id %>% unique(),
                    latent = c("nCount_RNA","nFeature_RNA","percent.mt","donor_organism.biomaterial_core.biomaterial_id"),
                    cell_type = c("Microglia"), nUMI = 2000, nFeature = 1000, logfc.threshold = 0, min.pct = 0.05)
DE1$ID <- "GSE157783 CN"
rlist::list.append(l,DE1) -> l
# oo@meta.data[oo@meta.data$cell_ident_MV%in%c("CAM","Microglia"),]$cell_ident_MV <- "Microglia/CAM"
saveRDS(l, file = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/FINAL/DE/SEPT_21/GSE157783_DE.RDS")


# GSE183068
readRDS(file = "E:/LOY/scLOY/processed_seurat/BRAIN/GSE183068_AD_PFC_introns_included.RDS") -> o
list( ) -> l

DE1 <- run_DE_multi(s_obj = subset(o),
                    sample = o@meta.data$donor_organism.biomaterial_core.biomaterial_id %>% unique(),
                    latent = c("nCount_RNA","nFeature_RNA","percent.mt","donor_organism.biomaterial_core.biomaterial_id"),
                    cell_type = c("Microglia"), nUMI = 2000, nFeature = 1000, logfc.threshold = 0, min.pct = 0.05)
DE1$ID <- "GSE183068 PFC"
rlist::list.append(l,DE1) -> l
# oo@meta.data[oo@meta.data$cell_ident_MV%in%c("CAM","Microglia"),]$cell_ident_MV <- "Microglia/CAM"
saveRDS(l, file = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/FINAL/DE/SEPT_21/GSE183068_DE.RDS")

# GSE167492
readRDS(file = "E:/LOY/scLOY/processed_seurat/BRAIN/GSE167492_AD_DLPFC_introns_included.RDS") -> o
readRDS(file = "E:/LOY/scLOY/processed_seurat/BRAIN/GSE167490_AD_DLPFC_introns_included.RDS") -> oo
merge(o,oo) -> o
list( ) -> l

subset(o, donor_organism.sex == "male") -> o

DE1 <- run_DE_multi(s_obj = o,
                    sample = o@meta.data$donor_organism.biomaterial_core.biomaterial_id %>% unique(),
                    latent = c("nCount_RNA","nFeature_RNA","percent.mt","donor_organism.biomaterial_core.biomaterial_id","GEO"),
                    cell_type = c("Microglia"), nUMI = 2000, nFeature = 1000, logfc.threshold = 0, min.pct = 0.05)
DE1$ID <- "GSE167494 DLPFC"
rlist::list.append(l,DE1) -> l
# oo@meta.data[oo@meta.data$cell_ident_MV%in%c("CAM","Microglia"),]$cell_ident_MV <- "Microglia/CAM"
saveRDS(l, file = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/FINAL/DE/SEPT_21/GSE167494_DE.RDS")


# GSE178175
# GSE173731
# GSE179590
# GSE186538
# GSE189202
readRDS(file = "E:/LOY/scLOY/processed_seurat/BRAIN/GSE183068_AD_PFC_introns_included.RDS") -> o
list( ) -> l

DE1 <- run_DE_multi(s_obj = subset(o),
                    sample = o@meta.data$donor_organism.biomaterial_core.biomaterial_id %>% unique(),
                    latent = c("nCount_RNA","nFeature_RNA","percent.mt","donor_organism.biomaterial_core.biomaterial_id"),
                    cell_type = c("Microglia"), nUMI = 2000, nFeature = 1000, logfc.threshold = 0, min.pct = 0.05)
DE1$ID <- "GSE183068 PFC"
rlist::list.append(l,DE1) -> l
# oo@meta.data[oo@meta.data$cell_ident_MV%in%c("CAM","Microglia"),]$cell_ident_MV <- "Microglia/CAM"
saveRDS(l, file = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/FINAL/DE/SEPT_21/GSE183068_DE.RDS")




### BUILD HEATMAP

list.files(path = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/FINAL/DE/SEPT_21/", full.names = T) -> files
lapply(files, function(x){

  readRDS(file = x) -> tmp
  rownames(tmp) -> tmp$gene
  data.table::rbindlist(tmp) -> tmp
  return(tmp)
})-> l
do.call("rbind",l) -> l


l[l$ID!= "GSE106936 FULL",] -> l

l[l$LOY_cells > 100, ] -> l



#### pull genes first
l[l$avg_log2FC >= 0.2 & l$FDR < 0.2, ]$gene %>% table %>% as.data.frame() -> df 
df[df$Freq>2,]$. -> l1

l[l$avg_log2FC <= -0.2 & l$FDR < 0.2, ]$gene %>% table %>% as.data.frame() -> df 
df[df$Freq>2,]$. -> l2
c(l1,l2) -> df

genes <- df

genes <- c("UTY","RPS4Y1","NLGN4Y",
           "DHRSX","CSF2RA","CD99","IL3RA","AKAP17A","SLC25A6","CD99P1",
           "PARD3B","S100Z","NAV3","ST6GALNAC3","TMEM176B","PRKCA",
           "ROBO1","RBFOX2","APOE",
           "FCGBP","DOCK1","DPYD") %>% unique

cohort <- c("GSE106936 EC","GSE160936 SSC","syn12514624 DLFPC", 
            "GSE183068 PFC",
            "GSE178265 SN",
            "GSE174367 PFC", "GSE167494 DLPFC", "GSE148822 OTC", 
            "GSE148822 OC",
            "GSE174332 MC","GSE137444 TC")

l[l$gene %in% genes,] -> l.


l. %>% dplyr::select(p_val_adj, p_val, avg_log2FC, gene, ID) -> mat

reshape2::dcast(data = mat, formula = ID ~ gene, value.var = "avg_log2FC") -> mat
mat$ID -> rownames(mat)
mat[,-1] -> mat

mat <- mat[cohort,genes]



# make sig * matrix
l. %>% dplyr::select(p_val_adj, p_val, avg_log2FC, gene, ID) -> sig

reshape2::dcast(data = sig, formula = ID ~ gene, value.var = "p_val_adj") -> sig
sig$ID -> rownames(sig)
sig[,-1] -> sig

sig <- sig[cohort,
           genes]



ifelse(test = sig < 0.001, yes = "***", no =
         ifelse(sig < 0.01, yes = "**", no =
                  ifelse(sig < 0.1, yes = "*", no = ""))) -> sig

sig[is.na(sig)] <- ""

colorRampPalette(colors = c("#196889","#FFFFFF","#DE2126")) -> k
breaksList = seq(-0.75, 0.75, by = 0.01)
pheatmap(mat = mat, breaks = breaksList, display_numbers = sig, angle_col = 45,
         cluster_cols = F, na_col = "grey40", gaps_col = c(3,10,16), fontsize_number = c(12),
         cluster_rows = F, border_color = "grey40",
         color = k(length(breaksList))
         )





#### DE gene table (number of DE genes)


list.files(path = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/FINAL/DE/SEPT_21/", full.names = T) -> files
lapply(files, function(x){

  readRDS(file = x) -> tmp
  rownames(tmp) -> tmp$gene
  data.table::rbindlist(tmp) -> tmp
  return(tmp)
})-> l
do.call("rbind",l) -> l


l[l$ID!= "GSE106936 FULL",] -> l

l[l$LOY_cells > 100, ] -> l


l[l$p_val_adj < 0.1, ] -> l

Repitools::annoGR2DF(genes(edb)) -> GL

dplyr::select(GL, gene_name, chr) -> dat

dat[dat$chr %in% c(1:22,"X","Y","MT"),] -> dat
dplyr::distinct(dat, gene_name, .keep_all = T) -> dat


dplyr::left_join(x = l, y = dat, by = c("gene"="gene_name")) -> l



# total sig genes by cohort
l[l$p_val_adj < 0.1,]$ID %>% table() %>% unclass() %>% as.data.frame() -> total_sig
names(total_sig) <- "TOTAL_SIG"
rownames(total_sig) -> total_sig$cohort
total_sig[,c(2,1)] -> total_sig

l[l$p_val_adj < 0.1 & l$avg_log2FC > 0 & l$chr %in% as.character(c(1:22)),]$ID %>% table() %>% unclass() %>% as.data.frame() -> up
names(up) <- "UP_SIG"
rownames(up) -> up$cohort
left_join(x = total_sig, y = up, by = c("cohort")) -> total_sig

l[l$p_val_adj < 0.1 & l$avg_log2FC < 0 & l$chr %in% as.character(c(1:22)),]$ID %>% table() %>% unclass() %>% as.data.frame() -> down
names(down) <- "DOWN_SIG"
rownames(down) -> down$cohort
left_join(x = total_sig, y = down, by = c("cohort")) -> total_sig


l[l$p_val_adj < 0.1 & l$avg_log2FC < 0 &
    l$chr=="Y",]$ID %>% table() %>% unclass() %>% as.data.frame() -> Y
names(Y) <- "DOWN_Y"
rownames(Y) -> Y$cohort
left_join(x = total_sig, y = Y, by = c("cohort")) -> total_sig


l[l$p_val_adj < 0.1 & l$avg_log2FC < 0 &
    l$gene %in% PAR_scLOY_genes$gene_name,]$ID %>% table() %>% unclass() %>% as.data.frame() -> PAR
names(PAR) <- "DOWN_PAR"
rownames(PAR) -> PAR$cohort
left_join(x = total_sig, y = PAR, by = c("cohort")) -> total_sig

l[l$p_val_adj < 0.1 & l$avg_log2FC < 0 &
    l$gene == "MT",]$ID %>% table() %>% unclass() %>% as.data.frame() -> PAR
names(PAR) <- "DOWN_MT"
rownames(PAR) -> PAR$cohort
left_join(x = total_sig, y = PAR, by = c("cohort")) -> total_sig

data.table::fwrite(x = total_sig, file = "a:/Dropbox/LOY_microglia_paper/GR_RESUB/MANUSCRIPT/TABLES/microglia_DE_chr.csv")



###### HEATMAP DATA (BESIDE HEATMAP)

list.files(path = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/FINAL/DE/SEPT_21/", full.names = T) -> files
lapply(files, function(x){

  readRDS(file = x) -> tmp
  rownames(tmp) -> tmp$gene
  data.table::rbindlist(tmp) -> tmp
  return(tmp)
})-> l
do.call("rbind",l) -> l


l[l$ID!= "GSE106936 FULL",] -> l

l[l$LOY_cells > 100, ] -> l



#### pull genes first

genes <- c("UTY","RPS4Y1","NLGN4Y",
           "DHRSX","CSF2RA","CD99","IL3RA","AKAP17A","SLC25A6","CD99P1",
           "PARD3B","S100Z","NAV3","ST6GALNAC3","TMEM176B","PRKCA",
           "ROBO1","RBFOX2","APOE",
           "FCGBP","DOCK1","DPYD") %>% unique

cohort <- c("GSE106936 EC","GSE160936 SSC","syn12514624 DLFPC", 
            "GSE183068 PFC",
            "GSE178265 SN",
            "GSE174367 PFC", "GSE167494 DLPFC", "GSE148822 OTC", 
            "GSE148822 OC",
            "GSE174332 MC","GSE137444 TC")

dplyr::distinct(l, ID, .keep_all = T) %>% dplyr::mutate(logFC = log2(expr_PAR_LOY / expr_PAR) ) -> tmp

tmp %>% dplyr::select(ID, logFC) %>% tibble::column_to_rownames("ID") %>% as.matrix() -> tmp
tmp[cohort,] %>% as.matrix() -> tmp

colorRampPalette(colors = c("#196889","#FFFFFF","#DE2126")) -> k
breaksList = seq(-0.75, 0.75, by = 0.01)
pheatmap::pheatmap(mat = tmp, breaks = breaksList, angle_col = 45, show_rownames = F, display_numbers = T,
         cluster_cols = F, na_col = "grey50", fontsize_number = c(11),
         cluster_rows = F, border_color = "grey40", legend = F,
         color = k(length(breaksList))
)
pheatmap::pheatmap(mat = tmp, breaks = breaksList, angle_col = 45, show_rownames = F, display_numbers = F,
                   cluster_cols = F, na_col = "grey50", fontsize_number = c(11),
                   cluster_rows = F, border_color = "grey40", legend = F,
                   color = k(length(breaksList))
)

dplyr::distinct(l, ID, .keep_all = T) %>% dplyr::mutate(logFC = log2(expr_PAR_LOY / expr_PAR) ) -> tmp

tmp %>% dplyr::select(ID, logFC) %>% tibble::column_to_rownames("ID") %>% as.matrix() -> tmp
tmp[cohort,] %>% as.matrix() -> tmp

tmp[,1] <- 0

colorRampPalette(colors = c("#196889","#FFFFFF","#DE2126")) -> k
breaksList = seq(-0.75, 0.75, by = 0.01)
pheatmap(mat = tmp, breaks = breaksList, angle_col = 45, show_rownames = F, display_numbers = F,
         cluster_cols = F, na_col = "grey50", fontsize_number = c(9),
         cluster_rows = F, border_color = "grey40", legend = F,
         color = k(length(breaksList))
)







###### LOY DE for each subject


#### SMITH AM
readRDS(file = "e:/LOY/scLOY/processed_seurat/MERGED_COHORT/GSE160936_cleaned_filtered.RDS") -> o
o@meta.data$`brain region (somatosensory or entorhinal cortex):ch1` -> o@meta.data$brain_region
o@meta.data$donor_organism.biomaterial_core.biomaterial_id <- o@meta.data$title
gsub(x = o@meta.data$donor_organism.biomaterial_core.biomaterial_id , pattern = " SSC", replacement = "") -> o@meta.data$donor_organism.biomaterial_core.biomaterial_id
gsub(x = o@meta.data$donor_organism.biomaterial_core.biomaterial_id , pattern = " EC", replacement = "") -> o@meta.data$donor_organism.biomaterial_core.biomaterial_id

lapply(o@meta.data$donor %>% unique(), function(x) {

  message(x)
  DE1 <- run_DE_multi(s_obj = subset(o, donor == x),
                    sample = o@meta.data$donor_organism.biomaterial_core.biomaterial_id %>% unique(),
                    latent = c("nCount_RNA","nFeature_RNA","percent.rb","specific_tissue"),
                    cell_type = c("Microglia"), logfc.threshold = 0, min.pct = 0.05)
  DE1$ID <- x

  annoGR2DF(genes(edb)) -> GL
  dplyr::select(GL, gene_name, chr) -> dat

  dat[dat$chr %in% c(1:22,"X","Y","MT"),] -> dat
  dplyr::distinct(dat, gene_name, .keep_all = T) -> dat

  dplyr::left_join(x = DE1, y = dat, by = c("gene"="gene_name")) -> DE1

  DE1$chr %>% as.character() -> DE1$chr
  DE1[DE1$gene %in% PAR_scLOY_genes$gene_name, ]$chr <- "PAR"


  DE1$label <- "Autosomal"
  DE1[DE1$chr %in% as.character(c(1:22)), ]$label <- "Autosomal"
  DE1[which(DE1$chr=="PAR"), ]$label <- "PAR"
  DE1[which(DE1$chr=="Y"), ]$label <- "MSY"
  DE1[which(DE1$chr=="MT"), ]$label <- "MT"
  DE1[is.na(DE1$chr), ]$label <- "Autosomal"

  return(DE1)
}) -> l

saveRDS(data.table::rbindlist(l),
        file = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/MANUSCRIPT/MAIN/Fig3/SMITH_AM_DONOR_DE_LOY_MICROGLIA.RDS")

##### OLAH
readRDS(file = "e:/LOY/scLOY/processed_seurat/MERGED_COHORT/syn12514624_cleaned_filtered.RDS") -> o

lapply(o@meta.data$library_id %>% unique(), function(x) {

  message(x)
  DE1 <- run_DE_multi(s_obj = subset(o, library_id == x),
                      sample = o@meta.data$donor_organism.biomaterial_core.biomaterial_id %>% unique(),
                      latent = c("nCount_RNA","nFeature_RNA","percent.rb"),
                      cell_type = c("Microglia"), nUMI = 1000, nFeature = 800, logfc.threshold = 0, min.pct = 0.05)
  DE1$ID <- x

  annoGR2DF(genes(edb)) -> GL
  dplyr::select(GL, gene_name, chr) -> dat

  dat[dat$chr %in% c(1:22,"X","Y","MT"),] -> dat
  dplyr::distinct(dat, gene_name, .keep_all = T) -> dat

  dplyr::left_join(x = DE1, y = dat, by = c("gene"="gene_name")) -> DE1

  DE1$chr %>% as.character() -> DE1$chr
  DE1[DE1$gene %in% PAR_scLOY_genes$gene_name, ]$chr <- "PAR"


  DE1$label <- "Autosomal"
  DE1[DE1$chr %in% as.character(c(1:22)), ]$label <- "Autosomal"
  DE1[which(DE1$chr=="PAR"), ]$label <- "PAR"
  DE1[which(DE1$chr=="Y"), ]$label <- "MSY"
  DE1[which(DE1$chr=="MT"), ]$label <- "MT"
  DE1[is.na(DE1$chr), ]$label <- "Autosomal"

  return(DE1)
}) -> l

saveRDS(data.table::rbindlist(l),
        file = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/MANUSCRIPT/MAIN/Fig3/OLAH_DONOR_DE_LOY_MICROGLIA.RDS")

##### GSE178265
readRDS(file = "e:/LOY/scLOY/processed_seurat/MERGED_COHORT/GSE178265_cleaned_filtered.RDS") -> o

lapply(c("pPDsHSrSNxi1963","pPDsHSrSNxi4775","pPDsHSrSNxi3887"), function(x) {

  message(x)
  DE1 <- run_DE_multi(s_obj = subset(o, donor_organism.biomaterial_core.biomaterial_id == x),
                      sample = o@meta.data$donor_organism.biomaterial_core.biomaterial_id %>% unique(),
                      latent = c("nCount_RNA","nFeature_RNA","percent.rb"),
                      cell_type = c("Microglia"), nUMI = 1000, nFeature = 800, logfc.threshold = 0, min.pct = 0.05)
  DE1$ID <- x

  annoGR2DF(genes(edb)) -> GL
  dplyr::select(GL, gene_name, chr) -> dat

  dat[dat$chr %in% c(1:22,"X","Y","MT"),] -> dat
  dplyr::distinct(dat, gene_name, .keep_all = T) -> dat

  dplyr::left_join(x = DE1, y = dat, by = c("gene"="gene_name")) -> DE1

  DE1$chr %>% as.character() -> DE1$chr
  DE1[DE1$gene %in% PAR_scLOY_genes$gene_name, ]$chr <- "PAR"


  DE1$label <- "Autosomal"
  DE1[DE1$chr %in% as.character(c(1:22)), ]$label <- "Autosomal"
  DE1[which(DE1$chr=="PAR"), ]$label <- "PAR"
  DE1[which(DE1$chr=="Y"), ]$label <- "MSY"
  DE1[which(DE1$chr=="MT"), ]$label <- "MT"
  DE1[is.na(DE1$chr), ]$label <- "Autosomal"

  return(DE1)
}) -> l

saveRDS(data.table::rbindlist(l),
        file = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/MANUSCRIPT/MAIN/Fig3/GSE178265_DONOR_DE_LOY_MICROGLIA.RDS")

##### GSE148822
readRDS(file = "e:/LOY/scLOY/processed_seurat/MERGED_COHORT/GSE148822_cleaned_filtered.RDS") -> o

o@meta.data$donor -> o@meta.data$donor_organism.biomaterial_core.biomaterial_id

lapply(o@meta.data$donor %>% unique(), function(x) {

  message(x)
  DE1 <- run_DE_multi(s_obj = subset(o, donor_organism.biomaterial_core.biomaterial_id == x),
                      sample = o@meta.data$donor_organism.biomaterial_core.biomaterial_id %>% unique(),
                      latent = c("nCount_RNA","nFeature_RNA","percent.rb","specific_tissue"),
                      cell_type = c("Microglia"), nUMI = 1000, nFeature = 800, logfc.threshold = 0, min.pct = 0.05)
  DE1$ID <- x

  annoGR2DF(genes(edb)) -> GL
  dplyr::select(GL, gene_name, chr) -> dat

  dat[dat$chr %in% c(1:22,"X","Y","MT"),] -> dat
  dplyr::distinct(dat, gene_name, .keep_all = T) -> dat

  dplyr::left_join(x = DE1, y = dat, by = c("gene"="gene_name")) -> DE1

  DE1$chr %>% as.character() -> DE1$chr
  DE1[DE1$gene %in% PAR_scLOY_genes$gene_name, ]$chr <- "PAR"


  DE1$label <- "Autosomal"
  DE1[DE1$chr %in% as.character(c(1:22)), ]$label <- "Autosomal"
  DE1[which(DE1$chr=="PAR"), ]$label <- "PAR"
  DE1[which(DE1$chr=="Y"), ]$label <- "MSY"
  DE1[which(DE1$chr=="MT"), ]$label <- "MT"
  DE1[is.na(DE1$chr), ]$label <- "Autosomal"

  return(DE1)
}) -> l

saveRDS(data.table::rbindlist(l),
        file = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/MANUSCRIPT/MAIN/Fig3/GSE148822_DONOR_DE_LOY_MICROGLIA.RDS")


##### GSE183068 PFC

readRDS(file = "e:/LOY/scLOY/processed_seurat/MERGED_COHORT/GSE183068_cleaned_filtered.RDS") -> o
readRDS(file = "E:/LOY/scLOY/processed_seurat/BRAIN/GSE183068_AD_PFC_introns_included.RDS") -> o

o@meta.data$donor <- o@meta.data$orig.ident
o@meta.data$donor_organism.biomaterial_core.biomaterial_id <- o@meta.data$orig.ident

subset(o, donor_organism.sex == "male") -> o
NormalizeData(o) -> o
o@meta.data$cell_ident_MV <- o@meta.data$TYPE

lapply(c("GSM5550475","GSM5550481","GSM5550501","GSM5550451","GSM5550491","GSM5550495","GSM5550503"), function(x) {
 
  message(x)
  
  DE1 <- run_DE_multi(s_obj = o,
                      sample = c(x),
                      latent = c("nCount_RNA","nFeature_RNA","percent.rb"),
                      cell_type = c("Microglia"), nUMI = 1000, nFeature = 800, logfc.threshold = 0, min.pct = 0.05)
  DE1$ID <- x
 
  annoGR2DF(genes(edb)) -> GL
  dplyr::select(GL, gene_name, chr) -> dat
  
  dat[dat$chr %in% c(1:22,"X","Y","MT"),] -> dat
  dplyr::distinct(dat, gene_name, .keep_all = T) -> dat
  
  dplyr::left_join(x = DE1, y = dat, by = c("gene"="gene_name")) -> DE1
  
  DE1$chr %>% as.character() -> DE1$chr
  DE1[DE1$gene %in% PAR_scLOY_genes$gene_name, ]$chr <- "PAR"
  
  
  DE1$label <- "Autosomal"
  DE1[DE1$chr %in% as.character(c(1:22)), ]$label <- "Autosomal"
  DE1[which(DE1$chr=="PAR"), ]$label <- "PAR"
  DE1[which(DE1$chr=="Y"), ]$label <- "MSY"
  DE1[which(DE1$chr=="MT"), ]$label <- "MT"
  DE1[is.na(DE1$chr), ]$label <- "Autosomal"
  
  return(DE1)
}) -> l

saveRDS(data.table::rbindlist(l),
        file = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/MANUSCRIPT/MAIN/Fig3/GSE183068_DONOR_DE_LOY_MICROGLIA.RDS")



##### GSE167490 PFC
readRDS(file = "E:/LOY/scLOY/processed_seurat/BRAIN/GSE167492_AD_DLPFC_introns_included.RDS") -> o
readRDS(file = "E:/LOY/scLOY/processed_seurat/BRAIN/GSE167490_AD_DLPFC_introns_included.RDS") -> oo
merge(o,oo) -> o

subset(o, donor_organism.sex == "male") -> o

lapply(c("GSM5106120","GSM5106119","GSM5106101","GSM5106122"), function(x) {
  
  message(x)
  DE1 <- run_DE_multi(s_obj = o,
                      sample = c(x),
                      latent = c("nCount_RNA","nFeature_RNA","percent.rb"),
                      cell_type = c("Microglia"), nUMI = 1000, nFeature = 800, logfc.threshold = 0, min.pct = 0.05)
  DE1$ID <- x
  
  annoGR2DF(genes(edb)) -> GL
  dplyr::select(GL, gene_name, chr) -> dat
  
  dat[dat$chr %in% c(1:22,"X","Y","MT"),] -> dat
  dplyr::distinct(dat, gene_name, .keep_all = T) -> dat
  
  dplyr::left_join(x = DE1, y = dat, by = c("gene"="gene_name")) -> DE1
  
  DE1$chr %>% as.character() -> DE1$chr
  DE1[DE1$gene %in% PAR_scLOY_genes$gene_name, ]$chr <- "PAR"
  
  
  DE1$label <- "Autosomal"
  DE1[DE1$chr %in% as.character(c(1:22)), ]$label <- "Autosomal"
  DE1[which(DE1$chr=="PAR"), ]$label <- "PAR"
  DE1[which(DE1$chr=="Y"), ]$label <- "MSY"
  DE1[which(DE1$chr=="MT"), ]$label <- "MT"
  DE1[is.na(DE1$chr), ]$label <- "Autosomal"
  
  return(DE1)
}) -> l

saveRDS(data.table::rbindlist(l),
        file = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/MANUSCRIPT/MAIN/Fig3/DONOR_DE/GSE167494_DONOR_DE_LOY_MICROGLIA.RDS")






### HEATMAP FOR DONORS

# load
list.files("h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/MANUSCRIPT/MAIN/Fig3/DONOR_DE",
           full.names = T, pattern = "RDS") -> files

#files[-1] -> files


lapply(files, function(x){
  readRDS(x) -> tmp
  return(tmp)
}) -> l
data.table::rbindlist(l) -> l

genes <- c("UTY","RPS4Y1","NLGN4Y",
           "DHRSX","CSF2RA","CD99","IL3RA","AKAP17A","SLC25A6","CD99P1",
           "PARD3B","S100Z","NAV3","ST6GALNAC3","TMEM176B","PRKCA",
           "ROBO1","RBFOX2","APOE","TPRG1",
           "CD163","DOCK1","DPYD") %>% unique


#genes <- c("UTY","RPS4Y1","NLGN4Y",
#           "DHRSX","CSF2RA","CD99",
#           "IL3RA","AKAP17A",
#           "SLC25A6","CD99P1","LINC02232",
#           "LNCAROD","NAV3","KCNIP1","CYFIP1","CD226",
#           "ROBO1","APOC1","APOE","CD163",,"DPYD",
#            "RBFOX2","DOCK1","CERNA2")
dplyr::distinct(l, ID, .keep_all = T) -> dist

dist[order(dist$LOY_cells),]$ID %>% rev() -> cohort


l[l$gene %in% genes,] -> l.


l. %>% dplyr::select(p_val_adj, p_val, avg_log2FC, gene, ID) -> mat

reshape2::dcast(data = mat, formula = ID ~ gene, value.var = "avg_log2FC") -> mat
mat$ID -> rownames(mat)
mat[,-1] -> mat



cohort <- c("A163/17","GSM5550503",
            "Microglia_MO_MCI3","A277/12","pPDsHSrSNxi3887","A096/14","pPDsHSrSNxi4775",
            "GSM5550475","GSM5550481","GSM5550491","GSM5106120",
            "A402/14","pPDsHSrSNxi1963","A053/11","A319/11","Donor11","A127/11","A297/16",
            "Donor8","Donor6","GSM5106119","GSM5106101","GSM5550495",
            "Microglia_MO_AD71",
            "Donor1","Donor18","Donor4","Donor2")
mat <- mat[cohort,genes]



# make sig * matrix
l. %>% dplyr::select(p_val_adj, p_val, avg_log2FC, gene, ID, FDR) -> sig

reshape2::dcast(data = sig, formula = ID ~ gene, value.var = "FDR") -> sig
sig$ID -> rownames(sig)
sig[,-1] -> sig

sig <- sig[cohort,
           genes]



ifelse(test = sig < 0.001, yes = "***", no =
         ifelse(sig < 0.01, yes = "**", no =
                  ifelse(sig < 0.05, yes = "*", no = ""))) -> sig

sig[is.na(sig)] <- ""

colorRampPalette(colors = c("#196889","#FFFFFF","#DE2126")) -> k
breaksList = seq(-0.5, 0.5, by = 0.01)
pheatmap(mat = mat, breaks = breaksList, display_numbers = sig, angle_col = 45,
         cluster_cols = F, na_col = "grey50", gaps_col = c(3,10,16), fontsize_number = c(12),
         cluster_rows = F, border_color = "grey50",
         color = k(length(breaksList))
)


#### SIDE DATA


# load
list.files("h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/MANUSCRIPT/MAIN/Fig3/DONOR_DE",
           pattern = "RDS", full.names = T) -> files

# files[-1] -> files


lapply(files, function(x){
  readRDS(x) -> tmp
  return(tmp)
}) -> l
data.table::rbindlist(l) -> l


genes <- c("UTY","RPS4Y1","NLGN4Y",
           "DHRSX","CSF2RA","CD99","IL3RA","AKAP17A","SLC25A6","CD99P1",
           "PARD3B","S100Z","NAV3","ST6GALNAC3","TMEM176B","PRKCA",
           "ROBO1","RBFOX2","APOE","TPRG1",
           "CD163","DOCK1","DPYD") %>% unique
dplyr::distinct(l, ID, .keep_all = T) -> dist

cohort <- c("A163/17","GSM5550503",
            "Microglia_MO_MCI3","A277/12","pPDsHSrSNxi3887","A096/14","pPDsHSrSNxi4775",
            "GSM5550475","GSM5550481","GSM5550491","GSM5106120",
            "A402/14","pPDsHSrSNxi1963","A053/11","A319/11","Donor11","A127/11","A297/16",
            "Donor8","Donor6","GSM5106119","GSM5106101","GSM5550495",
            "Microglia_MO_AD71",
            "Donor1","Donor18","Donor4","Donor2")

dplyr::distinct(l, ID, .keep_all = T) %>% dplyr::mutate(logFC = log2(expr_PAR_LOY / expr_PAR) ) -> tmp

tmp %>% dplyr::select(ID, logFC) %>% tibble::column_to_rownames("ID") %>% as.matrix() -> tmp
tmp[cohort,] %>% as.matrix() -> tmp

colorRampPalette(colors = c("#196889","#FFFFFF","#DE2126")) -> k
breaksList = seq(-0.5, 0.5, by = 0.01)
pheatmap(mat = tmp, breaks = breaksList, angle_col = 45, show_rownames = F, display_numbers = F,
         cluster_cols = F, na_col = "grey50", fontsize_number = c(11),
         cluster_rows = F, border_color = "grey40", legend = F,
         color = k(length(breaksList))
)

dplyr::distinct(l, ID, .keep_all = T) %>% dplyr::mutate(logFC = log2(expr_PAR_LOY / expr_PAR) ) -> tmp

tmp %>% dplyr::select(ID, logFC) %>% tibble::column_to_rownames("ID") %>% as.matrix() -> tmp
tmp[cohort,] %>% as.matrix() -> tmp

tmp[,1] <- 0

colorRampPalette(colors = c("#196889","#FFFFFF","#DE2126")) -> k
breaksList = seq(-0.75, 0.75, by = 0.01)
pheatmap(mat = tmp, breaks = breaksList, angle_col = 45, show_rownames = F, display_numbers = F,
         cluster_cols = F, na_col = "grey50", fontsize_number = c(9),
         cluster_rows = F, border_color = "grey40", legend = F,
         color = k(length(breaksList))
)


#### DE GENES PER DONOR

#### DE gene table (number of DE genes)


list.files(path = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/MANUSCRIPT/MAIN/Fig3/DONOR_DE", pattern = "RDS", full.names = T) -> files
files[-1] -> files
lapply(files, function(x){

  readRDS(file = x) -> tmp
  rownames(tmp) -> tmp$gene
  return(tmp)
})-> l
do.call("rbind",l) -> l

l[l$p_val_adj < 0.1, ] -> l

# total sig genes by subject
l[l$p_val_adj < 0.1,]$ID %>% table() %>% unclass() %>% as.data.frame() -> total_sig
names(total_sig) <- "TOTAL_SIG"
rownames(total_sig) -> total_sig$cohort
total_sig[,c(2,1)] -> total_sig

l[l$p_val_adj < 0.1 & l$avg_log2FC > 0 & l$label=="Autosomal",]$ID %>% table() %>% unclass() %>% as.data.frame() -> up
names(up) <- "UP_SIG"
rownames(up) -> up$cohort
left_join(x = total_sig, y = up, by = c("cohort")) -> total_sig

l[l$p_val_adj < 0.1 & l$avg_log2FC < 0 & l$label == "Autosomal",]$ID %>% table() %>% unclass() %>% as.data.frame() -> down
names(down) <- "DOWN_SIG"
rownames(down) -> down$cohort
left_join(x = total_sig, y = down, by = c("cohort")) -> total_sig


l[l$p_val_adj < 0.1 & l$avg_log2FC < 0 &
    l$label=="MSY",]$ID %>% table() %>% unclass() %>% as.data.frame() -> Y
names(Y) <- "DOWN_Y"
rownames(Y) -> Y$cohort
left_join(x = total_sig, y = Y, by = c("cohort")) -> total_sig


l[l$p_val_adj < 0.1 & l$avg_log2FC < 0 &
    l$label=="PAR",]$ID %>% table() %>% unclass() %>% as.data.frame() -> PAR
names(PAR) <- "DOWN_PAR"
rownames(PAR) -> PAR$cohort
left_join(x = total_sig, y = PAR, by = c("cohort")) -> total_sig

l[l$p_val_adj < 0.1 & l$avg_log2FC > 0 &
    l$label == "MT",]$ID %>% table() %>% unclass() %>% as.data.frame() -> MT
names(MT) <- "UP_MT"
rownames(MT) -> MT$cohort
left_join(x = total_sig, y = MT, by = c("cohort")) -> total_sig


data.table::fwrite(x = total_sig, file = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/MANUSCRIPT/TABLES/microglia_DE_donor_.csv")






### BARPLOTS

list.files(path = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/MANUSCRIPT/MAIN/Fig3/DONOR_DE", pattern = "RDS", full.names = T) -> files
files[-1] -> files
lapply(files, function(x){

  readRDS(file = x) -> tmp
  return(tmp)
})-> l
do.call("rbind",l) -> l


output <- "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/MANUSCRIPT/MAIN/Fig3/DONOR_DE"
## 163
l[l$ID == "A163/17", ] -> l.
l.[order(l.$avg_log2FC),] -> l.
ggbarplot(data = l.[l.$FDR < 0.05,], x = "gene", y = "avg_log2FC", fill = "label", color = NA) %>%
  ggpar(orientation = "horizontal", ylab = "", xlab = "", legend = "none", font.xtickslab = c(10, "bold"),
        palette = c("grey",get_palette(palette = "aaas", 10)[c(2,4,3)])) +
  rremove("y.text") + rremove("y.axis") + rremove("y.ticks") + geom_hline(yintercept = 0, size = 0.4) -> p0

ggsave(filename = paste0(output,"/A163_17_big.pdf"), units = "in", width = 4, height = 5,
       plot = p0, dpi = "retina", device = "pdf")

# MCI3
l[l$ID == "Microglia_MO_MCI3", ] -> l.
l.[order(l.$avg_log2FC),] -> l.
ggbarplot(data = l.[l.$FDR < 0.05,], x = "gene", y = "avg_log2FC", fill = "label", color = NA) %>%
  ggpar(orientation = "horizontal", ylab = "", xlab = "", legend = "none", font.xtickslab = c(10, "bold"),
        palette = c("grey",get_palette(palette = "aaas", 10)[c(2,3)])) +
  rremove("y.text") + rremove("y.axis") + rremove("y.ticks") + geom_hline(yintercept = 0, size = 0.4) + ylim(-3,1) -> p1

ggsave(filename = paste0(output,"/MCI3_small.pdf"), units = "in", width = 3, height = 3,
       plot = p1, dpi = "retina", device = "pdf")

## A277/12
l[l$ID == "A277/12", ] -> l.
l.[order(l.$avg_log2FC),] -> l.
ggbarplot(data = l.[l.$FDR < 0.05,], x = "gene", y = "avg_log2FC", fill = "label", color = NA) %>%
  ggpar(orientation = "horizontal", ylab = "", xlab = "", legend = "none", font.xtickslab = c(10, "bold"),
        palette = c("grey",get_palette(palette = "aaas", 10)[c(2,3)])) +
  rremove("y.text") + rremove("y.axis") + rremove("y.ticks") + geom_hline(yintercept = 0, size = 0.4) + ylim(-3.2,1) -> p2

ggsave(filename = paste0(output,"/A277_12_small.pdf"), units = "in", width = 3, height = 3,
       plot = p2, dpi = "retina", device = "pdf")

## A096/14
l[l$ID == "A096/14", ] -> l.
l.[order(l.$avg_log2FC),] -> l.
ggbarplot(data = l.[l.$FDR < 0.05,], x = "gene", y = "avg_log2FC", fill = "label", color = NA) %>%
  ggpar(orientation = "horizontal", ylab = "", xlab = "", legend = "none", font.xtickslab = c(10, "bold"),
        palette = c("grey",get_palette(palette = "aaas", 10)[c(2,3)])) +
  rremove("y.text") + rremove("y.axis") + rremove("y.ticks") + geom_hline(yintercept = 0, size = 0.4) + ylim(-3,1) -> p3

## AD71
l[l$ID == "Microglia_MO_AD71", ] -> l.
l.[order(l.$avg_log2FC),] -> l.
ggbarplot(data = l.[l.$FDR < 0.05,], x = "gene", y = "avg_log2FC", fill = "label", color = NA) %>%
  ggpar(orientation = "horizontal", ylab = "", xlab = "", legend = "none", font.xtickslab = c(10, "bold"),
        palette = c("grey",get_palette(palette = "aaas", 10)[c(2,3)])) +
  rremove("y.text") + rremove("y.axis") + rremove("y.ticks") + geom_hline(yintercept = 0, size = 0.4) + ylim(-3,1) -> p4

## pPDsHSrSNxi4775
l[l$ID == "pPDsHSrSNxi4775", ] -> l.
l.[order(l.$avg_log2FC),] -> l.
ggbarplot(data = l.[l.$FDR < 0.05,], x = "gene", y = "avg_log2FC", fill = "label", color = NA) %>%
  ggpar(orientation = "horizontal", ylab = "", xlab = "", legend = "none", font.xtickslab = c(10, "bold"),
        palette = c("grey",get_palette(palette = "aaas", 10)[c(2,3)])) +
  rremove("y.text") + rremove("y.axis") + rremove("y.ticks") + geom_hline(yintercept = 0, size = 0.4) + ylim(-3,1) -> p5

## pPDsHSrSNxi3887
l[l$ID == "pPDsHSrSNxi3887", ] -> l.
l.[order(l.$avg_log2FC),] -> l.
ggbarplot(data = l.[l.$FDR < 0.05,], x = "gene", y = "avg_log2FC", fill = "label", color = NA) %>%
  ggpar(orientation = "horizontal", ylab = "", xlab = "", legend = "none", font.xtickslab = c(10, "bold"),
        palette = c("grey",get_palette(palette = "aaas", 10)[c(2,3)])) +
  rremove("y.text") + rremove("y.axis") + rremove("y.ticks") + geom_hline(yintercept = 0, size = 0.4) + ylim(-3,1) -> p6

## A319
l[l$ID == "A319/11", ] -> l.
l.[order(l.$avg_log2FC),] -> l.
ggbarplot(data = l.[l.$FDR < 0.05,], x = "gene", y = "avg_log2FC", fill = "label", color = NA) %>%
  ggpar(orientation = "horizontal", ylab = "", xlab = "", legend = "none", font.xtickslab = c(10, "bold"),
        palette = c("grey",get_palette(palette = "aaas", 10)[c(2,3)])) +
  rremove("y.text") + rremove("y.axis") + rremove("y.ticks") + geom_hline(yintercept = 0, size = 0.4)+ ylim(-3,1) -> p7

## A297/16
l[l$ID == "A297/16", ] -> l.
l.[order(l.$avg_log2FC),] -> l.
ggbarplot(data = l.[l.$FDR < 0.05,], x = "gene", y = "avg_log2FC", fill = "label", color = NA) %>%
  ggpar(orientation = "horizontal", ylab = "", xlab = "", legend = "none", font.xtickslab = c(10, "bold"),
        palette = c("grey",get_palette(palette = "aaas", 10)[c(2,3)])) +
  rremove("y.text") + rremove("y.axis") + rremove("y.ticks") + geom_hline(yintercept = 0, size = 0.4)+ ylim(-3,1) -> p8


## A127/11
l[l$ID == "A127/11", ] -> l.
l.[order(l.$avg_log2FC),] -> l.
ggbarplot(data = l.[l.$FDR < 0.05,], x = "gene", y = "avg_log2FC", fill = "label", color = "label") %>%
  ggpar(orientation = "horizontal", ylab = "", xlab = "", legend = "none", font.xtickslab = c(10, "bold"),
        palette = c("grey",get_palette(palette = "aaas", 10)[c(2,3)])) +
  rremove("y.text") + rremove("y.axis") + rremove("y.ticks") + geom_hline(yintercept = 0, size = 0.4)+ ylim(-3,1) -> p9

## A053/11
l[l$ID == "A053/11", ] -> l.
l.[order(l.$avg_log2FC),] -> l.
ggbarplot(data = l.[l.$FDR < 0.05,], x = "gene", y = "avg_log2FC", fill = "label", color = "label") %>%
  ggpar(orientation = "horizontal", ylab = "", xlab = "", legend = "none", font.xtickslab = c(10, "bold"),
        palette = c("grey",get_palette(palette = "aaas", 10)[c(2,3)])) +
  rremove("y.text") + rremove("y.axis") + rremove("y.ticks") + geom_hline(yintercept = 0, size = 0.4)+ ylim(-3,1) -> p10


## A402/14
l[l$ID == "A402/14", ] -> l.
l.[order(l.$avg_log2FC),] -> l.
ggbarplot(data = l.[l.$FDR < 0.05,], x = "gene", y = "avg_log2FC", fill = "label", color = "label") %>%
  ggpar(orientation = "horizontal", ylab = "", xlab = "", legend = "none", font.xtickslab = c(10, "bold"),
        palette = c("grey",get_palette(palette = "aaas", 10)[c(2,3)])) +
  rremove("y.text") + rremove("y.axis") + rremove("y.ticks") + geom_hline(yintercept = 0, size = 0.4)+ ylim(-3,1) -> p11

ggarrange(plotlist = list(p1,p2,p3,p5,p6, p9)) -> out

ggsave(filename = paste0(output,"/ggarrange_LOY_DE.pdf"), units = "in", width = 8, height = 6,
       plot = out, dpi = "retina", device = "pdf")




###################### HYPER R  (GENE SET ENRICHMENT) ###########


## ADD PAR REGION HYPER

qusage::read.gmt("h:/LOY/scLOY/scLOY_NOV18/GMT/c1.all.v7.4.symbols_CUSTOM_PAR.gmt") -> gs
names(gs)[265] <- "chrXp22"
names(gs)[266] <- "chrXq11"

### order the list based on genomic coordinates

# pull first gene from each list element, find its location

annoGR2DF(genes(edb)) -> GL

lapply(names(gs),function(x){


  message(x)
  dplyr::select(GL, gene_name, chr, start, end) -> dat
  dat[dat$chr %in% c(1:22,"X","Y","MT"),] -> dat
  dplyr::distinct(dat, gene_name, .keep_all = T) -> dat


  if(dat[dat$gene_name==gs[[x]][1],] %>% nrow() > 0){dat[dat$gene_name==gs[[x]][1],]-> gene} else{

    dat[dat$gene_name==gs[[x]][5],] -> gene
  } el


  gene$band <- x

  return(gene)

}) -> order_chr



list.files(path = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/MANUSCRIPT/MAIN/Fig3/DONOR_DE/", pattern = "RDS", full.names = T) -> files
files[-1] -> files
lapply(files, function(x){

  readRDS(file = x) -> tmp
  return(tmp)
})-> l
do.call("rbind",l) -> l


output <- "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/MANUSCRIPT/MAIN/Fig3/DONOR_DE"
## 163
l[l$ID == "A163/17", ] -> l.
l.[order(l.$avg_log2FC),] -> l.
l.[which(l.$FDR < 0.05 & l.$label!="MSY"),]$gene -> sig

hyp_obj <- hypeR::hypeR(signature = sig, genesets = gs, background = l[l$ID == "A163/17", ]$gene)
hyp_show(hyp_obj)
hyp_dots(hyp_obj) -> p; ggpar(p,ggtheme = theme_bw(),font.tickslab = c(8,"bold"),
                                xlab = "c1.all.v7.4.symbols (GSEA)") +
  scale_color_viridis_c(option = "A", direction = -1, end = 0.7)
#hyp_to_table(hyp_obj, file_path=paste0("hypeR.txt")

ggsave(filename = paste0(output,"/ggarrange_LOY_DE.pdf"), units = "in", width = 8, height = 6,
       plot = out, dpi = "retina", device = "pdf")



list.files(path = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/MANUSCRIPT/MAIN/Fig3/DONOR_DE", pattern = "RDS", full.names = T) -> files

lapply(files, function(x){

  readRDS(file = x) -> tmp
  return(tmp)
})-> l
do.call("rbind",l) -> l

lapply(l$ID %>% unique(), function(x){
  message(x)
  l[l$ID == x, ] -> l.
  #l.[order(l.$avg_log2FC),] -> l.
  l.[which(l.$FDR < 0.1 & l.$label!="MSY"),]$gene -> sig

  if(length(sig)==0){return(NULL)}

  hyp_obj <- hypeR::hypeR(signature = sig, genesets = gs, background = l[l$ID == x, ]$gene)
  hyp_to_table(hyp_obj, file_path=paste0("hypeR.txt"))
  hyp_obj$as.data.frame() -> dat
  dat$subject <- x
  return(dat)
}) -> out

names(out) <- l$ID %>% unique()
plyr::compact(out) -> out

data.table::rbindlist(out) -> out

#out[out$pval < 0.1,] -> out.
dplyr::mutate(out, logFDR = -log10(fdr)) -> out

data.table::fwrite(x = out, file = "a:/Dropbox/LOY_microglia_paper/GR_RESUB/MANUSCRIPT/TABLES/DONOR_REGION_ENRICHMENT.csv")


ggdotplot(out, x = "label", y = "logFDR", fill = "overlap", binwidth = 0.2,
          size = 1) %>% ggpar(x.text.angle = 45, font.xtickslab = c(9)) +
          geom_hline(yintercept = 6.64, linetype = 'dotted')





#### DE gene table (number of DE genes) FDR


list.files(path = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/MANUSCRIPT/MAIN/Fig3/DONOR_DE", pattern = "*.RDS", full.names = T) -> files

lapply(files, function(x){

  readRDS(file = x) -> tmp
  rownames(tmp) -> tmp$gene
  return(tmp)
})-> l
do.call("rbind",l) -> l

l[l$FDR < 0.05, ] -> l

# total sig genes by subject
l[l$FDR < 0.05,]$ID %>% table() %>% unclass() %>% as.data.frame() -> total_sig
names(total_sig) <- "TOTAL_SIG"
rownames(total_sig) -> total_sig$cohort
total_sig[,c(2,1)] -> total_sig

l[l$FDR < 0.05 & l$avg_log2FC > 0 & l$label=="Autosomal",]$ID %>% table() %>% unclass() %>% as.data.frame() -> up
names(up) <- "UP_SIG"
rownames(up) -> up$cohort
left_join(x = total_sig, y = up, by = c("cohort")) -> total_sig

l[l$FDR < 0.05 & l$avg_log2FC < 0 & l$label == "Autosomal",]$ID %>% table() %>% unclass() %>% as.data.frame() -> down
names(down) <- "DOWN_SIG"
rownames(down) -> down$cohort
left_join(x = total_sig, y = down, by = c("cohort")) -> total_sig


l[l$FDR < 0.05 & l$avg_log2FC < 0 &
    l$label=="MSY",]$ID %>% table() %>% unclass() %>% as.data.frame() -> Y
names(Y) <- "DOWN_Y"
rownames(Y) -> Y$cohort
left_join(x = total_sig, y = Y, by = c("cohort")) -> total_sig


l[l$FDR < 0.05 & l$avg_log2FC < 0 &
    l$label=="PAR",]$ID %>% table() %>% unclass() %>% as.data.frame() -> PAR
names(PAR) <- "DOWN_PAR"
rownames(PAR) -> PAR$cohort
left_join(x = total_sig, y = PAR, by = c("cohort")) -> total_sig

l[l$FDR < 0.05 & l$avg_log2FC > 0 &
    l$label == "MT",]$ID %>% table() %>% unclass() %>% as.data.frame() -> MT
names(MT) <- "UP_MT"
rownames(MT) -> MT$cohort
left_join(x = total_sig, y = MT, by = c("cohort")) -> total_sig


data.table::fwrite(x = total_sig, file = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/MANUSCRIPT/TABLES/microglia_DE_donor_FDR.csv")




##### ENRICHMENT of LOY DE genes
data.table::fread(file = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/MANUSCRIPT/TABLES/dumanki_LATE_genes.csv", header = T) -> LATE
qusage::read.gmt("h:/LOY/scLOY/scLOY_NOV18/GMT/c2.cp.v7.4.symbols.gmt") -> gs1
qusage::read.gmt("h:/LOY/scLOY/scLOY_NOV18/GMT/c5.go.v7.4.symbols.gmt") -> gs2

rlist::list.append(gs1, gs2) -> gs

list.files(path = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/MANUSCRIPT/MAIN/Fig3/DONOR_DE", pattern = "*.RDS", full.names = T) -> files

lapply(files, function(x){

  readRDS(file = x) -> tmp

  return(tmp)
})-> l
do.call("rbind",l) -> l


#LATE[LATE$Location%in%c("Autosomal gene", "PAR gene", "X gene"),]$Gene -> sig
#sig[!(grepl(x = sig, pattern = "[RPSL]"))] -> sig

l[l$FDR < 0.05 & l$label!="MSY",] -> l.
# [l$FDR < 0.05 & l$label!="MSY" & l$label!="MT",] -> l.

l.$gene %>% unique() -> sig

# hyp_obj <- hypeR::hypeR(signature = sig, genesets = gs, background = 15000)
hyp_obj <- hypeR::hypeR(signature = sig, genesets = gs1, background = l[l$label!="MSY" , ]$gene %>% unique())
hyp_dots(hyp_obj) -> p; ggpar(p,ggtheme = theme_bw(),xlab = "Pathways", xlim = c(0,15), orientation = "horizontal") +
  scale_color_viridis_c(option = "A", direction = -1, end = 1, limits = c(0,0.15)) + rremove("grid") -> plot
hyp_obj$as.data.frame() -> dat
data.table::fwrite(x = dat,
                   file = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/MANUSCRIPT/TABLES/ALL_MICROGLIA_LATE_REACTOME_BIOCARTA_HYPER.csv")
ggsave(filename = paste0("h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/MANUSCRIPT/SUPP/GENE_SET/ALL_MICROGLIA_LATE_REACTOME_BIOCARTA_HYPER.pdf"),
       units = "in", width = 10, height = 4,
       plot = plot, dpi = "retina", device = "pdf")


hyp_obj <- hypeR::hypeR(signature = sig, genesets = gs2, background = l[l$label!="MSY" , ]$gene %>% unique())
hyp_obj$as.data.frame() -> dat2

hyp_dots(hyp_obj) -> p; ggpar(p,ggtheme = theme_bw(),xlab = "GO Terms", xlim = c(0.15,0), orientation = "horizontal") +
 scale_color_viridis_c(option = "A", direction = -1, end = 1, limits = c(0,0.15)) + rremove("grid") -> plot
data.table::fwrite(x = dat2, file = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/MANUSCRIPT/TABLES/ALL_MICROGLIA_LATE_GO_HYPER.csv")
ggsave(filename = paste0("h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/MANUSCRIPT/SUPP/GENE_SET/ALL_MICROGLIA_LATE_GO_HYPER.pdf"),
       units = "in", width = 10, height = 4,
       plot = plot, dpi = "retina", device = "pdf")



### HEATMAP FOR DONOR (GENES OF INTEREST IN GENE SET ENRICHMENT HITS)

### HEATMAP FOR DONORS

# load
list.files("h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/MANUSCRIPT/MAIN/Fig3/DONOR_DE",
           pattern = "*.RDS", full.names = T) -> files

lapply(files, function(x){
  readRDS(x) -> tmp
  return(tmp)
}) -> l
data.table::rbindlist(l) -> l


genes <- c("B2M","APOE","HLA-C","HLA-B","CD163","APOC1","HDLBP","IFITM3")
dplyr::distinct(l, ID, .keep_all = T) -> dist

dist[order(dist$LOY_cells),]$ID %>% rev() -> cohort


l[l$gene %in% genes,] -> l.


l. %>% dplyr::select(p_val_adj, p_val, avg_log2FC, gene, ID) -> mat

reshape2::dcast(data = mat, formula = ID ~ gene, value.var = "avg_log2FC") -> mat
mat$ID -> rownames(mat)
mat[,-1] -> mat

mat <- mat[cohort,genes]



# make sig * matrix
l. %>% dplyr::select(p_val_adj, p_val, avg_log2FC, gene, ID, FDR) -> sig

reshape2::dcast(data = sig, formula = ID ~ gene, value.var = "FDR") -> sig
sig$ID -> rownames(sig)
sig[,-1] -> sig

sig <- sig[cohort,
           genes]



ifelse(test = sig < 0.001, yes = "***", no =
         ifelse(sig < 0.01, yes = "**", no =
                  ifelse(sig < 0.05, yes = "*", no = ""))) -> sig

sig[is.na(sig)] <- ""

colorRampPalette(colors = c("#196889","#FFFFFF","#DE2126")) -> k
breaksList = seq(-0.5, 0.5, by = 0.01)
pheatmap(mat = mat, breaks = breaksList, display_numbers = sig, angle_col = 45,
         cluster_cols = F, na_col = "grey50",  fontsize_number = c(12),
         cluster_rows = F, border_color = "grey40",
         color = k(length(breaksList))
)




## HYPER ON FULL COHORT LOY DE ##############

list.files(path = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/FINAL/DE/SEPT_21/", pattern = "*.RDS", full.names = T) -> files
lapply(files, function(x){

  readRDS(file = x) -> tmp
  rownames(tmp) -> tmp$gene
  data.table::rbindlist(tmp) -> tmp
  return(tmp)
})-> l
do.call("rbind",l) -> l


l[l$ID!= "GSE106936 FULL",] -> l

l[l$LOY_cells > 100, ] -> l

annoGR2DF(genes(edb)) -> GL
dplyr::select(GL, gene_name, chr) -> dat

dat[dat$chr %in% c(1:22,"X","Y","MT"),] -> dat
dplyr::distinct(dat, gene_name, .keep_all = T) -> dat

dplyr::left_join(x = l, y = dat, by = c("gene"="gene_name")) -> l

l$chr %>% as.character() -> l$chr
l[l$gene %in% PAR_scLOY_genes$gene_name, ]$chr <- "PAR"


l$label <- "Autosomal"
l[l$chr %in% as.character(c(1:22)), ]$label <- "Autosomal"
l[which(l$chr=="PAR"), ]$label <- "PAR"
l[which(l$chr=="Y"), ]$label <- "MSY"
l[which(l$chr=="MT"), ]$label <- "MT"
l[is.na(l$chr), ]$label <- "Autosomal"

qusage::read.gmt("h:/LOY/scLOY/scLOY_NOV18/GMT/c1.all.v7.4.symbols_CUSTOM_PAR.gmt") -> gs
names(gs)[265] <- "chrXp22"
names(gs)[266] <- "chrXq11"

lapply(l$ID %>% unique(), function(x){
  message(x)
  l[l$ID == x, ] -> l.
  #l.[order(l.$avg_log2FC),] -> l.
  l.[which(l.$FDR < 0.05 & l.$label!="MSY"),]$gene -> sig

  if(length(sig)==0){return(NULL)}

  hyp_obj <- hypeR::hypeR(signature = sig, genesets = gs, background = l[l$ID == x, ]$gene)
  hyp_show(hyp_obj)
  #hyp_dots(hyp_obj) -> p; ggpar(p,ggtheme = theme_bw(),xlab = "c1.all.v7.4.symbols (GSEA)") +
  # scale_color_viridis_c(option = "A", direction = -1, end = 0.7)
  #hyp_to_table(hyp_obj, file_path=paste0("hypeR.txt")

  # ggsave(filename = paste0(output,"/ggarrange_LOY_DE.pdf"), units = "in", width = 8, height = 6,
  #        plot = out, dpi = "retina", device = "pdf")
  #
  #
  hyp_obj$as.data.frame() -> dat
  dat$subject <- x
  return(dat)
}) -> out

names(out) <- l$ID %>% unique()
plyr::compact(out) -> out

data.table::rbindlist(out) -> out

#out[out$pval < 0.1,] -> out.
dplyr::mutate(out, logFDR = -log2(fdr)) -> out






### LATE GENE TABLE
list.files("h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/MANUSCRIPT/MAIN/Fig3/DONOR_DE",
           pattern = "*.RDS", full.names = T) -> files

# files[-1] -> files


lapply(files, function(x){
  readRDS(x) -> tmp
  return(tmp)
}) -> l
data.table::rbindlist(l) -> l



l[l$FDR < 0.05,] %>% dplyr::group_by(gene) %>%
  dplyr::summarise(n = dplyr::n(), chr = chr, label = label , mean_logFC = mean(avg_log2FC)) %>%
  dplyr::distinct(gene, .keep_all = T)-> out

data.table::fwrite(x = out, file = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/MANUSCRIPT/TABLES/DE_GENES_SUBJECTS_MICROGLIA.csv")





#### COMPARING LATE GENES WITH FORSBERG AND GWAS LOCI


data.table::fread(file = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/MANUSCRIPT/TABLES/DE_GENES_SUBJECTS_MICROGLIA.csv") -> out
data.table::fread(file = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/MANUSCRIPT/TABLES/GWAS_loci.csv", header = F) -> gwas
data.table::fread(file = "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/MANUSCRIPT/TABLES/dumanki_LATE_genes.csv", header = T) -> LATE

GeneOverlap::newGeneOverlap(listA = out$gene %>% unique(), listB = gwas$V1, spec = "hg19.gene") -> overlap
GeneOverlap::testGeneOverlap(overlap) -> overlap



GeneOverlap::newGeneOverlap(listA = out[out$label == "Autosomal",]$gene %>% unique(),
                            listB = LATE[LATE$Location=="Autosomal gene",]$Gene, spec = "hg19.gene") -> overlap
GeneOverlap::testGeneOverlap(overlap) -> overlap
overlap@intersection -> int


### LATE GENE TABLE
list.files("h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/MANUSCRIPT/MAIN/Fig3/DONOR_DE",
           pattern = "*.RDS", full.names = T) -> files

# files[-1] -> files


lapply(files, function(x){
  readRDS(x) -> tmp
  return(tmp)
}) -> l
data.table::rbindlist(l) -> l

l[l$FDR < 0.05 & l$gene %in% int,]






