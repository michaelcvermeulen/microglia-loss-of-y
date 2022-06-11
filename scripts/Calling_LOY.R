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
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
edb <- EnsDb.Hsapiens.v86


source("E:/LOY/scLOY/scripts/annotation_functions.R")
devtools::source_url("https://raw.githubusercontent.com/michaelcvermeulen/microglia-loss-of-y/main/scripts/LOYtools/functions.R")

### call LOY stats across all merged cohorts
### Seurat objects used for this LOY analysis are available upon request. 


## build stats files
## Seurat objects are in MERGED COHORT 
files <- list.files(path = "e:/LOY/scLOY/processed_seurat/MERGED_COHORT/", full.names = T, recursive = T)
nFeature <- 1000
meta_output <- "e:/LOY/scLOY/results/LOY_DE/APRIL28/"

for(i in files){
  readRDS(i) -> oo
  check_object_meta(oo@meta.data)
  
  oo@meta.data$TYPE -> oo@meta.data$cell_ident_MV
  DefaultAssay(oo) <- "RNA"
  print(paste0(oo@meta.data$cohort %>% unique()))
  if(oo@meta.data[oo@meta.data$cell_ident_MV %>% is.na(), ] %>% nrow() > 1){message("error in cell ident"); next}
  
  for(j in oo@meta.data$cell_ident_MV %>% unique()){
    flag <- T
    message(i)
    message(j)
    
    gsub(x = basename(i), replacement = "", pattern = ".RDS") -> name
    paste0(meta_output,"/",name,"/",j) -> dir_name
    dir.create(path = dir_name, showWarnings = FALSE, recursive = T)
    
    subset(oo, cell_ident_MV == j) -> obj
    
    
    unique(obj@meta.data$donor_organism.biomaterial_core.biomaterial_id) -> samples
    for(sam in samples){
      tryCatch({
        
        for(it in seq(1500, 5000, 500)){
          call_stats(s_obj = obj,
                     cell_type = j,
                     nUMI = it,
                     nFeature = nFeature,
                     sample = sam) -> stats
          
          data.table::fwrite(file = paste0(dir_name, "/",sam,
                                           "_samples_stat_",
                                           paste0(obj@meta.data$GEO %>% unique()),
                                           "_",it,"_",nFeature,"_",j,
                                           ".csv"),
                             sep = ",", x = stats)
        }}, error = function(e){message("Stat sample failed ",j); flag <- NULL})
      if(is.null(flag)){flag <- T; next}
    }
  }
}

## build stats files for samples with multiple samples from 1 individual
files <- list.files(path = "e:/LOY/scLOY/processed_seurat/MERGED_COHORT/", full.names = T, recursive = T)

c("e:/LOY/scLOY/processed_seurat/MERGED_COHORT/GSE148822_cleaned_filtered.RDS",
  "e:/LOY/scLOY/processed_seurat/MERGED_COHORT/GSE160936_cleaned_filtered.RDS") -> files

nFeature <- 1000
meta_output <- "e:/LOY/scLOY/results/LOY_DE/STATS_JULY21/"
for(i in files){
  readRDS(i) -> oo
  check_object_meta(oo@meta.data)
  
  oo@meta.data$TYPE -> oo@meta.data$cell_ident_MV
  DefaultAssay(oo) <- "RNA"
  print(paste0(oo@meta.data$cohort %>% unique()))
  if(oo@meta.data[oo@meta.data$cell_ident_MV %>% is.na(), ] %>% nrow() > 1){message("error in cell ident"); next}
  
  for(j in oo@meta.data$cell_ident_MV %>% unique()){
    flag <- T
    message(i)
    message(j)
    
    gsub(x = basename(i), replacement = "", pattern = ".RDS") -> name
    paste0(meta_output,"/",name,"/",j) -> dir_name
    dir.create(path = dir_name, showWarnings = FALSE, recursive = T)
    
    subset(oo, cell_ident_MV == j) -> obj
    
    obj@meta.data$donor -> obj@meta.data$donor_organism.biomaterial_core.biomaterial_id
    gsub(x= obj@meta.data$donor_organism.biomaterial_core.biomaterial_id, pattern = "/", replacement = "_") -> obj@meta.data$donor_organism.biomaterial_core.biomaterial_id
    unique(obj@meta.data$donor_organism.biomaterial_core.biomaterial_id) -> samples
    for(sam in samples){
      tryCatch({
        
        
        
        for(it in seq(1500, 5000, 500)){
          call_stats(s_obj = obj,
                     cell_type = j,
                     nUMI = it,
                     nFeature = nFeature,
                     sample = sam) -> stats
          
          data.table::fwrite(file = paste0(dir_name, "/",sam,
                                           "_samples_stat_",
                                           obj@meta.data$GEO %>% unique(),
                                           "_",it,"_",nFeature,"_",j,
                                           ".csv"),
                             sep = ",", x = stats)
        }}, error = function(e){message("Stat sample failed ",j); flag <- NULL})
      if(is.null(flag)){flag <- T; next}
    }
  }
}

output <- "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/FINAL/DE/COMPARE_CELL_TYPES_3"
for(d in c(1500,2000,2500,3000,3500,4000,4500,5000)){
  message(d)
  # read in stats files
  list.files(path = "e:/LOY/scLOY/results/LOY_DE/APRIL5/",
             recursive = T,
             pattern = paste0(d,"_1000"), full.names = T) -> files
  
  lapply(X = files, FUN = function(x){
    
    message(x)
    data.table::fread(file = x) -> tmp
    tmp$file <- basename(x)
    dirname(x) %>% dirname() %>% basename() -> tmp$cohort
    return(tmp)
    
  }) -> l
  data.table::rbindlist(l , fill = T) -> l
  l -> save
  names(l)[1] <- "cell_type"
  
  
  dplyr::mutate(l, TOTAL_CELLS = LOY_cells + NORMAL_cells, LOY_prop = LOY_cells / TOTAL_CELLS) -> l
  l[!(l$cell_type %>% is.na()),] -> l
  
  l$diagnosis -> l$CLEAN
  
  l[l$diagnosis=="PN",]$CLEAN <- "Control"
  
  l[l$CLEAN %in% c("Non-disease control","Non-symptomatic",
                   "Control","COVID-19_Moderate","MDD","Cancer","Epilepsy",
                   "Control epilepsy", "Suicide","TLE","Donation after circulatory death",
                   "Healthy","Healthy control","Influenza patient","Bronchitis","Intracranial haemorrage",
                   "healthy control","none","control","Ctrl","Healthy_Healthy",
                   "normal","ASD","Mild COVID")]$CLEAN <- "Control"
  
  l[l$CLEAN %in% c("Unknown","Not listed","unknown")]$CLEAN <- "Unknown"
  l[is.na(l$CLEAN),]$CLEAN <- "Unknown"
  
  l[l$CLEAN %in% c("Alzheimer's disease")]$CLEAN <- "AD"
  
  
  l$donor_organism.age %>% as.numeric() -> l$donor_organism.age
  
  
  
  
  ### clean up cancer diagnosis / tumor tissue
  l$diagnosis %>% unique
  l$diagnosis -> l$cancer_diagnosis
  l[l$cancer_diagnosis %in% c("NSCLC","LUAD","Brain cancer","Esophageal cancer with liver mets",
                              "Epilepsy caused by brain tumor","squamous cell carcinoma of buccal mucosa",
                              "hypertension, diabetes, prostate cancer","Colon cancer, cholangiocarcinoma",
                              "tumor and gastritis","oral cavity carcinoma",
                              "tongue cancer","carcinoma of supraglottis","human papilloma virus infection||tongue cancer",
                              "human papilloma virus infection||tonsil carcinoma","laryngeal carcinoma",
                              "lower gum cancer","endobronchial carcinoid")]$cancer_diagnosis <- "Cancer diagnosis"
  l[!(l$cancer_diagnosis %in% c("Cancer diagnosis")),]$cancer_diagnosis <- "No cancer diagnosis"
  
  ### clean up neuro
  
  
  
  l$neuro_degen_type -> l$neuro_degen_diagnosis
  l[l$neuro_degen_diagnosis %in% c(F,NA,"None","Control",
                                   "Control epilepsy",
                                   "Non-disease control",
                                   "Epilepsy control","control"),]$neuro_degen_diagnosis <- "Control"
  l[l$neuro_degen_diagnosis=="MCI",]$neuro_degen_diagnosis <- "MCI"
  l[l$neuro_degen_diagnosis=="LBD",]$neuro_degen_diagnosis <- "LBD/FTLD"
  l[l$neuro_degen_diagnosis=="IPD",]$neuro_degen_diagnosis <- "PD"
  l[l$neuro_degen_diagnosis=="FTLD",]$neuro_degen_diagnosis <- "LBD/FTLD"
  l[l$neuro_degen_diagnosis=="AD/CAA",]$neuro_degen_diagnosis <- "AD"
  l[l$neuro_degen_diagnosis=="AD/TLBD",]$neuro_degen_diagnosis <- "AD"
  
  data.table::fwrite(x= l, file = paste0(output,"/LOY_prop_table_",d,"_1000.txt"), sep = "\t")
  
}

