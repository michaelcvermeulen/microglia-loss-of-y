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
files <- list.files(path = "e:/LOY/scLOY/processed_seurat/MERGED_COHORT/", full.names = T, recursive = T)
nFeature <- 1000
meta_output <- "e:/LOY/scLOY/results/LOY_DE/APRIL5/"

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

output <- "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/FINAL/DE/COMPARE_CELL_TYPES_2"
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

### start
data.table::fread(paste0(output,"/LOY_prop_table_3000_1000.txt")) -> l

l[l$cohort=="GSE174332_cleaned_filtered",] %>% View()


total <- 75
Y <- 250
l[l$cell_type%in%c("Oligodendrocyte","Microglia","Astrocyte","Neuron","OPC"),] -> dat
l[l$cell_type%in%c("Microglia"),] -> dat
dat[!grepl(dat$file, pattern = "MD5893"),] -> dat
scales::rescale(x = scale(dat$sum_Y_exp), to = c(0,1)) -> dat$scale_sum_Y_exp
scales::rescale(x = scale(dat$TOTAL_CELLS), to = c(0,1)) -> dat$scale_TOTAL_CELLS
scales::rescale(x = scale(dat$median_nFeature_NORMAL), to = c(0,1)) -> dat$scale_median_nFeature_NORMAL

dat$scale_sum_Y_exp * 1.5 + dat$scale_TOTAL_CELLS + dat$scale_median_nFeature_NORMAL -> dat$Y_score
scales::rescale(dat$Y_score, to = c(0,1)) -> dat$Y_score_scale

dat[dat$sum_Y_exp > Y & dat$TOTAL_CELLS > total,] -> dat

ggscatter(dat[dat$sum_Y_exp > Y & dat$TOTAL_CELLS > total,],  x = "sum_Y_exp", y = "LOY_prop", size = 1, facet.by = "cell_type")
ggscatter(dat[dat$sum_Y_exp > Y & dat$TOTAL_CELLS > total,],  x = "Y_genes", y = "LOY_prop", size = 1, facet.by = "cell_type")
ggscatter(dat[dat$sum_Y_exp > Y & dat$TOTAL_CELLS > total,],  x = "Y_UMI", y = "LOY_prop", size = 1, facet.by = "cell_type")
ggscatter(dat[dat$sum_Y_exp > Y & dat$TOTAL_CELLS > total,],  x = "TOTAL_CELLS", y = "LOY_prop", size = 1, facet.by = "cell_type")
ggscatter(dat[dat$sum_Y_exp > Y & dat$TOTAL_CELLS > total,],  x = "Y_score_scale", y = "LOY_prop", size = 1, facet.by = "cell_type")


ggscatter(dat, y = "LOY_prop", x = "sum_Y_exp", alpha = 0.8, size = "TOTAL_CELLS",
          color = "neuro_degen_diagnosis", label = "donor_organism.age",
          add.params = list(color = "red")) + geom_smooth(se = T , method = "lm")


ggscatter(dat, y = "LOY_prop", x = "donor_organism.age", alpha = 0.8, size = "TOTAL_CELLS",
          color = "neuro_degen_diagnosis", label = "donor_organism.age",
          add.params = list(color = "red"))



gghistogram(dat, x = "sum_Y_exp", fill = "cell_type", density = T, add = "median", palette = "aaas") +
  facet_grid(~ cell_type) -> p;
ggpar(p,legend = "none",x.text.angle = 45, xlab = "Mean sum of Y gene percent expressed")

gghistogram(dat, x = "median_nFeature_NORMAL", fill = "cell_type", density = T, add = "median") +
  facet_grid(~ cell_type) -> p;
ggpar(p,legend = "none",x.text.angle = 45)

ggviolin(dat[dat$sum_Y_exp > Y & dat$TOTAL_CELLS > total,],  x = "cell_type", y = "LOY_prop", size = 1, add = "boxplot")

## LOY score
scales::rescale(x = scale(dat$sum_Y_exp), to = c(0,1)) -> dat$scale_sum_Y_exp
scales::rescale(x = scale(dat$TOTAL_CELLS), to = c(0,1)) -> dat$scale_TOTAL_CELLS
scales::rescale(x = scale(dat$median_nFeature_NORMAL), to = c(0,1)) -> dat$scale_median_nFeature_NORMAL

### proportion FIGURE 1 #######
data.table::fread(paste0(output,"/LOY_prop_table_3500_1000.txt")) -> l

l[l$cell_type=="CAM",]$cell_type <- "Microglia"

l[l$TOTAL_CELLS > 50 & l$sum_Y_exp > 250, ] -> l
l[l$cell_type%in%c("Oligodendrocyte","Microglia","Astrocyte","Neuron","OPC","Pericyte"),] -> dat
dat[!grepl(dat$file, pattern = "MD5893"),] -> dat

ggscatter(l, y = "LOY_prop", x = "sum_Y_exp", facet.by = "cell_type")


ggscatter(dat, x = "donor_organism.age", y = "LOY_prop", add = "reg.line", add.params = list(color = "red")) + stat_cor()

dat[dat$cell_type == "Microglia",] -> dat.
rescale(x = scale(dat.$TOTAL_CELLS), to = c(0,1)) -> out
rescale(x = scale(dat.$sum_Y_exp), to = c(0,1)) -> out2
out*out2 -> out
lm(formula = LOY_prop ~  donor_organism.age, dat = dat., weights = out) %>% summary()


ggscatter(dat., x = "donor_organism.age", y = "LOY_prop", add = "reg.line", add.params = list(color = "red")) + stat_cor()



ggboxplot(dat[dat$cell_type %in% c("Microglia","Astrocyte","OPC","Oligodendrocyte","Neuron")],
          x = "cell_type", y = "LOY_prop")
ggboxplot(dat, x = "cancer_diagnosis", y = "LOY_prop") + stat_compare_means()
lm(dat, formula = LOY_prop ~ CLEAN + donor_organism.age ) %>% summary()

ggboxplot(dat, x = "cancer_diagnosis", y = "LOY_prop") + stat_compare_means()

factor(x = dat$cell_type, levels = c("Microglia","Astrocyte","Neuron","Oligodendrocyte","OPC","Endothelial"),
       labels = c("Microglia","Astrocyte","Neuron","Oligo","OPC","Endothelial")) -> dat$cell_type
factor(dat$neuro_degen_diagnosis, levels = c("Control","AD","HD","PD","LBD/FTLD","ALS","Brain cancer","Seizure disorder")) -> dat$neuro_degen_diagnosis
dat <- dat[dat$neuro_degen_diagnosis %in% c("AD", "Control","HD","PD","LBD/FTLD","ALS","Seizure disorder"),]
ggboxplot(dat[dat$cell_type %in% c("Microglia","Astrocyte","Neuron","OPC","Oligo"),],
          x = "neuro_degen_diagnosis", font.xtickslab = c(8), size = 0.3,
          xlab = "", ylab = "Loss of Y proportion", fill = "neuro_degen_diagnosis", palette = "aaas", add = "jitter", add.params = list(size = 0.9),
          y = "LOY_prop") %>% ggpar(x.text.angle = 45, ylim = c(0,0.7), legend = "none") -> tmp
tmp + facet_grid(~ cell_type) -> tmp

ggboxplot(dat[dat$cell_type %in% c("Microglia","Astrocyte","Neuron","OPC","Oligo"),],
          x = "neuro_degen_diagnosis", font.xtickslab = c(8), size = 0.3,
          xlab = "", ylab = "Loss of Y proportion", fill = "neuro_degen_diagnosis", palette = "aaas",
          y = "LOY_prop") %>% ggpar(x.text.angle = 45, ylim = c(0,0.7), legend = "none") -> tmp
tmp + facet_grid(~ cell_type) -> tmp

ggsave(plot = tmp, filename = paste0(output,"/NEURO_LOY_BRAIN_CELL_TYPES_3000_1000_no_jitter.tiff"), dpi = "retina", units = "cm", height = 9, width = 20)

## LOY by sample




### bin by age

for(t in c(50,100,150,200,250,300)){
  for(y in c(200,225,250,275,300)){
    for(dep in c(1500,2000,2500,3000,3500,4000,4500,5000)){
      
      
      message(paste0(t,"_",y,"_",dep,".pdf"))
      output <- "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/FINAL/DE/COMPARE_CELL_TYPES"
      data.table::fread(paste0(output,"/LOY_prop_table_",dep,"_1000.txt")) -> l
      #data.table::fread(paste0(output,"/LOY_prop_table.txt")) -> l
      l[l$TOTAL_CELLS > t & l$sum_Y_exp > y, ] -> dat
      dat[!grepl(dat$file, pattern = "MD5893"),] -> dat
      
      cut(x = dat$donor_organism.age, breaks = c(-20,0,30,60,80,100),
          labels = c("Prenatal","0-29","30-59","60-79","80-100")) -> dat$quantile_age
      dat[!is.na(dat$donor_organism.age),] -> tmp
      
      
      factor(x = dat$cell_type, levels = c("Microglia","Astrocyte","Neuron","Oligodendrocyte","OPC","Monocyte","Macrophage","T cell","NK"),
             labels = c("Microglia","Astrocyte","Neuron","OL","OPC","Monocyte","Macrophage","T cell", "NK")) -> dat$cell_type
      
      dat[!is.na(dat$donor_organism.age),] -> tmp
      tmp[tmp$cell_type %in% c("Microglia","Astrocyte","Neuron","Oligodendrocyte","OPC","Monocyte","Macrophage","T cell","NK","B cell","cDC") &
            !(tmp$quantile_age %in% c("Prenatal")),] -> tmp
      ggboxplot(tmp, x = "quantile_age", fill = "quantile_age", size = 0.2,
                xlab = "", ylab = "Loss of Y proportion", add = "jitter", add.params = list(size = 0.3),
                y = "LOY_prop") %>% ggpar(x.text.angle = 45, ylim = c(0,0.7), legend = "none", font.xtickslab = c(8)) -> e
      e + facet_grid(~ cell_type) + scale_fill_brewer(palette = "Greens") -> e
      ggsave(plot = e, filename = paste0(output,"/AGE_BRAIN_CELL_TYPES_jitter",t,"_",y,"_",dep,".tiff"), dpi = "retina", units = "in",
             height = 4.76, width = 8.95)
      
      tmp[tmp$cell_type %in% c("Monocyte") & !(tmp$CLEAN %in% c("Unknown")) &
            !(tmp$quantile_age %in% c("Prenatal")),] -> tmp.
      
      lm(formula = LOY_prop ~ donor_organism.age + CLEAN, data = tmp.) %>% summary()
      
      
      hist(tmp.$LOY_prop)
      
      
      tmp[tmp$cell_type == "Monocyte" & CLEAN == "Control" & tissue %in% c("PBMC","peripheral blood"),] %>%
        ggscatter(x = "donor_organism.age", y = "LOY_prop")
      lm(data = tmp[tmp$cell_type == "Microglia" & CLEAN == "Control",], formula = LOY_prop ~ donor_organism.age ) %>% summary()
      lm(data = tmp[tmp$cell_type == "Astrocyte",], formula = LOY_prop ~ donor_organism.age ) %>% summary()
      lm(data = tmp[tmp$cell_type == "Neuron",], formula = LOY_prop ~ donor_organism.age ) %>% summary()
      
      
      tmp[tmp$cell_type %in% c("Microglia","Astrocyte","Neuron","Oligo","OPC") &
            !(tmp$quantile_age %in% c("Prenatal")),] -> tmp
      ggboxplot(tmp, x = "quantile_age", fill = "quantile_age",
                xlab = "", ylab = "Loss of Y proportion",
                y = "LOY_prop") %>% ggpar(x.text.angle = 45, ylim = c(0,0.7), legend = "none", font.xtickslab = c(8)) -> e
      e + facet_grid(~ cell_type) + scale_fill_brewer(palette = "Greens") -> e
      ggsave(plot = e, filename = paste0(output,"/AGE_BRAIN_CELL_TYPES_",t,"_",y,"_",dep,".pdf"), dpi = "retina", units = "cm", height = 9, width = 13)
      
      
      factor(x = tmp$cell_type, levels = c("T cell","Monocyte","Macrophage","Microglia","Neuron","OPC","Oligodendrocyte","Astrocyte"),
      ) -> tmp$cell_type
      tmp[tmp$cell_type %in% c("T cell","Monocyte","Macrophage","Microglia","Neuron","OPC","Oligodendrocyte","Astrocyte") &
            !(tmp$quantile_age %in% c("Prenatal")),] %>%
        ggboxplot(x = "quantile_age", fill = "quantile_age", size = 0.2,
                  xlab = "", ylab = "Loss of Y proportion", add = "jitter", add.params = list(size = 0.3),
                  y = "LOY_prop") %>% ggpar(x.text.angle = 45, ylim = c(0,0.7), legend = "none", font.xtickslab = c(8)) -> e
      e + facet_grid(~ cell_type) + scale_fill_brewer(palette = "Greens") +
        theme(strip.text.x = element_text(size = 7)) -> e
      e + ylim(0,0.6) -> e
      ggsave(plot = e, filename = paste0(output,"/AGE_BRAIN_IMMUNE_CELL_TYPES_jitter",t,"_",y,"_",dep,".pdf"), dpi = "retina", units = "cm", height = 9, width = 13)
      
      
      tmp[tmp$cell_type %in% c("Microglia") &
            !(tmp$quantile_age %in% c("Prenatal")),] %>%
        ggscatter(x = "donor_organism.age", y = "LOY_prop") + stat_cor(method = "pearson") -> p
      ggsave(plot = p, filename = paste0(output,"/MICROGLIA_age_SCATTER_jitter",t,"_",y,"_",dep,".pdf"), dpi = "retina", units = "cm", height = 9, width = 13)
      
    }}}

lm(data = tmp[tmp$cell_type == "Microglia",], formula = LOY_prop ~ donor_organism.age ) %>% summary()
lm(data = tmp[tmp$cell_type == "Astrocyte",], formula = LOY_prop ~ donor_organism.age ) %>% summary()
lm(data = tmp[tmp$cell_type == "Neuron",], formula = LOY_prop ~ donor_organism.age ) %>% summary()




# how many samples of each cell-type have available age data ######



## LOY prop by cell type ###
output <- "h:/LOY/scLOY/scLOY_NOV18/plots/MARCH_12/FINAL/DE/COMPARE_CELL_TYPES"
data.table::fread(paste0(output,"/LOY_prop_table_3000_1000.txt")) -> l
l[l$TOTAL_CELLS > 75 & l$sum_Y_exp > 250, ] -> dat
dat[!grepl(dat$file, pattern = "MD5893"),] -> dat

dat[dat$cell_type %in% c("Microglia","Astrocyte","Neuron","Oligodendrocyte","OPC"),] -> dat

stderr <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}
cut(x = dat$donor_organism.age, breaks = c(-20,0,30,60,80,100),
    labels = c("Prenatal","0-29","30-59","60-79","80-100")) -> dat$quantile_age

dat %>% dplyr::group_by(cell_type,neuro_degen_diagnosis) %>%
  dplyr::summarise(median_LOY = median(LOY_prop), se_LOY = stderr(LOY_prop, na.rm =T), count = dplyr::n()) -> tmp
reshape2::dcast(tmp, cell_type ~ neuro_degen_diagnosis, value.var = "median_LOY")

dat %>% dplyr::group_by(cell_type,quantile_age) %>%
  dplyr::summarise(median_LOY = median(LOY_prop), se_LOY = stderr(LOY_prop, na.rm =T), count = dplyr::n()) -> tmp
reshape2::dcast(tmp, cell_type ~ quantile_age, value.var = "median_LOY")


factor(x = dat$cell_type, levels = c("Microglia","Astrocyte","Neuron","Oligodendrocyte","OPC","Endothelial"),
       labels = c("Microglia","Astrocyte","Neuron","Oligo","OPC","Endothelial")) -> dat$cell_type

ggboxplot(dat, x = "cell_type", y = "LOY_prop")



### STATS MICROGLIA BETWEEN CONTROL AND AD
dat[dat$cell_type == "Microglia" & is.na(dat$sequencing_method2),]$sequencing_method2 <- "single-cell"
dat[dat$cell_type == "Microglia" & is.na(dat$sequencing_method),]$sequencing_method <- "10X 3' v2"
dat[dat$cell_type == "Microglia",]$sequencing_method2 %>% toupper() -> dat[dat$cell_type == "Microglia",]$sequencing_method2
dat[dat$cell_type == "Microglia",]$sequencing_method %>% toupper() -> dat[dat$cell_type == "Microglia",]$sequencing_method

dat[dat$cell_type == "Microglia" & neuro_degen_diagnosis %in% c("AD","Control"),] -> tmp
stringr::str_split(string = tmp$file, pattern = "_GSE", simplify = T)[,2] -> tmp$cohort



tmp$neuro_degen_diagnosis %>% as.character() -> tmp$neuro_degen_diagnosis
tmp[!is.na(tmp$donor_organism.age) ,] -> tmp

wilcox.test(LOY_prop ~ neuro_degen_diagnosis + donor_organism.age, data = tmp)

lm(dat = tmp, formula = neuro_degen_diagnosis ~ LOY_prop)

anova_test(data = tmp, formula = LOY_prop ~ donor_organism.age + cohort + neuro_degen_diagnosis)


#### STATS LOY and QC METRICS #############

l[l$TOTAL_CELLS > 50 & l$sum_Y_exp > 0, ] -> dat
dat[!grepl(dat$file, pattern = "MD5893"),] -> dat

ggscatter(data = dat[dat$cell_type %in% c("Neuron","Microglia","Astrocyte")], palette = "lancet", size = 1.1,
          x = "sum_Y_exp", y = "LOY_prop", color = "neuro_degen_diagnosis", facet.by = "cell_type") %>% ggpar(legend.title = "")

ggscatter(data = dat[dat$cell_type %in% c("Macrophage","Monocyte")], palette = "lancet", size = 1.1,
          x = "sum_Y_exp", y = "LOY_prop", color = "cancer_diagnosis", facet.by = "cell_type") %>% ggpar(legend.title = "")


l[l$TOTAL_CELLS > 50 & l$sum_Y_exp > 200, ] -> dat
ggscatter(data = dat[dat$cell_type %in% c("Neuron","Microglia","Astrocyte")], palette = "lancet", size = 1.1,
          x = "median_nUMI_NORMAL", y = "LOY_prop", color = "neuro_degen_diagnosis", facet.by = "cell_type") %>% ggpar(legend.title = "")


################ load Single-cell single sample #########

files <- list.files(path = "e:/LOY/scLOY/results/LOY_DE/STATS_SAMPLE/", pattern = "*.csv",
                    full.names = T, recursive = T)
files <- files[!grepl(x = files , pattern = "All")]
files <- files[!grepl(x = files , pattern = "stat")]


lapply(X = files, FUN = function(x){
  
  data.table::fread(x, data.table = F) -> dat
  dat[dat$LOY_cells >= 25, ] -> dat
  
  if(nrow(dat)==0){return(NULL)} else{
    message(x)
    
    dat$sample <- basename(x)
    dat$NORMAL_pct_exp <- dat$pct.2
    dat$DIR <- ifelse(test = dat$avg_log2FC > 0, yes = "UP", no = "DOWN")
    dat$FDR -> dat$adj_pval
    dat$file <- x
    dat[dat$p_val < 0.05,] -> dat
    
    return(dat)
  }
  
}) -> l
data.table::rbindlist(l, fill = T) -> l




#### BAR PLOT | NUMBER OF EACH CELL TYPE
data.table::fread(paste0(output,"/LOY_prop_table_3000_1000.txt")) -> l

#l[l$cell_type=="CAM",]$cell_type <- "Microglia"

l[l$TOTAL_CELLS > 75 & l$sum_Y_exp > 250, ] -> l
l[l$cell_type%in%c("Oligodendrocyte","Microglia","Astrocyte","Neuron","OPC","Pericyte"),] -> dat
dat[!grepl(dat$file, pattern = "MD5893"),] -> dat

dat %>% group_by(cell_type) %>% summarise(LOY = sum(LOY_cells), NORMAL = sum(NORMAL_cells)) %>%
  dplyr::mutate(total = LOY + NORMAL, prop = LOY / total) -> out

ggboxplot(dat, x = "cell_type", y = "LOY_prop")


