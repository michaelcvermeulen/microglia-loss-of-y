########################################## 
run_DE <- function(s_obj, sample, cell_type, mode = "LOY",
                   nUMI = 2500, nFeature = 1000, logfc.threshold = 0, min.pct = 0.05){
  
  message(sample," | ", cell_type, " ")
  cells <- s_obj@meta.data[s_obj@meta.data$donor_organism.biomaterial_core.biomaterial_id == sample &
                             s_obj@meta.data$cell_ident_MV == cell_type &
                             s_obj@meta.data$nCount_RNA >= nUMI &
                             s_obj@meta.data$nFeature_RNA >= nFeature,] %>% rownames()
  
  subset(s_obj, cells = cells) -> s_obj
  
  s_obj@meta.data$cohort %>% unique() -> cohort
  
  if(mode == "LOY"){
    DE <- FindMarkers(s_obj, 
                      ident.1 = "LOY", ident.2 = "NORMAL", 
                      group.by = "LOY", 
                      test.use = "MAST", 
                      latent.vars = c("nCount_RNA","percent.rb","nFeature_RNA","donor_organism.biomaterial_core.biomaterial_id"), 
                      logfc.threshold = logfc.threshold)
 } else {
    DE <- FindMarkers(s_obj, 
                      ident.1 = "0", ident.2 = "4", 
                      group.by = "LOY_levels", 
                      test.use = "MAST", 
                      latent.vars = c("nCount_RNA","percent.rb","nFeature_RNA","donor_organism.biomaterial_core.biomaterial_id"), 
                      logfc.threshold = logfc.threshold)
 }
  
  rownames(DE) -> DE$gene
  DE$nUMI <- nUMI
  DE$nFeature <- nFeature
  DE$cell_type <- cell_type
  DE$mode <- mode
  DE$sample <- sample
  DE$cohort <- cohort
  DE$LOY_cells <- s_obj@meta.data[s_obj@meta.data$LOY=="LOY",] %>% nrow()
  DE$NORMAL_cells <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",] %>% nrow()
  DE$LEVEL2_cells <- s_obj@meta.data[s_obj@meta.data$LOY_levels=="2",] %>% nrow()
  DE$LEVEL3_cells <- s_obj@meta.data[s_obj@meta.data$LOY_levels=="3",] %>% nrow()
  DE$LEVEL4_cells <- s_obj@meta.data[s_obj@meta.data$LOY_levels=="4",] %>% nrow()
  
  
  DE$median_nUMI_LOY4 <- s_obj@meta.data[s_obj@meta.data$LOY_levels=="4",]$nCount_RNA %>% median()
  DE$median_nUMI_LOY <- s_obj@meta.data[s_obj@meta.data$LOY_levels=="0",]$nCount_RNA %>% median()
  DE$median_nUMI_NORMAL <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$nCount_RNA %>% median()
  
  DE$median_nFeature_LOY4 <- s_obj@meta.data[s_obj@meta.data$LOY_levels=="4",]$nFeature_RNA %>% median()
  DE$median_nFeature_LOY <- s_obj@meta.data[s_obj@meta.data$LOY_levels=="0",]$nFeature_RNA %>% median()
  DE$median_nFeature_NORMAL <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$nFeature_RNA %>% median()
  
  DE$expr_Y <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$normY %>% median()
  DE$expr_PAR <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$normPAR %>% median()
  DE$expr_PAR_LOY <- s_obj@meta.data[s_obj@meta.data$LOY=="LOY",]$normPAR %>% median()
  
  DE$PAR_UMI <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$PAR_UMI %>% mean()
  DE$PAR_genes <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$PAR_genes %>% mean()
  DE$PAR_score <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$PAR_score1 %>% mean()
  DE$PAR_score_LOY <- s_obj@meta.data[s_obj@meta.data$LOY=="LOY",]$PAR_score1 %>% mean()
  
  DE$Y_UMI <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$chrY_UMI %>% mean()
  DE$Y_genes <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$chrY_genes %>% mean()
  DE$Y_score <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$Y_score1 %>% mean()
  
  Y <- c("RPS4Y1","ZFY","LINC00278","PRKY","USP9Y","DDX3Y","UTY","NLGN4Y","TTTY14","EIF1AY","KDM5D")
  PAR <- c("GTPBP6","CSF2RA","IL3RA","SLC25A6","ASMTL","P2RY8","AKAP17A","DHRSX","CD99","ZBED1")
  OTHER <- c("APOE","ROBO1","OXR1","CACNA1A","FCGBP","RORA","SLIT1","SLIT2","ROBO2","ROBO3","ROBO4")
  
  ### pull data from dotplots for pct exp info + avg exp
  # add to the DE table for genes of interest 
  dat1 <- Seurat::DotPlot(s_obj, group.by = "LOY",
                          features = expressed_Y_genes(s_obj))
  
  dat2 <- Seurat::DotPlot(s_obj, group.by = "LOY",
                          features = expressed_PAR_genes(s_obj))
  
  dat3 <- Seurat::DotPlot(s_obj, group.by = "LOY",
                          features = unique(OTHER[OTHER %in% rownames(s_obj)]))
  
  tmp <- dat1$data[dat1$data$id == "NORMAL",]
  sum(tmp$pct.exp) -> DE$sum_Y_exp
  tmp <- dat2$data[dat2$data$id == "NORMAL",]
  sum(tmp$pct.exp) -> DE$sum_PAR_exp
  tmp <- dat3$data[dat3$data$id == "NORMAL",]
  sum(tmp$pct.exp) -> DE$sum_OTHER_exp
  
  dat1 <- dat1$data[dat1$data$id == "NORMAL" & dat1$data$features.plot %in% c(Y),]
  dat2 <- dat2$data[dat2$data$id == "NORMAL" & dat2$data$features.plot %in% c(PAR),]
  dat3 <- dat3$data[dat3$data$id == "NORMAL" & dat3$data$features.plot %in% c(OTHER),]
  
  rbind(dat1, dat2, dat3) -> dat
  
  
  lapply(c(Y,PAR,OTHER), function(x) {
    
    
    if(nrow(dat[dat$features.plot==x,]) < 1){ 
      data.frame(avg.exp = NA, pct.exp = NA) -> out 
      rownames(out) <- x
      names(out) <- paste0(names(out),"_",x)
    } else {
      dat[dat$features.plot==x,] %>% dplyr::select(avg.exp,pct.exp) %>% as.data.frame() -> out
      names(out) <- paste0(names(out),"_",x)
    }
    
    return(out)
    
  }) -> l
  do.call("cbind",l) -> l
  
  
  cbind(DE,l) -> DE
  
  p.adjust(DE$p_val, method = "fdr") -> DE$FDR
  p.adjust(DE$p_val, method = "BH") -> DE$bh
  p.adjust(DE$p_val, method = "hochberg") -> DE$hoch
  p.adjust(DE$p_val, method = "bonferroni") -> DE$bonf
  
  return(DE)
  
} 

#######################################
loop_through_DE <- function(o, nUMI, nFeature, logFC, min_cells, output, type, min.pct = 0.05, mode = "LOY"){
  
    cells <- o@meta.data[o@meta.data$TYPE == type & 
                  o@meta.data$donor_organism.sex == "male" & 
                  o@meta.data$scDblFinder.class == "singlet",] %>% rownames()
    subset(o, cells = cells) -> o
    o@meta.data$cell_ident_MV <- o@meta.data$TYPE
    
    DefaultAssay(o) <- "RNA"
    
    for(i in o@meta.data$donor_organism.biomaterial_core.biomaterial_id %>% unique()){
        sample <- i
        
        
        o@meta.data[o@meta.data$donor_organism.biomaterial_core.biomaterial_id == sample & 
                      o@meta.data$cell_ident_MV == type & 
                      o@meta.data$nCount_RNA >= nUMI & 
                      o@meta.data$nFeature_RNA >= nFeature,] -> dat
        
        cohort <- dat$cohort %>% unique()
        message(paste0(i," : ",cohort))
        
        data.frame(LOY = dat[dat$LOY=="LOY",] %>% nrow(), 
                   NORMAL = dat[dat$LOY=="NORMAL",] %>% nrow(),
                   LOY4 = dat[dat$LOY_levels=="4",] %>% nrow()) -> df
        
        if(mode == "LOY"){
          if(df[,1] < min_cells | df[,2] < min_cells ){ message("not enough cells"); next }
          
          DE <- run_DE(s_obj = o, sample = sample, cell_type = type, 
                       mode = "LOY", nUMI = nUMI, nFeature = nFeature, logfc.threshold = logFC, min.pct = min.pct)
          gsub(x = type, pattern = "[\\/\\+]", replacement = "_") -> cell_type
          gsub(x = sample, pattern = "[\\/\\+]", replacement = "_") -> sample
          gsub(x = sample, pattern = " ", replacement = "_") -> sample
          
          #LOY_dataset_summary(df = dat, sample_size_min = min_cells) -> LOY_meta
          #DE <- dplyr::left_join(x = DE, y = LOY_meta , 
          #                       by = c("sample"="donor_organism.biomaterial_core.biomaterial_id",
          #                              "cell_type"="cell_ident_MV"))
          
          data.table::fwrite(x = DE, sep = ",", 
                             file = gsub(x = paste0(output,"/",cohort,"_",sample,"_",cell_type,"_",
                                                    nUMI,"_",nFeature,"_",mode,".csv"), pattern = " ", replacement = ""))
          
        } else{
          if(df[,1] < min_cells | df[,3] < min_cells ){ message("not enough cells"); next  }
          
          
          DE <- run_DE(s_obj = o, sample = sample, cell_type = type, 
                       mode = "LOY_levels", nUMI = nUMI, nFeature = nFeature, 
                       logfc.threshold = logFC, min.pct = 0.05)
          
          #LOY_dataset_summary(df = dat, sample_size_min = min_cells) -> LOY_meta
          #DE <- dplyr::left_join(x = DE, y = LOY_meta , 
          #                       by = c("sample"="donor_organism.biomaterial_core.biomaterial_id",
          #                              "cell_type"="cell_ident_MV"))
          
          gsub(x = type, pattern = "[\\/\\+]", replacement = "_") -> cell_type
          gsub(x = sample, pattern = "[\\/\\+]", replacement = "_") -> sample
          gsub(x = sample, pattern = " ", replacement = "_") -> sample
          data.table::fwrite(x = DE, sep = ",", 
                             file = gsub(x = paste0(output,"/",cohort,"_",sample,"_",type,"_",
                                                    nUMI,"_",nFeature,"_",mode,".csv"), pattern = "[ ]", replacement = ""))
        }
      } # end loop 2 
    } # end loop 1
  


#######################################
loop_through_DE_mic <- function(list, nUMI, nFeature, logFC, min_cells, output,min.pct = 0.05, mode = "LOY"){
  
  for(k in list){
    
    message(k)
    readRDS(k) -> o
    subset(o, donor_organism.sex == "male" & cell_subtype_MV %in% c("CAM","Microglia")) -> o
    o@meta.data$cell_ident_MV <- "Microglia/CNS-Macrophage"
    #o@meta.data$donor_organism.biomaterial_core.biomaterial_id <- stringr::str_split(string = o@meta.data$title, pattern = " ", simplify = T)[,1]
    #o@meta.data$donor_organism.biomaterial_core.biomaterial_id <- o@meta.data$donor
    DefaultAssay(o) <- "RNA"
    
    for(i in o@meta.data$donor_organism.biomaterial_core.biomaterial_id %>% unique()){
      for(j in o@meta.data$cell_ident_MV %>% unique()){
        
        sample <- i
        cell_type <- j
        cohort <- o@meta.data$cohort %>% unique()
        
        
        o@meta.data[o@meta.data$donor_organism.biomaterial_core.biomaterial_id == sample & 
                      o@meta.data$cell_ident_MV == cell_type & 
                      o@meta.data$nCount_RNA >= nUMI & 
                      o@meta.data$nFeature_RNA >= nFeature,] -> dat
        
        data.frame(LOY = dat[dat$LOY=="LOY",] %>% nrow(), 
                   NORMAL = dat[dat$LOY=="NORMAL",] %>% nrow(),
                   LOY4 = dat[dat$LOY_levels=="4",] %>% nrow()) -> df
        
        if(mode == "LOY"){
          if(df[,1] < min_cells | df[,2] < min_cells ){ message("not enough cells"); next }
          
          DE <- run_DE(s_obj = o, sample = sample, cell_type = cell_type, 
                       mode = "LOY", nUMI = nUMI, nFeature = nFeature, logfc.threshold = logFC, min.pct = min.pct)
          gsub(x = cell_type, pattern = "[\\/\\+]", replacement = "_") -> cell_type
          gsub(x = sample, pattern = "[\\/\\+]", replacement = "_") -> sample
          gsub(x = sample, pattern = " ", replacement = "_") -> sample
          
          #LOY_dataset_summary(df = dat, sample_size_min = min_cells) -> LOY_meta
          #DE <- dplyr::left_join(x = DE, y = LOY_meta , 
          #                       by = c("sample"="donor_organism.biomaterial_core.biomaterial_id",
          #                              "cell_type"="cell_ident_MV"))
          
          data.table::fwrite(x = DE, sep = ",", 
                             file = gsub(x = paste0(output,"/",cohort,"_",sample,"_",cell_type,"_",
                                                    nUMI,"_",nFeature,"_",mode,".csv"), pattern = " ", replacement = ""))
          
        } else{
          if(df[,1] < min_cells | df[,3] < min_cells ){ message("not enough cells"); next  }
          
          
          DE <- run_DE(s_obj = o, sample = sample, cell_type = cell_type, 
                       mode = "LOY_levels", nUMI = nUMI, nFeature = nFeature, 
                       logfc.threshold = logFC, min.pct = 0.05)
          
          #LOY_dataset_summary(df = dat, sample_size_min = min_cells) -> LOY_meta
          #DE <- dplyr::left_join(x = DE, y = LOY_meta , 
          #                       by = c("sample"="donor_organism.biomaterial_core.biomaterial_id",
          #                              "cell_type"="cell_ident_MV"))
          
          gsub(x = cell_type, pattern = "[\\/\\+]", replacement = "_") -> cell_type
          gsub(x = sample, pattern = "[\\/\\+]", replacement = "_") -> sample
          gsub(x = sample, pattern = " ", replacement = "_") -> sample
          data.table::fwrite(x = DE, sep = ",", 
                             file = gsub(x = paste0(output,"/",cohort,"_",sample,"_",cell_type,"_",
                                                    nUMI,"_",nFeature,"_",mode,".csv"), pattern = "[ ]", replacement = ""))
        }
      } # end loop 2 
    } # end loop 1
  }
}




##################################
test_DE <- function(s_obj, sample, cell_type, mode = "LOY", logfc.threshold = 0.1, 
                    min = 1000, max = 5000, iter = 100){
  
  s_obj@meta.data[s_obj@meta.data$donor_organism.biomaterial_core.biomaterial_id == sample & 
                    s_obj@meta.data$cell_ident_MV == cell_type,] -> dat
  
  data.frame(LOY = dat[dat$LOY=="LOY",] %>% nrow(), 
             NORMAL = dat[dat$LOY=="NORMAL",] %>% nrow(),
             LOY4 = dat[dat$LOY_levels=="4",] %>% nrow()) -> df
  
  
  lapply(seq(min,max,iter), function(x){
    
    DE <- run_DE(s_obj = s_obj, sample = sample, cell_type = cell_type, 
                 mode = mode, nUMI = x, nFeature = 0, logfc.threshold = logfc.threshold)
    
    return(DE)
    
  }) -> l
  do.call("rbind",l) -> l
  
  return(l)
}



#############################################
plot_test_DE <- function(l, p_cutoff = 0.05){
  
  l[l$p_val_adj < p_cutoff,] -> l
  l$logP <- -log(l$p_val_adj)
  l[l$logP == Inf,]$logP <- max(l[l$logP != Inf,]$logP) 
  
  table(l[l$avg_log2FC > 0,]$gene, l[l$avg_log2FC > 0,]$nUMI) -> Up
  colSums(Up) -> Up
  
  table(l[l$avg_log2FC < 0,]$gene, l[l$avg_log2FC < 0,]$nUMI) -> Down
  colSums(Down) -> Down
  
  rbind(Up,Down) -> dat
  t(dat) %>% as.data.frame() -> dat
  rownames(dat) -> dat$nUMI
  dat$Down <- -(dat$Down)
  reshape2::melt(dat) -> dat
  
  ggbarplot(data = dat, x = "nUMI", y = "value", fill = "variable", 
            
            color = "variable") -> plot
  plot <- ggpar(plot, x.text.angle = 45, palette = get_palette("aaas",2)[c(2,1)], 
                xlab = "nUMI Threshold", ylab = "DE genes" ,
                legend.title = "")
  return(plot)
}

#############################################
plot_genes <- function(l, genes, mode = "logFC"){
  
  l$logP <- -log(l$p_val_adj)
  l[l$logP == Inf,]$logP <- max(l[l$logP != Inf,]$logP)
  l[l$gene %in% genes,] -> l
  
  
  if(mode == "logFC"){
  ggline(data = l, x = "nUMI", y = "avg_log2FC", fill = "gene", 
            color = "gene") -> plot
  plot <- ggpar(plot, x.text.angle = 45,  
                xlab = "nUMI Threshold", ylab = "avg log2FC" ,
                legend.title = "")
  return(plot)}
  
  if(mode == "p"){
    ggline(data = l, x = "nUMI", y = "logP", fill = "gene", 
           color = "gene") -> plot
    plot <- ggpar(plot, x.text.angle = 45,  
                  xlab = "nUMI Threshold", ylab = "-logP" ,
                  legend.title = "")
    return(plot)}
  
}


##################################



call_stats <- function(s_obj, cell_type, sample, nUMI, nFeature){
  
  s_obj@meta.data[s_obj@meta.data$donor_organism.biomaterial_core.biomaterial_id %in% sample & 
                  s_obj@meta.data$cell_ident_MV %in% cell_type & 
                  s_obj@meta.data$nCount_RNA >= nUMI & 
                  s_obj@meta.data$nFeature_RNA >= nFeature,] %>% rownames() -> cells
                  
  subset(s_obj, cells = cells ) -> s_obj
  
  DE <- cell_type %>% as.data.frame()
  cohort <- s_obj@meta.data$GEO %>% unique
  DE$nUMI_cut <- nUMI
  DE$nFeature_cut <- nFeature
  DE$sample <- paste(sample,sep = ", ", collapse = " ")
  DE$LOY_cells <- s_obj@meta.data[s_obj@meta.data$LOY=="LOY",] %>% nrow()
  DE$NORMAL_cells <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",] %>% nrow()
  DE$LEVEL2_cells <- s_obj@meta.data[s_obj@meta.data$LOY_levels=="2",] %>% nrow()
  DE$LEVEL3_cells <- s_obj@meta.data[s_obj@meta.data$LOY_levels=="3",] %>% nrow()
  DE$LEVEL4_cells <- s_obj@meta.data[s_obj@meta.data$LOY_levels=="4",] %>% nrow()
  
  DE$median_nUMI_LOY4 <- s_obj@meta.data[s_obj@meta.data$LOY_levels=="4",]$nCount_RNA %>% median()
  DE$median_nUMI_LOY <- s_obj@meta.data[s_obj@meta.data$LOY_levels=="0",]$nCount_RNA %>% median()
  DE$median_nUMI_NORMAL <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$nCount_RNA %>% median()
  
  DE$median_nFeature_LOY4 <- s_obj@meta.data[s_obj@meta.data$LOY_levels=="4",]$nFeature_RNA %>% median()
  DE$median_nFeature_LOY <- s_obj@meta.data[s_obj@meta.data$LOY_levels=="0",]$nFeature_RNA %>% median()
  DE$median_nFeature_NORMAL <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$nFeature_RNA %>% median()
  
  DE$expr_Y <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$normY %>% median()
  DE$expr_PAR <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$normPAR %>% median()
  DE$expr_PAR_LOY <- s_obj@meta.data[s_obj@meta.data$LOY=="LOY",]$normPAR %>% median()
  
  DE$PAR_UMI <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$PAR_UMI %>% mean()
  DE$PAR_genes <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$PAR_genes %>% mean()
  DE$PAR_score <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$PAR_score1 %>% mean()
  DE$PAR_score_LOY <- s_obj@meta.data[s_obj@meta.data$LOY=="LOY",]$PAR_score1 %>% mean()
  
  DE$Y_UMI <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$chrY_UMI %>% mean()
  DE$Y_genes <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$chrY_genes %>% mean()
  DE$Y_score <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$Y_score1 %>% mean()
  
  DE$sequencing_method2 <- paste(s_obj@meta.data$sequencing_method2 %>% unique(),sep = " ", collapse = " ")
  DE$tissue <- paste(s_obj@meta.data$tissue %>% unique(),sep = " ", collapse = " ")
  DE$specific_tissue <- paste(s_obj@meta.data$specific_tissue %>% unique(),sep = " ", collapse = " ")
  DE$diagnosis <- paste(s_obj@meta.data$diagnosis %>% unique(),sep = " ", collapse = " ")
  DE$neuro_degen <- paste(s_obj@meta.data$neuro_degen %>% unique(),sep = " ", collapse = " ")
  DE$neuro_degen_type <- paste(s_obj@meta.data$neuro_degen_type %>% unique(),sep = " ", collapse = " ")
  DE$smoking <- paste(s_obj@meta.data$smoking %>% unique(),sep = " ", collapse = " ")
  DE$donor_organism.age <- paste(s_obj@meta.data$donor_organism.age %>% unique(),sep = " ", collapse = " ")
  DE$sequencing_method <- paste(s_obj@meta.data$sequencing_method %>% unique(),sep = " ", collapse = " ")
  DE$cancer_diagnosis <- paste(s_obj@meta.data$diagnosis %>% unique(),sep = " ", collapse = " ")
  DE$donor_organism.sex <- paste(s_obj@meta.data$donor_organism.sex %>% unique(),sep = " ", collapse = " ")

  Y <- c("RPS4Y1","ZFY","LINC00278","PRKY","USP9Y","DDX3Y","UTY","NLGN4Y","TTTY14","EIF1AY","KDM5D")
  PAR <- c("GTPBP6","CSF2RA","IL3RA","SLC25A6","ASMTL","P2RY8","AKAP17A","DHRSX","CD99","ZBED1")
  OTHER <- c("ROBO1","OXR1","CACNA1A","FCGBP","RORA","SLIT1","SLIT2","ROBO2","ROBO3","ROBO4")
  
  ### pull data from dotplots for pct exp info + avg exp
  # add to the DE table for genes of interest 
  
  dat1 <- Seurat::DotPlot(s_obj, group.by = "LOY",
                          features = expressed_Y_genes(s_obj))
  
  dat2 <- Seurat::DotPlot(s_obj, group.by = "LOY",
                          features = expressed_PAR_genes(s_obj))
  
  dat3 <- Seurat::DotPlot(s_obj, group.by = "LOY",
                          features = unique(OTHER[OTHER %in% rownames(s_obj)]))
  
  
  tmp <- dat1$data[dat1$data$id == "NORMAL",]
  sum(tmp$pct.exp) -> DE$sum_Y_exp
  tmp <- dat2$data[dat2$data$id == "NORMAL",]
  sum(tmp$pct.exp) -> DE$sum_PAR_exp
  tmp <- dat3$data[dat3$data$id == "NORMAL",]
  sum(tmp$pct.exp) -> DE$sum_OTHER_exp
  
  dat1 <- dat1$data[dat1$data$id == "NORMAL" & dat1$data$features.plot %in% c(Y),]
  dat2 <- dat2$data[dat2$data$id == "NORMAL" & dat2$data$features.plot %in% c(PAR),]
  dat3 <- dat3$data[dat3$data$id == "NORMAL" & dat3$data$features.plot %in% c(OTHER),]
  
  rbind(dat1, dat2, dat3) -> dat
  
  
  lapply(c(Y,PAR,OTHER), function(x) {
    
    
    if(nrow(dat[dat$features.plot==x,]) < 1){ 
      data.frame(avg.exp = NA, pct.exp = NA) -> out 
      rownames(out) <- x
      names(out) <- paste0(names(out),"_",x)
    } else {
      dat[dat$features.plot==x,] %>% dplyr::select(avg.exp,pct.exp) %>% as.data.frame() -> out
      names(out) <- paste0(names(out),"_",x)
    }
    
    return(out)
    
  }) -> l
  do.call("cbind",l) -> l
  
  
  cbind(DE,l) -> DE
  
  
  return(DE)
}



########################################## 
run_DE_multi <- function(s_obj, sample, cell_type, mode = "LOY",
                   latent = c("nCount_RNA","nFeature_RNA","percent.rb"),
                   nUMI = 2500, nFeature = 1000, logfc.threshold = 0, min.pct = 0.05){
  tryCatch({
    message("sample: ", paste(sample, sep = ", ", collapse = " "), " ", cell_type, " ")
    vars <- paste0("nCount_RNA","percent.rb","nFeature_RNA")
    
    cells <- s_obj@meta.data[s_obj@meta.data$donor_organism.biomaterial_core.biomaterial_id %in% sample &
                               s_obj@meta.data$cell_ident_MV %in% cell_type &
                               s_obj@meta.data$nCount_RNA >= nUMI &
                               s_obj@meta.data$nFeature_RNA >= nFeature,] %>% rownames()
    
    subset(s_obj, cells = cells) -> s_obj
    
    
    s_obj@meta.data[s_obj@meta.data$LOY == "LOY",] %>% nrow() -> LOY
    s_obj@meta.data[s_obj@meta.data$LOY == "NORMAL",] %>% nrow() -> NORMAL
    s_obj@meta.data$cohort %>% unique() -> cohort
    
    if(LOY < 3 | NORMAL < 3){return(NULL)}
    
    if(mode == "LOY"){
      DE <- FindMarkers(s_obj, 
                        ident.1 = "LOY", ident.2 = "NORMAL", 
                        group.by = "LOY", 
                        test.use = "MAST", 
                        latent.vars = latent, 
                        logfc.threshold = logfc.threshold)
    } else {
      DE <- FindMarkers(s_obj, 
                        ident.1 = "0", ident.2 = "4", 
                        group.by = "LOY_levels", 
                        test.use = "MAST", 
                        latent.vars = latent, 
                        logfc.threshold = logfc.threshold)
    }
    
    rownames(DE) -> DE$gene
    DE$nUMI <- nUMI
    DE$nFeature <- nFeature
    DE$cell_type <- cell_type
    DE$mode <- mode
    DE$sample <- "all"
    DE$cohort <- cohort
    DE$LOY_cells <- s_obj@meta.data[s_obj@meta.data$LOY=="LOY",] %>% nrow()
    DE$NORMAL_cells <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",] %>% nrow()
    DE$LEVEL2_cells <- s_obj@meta.data[s_obj@meta.data$LOY_levels=="2",] %>% nrow()
    DE$LEVEL3_cells <- s_obj@meta.data[s_obj@meta.data$LOY_levels=="3",] %>% nrow()
    DE$LEVEL4_cells <- s_obj@meta.data[s_obj@meta.data$LOY_levels=="4",] %>% nrow()
    
    
    DE$median_nUMI_LOY4 <- s_obj@meta.data[s_obj@meta.data$LOY_levels=="4",]$nCount_RNA %>% median()
    DE$median_nUMI_LOY <- s_obj@meta.data[s_obj@meta.data$LOY_levels=="0",]$nCount_RNA %>% median()
    DE$median_nUMI_NORMAL <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$nCount_RNA %>% median()
    
    DE$median_nFeature_LOY4 <- s_obj@meta.data[s_obj@meta.data$LOY_levels=="4",]$nFeature_RNA %>% median()
    DE$median_nFeature_LOY <- s_obj@meta.data[s_obj@meta.data$LOY_levels=="0",]$nFeature_RNA %>% median()
    DE$median_nFeature_NORMAL <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$nFeature_RNA %>% median()
    
    DE$expr_Y <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$normY %>% median()
    DE$expr_PAR <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$normPAR %>% median()
    DE$expr_PAR_LOY <- s_obj@meta.data[s_obj@meta.data$LOY=="LOY",]$normPAR %>% median()
    
    DE$PAR_UMI <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$PAR_UMI %>% mean()
    DE$PAR_genes <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$PAR_genes %>% mean()
    DE$PAR_score <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$PAR_score1 %>% mean()
    DE$PAR_score_LOY <- s_obj@meta.data[s_obj@meta.data$LOY=="LOY",]$PAR_score1 %>% mean()
    
    DE$Y_UMI <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$chrY_UMI %>% mean()
    DE$Y_genes <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$chrY_genes %>% mean()
    DE$Y_score <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$Y_score1 %>% mean()
    
    Y <- c("RPS4Y1","ZFY","LINC00278","PRKY","USP9Y","DDX3Y","UTY","NLGN4Y","TTTY14","EIF1AY","KDM5D")
    PAR <- c("GTPBP6","CSF2RA","IL3RA","SLC25A6","ASMTL","P2RY8","AKAP17A","DHRSX","CD99","ZBED1")
    OTHER <- c("ROBO1","OXR1","CACNA1A","FCGBP","RORA","SLIT1","SLIT2","ROBO2","ROBO3","ROBO4")
    
    ### pull data from dotplots for pct exp info + avg exp
    # add to the DE table for genes of interest 
    
    
      
      dat1 <- Seurat::DotPlot(s_obj, group.by = "LOY",
                            features = expressed_Y_genes(s_obj))
    
      dat2 <- Seurat::DotPlot(s_obj, group.by = "LOY",
                            features = expressed_PAR_genes(s_obj))
    
      dat3 <- Seurat::DotPlot(s_obj, group.by = "LOY",
                            features = unique(OTHER[OTHER %in% rownames(s_obj)]))
    
    
    tmp <- dat1$data[dat1$data$id == "NORMAL",]
    sum(tmp$pct.exp) -> DE$sum_Y_exp
    tmp <- dat2$data[dat2$data$id == "NORMAL",]
    sum(tmp$pct.exp) -> DE$sum_PAR_exp
    tmp <- dat3$data[dat3$data$id == "NORMAL",]
    sum(tmp$pct.exp) -> DE$sum_OTHER_exp
    
    dat1 <- dat1$data[dat1$data$id == "NORMAL" & dat1$data$features.plot %in% c(Y),]
    dat2 <- dat2$data[dat2$data$id == "NORMAL" & dat2$data$features.plot %in% c(PAR),]
    dat3 <- dat3$data[dat3$data$id == "NORMAL" & dat3$data$features.plot %in% c(OTHER),]
    
    rbind(dat1, dat2, dat3) -> dat
    
    
    lapply(c(Y,PAR,OTHER), function(x) {
      
      
      if(nrow(dat[dat$features.plot==x,]) < 1){ 
        data.frame(avg.exp = NA, pct.exp = NA) -> out 
        rownames(out) <- x
        names(out) <- paste0(names(out),"_",x)
      } else {
        dat[dat$features.plot==x,] %>% dplyr::select(avg.exp,pct.exp) %>% as.data.frame() -> out
        names(out) <- paste0(names(out),"_",x)
      }
      
      return(out)
      
    }) -> l
    do.call("cbind",l) -> l
    
    
    cbind(DE,l) -> DE
    
    p.adjust(DE$p_val, method = "fdr") -> DE$FDR
    p.adjust(DE$p_val, method = "BH") -> DE$bh
    p.adjust(DE$p_val, method = "hochberg") -> DE$hoch
    p.adjust(DE$p_val, method = "bonferroni") -> DE$bonf
    
    return(DE)
  }, error=function(e){message("Error"); return(NULL)})
} 


run_DE_multi_simple <- function(s_obj, sample, cell_type, mode = "LOY",
                         nUMI = 2500, nFeature = 1000, logfc.threshold = 0, min.pct = 0.05){
  tryCatch({
    message(paste(sample, sep = " ,"), " | ", cell_type, " ")
    cells <- s_obj@meta.data[s_obj@meta.data$donor_organism.biomaterial_core.biomaterial_id %in% sample &
                               s_obj@meta.data$cell_ident_MV %in% cell_type &
                               s_obj@meta.data$nCount_RNA >= nUMI &
                               s_obj@meta.data$nFeature_RNA >= nFeature,] %>% rownames()
    
    subset(s_obj, cells = cells) -> s_obj
    
    s_obj@meta.data$cohort %>% unique() -> cohort
    
    if(mode == "LOY"){
      DE <- FindMarkers(s_obj, 
                        ident.1 = "LOY", ident.2 = "NORMAL", 
                        group.by = "LOY", 
                        test.use = "MAST", 
                        latent.vars = c("nCount_RNA","percent.rb","nFeature_RNA"), 
                        logfc.threshold = logfc.threshold)
    } else {
      DE <- FindMarkers(s_obj, 
                        ident.1 = "0", ident.2 = "4", 
                        group.by = "LOY_levels", 
                        test.use = "MAST", 
                        latent.vars = c("nCount_RNA","percent.rb","nFeature_RNA"), 
                        logfc.threshold = logfc.threshold)
    }
    
    rownames(DE) -> DE$gene
    DE$nUMI <- nUMI
    DE$nFeature <- nFeature
    DE$cell_type <- cell_type
    DE$mode <- mode
    DE$sample <- "all"
    DE$cohort <- cohort
    DE$LOY_cells <- s_obj@meta.data[s_obj@meta.data$LOY=="LOY",] %>% nrow()
    DE$NORMAL_cells <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",] %>% nrow()
    DE$LEVEL2_cells <- s_obj@meta.data[s_obj@meta.data$LOY_levels=="2",] %>% nrow()
    DE$LEVEL3_cells <- s_obj@meta.data[s_obj@meta.data$LOY_levels=="3",] %>% nrow()
    DE$LEVEL4_cells <- s_obj@meta.data[s_obj@meta.data$LOY_levels=="4",] %>% nrow()
    
    
    DE$median_nUMI_LOY4 <- s_obj@meta.data[s_obj@meta.data$LOY_levels=="4",]$nCount_RNA %>% median()
    DE$median_nUMI_LOY <- s_obj@meta.data[s_obj@meta.data$LOY_levels=="0",]$nCount_RNA %>% median()
    DE$median_nUMI_NORMAL <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$nCount_RNA %>% median()
    
    DE$median_nFeature_LOY4 <- s_obj@meta.data[s_obj@meta.data$LOY_levels=="4",]$nFeature_RNA %>% median()
    DE$median_nFeature_LOY <- s_obj@meta.data[s_obj@meta.data$LOY_levels=="0",]$nFeature_RNA %>% median()
    DE$median_nFeature_NORMAL <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$nFeature_RNA %>% median()
    
    DE$expr_Y <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$normY %>% median()
    DE$expr_PAR <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$normPAR %>% median()
    DE$expr_PAR_LOY <- s_obj@meta.data[s_obj@meta.data$LOY=="LOY",]$normPAR %>% median()
    
    DE$PAR_UMI <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$PAR_UMI %>% mean()
    DE$PAR_genes <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$PAR_genes %>% mean()
    DE$PAR_score <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$PAR_score1 %>% mean()
    DE$PAR_score_LOY <- s_obj@meta.data[s_obj@meta.data$LOY=="LOY",]$PAR_score1 %>% mean()
    
    DE$Y_UMI <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$chrY_UMI %>% mean()
    DE$Y_genes <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$chrY_genes %>% mean()
    DE$Y_score <- s_obj@meta.data[s_obj@meta.data$LOY=="NORMAL",]$Y_score1 %>% mean()
    
    p.adjust(DE$p_val, method = "fdr") -> DE$FDR
    p.adjust(DE$p_val, method = "BH") -> DE$bh
    p.adjust(DE$p_val, method = "hochberg") -> DE$hoch
    p.adjust(DE$p_val, method = "bonferroni") -> DE$bonf
    
    return(DE)
  }, error=function(e){message("Error"); return(NULL)})
} 