

################### 
PAR_scLOY_genes <- readRDS(url("https://github.com/michaelcvermeulen/microglia-loss-of-y/blob/main/data/PAR_scLOY_genes.RDS?raw=true"))
Y_scLOY_genes <- readRDS(url("https://github.com/michaelcvermeulen/microglia-loss-of-y/blob/main/data/Y_scLOY_genes.RDS?raw=true"))
X_scLOY_genes <- readRDS(url("https://github.com/michaelcvermeulen/microglia-loss-of-y/blob/main/data/X_scLOY_genes.RDS?raw=true"))




##########################
call_LOY <- function(s_obj,
                     genes_to_exclude = c("LINC00266-4P","PRORY","AC022486.1",
                                                 "TTTY4C","AC012078.2",
                                                 "TTTY4B","TBL1Y","AMELY",
                                                 "TTTY10","TTTY2B",
                                                 "HSFY2","TTTY22",
                                                 "HSFY1",
                                                 "NLGN4Y-AS1",
                                                 "ZFY-AS1",
                                                 "PCDH11Y",
                                                 "MAFIP")){


  suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
  edb <- EnsDb.Hsapiens.v86


  ## Query PAR genes
  grf <- GRangesFilter(GRanges("X", ranges = IRanges(10001, 2781479)), type = "any")
  gn <- as.data.frame(genes(edb, filter = grf))
  dplyr::select(gn, gene_id, strand, gene_name, gene_biotype) -> PAR_scLOY_genes

  grf2 <- GRangesFilter(GRanges("X", ranges = IRanges(155701383, 156030895)), type = "any")
  gn2 <- as.data.frame(genes(edb, filter = grf2))
  dplyr::select(gn, gene_id, strand, gene_name, gene_biotype) -> PAR2_scLOY_genes

  rbind(PAR_scLOY_genes,PAR2_scLOY_genes) -> PAR_scLOY_genes
  PAR_scLOY_genes[!(PAR_scLOY_genes$gene_biotype %in% c("LRG_gene")),] -> PAR_scLOY_genes

  ## Query chrY genes
  grf <- GRangesFilter(GRanges("Y", ranges = IRanges(1, 59373566)), type = "any")
  gn <- as.data.frame(genes(edb, filter = grf))
  dplyr::select(gn, gene_id, strand, gene_name, gene_biotype) -> Y_scLOY_genes

  message("calling LOY")
  PAR_scLOY_genes$gene_name -> PAR_genes
  unique(PAR_genes[PAR_genes %in% rownames(s_obj)]) -> PAR_genes

  Y_scLOY_genes$gene_name -> Y_genes
  Y_genes[Y_genes %in% rownames(s_obj)] -> Y_genes
  if(!(is.null(genes_to_exclude))){Y_genes <- Y_genes[!(Y_genes %in% genes_to_exclude)]; print(paste0("removed gene: ", genes_to_exclude))}

  print(paste0(length(Y_genes)," chromosome Y genes detected!"))
  print(paste0(Y_genes))

  s_obj@meta.data -> df
  s_obj@assays$RNA@counts -> counts
  counts[Y_genes,,drop=F] -> counts_chrY
  counts[PAR_genes,,drop=F] -> counts_PAR
  s_obj@assays$RNA@data -> norm
  norm[Y_genes,,drop=F] -> norm_chrY

  as.data.frame(colMeans(as.matrix(norm_chrY))) -> chrY_norm; names(chrY_norm) <- "chrY_norm"
  as.data.frame(colSums(as.matrix(counts_chrY))) -> chrY_UMI; names(chrY_UMI) <- "chrY_UMI"
  as.data.frame(colSums(as.matrix(counts_chrY > 0))) ->  chrY_genes; names(chrY_genes) <- "chrY_genes"
  as.data.frame(colSums(as.matrix(counts_PAR))) -> PAR_UMI;  names(PAR_UMI) <- "PAR_UMI"
  as.data.frame(colSums(as.matrix(counts_chrY > 0))) ->  PAR_genes; names(PAR_genes) <- "PAR_genes"

  df$chrY_genes <- chrY_genes$chrY_genes
  df$chrY_UMI <- chrY_UMI$chrY_UMI
  df$chrY_norm <- chrY_norm$chrY_norm
  df$PAR_UMI <- PAR_UMI$PAR_UMI
  df$PAR_genes <- PAR_genes$PAR_genes


  tibble::rownames_to_column(df, "CB") -> df
  dplyr::mutate(df, LOY = ifelse(test = chrY_UMI == 0 , yes = "LOY", no = "NORMAL"))  -> df
  dplyr::mutate(df, LOY_levels = ifelse(test = chrY_UMI > 2 & chrY_genes > 2 , yes = 4,
                                        no = ifelse(test = chrY_UMI > 1 & chrY_genes > 1, yes = 3,
                                                    no = ifelse(test = chrY_UMI > 1 | chrY_genes > 1, yes = 2,
                                                                no = ifelse(test = chrY_UMI >0, yes = 1, no = 0)))))  -> df
  ifelse(df$LOY_levels == 0, yes = 0, no = ifelse(df$LOY_levels %in% c(1,2), yes = 1, no = 2) ) -> df$LOY_levels_condensed

  dplyr::mutate(df, normY = log1p((chrY_UMI / nCount_RNA) * 10000)) -> df
  dplyr::mutate(df, normPAR = log1p((PAR_UMI / nCount_RNA) * 10000)) -> df


  tibble::column_to_rownames(df, "CB") -> df
  s_obj@meta.data <- df

  return(s_obj)

}




###############################
LOY_genes_to_exclude_by_sex <- function(s_obj, plot=T, ratio=5){
  
  # pull normalized matrix
  s_obj[["RNA"]]@data -> mat
  Y <- unique(Y_scLOY_genes[Y_scLOY_genes$gene_name %in% rownames(mat),]$gene_name)
  mat[Y,] -> mat
  mat %>% as.matrix() %>% as.data.frame() -> mat
  mat %>% t() -> mat
  reshape2::melt(mat) -> mat.melt
  names(mat.melt) <- c("CB","Gene","nUMI") 
  
  dplyr::left_join(x = mat.melt , y = tibble::rownames_to_column(s_obj@meta.data,"CB"), by = "CB") -> t
  t %>% dplyr::select(CB, Gene, nUMI,
                      donor_organism.biomaterial_core.biomaterial_id,
                      donor_organism.sex) -> t
  
  t %>% dplyr::group_by(Gene,donor_organism.biomaterial_core.biomaterial_id) %>% dplyr::summarize(mean_UMI = mean(nUMI)) -> tmp
  reshape2::dcast(tmp, Gene ~ donor_organism.biomaterial_core.biomaterial_id, value.var = "mean_UMI") -> tmp
  tmp %>% tibble::column_to_rownames("Gene") -> tmp
  
  annotation_c <- data.frame(names = names(tmp))
  merge.data.frame(x = annotation_c, y = dplyr::distinct(s_obj@meta.data,donor_organism.biomaterial_core.biomaterial_id, .keep_all = T), 
                   by.x = "names", by.y = "donor_organism.biomaterial_core.biomaterial_id") -> annotation_c
  annotation_c[order(annotation_c$donor_organism.sex),] -> annotation_c
  rownames(annotation_c) <- NULL
  as.data.frame(annotation_c) %>% tibble::column_to_rownames("names") -> annotation_c
  
  dplyr::select(annotation_c, donor_organism.sex) -> annotation_c
  names(annotation_c) <- "sex"
  
  
  my_colour = list(
    sex = c(male = "#5977ff", female = "#f74747")
  )
  
  
  if(plot == T){
    pheatmap::pheatmap(mat = tmp, main = "Mean UMIs detected per cell for each chromosome Y gene",
                       display_numbers = T,
                       cluster_cols = F, 
                       cluster_rows = F, 
                       color = viridis(15,end = 0.8), 
                       annotation_col = annotation_c, 
                       annotation_colors = my_colour,
                       show_colnames = F, number_color = "white", 
                       
    ) 
  }
  
  
  # sex based summary
  
  t %>% dplyr::group_by(Gene,donor_organism.sex) %>% dplyr::summarize(mean_UMI = mean(nUMI)) -> tmp
  reshape2::dcast(tmp, Gene ~ donor_organism.sex, value.var = "mean_UMI") -> tmp
  dplyr::mutate(tmp, ratio = (male / female)) -> df
  
  df[df$ratio < ratio,]$Gene %>% as.character() %>% na.omit() -> genes
  
  
  print(paste0("Unreliable genes for Y detection: ", genes))
  
  return(as.vector(genes))
  
  
}


########################
LOY_genes_to_exclude_expr <- function(s_obj, pct.exp = 0.05, avg.exp = 0.05){
  # pull normalized matrix
  s_obj[["RNA"]]@data -> mat
  Y <- unique(Y_scLOY_genes[Y_scLOY_genes$gene_name %in% rownames(mat),]$gene_name)
  
  Idents(s_obj) <- "donor_organism.sex"
  DotPlot(s_obj, features = Y) -> p
  
  p$data[p$data$id == "male", ] -> dat
  
  dat[dat$pct.exp < pct.exp | dat$avg.exp < avg.exp,]$features.plot %>% as.character() %>% na.omit() -> genes
  
  print(paste0("Unreliable genes for Y detection: ", genes))
  
  return(as.vector(genes))
  
}


########################
check_object <- function(s_obj, return = F){
  error<-0
  
  c("cohort", "nCount_RNA", "nFeature_RNA", "percent.mt", 
    "percent.rb", "PAR_score1", "Y_score1", "cell_ident_MV", "cell_broad_MV", 
    "donor_organism.biomaterial_core.biomaterial_id", "donor_organism.sex",
    "donor_organism.age", "tissue","specific_tissue","GEO",
    "chrY_genes", "chrY_UMI", "PAR_UMI", "PAR_genes", "LOY", "normY", "normPAR", 
    "PAR_score_cell_ident_MV", "Y_score_cell_ident_MV") -> required_columns
  
  tryCatch(
    expr = {
      df <- s_obj@meta.data %>% dplyr::select(required_columns)
    },
    error = function(e){ 
      print("ERROR: structure of meta.data is incorrect")
      
      required_columns[!(required_columns %in% names(s_obj@meta.data) )] -> missing
      required_columns[(required_columns %in% names(s_obj@meta.data) )] -> present
      
      print(paste0("missing: ", missing))
      print(paste0("present: ", present))
      
      error<-1
      return(NULL)
    },
    warning=function(cond) {
      
      # Choose a return value in case of warning
      return(NULL)
    },
    finally={
      # NOTE:
      # Here goes everything that should be executed at the end,
      # regardless of success or error.
      # If you want more than one expression to be executed, then you 
      # need to wrap them in curly brackets ({...}); otherwise you could
      # just have written 'finally=<expression>' 
      
    }
    
  )
  
  if(error==0){print("Object looks good")
    if(isTRUE(return)){
      print("Returning formatted data.frame")
      return(df %>% as.data.frame())
    }}
  
}

########################
check_object_meta <- function(meta, return = F){
  error<-0
  
  c("cohort", "nCount_RNA", "nFeature_RNA", "percent.mt", 
    "percent.rb", "PAR_score1", "Y_score1", "cell_ident_MV", "cell_broad_MV", 
    "donor_organism.biomaterial_core.biomaterial_id", "donor_organism.sex",
    "donor_organism.age","tissue","specific_tissue","GEO", "sequencing_method2",
    "chrY_genes", "chrY_UMI", "PAR_UMI", "PAR_genes", "LOY", "normY", "normPAR", 
    "PAR_score_cell_ident_MV", "Y_score_cell_ident_MV",
    "diagnosis","cancer_diagnosis","neuro_degen","neuro_degen_type","smoking",
    "sequencing_method","introns") -> required_columns
  
  tryCatch(
    expr = {
      df <- meta %>% dplyr::select(required_columns)
    },
    error = function(e){ 
      print("ERROR: structure of meta.data is incorrect")
      
      required_columns[!(required_columns %in% names(meta) )] -> missing
      required_columns[(required_columns %in% names(meta) )] -> present
      
      print(paste0("missing: ", missing))
      print(paste0("present: ", present))
      
      error<-1
      return(NULL)
    },
    warning=function(cond) {
      
      # Choose a return value in case of warning
      return(NULL)
    },
    finally={
      # NOTE:
      # Here goes everything that should be executed at the end,
      # regardless of success or error.
      # If you want more than one expression to be executed, then you 
      # need to wrap them in curly brackets ({...}); otherwise you could
      # just have written 'finally=<expression>' 
      
    }
    
  )
  
  if(error==0){print("Object looks good")
    if(isTRUE(return)){
      print("Returning formatted data.frame")
      return(df %>% as.data.frame())
    }}
  
}

########################
pull_genes_of_interest <- function(s_obj, features){
  GetAssayData(o) -> mat
  mat[features,] %>% as.data.frame() -> mat
  return(mat)
}


########################
module_score_by_cluster <- function(s_obj, ident){
  
  if(!(ident %in% names(s_obj@meta.data))){print("Ident doesnt exist"); return(s_obj)}
  if(!(c("donor_organism.sex") %in% names(s_obj@meta.data))){print("Donor_organism.sex doesnt exist. Run call_sex() first"); return(s_obj)}
  
  s_obj@meta.data %>% dplyr::select(ident) -> idents
  rownames(idents) -> idents$CB
  
  # run male, then female
  cells <- s_obj@meta.data[s_obj@meta.data$donor_organism.sex == "male",] %>% rownames()
  if(length(cells) == 0){s_obj.m <- NULL} else(s_obj.m <- subset(s_obj, cells = cells))
  
  cells <- s_obj@meta.data[s_obj@meta.data$donor_organism.sex == "female",] %>% rownames()
  if(length(cells) == 0){s_obj.f <- NULL} else(s_obj.f <- subset(s_obj, cells = cells))
  
  if(is.null(s_obj.m)){lm <- NULL} else {
    lapply(idents[,1] %>% unique(), function(x) {
      
      message("Working on ", x)
      cells <- idents[idents[,1]==x,]$CB
      cells <-  cells[cells %in% rownames(s_obj.m@meta.data)]
      if(length(cells) == 0){return()}
      s_obj. <- subset(s_obj.m, cells = cells)
      
      list(PAR_scLOY_genes$gene_name[PAR_scLOY_genes$gene_name %in% rownames(s_obj.)]) -> P
      Seurat::AddModuleScore(object = s_obj., features = P, name = "PAR_score") -> s_obj.
      list(Y_scLOY_genes$gene_name[Y_scLOY_genes$gene_name %in% rownames(s_obj.)]) -> Y
      Seurat::AddModuleScore(object = s_obj., features = Y, name = "Y_score", ) -> s_obj.
      
      ## return cell_id and scores 
      out <- s_obj.@meta.data %>% tibble::rownames_to_column("CB") %>% dplyr::select(CB,PAR_score1,Y_score1)
      names(out) <- c("CB",paste0("PAR_score_",ident),paste0("Y_score_",ident))
      
      out <- out %>% as.data.frame()
      return(out)
    }) -> l
    do.call("rbind",l) -> lm }
  
  if(is.null(s_obj.f)){lf <- NULL} else {
    lapply(idents[,1] %>% unique(), function(x) {
      
      message("Working on ", x)
      cells <- idents[idents[,1]==x,]$CB
      cells <-  cells[cells %in% rownames(s_obj.f@meta.data)]
      if(length(cells) == 0){return()}
      s_obj. <- subset(s_obj.f, cells = cells)
      
      list(PAR_scLOY_genes$gene_name[PAR_scLOY_genes$gene_name %in% rownames(s_obj.)]) -> P
      Seurat::AddModuleScore(object = s_obj., features = P, name = "PAR_score") -> s_obj.
      list(Y_scLOY_genes$gene_name[Y_scLOY_genes$gene_name %in% rownames(s_obj.)]) -> Y
      Seurat::AddModuleScore(object = s_obj., features = Y, name = "Y_score", ) -> s_obj.
      
      ## return cell_id and scores 
      out <- s_obj.@meta.data %>% tibble::rownames_to_column("CB") %>% dplyr::select(CB,PAR_score1,Y_score1)
      names(out) <- c("CB",paste0("PAR_score_",ident),paste0("Y_score_",ident))
      
      out <- out %>% as.data.frame()
      return(out)
    }) -> l
    do.call("rbind",l) -> lf }
  
  rbind(lm,lf) -> l
  
  ## attach scores to s_obj and return
  dplyr::left_join(x = tibble::rownames_to_column(s_obj@meta.data,"CB"), y = l, by = "CB") %>% 
    tibble::column_to_rownames("CB") -> s_obj@meta.data
  
  message("Added ",paste0(names(l)[-1], collapse = " - ")," to Seurat object")
  return(s_obj)
}


########################
call_sample_sex <- function(s_obj, overwrite = F, ratio = 0.4){
  Idents(s_obj) <- "donor_organism.biomaterial_core.biomaterial_id"
  
  
  if(c("donor_organism.sex") %in% names(s_obj@meta.data) & overwrite == F){print("sample sex meta is already available in donor_organism.sex")
    AverageExpression(s_obj, features = c("RPS4Y1","XIST","UTY")) -> expr
    expr$RNA %>% t() %>% as.data.frame() -> expr
    
    ifelse(test = ((expr$RPS4Y1 / expr$XIST) > ratio), "male", "female") -> expr$donor_organism.sex 
    print(expr)
    
    return(s_obj)
  }
  
  
  AverageExpression(s_obj, features = c("RPS4Y1","XIST","UTY"), group.by = "donor_organism.biomaterial_core.biomaterial_id") -> expr
  expr$RNA %>% t() %>% as.data.frame() -> expr
  
  ifelse(test = ((expr$RPS4Y1 + expr$UTY) / expr$XIST) > ratio, "male", "female") -> expr$donor_organism.sex 
  
  tibble::rownames_to_column(expr,"sample") %>% dplyr::select(sample,donor_organism.sex) -> expr
  
  s_obj@meta.data -> df
  
  
  dplyr::left_join(x = tibble::rownames_to_column(df,"CB"),
                   y = expr, by = c("donor_organism.biomaterial_core.biomaterial_id"="sample")
  ) -> df
  
  tibble::column_to_rownames(df, "CB") -> s_obj@meta.data
  print(expr)
  return(s_obj)
  
}


########################
volcano <- function(table, include_Y=F){
  
  PAR_scLOY_genes[PAR_scLOY_genes$gene_name %in% rownames(table),]$gene_name -> n
  rownames(table) -> table$gene
  
  if(isFALSE(include_Y)){
    table[!(table$gene %in% Y_scLOY_genes$gene_name),] -> table
  }
  
  EnhancedVolcano::EnhancedVolcano(toptable = table, gridlines.major = F, gridlines.minor = F, 
                                   lab = rownames(table), FCcutoff = 0.25, labSize = 4, border = 'full',
                                   borderWidth = 0, col=c('grey', 'grey', 'grey', 'red'),
                                   borderColour = 'black', colAlpha = 1, legendLabels = "",
                                   x = "avg_log2FC", xlim = c(min(table$avg_log2FC)-0.1,max(table$avg_log2FC)+0.1),
                                   title = "", subtitle = "", 
                                   y = "p_val") -> p
  p
  return(p)
  
}


########################
dataset_summary <- function(s_obj, sample_size_min = 20){
  
  # plot the number of LOY cells in each cell type for each sample 
  s_obj@meta.data -> df
  
  df %>% dplyr::group_by(donor_organism.biomaterial_core.biomaterial_id, LOY, cell_ident_MV, donor_organism.sex) %>% 
    dplyr::summarize(counts = dplyr::n()) %>% 
    reshape2::dcast(donor_organism.biomaterial_core.biomaterial_id + cell_ident_MV + donor_organism.sex ~ LOY, value.var = "counts") -> tmp
  
  if(tmp[is.na(tmp$LOY),]$LOY %>% length() > 0){tmp[is.na(tmp$LOY),]$LOY <- 0} 
  if(tmp[is.na(tmp$NORMAL),]$NORMAL %>% length() > 0){tmp[is.na(tmp$NORMAL),]$NORMAL <- 0} 
  
  dplyr::mutate(tmp, LOY_prop = LOY / (LOY + NORMAL), TOTAL = LOY + NORMAL) -> tmp
  
  ## remove female samples
  tmp[tmp$donor_organism.sex=="male",] -> tmp2
  
  ## plot LOY for each cell type across all individuals
  # order celltypes by prop 
  factor(tmp2$cell_ident_MV, levels = tmp2[order(-tmp2$LOY_prop),]$cell_ident_MV %>% unique()) -> tmp2$cell_ident_MV
  tmp2[tmp2$TOTAL<=sample_size_min,]$LOY_prop <- 0
  
  ggbarplot(data = tmp2, x = "cell_ident_MV",
            y = "LOY_prop", position = position_dodge(0.6), palette = "aaas",
            fill = "donor_organism.biomaterial_core.biomaterial_id") %>% ggpar(x.text.angle = 45, xlab = "", ylab = "Loss of Y proportion",
                                                                               legend.title = "Sample") -> p1
  
  # plot LOY for each sample
  tmp[tmp$donor_organism.sex=="male",] -> tmp2
  tmp2 %>% group_by(donor_organism.biomaterial_core.biomaterial_id) %>% summarize(LOY = sum(LOY) / sum(TOTAL)) -> total
  ggbarplot(data = total, x = "donor_organism.biomaterial_core.biomaterial_id", fill = "donor_organism.biomaterial_core.biomaterial_id",
            y = "LOY", palette = "aaas") %>% ggpar(x.text.angle = 45, xlab = "", ylab = "Loss of Y proportion",legend = "none") -> p2
  
  
  ## cell_type proportions by sample 
  df %>% group_by(donor_organism.biomaterial_core.biomaterial_id, cell_ident_MV) %>% 
    summarise(total = dplyr::n()) -> a
  
  df %>% group_by(donor_organism.biomaterial_core.biomaterial_id) %>% summarise(total_cells = dplyr::n()) -> totals
  
  dplyr::left_join(x =a , y = totals, by = "donor_organism.biomaterial_core.biomaterial_id") ->a
  dplyr::mutate(a, prop = (total / total_cells)) -> a
  
  ggbarplot(data = a, x = "donor_organism.biomaterial_core.biomaterial_id", y = "prop", fill = "cell_ident_MV", 
            palette = Seurat::DiscretePalette(palette = "stepped", n = length(a$cell_ident_MV))) %>% 
    ggpar(legend = "bottom", xlab = "", ylab = "Cell-type proportion", x.text.angle = 45, legend.title = "Celltype annotations") -> p3
  
  ggarrange(plotlist = list(p1,p2,p3), widths = c(2,1)) -> p
  
  return(p) 
  
}

########################
run_DE <- function(s_obj){
  
  s_obj@meta.data -> df
  
  obj <- list()
  n <- NULL
  for(i in df[df$donor_organism.sex=="male",]$donor_organism.biomaterial_core.biomaterial_id %>% unique()){
    for(j in df$cell_ident_MV %>% unique()){
      
      print(i)
      print(j)
      LOY_DE(s_obj = s_obj, cell_type = j, subject = i, top = "95%", bottom = "5%") -> d
      
      if(is.null(d)){next}
      print(head(d,10))
      
      rlist::list.append(obj,d) -> obj
      n <- c(n,paste0(i,"_",j))
      
    }
  }
  obj <- list()
  n <- NULL
  for(i in df[df$donor_organism.sex=="male",]$donor_organism.biomaterial_core.biomaterial_id %>% unique()){
    for(j in df$cell_ident_MV %>% unique()){
      
      print(i)
      print(j)
      LOY_DE(s_obj = s_obj, cell_type = j, subject = i, top = "95%", bottom = "5%") -> d
      
      if(is.null(d)){next}
      print(head(d,10))
      
      rlist::list.append(obj,d) -> obj
      n <- c(n,paste0(i,"_",j))
      
    }
  }
  names(obj) <- n
  return(obj)
  
}


########################
plot_expression <- function(s_obj, plot_Y=F, plot_PAR=F){
  
  # pull normalized matrix
  s_obj[["RNA"]]@data -> mat
  
  if(plot_Y==T & plot_PAR==T){
    Y <- unique(Y_scLOY_genes[Y_scLOY_genes$gene_name %in% rownames(mat),]$gene_name)
    P <- unique(PAR_scLOY_genes[PAR_scLOY_genes$gene_name %in% rownames(mat),]$gene_name)
    mat[c(Y,P),] -> mat
  }
  if(plot_Y==T & plot_PAR==F){
    Y <- unique(Y_scLOY_genes[Y_scLOY_genes$gene_name %in% rownames(mat),]$gene_name)
    mat[Y,] -> mat
  }
  if(plot_Y==F & plot_PAR==T){
    P <- unique(PAR_scLOY_genes[PAR_scLOY_genes$gene_name %in% rownames(mat),]$gene_name)
    mat[P,] -> mat
  }
  if(plot_Y==F & plot_PAR==F){print("must pick genes to plot"); return(NULL)}
  
  mat %>% as.matrix() %>% as.data.frame() -> mat
  mat %>% t() -> mat
  reshape2::melt(mat) -> mat.melt
  names(mat.melt) <- c("CB","Gene","nUMI") 
  
  dplyr::left_join(x = mat.melt , y = tibble::rownames_to_column(s_obj@meta.data,"CB"), by = "CB") -> t
  t %>% dplyr::select(CB, Gene, nUMI,
                      donor_organism.biomaterial_core.biomaterial_id,
                      donor_organism.sex) -> t
  
  t %>% group_by(Gene,donor_organism.biomaterial_core.biomaterial_id) %>% summarize(mean_UMI = mean(nUMI)) -> tmp
  reshape2::dcast(tmp, Gene ~ donor_organism.biomaterial_core.biomaterial_id, value.var = "mean_UMI") -> tmp
  tmp %>% tibble::column_to_rownames("Gene") -> tmp
  
  annotation_c <- data.frame(names = names(tmp))
  merge.data.frame(x = annotation_c, y = dplyr::distinct(s_obj@meta.data,donor_organism.biomaterial_core.biomaterial_id, .keep_all = T), 
                   by.x = "names", by.y = "donor_organism.biomaterial_core.biomaterial_id") -> annotation_c
  annotation_c[order(annotation_c$donor_organism.sex),] -> annotation_c
  rownames(annotation_c) <- NULL
  as.data.frame(annotation_c) %>% tibble::column_to_rownames("names") -> annotation_c
  
  dplyr::select(annotation_c, donor_organism.sex) -> annotation_c
  names(annotation_c) <- "sex"
  
  
  my_colour = list(
    sex = c(male = "#5977ff", female = "#f74747")
  )
  
  
  
  pheatmap::pheatmap(mat = tmp, main = "",
                     display_numbers = T,
                     cluster_cols = F, 
                     cluster_rows = T, 
                     color = viridis(15,end = 0.8), 
                     annotation_col = annotation_c, 
                     annotation_colors = my_colour,
                     show_colnames = T, number_color = "white")
  
}



########################
Marker_heatmap2 <- function(expr, gene, order){
  
  gene <- gene[gene %in% names(expr)]
  expr <- expr[, c(gene, "label")]
  type_expr <- expr %>% tidyr::nest(-label) %>% dplyr::rename(expr = data) %>% 
    dplyr::mutate(colmeans = purrr::map(.x = expr, .f = function(.x) {
      colMeans(.x)
    }))
  type_mean_expr <- type_expr$colmeans %>% as.data.frame() %>% 
    tibble::remove_rownames() %>% t() %>% as.data.frame() %>% 
    tibble::remove_rownames()
  rownames(type_mean_expr) <- type_expr$label
  colnames(type_mean_expr) <- colnames(expr)[-ncol(expr)]
  sub_expr <- type_mean_expr
  sub_expr <- sub_expr %>% as.tibble() %>% dplyr::mutate_all(funs((. - 
                                                                     mean(.))/sd(.))) %>% t()
  colnames(sub_expr) <- type_expr$label
  get_label <- function(num) {
    v <- sub_expr[num, ]
    colnames(sub_expr)[which(v == max(v))]
  }
  sub_expr <- sub_expr %>% tibble::as.tibble() %>% dplyr::mutate(group = purrr::map_chr(1:length(gene), 
                                                                                        get_label))
  sub_expr <- as.data.frame(sub_expr)
  rownames(sub_expr) <- gene
  sub_expr <- sub_expr %>% dplyr::mutate(gene = gene) %>% tidyr::gather(key = "cell_type", 
                                                                        value = "zscore", -group, -gene) 
  
  
  factor(sub_expr$group, levels = ord, ordered = T) -> sub_expr$group
  sub_expr[order(sub_expr$group,sub_expr$zscore),] -> sub_expr
  
  
  p <- sub_expr %>% ggplot(aes(factor(gene, levels = unique(sub_expr$gene)), 
                               factor(cell_type, levels = ord))) + 
    geom_point(aes(size = zscore, colour = zscore)) + theme(strip.text.x = element_blank(),
                                                            axis.title = element_text(size = 15),
                                                            axis.text = element_text(size = 13),
                                                            legend.title = element_text(size = 13),
                                                            legend.text = element_text(size = 13),
                                                            axis.text.y = element_text(color = "black"), 
                                                            axis.text.x = element_text(color = "black", 
                                                                                       angle = -90, hjust = 0), 
                                                            panel.background = element_rect(colour = "black", fill = "white"), 
                                                            panel.grid = element_line(colour = "grey", linetype = "dashed"), 
                                                            panel.grid.major = element_line(colour = "grey", 
                                                                                            linetype = "dashed", size = 0.2)) + 
    facet_grid(. ~ group, scales = "free", space = "free") + 
    scale_colour_distiller(palette = "RdYlBu") + labs(x = "", y = "")
  
  return(p)
}

########################
LOY_dataset_summary <- function(df, sample_size_min = 20){
  # does the object have donor_organism.sex cell_ident_MV donor_organism.biomaterial_core.biomaterial_id?
  c("diagnosis","donor_organism.age","tissue") %in% names(df) -> tmp
  
  if(tmp[1] == FALSE){df$diagnosis <- NA}
  if(tmp[2] == FALSE){df$donor_organism.age <- NA}
  if(tmp[3] == FALSE){df$tissue <- NA}
  
  df %>% dplyr::group_by(donor_organism.biomaterial_core.biomaterial_id, LOY, cell_ident_MV, donor_organism.sex, diagnosis,
                         cohort, donor_organism.age, cell_broad_MV ,tissue) %>% 
    dplyr::summarize(counts = dplyr::n()) %>% 
    reshape2::dcast(donor_organism.biomaterial_core.biomaterial_id + cell_ident_MV + cell_broad_MV + tissue + 
                      donor_organism.sex + donor_organism.age + cohort + diagnosis ~ LOY, value.var = "counts") -> tmp
  
  if(tmp[is.na(tmp$LOY),]$LOY %>% length() > 0){tmp[is.na(tmp$LOY),]$LOY <- 0} 
  if(tmp[is.na(tmp$NORMAL),]$NORMAL %>% length() > 0){tmp[is.na(tmp$NORMAL),]$NORMAL <- 0} 
  
  dplyr::mutate(tmp, LOY_prop = LOY / (LOY + NORMAL), TOTAL = LOY + NORMAL) -> tmp
  
  ## remove female samples
  tmp[tmp$donor_organism.sex=="male",] -> tmp2
  
  ## plot LOY for each cell type across all individuals
  # order celltypes by prop 
  factor(tmp2$cell_ident_MV, levels = tmp2[order(-tmp2$LOY_prop),]$cell_ident_MV %>% unique()) -> tmp2$cell_ident_MV
  
  tmp2[tmp2$TOTAL > sample_size_min,] -> tmp2
  
  return(tmp2)
}


########################
expressed_Y_genes <- function(o){
  return(Y_scLOY_genes$gene_name[Y_scLOY_genes$gene_name %in% rownames(o)] %>% unique())
}

########################
expressed_PAR_genes <- function(o){
  return(PAR_scLOY_genes$gene_name[PAR_scLOY_genes$gene_name %in% rownames(o)] %>% unique())
}

########################
expressed_X_genes <- function(o){
  return(X_scLOY_genes$gene_name[X_scLOY_genes$gene_name %in% rownames(o)] %>% unique())
}

























