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
