### LOAD FUNCTIONS





load_Seurat_object <- function(o, min_reads = 500, max_reads = 10000000, mt.filter = 20){

  message("finding Y chromosome genes")
  Y[Y$chrY %in% rownames(o),]$chrY -> chrY


  o[["percent.mt"]] <- Seurat::PercentageFeatureSet(o, pattern = "^MT-")
  o[["percent.rb"]] <- Seurat::PercentageFeatureSet(o, pattern = "^RP[SL]")

  message("filtering cells")
  o <- subset(o, subset = nFeature_RNA > min_reads & nFeature_RNA < max_reads & percent.mt < mt.filter)

  message("normalizing counts")
  o <- Seurat::NormalizeData(o, scale.factor = 10000)

  message("finding variable features")
  o <- FindVariableFeatures(o, selection.method = "vst", nfeatures = 2000)

  message("scaling data")
  o <- ScaleData(o)

  #message("dimension reduction")
  #o <- RunPCA(o)
  #o <- RunUMAP(o, dims = 1:10)

  #message("clustering")
  #o <- FindNeighbors(o, dims = 1:10)
  #o <- FindClusters(o, resolution = 0.5)

  #message("calling doublets")
  #sweep.res.list_o <- DoubletFinder::paramSweep_v3(o, PCs = 1:10, sct = FALSE)
  #sweep.stats_o <- DoubletFinder::summarizeSweep(sweep.res.list_o, GT = FALSE)
  #bcmvn_o <- DoubletFinder::find.pK(sweep.stats_o)
  #bcmvn_o[which(bcmvn_o$BCmetric==min(bcmvn_o$BCmetric)),]$pK %>% as.numeric() -> pk
  #homotypic.prop <- modelHomotypic(o@meta.data$seurat_clusters)           ## ex: annotations <- o@meta.data$ClusteringResults
  #nExp_poi <- round(0.05*nrow(o@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  #nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

  #o <- doubletFinder_v3(o, PCs = 1:10, pN = 0.25, pK = 0.14, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

  message("calling LOY")
  o@meta.data -> df
  o@assays$RNA@counts -> counts
  counts[chrY,] -> counts_chrY

  o@assays$RNA@data -> norm
  norm[chrY,] -> norm_chrY

  as.data.frame(colMeans(as.matrix(norm_chrY))) -> chrY_norm; names(chrY_norm) <- "chrY_norm"
  as.data.frame(colSums(as.matrix(counts_chrY))) -> chrY_UMI; names(chrY_UMI) <- "chrY_UMI"
  as.data.frame(colSums(as.matrix(counts_chrY > 0))) ->  chrY_genes; names(chrY_genes) <- "chrY_genes"



  df$chrY_genes <- chrY_genes$chrY_genes
  df$chrY_UMI <- chrY_UMI$chrY_UMI
  df$chrY_norm <- chrY_norm$chrY_norm
  tibble::rownames_to_column(df, "CB") -> df
  dplyr::mutate(df, LOY = ifelse(test = chrY_UMI == 0 , yes = "LOY", no = "NORMAL"))  -> df
  dplyr::mutate(df, normY = log1p((chrY_UMI / nCount_RNA) * 10000)) -> df
  dplyr::mutate(df, normY = log1p((chrY_UMI / nCount_RNA) * 10000)) -> df

  as.data.frame(quantile(df$nFeature_RNA, probs = seq(0,1,0.2))) -> q
  df$gene_cut <- cut(df$nFeature_RNA, breaks = as.vector(q[,1]))

  as.data.frame(quantile(df$nCount_RNA, probs = seq(0,1,0.2))) -> q
  cut(df$nCount_RNA, breaks = as.vector(q[,1])) -> df$read_cut
  tibble::column_to_rownames(df, "CB") -> df
  o@meta.data <- df



}



#' load_matrix_h5
#'
#' loads the h5 matrix, outputs formatted SEURAT object
#'
#' @param h5
#' @param project_name
#' @param min.cells
#' @param min.features
#' @param min_reads
#' @param max_reads
#' @param mt.filter
#'
#' @return formatted Seurat object
#' @export
#'
#' @examples
load_matrix_h5 <- function(h5, aggr, project_name, min.cells = 3, min.features = 200, regress_counts = F,
                           min_reads = 200, max_reads = 10000000, mt.filter = 20, output = NULL, genes_to_remove = NULL){

  message("loading expression matrix")
  Seurat::Read10X_h5(filename = h5) -> o

  message("making Seurat object")
  o <- Seurat::CreateSeuratObject(o, names.field = 2, names.delim = "-",
                                  project = project_name,
                                  min.cells = min.cells,
                                  min.features = min.features)

  read.csv(aggr) -> a
  a$seurat_id <- as.character(1:nrow(a))
  dplyr::select(a, library_id, seurat_id, sex) -> a

  as.data.frame(o@meta.data) -> df
  tibble::rownames_to_column(df,"CB") -> df
  dplyr::left_join(x = df, y = a, by = c("orig.ident"="seurat_id")) -> df
  tibble::column_to_rownames(df, 'CB') -> df
  df -> o@meta.data

  o[["percent.mt"]] <- Seurat::PercentageFeatureSet(o, pattern = "^MT-")
  o[["percent.rb"]] <- Seurat::PercentageFeatureSet(o, pattern = "^RP[SL]{1}[[:digit:]]+[^K]")

  message("filtering cells")
  # o <- subset(o, subset = nFeature_RNA > min_reads & nFeature_RNA < max_reads)

  message("normalizing counts")
  o <- Seurat::NormalizeData(o, scale.factor = 10000)

  message("finding variable features")
  o <- Seurat::FindVariableFeatures(o, selection.method = "vst", nfeatures = 2000)

  if(isTRUE(regress_counts)){
    message("scaling data and regressing out nCount_RNA")
    o <- Seurat::ScaleData(o, vars.to.regress = "nCount_RNA")} else {o <- Seurat::ScaleData(o)
    message("scaling data")}


  #message("dimension reduction")
  #o <- RunPCA(o)
  #o <- RunUMAP(o, dims = 1:10)

  #message("clustering")
  #o <- FindNeighbors(o, dims = 1:10)
  #o <- FindClusters(o, resolution = 0.5)

  #message("calling doublets")
  #sweep.res.list_o <- DoubletFinder::paramSweep_v3(o, PCs = 1:10, sct = FALSE)
  #sweep.stats_o <- DoubletFinder::summarizeSweep(sweep.res.list_o, GT = FALSE)
  #bcmvn_o <- DoubletFinder::find.pK(sweep.stats_o)
  #bcmvn_o[which(bcmvn_o$BCmetric==min(bcmvn_o$BCmetric)),]$pK %>% as.numeric() -> pk
  #homotypic.prop <- modelHomotypic(o@meta.data$seurat_clusters)           ## ex: annotations <- o@meta.data$ClusteringResults
  #nExp_poi <- round(0.05*nrow(o@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  #nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

  #o <- doubletFinder_v3(o, PCs = 1:10, pN = 0.25, pK = 0.14, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

  message("calling LOY")
  PAR_scLOY_genes$gene_name -> PAR_genes
  PAR_genes[PAR_genes %in% rownames(o)] -> PAR_genes

  Y_scLOY_genes$gene_name -> Y_genes
  Y_genes[Y_genes %in% rownames(o)] -> Y_genes
  if(!(is.null(genes_to_remove))){Y_genes <- Y_genes[!(Y_genes %in% genes_to_remove)]; print(paste0("removed gene: ", genes_to_remove))}

  print(paste0(length(Y_genes)," chromosome Y genes detected!"))
  print(paste0(Y_genes))

  o@meta.data -> df
  o@assays$RNA@counts -> counts
  counts[Y_genes,] -> counts_chrY
  counts[PAR_genes,] -> counts_PAR
  o@assays$RNA@data -> norm
  norm[Y_genes,] -> norm_chrY




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

  as.data.frame(quantile(df$nFeature_RNA, probs = seq(0,1,0.2))) -> q
  df$gene_cut <- cut(df$nFeature_RNA, breaks = as.vector(q[,1]))

  as.data.frame(quantile(df$nCount_RNA, probs = seq(0,1,0.2))) -> q
  cut(df$nCount_RNA, breaks = as.vector(q[,1])) -> df$read_cut
  tibble::column_to_rownames(df, "CB") -> df
  o@meta.data <- df

  if(!is.null(output)){
    message("saving Seurat object and metadata")
    saveRDS(object = o, file = paste0(output,"/",project_name,"_seurat_object.RDS"))
    saveRDS(object = o@meta.data, file = paste0(output,"/",project_name,"_meta_data.RDS"))
  }

  return(o)
}


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
  as.data.frame(colSums(as.matrix(counts_PAR > 0))) ->  PAR_genes; names(PAR_genes) <- "PAR_genes"

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



#' load_matrix_mtx
#'
#' @description needs to be tested
#'
#'
#' @param mtx
#' @param mtx_row
#' @param mtw_col
#' @param aggr
#' @param project_name
#' @param min.cells
#' @param min.features
#' @param min_reads
#' @param max_reads
#' @param mt.filter
#' @param output
#'
#' @return
#' @export
#'
#' @examples
load_matrix_mtx <- function(mtx, mtx_row, mtw_col, aggr, project_name, min.cells = 3, min.features = 300,
                            min_reads = 200, max_reads = 10000000, mt.filter =20, output= NULL){


  message("loading expression matrix")
  Matrix::readMM(mtx) -> o

  read.table(mtx_row) -> rows
  read.table(mtx_col) -> cols

  ## if rownames are ensembl ids switch them to hgnc


  dimnames(mtx) = list(as.character(rows$V1),as.character(cols$V1))

  Seurat::CreateSeuratObject(counts = o, project = project , min.cells = 3, min.features = 200, assay = "RNA") -> o

  message("finding Y chromosome genes")
  Y[Y$chrY %in% rownames(o),]$chrY -> chrY

  o[["percent.mt"]] <- Seurat::PercentageFeatureSet(o, pattern = "^MT-")
  o[["percent.rb"]] <- Seurat::PercentageFeatureSet(o, pattern = "^RP[SL]")

  message("filtering cells")
  o <- subset(o, subset = nFeature_RNA > min_reads & nFeature_RNA < max_reads & percent.mt < mt.filter)

  message("normalizing counts")
  Seurat::NormalizeData(o, scale.factor = 10000) -> o

  message("calling LOY")
  o@meta.data -> df
  o@assays$RNA@counts -> counts
  counts[chrY,] -> counts_chrY

  o@assays$RNA@data -> norm
  norm[chrY,] -> norm_chrY

  as.data.frame(colMeans(as.matrix(norm_chrY))) -> chrY_norm; names(chrY_norm) <- "chrY_norm"
  as.data.frame(colSums(as.matrix(counts_chrY))) -> chrY_UMI; names(chrY_UMI) <- "chrY_UMI"
  as.data.frame(colSums(as.matrix(counts_chrY > 0))) ->  chrY_genes; names(chrY_genes) <- "chrY_genes"

  df$chrY_genes <- chrY_genes$chrY_genes
  df$chrY_UMI <- chrY_UMI$chrY_UMI
  df$chrY_norm <- chrY_norm$chrY_norm
  df %>% rownames_to_column('CB') -> df
  dplyr::mutate(df, LOY = ifelse(test = chrY_UMI == 0 , yes = "LOY", no = "NORMAL"))  -> df
  dplyr::mutate(df, normY = log1p((chrY_UMI / nCount_RNA) * 10000)) -> df


  as.data.frame(quantile(df$nFeature_RNA, probs = seq(0,1,0.2))) -> q
  df$gene_cut <- cut(df$nFeature_RNA, breaks = as.vector(q[,1]))

  as.data.frame(quantile(df$nCount_RNA, probs = seq(0,1,0.2))) -> q
  cut(df$nCount_RNA, breaks = as.vector(q[,1])) -> df$read_cut

  gsub(x = df$CB, pattern = "[ACGT]+-", replacement = "") -> df$sample

  #df$cell_id -> row.names(df)
  read.csv(aggr) -> a
  a$seurat_id <- as.character(1:nrow(a))
  dplyr::select(a, library_id, seurat_id) -> a


  dplyr::left_join(x = df, y = a, by = c("sample"="seurat_id")) %>% column_to_rownames('CB') -> df
  df -> o@meta.data

  if(!is.null(output)){
    message("saving Seurat object and metadata")
    saveRDS(object = o, file = paste0(output,"/",project_name,"/seurat_object.RDS"))
    saveRDS(object = o@meta.data, file = paste0(output,"/",project_name,"/meta_data.RDS"))
  }

  return(o)
}


## add roxygen2 header
## this function takes in output from vartrix, filters for PAR snps, finds the LOY haplotype consensus, and adds info to meta data
load_vartrix <- function(o, consensus, snps, QUAL=200, UMI_cutoff=12, plots = F){

  read.csv(snps) -> snps
  readRDS(consensus) -> consensus
  o@meta.data -> meta

  ## process snps
  stringr::str_split(snps$locus, "_", simplify = T) -> tmp
  snps$rsid <- tmp[,3]
  snps$contig <- tmp[,1]
  snps$location <- tmp[,2] %>% as.numeric()

  # only keep known sites of variation (ie with rsid)
  snps[snps$rsid!=".",] -> snps

  dplyr::mutate(snps, total = LOY_total + NORMAL_total) -> snps


  # watch the chrX | X when using different reference + PAR coordinates
  # also only investigating PAR1 here
  snps[snps$contig=="chrX" & snps$QUAL >= QUAL & snps$LOY_total >= UMI_cutoff & snps$NORMAL_total >= UMI_cutoff & snps$location < 2781479,] -> snps.
  snps.$rsid -> quality_sites
  snps.[order(snps.$location),] -> snps.

  if(isTRUE(plots)){
    dplyr::select(snps., locus,gene,region,LOY_total,LOY_ref_ref_percent,LOY_alt_alt_percent,LOY_alt_ref_percent,rsid) -> L
    reshape2::melt(L) -> L.m
    L.m[L.m$variable %in% c("LOY_ref_ref_percent","LOY_alt_alt_percent","LOY_alt_ref_percent"),] -> L.m
    L.m$variable %>% factor(labels = c("REF","ALT","HET")) -> L.m$variable
    ggbarplot(data = L.m, x = "rsid", y = "value", fill = "variable", palette = "aaas") %>%
    ggpar(x.text.angle = 0, xlab = FALSE, ylab = "Proportion of cells (LOY)", legend.title = "", tickslab = T) -> a
    a + rremove("x.text") -> a

    snps. %>% dplyr::select(locus,gene,region, NORMAL_total,NORMAL_ref_ref_percent,NORMAL_alt_alt_percent,NORMAL_alt_ref_percent,rsid) -> N
    reshape2::melt(N) -> N.m
    N.m[N.m$variable %in% c("NORMAL_ref_ref_percent","NORMAL_alt_alt_percent","NORMAL_alt_ref_percent"),] -> N.m
    N.m$variable %>% factor(labels = c("REF","ALT","HET")) -> N.m$variable
    ggbarplot(data = N.m, x = "rsid", y = "value", fill = "variable", palette = "aaas") %>%
    ggpar(x.text.angle = 45, xlab = "RSID", ylab = "Proportion of cells (Normal)", legend.title = "") -> b

    ggarrange(plotlist = list(a,b), ncol = 1, nrow = 2,  align = "hv", common.legend = T ) -> pp
    print(pp)

  }

  # determine consensus Y haplotype
  ifelse(test = snps.$LOY_ref_ref_percent > snps.$LOY_alt_alt_percent, yes = "REF", no = "ALT") -> snps.$Y_haplo

  consensus[consensus$UMI_consensus!=0 & consensus$rsid%in%quality_sites,] -> consensus.
  dplyr::select(snps., -contig, -location, -locus) -> snps.
  merge.data.frame(x = consensus., y = snps., by.x = "rsid", by.y = "rsid") -> df

  df[which(df$UMI_call=="ref/ref"),]$UMI_call <- "REF"
  df[which(df$UMI_call=="alt/alt"),]$UMI_call  <- "ALT"
  df[which(df$UMI_call=="alt/ref"),]$UMI_call  <- "HET"
  as.character(df$UMI_call) -> df$UMI_call
  as.character(df$Y_haplo) -> df$Y_haplo

  ifelse(df$UMI_call==df$Y_haplo, "AGREE","DISAGREE") -> df$agree
  with(df, table(CB,agree)) %>% as.data.frame.matrix()  -> tbl
  merge.data.frame(x = df, y = tbl, by.x = "CB", by.y = 0) -> tbl

  dplyr::distinct(tbl,CB, .keep_all = TRUE) %>% dplyr::select(CB,AGREE,DISAGREE) -> ASE_cell_attributes

  merge.data.frame(x = meta, y = ASE_cell_attributes, by.x = 0, by.y = "CB",all.x = T) -> df
  tibble::column_to_rownames(df,"Row.names") -> df


  return(list(snps,df,consensus))

}





#' load_ASE
#'
#' load ASEReadCounter files, format and output as data.frame
#'
#' @param ase_dir
#' @param sample
#' @param meta
#'
#' @return
#' @export
#'
#' @examples
load_ASE <- function(ase_dir, sample, meta){



  unique(list.files(path = paste0(ase_dir,"/",sample),
             pattern = paste0(sample), full.names = T)) -> files

  list() -> l
  for(i in files){
    message(i)
    data.table::fread(i) -> N
    N$cell <- gsub(pattern = paste0(sample,"_"), replacement = "", x = basename(i))
    N$sample <- sample
    N[N$totalCount >= 1 & N$variantID!="." , ] -> N

    dplyr::mutate(N, REF_PERCENT = refCount / totalCount,
           ALT_PERCENT = altCount /totalCount  ) -> N
    dplyr::mutate(N, BIAS = abs(REF_PERCENT - 0.5)) -> N



    rlist::list.append(l, N) -> l
  }
  do.call("rbind",l) -> df

  dplyr::mutate(df, expression_REF = ifelse(refCount > 0, TRUE, FALSE),
         expression_ALT = ifelse(altCount > 0, TRUE, FALSE)) -> df
  dplyr::mutate(df, expression_BOTH = ifelse(expression_REF==T & expression_ALT==T, TRUE, FALSE)) -> df

  df[df$expression_BOTH==T,]$expression_ALT <- NA
  df[df$expression_BOTH==T,]$expression_REF <- NA



  return(df)
}



LOY_summary_table <- function(s_obj, dataset, author, introns){

  s_obj@meta.data -> md
  dplyr::select(tibble::rownames_to_column(md,"CB"), CB, library_id, sex, nCount_RNA, nFeature_RNA, percent.mt, percent.rb, chrY_genes, chrY_UMI, PAR_UMI, PAR_genes, LOY, normY, normPAR,
                RNA_snn_res.0.3, RNA_snn_res.0.5, RNA_snn_res.0.7, RNA_snn_res.0.9, RNA_snn_res.1.2, broad_celltypes, cell_type) -> md

  md$data_set <- dataset
  md$author <- author
  md$introns <- introns

  ## add chrY genes and PAR genes if they arent in the dataset set NA
  rownames(s_obj)[rownames(s_obj) %in% PAR_scLOY_genes$gene_name] -> p
  rownames(s_obj)[rownames(s_obj) %in% Y_scLOY_genes$gene_name] -> y

  s_obj@assays$RNA@counts[y,] %>% as.matrix() %>% as.data.frame() %>% t() -> y_counts
  s_obj@assays$RNA@counts[p,] %>% as.matrix() %>% as.data.frame() %>% t() -> par_counts

  cbind(md, y_counts, par_counts) -> md

  return(md)
}



