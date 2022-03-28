suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
edb <- EnsDb.Hsapiens.v86


## Query PAR genes
grf <- AnnotationFilter::GRangesFilter(GRanges("X", ranges = IRanges(10001, 2781479)), type = "any")
gn <- as.data.frame(ensembldb::genes(edb, filter = grf))
dplyr::select(gn, gene_id, strand, gene_name, gene_biotype) -> PAR_scLOY_genes

grf2 <- AnnotationFilter::GRangesFilter(GRanges("X", ranges = IRanges(155701383, 156030895)), type = "any")
gn2 <- as.data.frame(ensembldb::genes(edb, filter = grf2))
dplyr::select(gn, gene_id, strand, gene_name, gene_biotype) -> PAR2_scLOY_genes

rbind(PAR_scLOY_genes,PAR2_scLOY_genes) -> PAR_scLOY_genes
PAR_scLOY_genes[!(PAR_scLOY_genes$gene_biotype %in% c("LRG_gene")),] -> PAR_scLOY_genes

## Query chrY genes
grf <- AnnotationFilter::GRangesFilter(GRanges("Y", ranges = IRanges(1, 59373566)), type = "any")

gn <- as.data.frame(ensembldb::genes(edb, filter = grf))
dplyr::select(gn, gene_id, strand, gene_name, gene_biotype) -> Y_scLOY_genes

usethis::use_data(Y_scLOY_genes,overwrite = TRUE)
usethis::use_data(PAR_scLOY_genes,overwrite = TRUE)

# Query X genes
grf <- AnnotationFilter::GRangesFilter(GRanges("X", ranges = IRanges(2781480, 155701382)), type = "any")
gn <- as.data.frame(ensembldb::genes(edb, filter = grf))
dplyr::select(gn, gene_id, strand, gene_name, gene_biotype) -> X_scLOY_genes

usethis::use_data(X_scLOY_genes,overwrite = TRUE)

###### MOUSE
suppressPackageStartupMessages(library(EnsDb.Mmusculus.v79))
edb <- EnsDb.Mmusculus.v79
## Query PAR genes
grf <- AnnotationFilter::GRangesFilter(GRanges("X", ranges = IRanges(169969759, 170931299)), type = "any")
gn <- as.data.frame(ensembldb::genes(edb, filter = grf))
dplyr::select(gn, gene_id, strand, gene_name, gene_biotype) -> PAR_scLOY_genes_mouse
PAR_scLOY_genes[!(PAR_scLOY_genes$gene_biotype %in% c("LRG_gene")),] -> PAR_scLOY_genes_mouse

## Query chrY genes
grf <- AnnotationFilter::GRangesFilter(GRanges("Y", ranges = IRanges(1, 91744698)), type = "any")
gn <- as.data.frame(ensembldb::genes(edb, filter = grf))
dplyr::select(gn, gene_id, strand, gene_name, gene_biotype) -> Y_scLOY_genes_mouse


usethis::use_data(Y_scLOY_genes_mouse,overwrite = TRUE)
usethis::use_data(PAR_scLOY_genes_mouse,overwrite = TRUE)




