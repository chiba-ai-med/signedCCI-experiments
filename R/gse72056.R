# gse72056.R
# Load GSE72056 (Tirosh et al. 2016 melanoma TME), assign cell types, and
# return a per-cell-type mean-expression matrix. Shared by the real-data scripts.

suppressMessages(library(data.table))

#' Read GSE72056 and build a (gene x cell type) mean-expression matrix.
#'
#' Cell types: malignant -> Tumor; non-malignant codes -> Bcell/Macrophage/
#' Endo/CAF/NK; T cells split into Treg (FOXP3+), CD8 (CD8A/CD8B+), else CD4conv.
#'
#' @param infile path to GSE72056_melanoma_single_cell_revised_v2.txt.gz.
#' @param min_cells keep only cell types with at least this many cells.
#' @return list(expr_ct = gene x cell type means, counts = table of all calls,
#'   celltype = per-cell labels, cell_names = kept cell types).
load_gse72056_celltype_expr <- function(infile, min_cells = 20L) {
  stopifnot(file.exists(infile))
  dt <- fread(cmd = paste("zcat", infile), header = TRUE, sep = "\t",
              showProgress = FALSE)
  cell_ids <- colnames(dt)[-1]
  labels   <- dt[[1]]
  anno_malig <- as.integer(unlist(dt[2, -1]))   # 1=no, 2=yes, 0=unresolved
  anno_ct    <- as.integer(unlist(dt[3, -1]))   # 1=T,2=B,3=Macro,4=Endo,5=CAF,6=NK

  genes <- labels[-(1:3)]
  expr_all <- as.matrix(dt[-(1:3), -1])
  storage.mode(expr_all) <- "double"
  rm(dt); gc()

  if (anyDuplicated(genes)) {
    ord <- order(genes, -matrixStats::rowMaxs(expr_all))
    expr_all <- expr_all[ord, ]; genes <- genes[ord]
    keep <- !duplicated(genes)
    expr_all <- expr_all[keep, ]; genes <- genes[keep]
  }
  rownames(expr_all) <- genes

  g <- function(sym) if (sym %in% rownames(expr_all)) expr_all[sym, ] else rep(0, ncol(expr_all))
  ct_code_map <- c("1" = "Tcell", "2" = "Bcell", "3" = "Macrophage",
                   "4" = "Endo", "5" = "CAF", "6" = "NK")

  celltype <- rep(NA_character_, length(cell_ids))
  celltype[anno_malig == 2] <- "Tumor"
  nonmalig <- which(anno_malig == 1)
  celltype[nonmalig] <- ct_code_map[as.character(anno_ct[nonmalig])]

  FOXP3 <- g("FOXP3"); CD8A <- g("CD8A"); CD8B <- g("CD8B")
  is_T <- which(celltype == "Tcell")
  for (i in is_T) {
    if (FOXP3[i] > 0) celltype[i] <- "Treg"
    else if (CD8A[i] > 0 || CD8B[i] > 0) celltype[i] <- "CD8"
    else celltype[i] <- "CD4conv"
  }

  counts <- sort(table(celltype, useNA = "no"), decreasing = TRUE)
  keep_types <- names(counts)[counts >= min_cells]

  expr_ct <- sapply(keep_types, function(ct)
    rowMeans(expr_all[, which(celltype == ct), drop = FALSE]))
  expr_ct <- as.matrix(expr_ct)

  list(expr_ct = expr_ct, counts = counts, celltype = celltype,
       cell_names = keep_types)
}
