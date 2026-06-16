# gse103322.R
# Load GSE103322 (Puram et al. 2017, head & neck squamous cell carcinoma; same
# lab/format as GSE72056), assign cell types, return a (gene x cell type) mean
# expression matrix. Mirrors R/gse72056.R so the downstream analysis is identical.

suppressMessages(library(data.table))

#' @param infile path to GSE103322_HNSCC_all_data.txt.gz
#' @param min_cells keep cell types with at least this many cells
#' @return list(expr_ct, counts, celltype, cell_names)
load_gse103322_celltype_expr <- function(infile, min_cells = 20L) {
  stopifnot(file.exists(infile))
  dt <- fread(cmd = paste("zcat", infile), header = TRUE, sep = "\t",
              showProgress = FALSE)
  labels <- dt[[1]]
  # Body rows: 1 Maxima, 2 Lymph node, 3 cancer(0/1), 4 non-cancer(0/1),
  # 5 non-cancer cell type (text), 6+ genes.
  cancer <- as.integer(unlist(dt[3, -1]))
  nctype <- as.character(unlist(dt[5, -1]))
  nctype <- gsub("^-", "", trimws(nctype))            # "-Fibroblast" -> "Fibroblast"

  genes <- gsub("'", "", labels[-(1:5)])              # strip single quotes
  expr_all <- as.matrix(dt[-(1:5), -1])
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
  type_map <- c("T cell" = "Tcell", "B cell" = "Bcell", "Macrophage" = "Macrophage",
                "Dendritic" = "DC", "Endothelial" = "Endo", "Fibroblast" = "CAF",
                "Mast" = "Mast", "myocyte" = "Myocyte")

  celltype <- rep(NA_character_, length(cancer))
  celltype[cancer == 1] <- "Tumor"
  noncancer <- which(cancer != 1)
  celltype[noncancer] <- type_map[nctype[noncancer]]

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
