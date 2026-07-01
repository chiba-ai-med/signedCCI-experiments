# gse120575.R
# Load GSE120575 (Sade-Feldman 2018): CD45+ immune cells from melanoma patients
# on checkpoint immunotherapy, with Responder / Non-responder and Pre/Post labels.
# No tumor cells. Cell types assigned from markers (no prior annotation).
#
# The full matrix is 55,737 genes x 16,291 cells -> reading it whole overflows
# R's 2^31-byte string limit in fread. We only need ~200 genes (markers, receiver
# panels, and signedLRBase LR genes), so we extract just those rows with awk.

suppressMessages(library(data.table))

.gse120575_markers <- list(
  lineage = list(
    Tcell     = c("CD3D","CD3E","CD3G","TRAC"),
    NK        = c("NKG7","GNLY","KLRD1","NCAM1","KLRF1"),
    Bcell     = c("MS4A1","CD79A","CD79B","CD19"),
    Plasma    = c("MZB1","IGHG1","JCHAIN","DERL3"),
    Macrophage= c("LYZ","CD68","CD14","C1QA","C1QB","CD163","FCGR3A"),
    DC        = c("FCER1A","CD1C","CLEC9A","LILRA4")),
  Tsub = c("FOXP3","CD8A","CD8B","CD4"),
  panel = c("GZMA","GZMB","GZMK","PRF1","NKG7","GNLY","KLRD1","IFNG",
            "PDCD1","HAVCR2","LAG3","TIGIT","CTLA4","TOX","ENTPD1","BTLA",
            "TNF","IL1B","CXCL9","CXCL10","CXCL11","CD86","IL12B","NOS2",
            "MRC1","CD163","MSR1","IL10","ARG1","MARCO","CCL22"))

#' @param extra_genes additional gene symbols to read (e.g. signedLRBase LR genes)
#' @return list(expr [log2(TPM/10+1), gene x cell], celltype, response,
#'   timepoint, patient, barcodes)
load_gse120575 <- function(tpm_file = "data/GSE120575_TPM.txt.gz",
                           meta_file = "data/GSE120575_patient_ID_single_cells.txt.gz",
                           extra_genes = character(0)) {
  need <- unique(c(unlist(.gse120575_markers$lineage), .gse120575_markers$Tsub,
                   .gse120575_markers$panel, extra_genes))
  gf <- tempfile(); writeLines(need, gf)
  of <- tempfile(fileext = ".txt")
  # keep header (row1) + patient/timepoint (row2) + rows whose gene (col1) is needed
  system(sprintf("zcat %s | awk -F'\\t' 'NR==FNR{g[$1]=1;next} FNR<=2 || ($1 in g)' %s - > %s",
                 shQuote(tpm_file), shQuote(gf), shQuote(of)))
  raw <- as.data.frame(fread(of, header = FALSE, sep = "\t", fill = TRUE, showProgress = FALSE))
  unlink(c(gf, of))
  # Gene rows carry a trailing empty field (16293 vs 16292); keep only the
  # non-empty header (barcode) columns to stay aligned.
  hdr <- as.character(raw[1, ])
  bcol <- which(nzchar(hdr) & !is.na(hdr))        # barcode columns
  barcodes <- hdr[bcol]
  pt_tp <- as.character(raw[2, bcol])
  genes <- as.character(raw[-(1:2), 1])
  expr <- as.matrix(raw[-(1:2), bcol]); storage.mode(expr) <- "double"
  rm(raw); gc()
  expr <- log2(expr / 10 + 1)
  if (anyDuplicated(genes)) { keep <- !duplicated(genes); expr <- expr[keep, ]; genes <- genes[keep] }
  rownames(expr) <- genes; colnames(expr) <- barcodes

  timepoint <- sub("_.*$", "", pt_tp)
  patient   <- sub("^[^_]*_", "", pt_tp)

  meta <- fread(cmd = paste("zcat", meta_file), header = FALSE, sep = "\t", fill = TRUE, showProgress = FALSE)
  smp <- meta[grepl("^Sample [0-9]", meta$V1, useBytes = TRUE)]
  response <- unname(setNames(smp$V6, smp$V2)[barcodes])

  # ---- Marker-based cell-type assignment ---------------------------------
  zrow <- function(g) { g <- intersect(g, rownames(expr))
    if (!length(g)) return(NULL); t(scale(t(expr[g, , drop = FALSE]))) }
  score <- sapply(.gse120575_markers$lineage, function(gs) { z <- zrow(gs)
    if (is.null(z)) rep(0, ncol(expr)) else { z[is.na(z)] <- 0; colMeans(z) } })
  lin <- colnames(score)[max.col(score, ties.method = "first")]

  gv <- function(s) if (s %in% rownames(expr)) expr[s, ] else rep(0, ncol(expr))
  FOXP3 <- gv("FOXP3"); CD8A <- gv("CD8A"); CD8B <- gv("CD8B")
  celltype <- lin
  isT <- lin == "Tcell"
  celltype[isT] <- ifelse(FOXP3[isT] > 0, "Treg",
                    ifelse(CD8A[isT] > 0 | CD8B[isT] > 0, "CD8", "CD4conv"))
  celltype[lin == "Plasma"] <- "Bcell"

  list(expr = expr, celltype = celltype, response = response,
       timepoint = timepoint, patient = patient, barcodes = barcodes)
}

#' Per-cell-type mean expression for a subset of cells.
celltype_mean_expr <- function(d, mask, min_cells = 20L) {
  ct <- d$celltype[mask]; E <- d$expr[, mask, drop = FALSE]
  tab <- table(ct); keep <- names(tab)[tab >= min_cells]
  m <- sapply(keep, function(k) rowMeans(E[, ct == k, drop = FALSE]))
  list(expr_ct = as.matrix(m), counts = sort(tab, decreasing = TRUE), cell_names = keep)
}
