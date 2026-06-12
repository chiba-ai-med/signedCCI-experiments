# run_real.R
# Signed cell-cell interaction (CCI) on real immune x tumor data:
# GSE72056 (Tirosh et al. 2016, melanoma tumor micro-environment).
#
# Pipeline: read counts -> assign cell types (incl. CD8/Treg from T cells) ->
# per-cell-type mean expression -> build_lr_matrices(expr, lr_table) ->
# signed_cci() -> report signed PageRank + indirect (-)x(-)=(+) paths.
#
# Run: Rscript experiments/real/run_real.R

suppressMessages({
  source("R/signed_lr_utils.R")
  library(data.table)
})

set.seed(1)
infile <- "data/GSE72056_melanoma_single_cell_revised_v2.txt.gz"
outdir <- "results/real"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(infile))

# ---- Read the matrix -----------------------------------------------------
# Layout: col 1 = label/gene; row1 tumor id, row2 malignant code,
# row3 non-malignant cell-type code, rows 4+ = genes. log2(TPM/10+1) values.
cat("Reading", infile, "...\n")
dt <- fread(cmd = paste("zcat", infile), header = TRUE, sep = "\t",
            showProgress = FALSE)
cell_ids <- colnames(dt)[-1]
labels   <- dt[[1]]

anno_malig <- as.integer(unlist(dt[2, -1]))   # 1=no, 2=yes, 0=unresolved
anno_ct    <- as.integer(unlist(dt[3, -1]))   # 1=T,2=B,3=Macro,4=Endo,5=CAF,6=NK

# Expression matrix (genes x cells); collapse duplicate gene symbols by max.
genes <- labels[-(1:3)]
expr_all <- as.matrix(dt[-(1:3), -1])
storage.mode(expr_all) <- "double"
rm(dt); gc()
if (anyDuplicated(genes)) {
  cat("Collapsing", sum(duplicated(genes)), "duplicate gene symbols (max)...\n")
  keep <- !duplicated(genes)
  ord <- order(genes, -matrixStats::rowMaxs(expr_all))
  expr_all <- expr_all[ord, ]; genes <- genes[ord]
  keep <- !duplicated(genes)
  expr_all <- expr_all[keep, ]; genes <- genes[keep]
}
rownames(expr_all) <- genes
cat("Matrix:", nrow(expr_all), "genes x", ncol(expr_all), "cells\n")

# ---- Cell-type assignment ------------------------------------------------
g <- function(sym) if (sym %in% rownames(expr_all)) expr_all[sym, ] else rep(0, ncol(expr_all))

ct_code_map <- c("1" = "Tcell", "2" = "Bcell", "3" = "Macrophage",
                 "4" = "Endo", "5" = "CAF", "6" = "NK")

celltype <- rep(NA_character_, length(cell_ids))
celltype[anno_malig == 2] <- "Tumor"
nonmalig <- which(anno_malig == 1)
celltype[nonmalig] <- ct_code_map[as.character(anno_ct[nonmalig])]

# Subdivide T cells into CD8 / Treg / CD4conv using markers.
FOXP3 <- g("FOXP3"); CD8A <- g("CD8A"); CD8B <- g("CD8B")
is_T <- which(celltype == "Tcell")
for (i in is_T) {
  if (FOXP3[i] > 0) {
    celltype[i] <- "Treg"
  } else if (CD8A[i] > 0 || CD8B[i] > 0) {
    celltype[i] <- "CD8"
  } else {
    celltype[i] <- "CD4conv"
  }
}

tab <- sort(table(celltype, useNA = "no"), decreasing = TRUE)
cat("\n=== Cell-type counts ===\n"); print(tab)

# Keep cell types with enough cells; drop NA / tiny groups.
min_cells <- 20
keep_types <- names(tab)[tab >= min_cells]
cat("\nKeeping cell types (>=", min_cells, "cells):",
    paste(keep_types, collapse = ", "), "\n")

# ---- Per-cell-type mean expression ---------------------------------------
expr_ct <- sapply(keep_types, function(ct) {
  rowMeans(expr_all[, which(celltype == ct), drop = FALSE])
})
expr_ct <- as.matrix(expr_ct)            # gene x cell type
cat("Per-cell-type expression matrix:", nrow(expr_ct), "x", ncol(expr_ct), "\n")

# ---- Signed LR table (unique pairs; one biological sign each) -------------
# In expression-driven (real) mode the sign is a property of the LR pair, not a
# directed edge, so each pair gets a single global sign.
lr_table <- data.frame(
  ligand   = c("TGFB1","IL10","CD274","FASLG","LGALS9","PVR",
               "CD80","CD86","IL12A","IFNG","CXCL9","CCL5"),
  receptor = c("TGFBR1","IL10RA","PDCD1","FAS","HAVCR2","TIGIT",
               "CD28","CD28","IL12RB1","IFNGR1","CXCR3","CCR5"),
  sign     = c(-1, -1, -1, -1, -1, -1,
                1,  1,  1,  1,  1,  1),
  stringsAsFactors = FALSE
)
lr_table$pair_name <- paste0(lr_table$ligand, "_", lr_table$receptor)

# Report which LR pairs are usable (both genes present).
present <- lr_table$ligand %in% rownames(expr_ct) &
           lr_table$receptor %in% rownames(expr_ct)
cat("\n=== Signed LR pairs (", sum(present), "of", nrow(lr_table), "usable) ===\n")
print(cbind(lr_table[, c("ligand","receptor","sign")], usable = present))
lr_table <- lr_table[present, ]

# ---- Build matrices and run signed CCI -----------------------------------
cell_names <- keep_types
m <- build_lr_matrices(lr_table, cell_names, expr = expr_ct)
res <- signed_cci(m$lig_expr, m$rec_expr, m$lr_sign, cell_names,
                  max_hop = 3L, damping = 0.85)

A_pos <- res$cci$A_pos; A_neg <- res$cci$A_neg   # sender x receiver
cat("\n=== A_pos (sender -> receiver), rounded ===\n"); print(round(A_pos, 2))
cat("\n=== A_neg (sender -> receiver), rounded ===\n"); print(round(A_neg, 2))

# ---- Signed PageRank -----------------------------------------------------
pr <- res$pagerank
cat("\n=== Signed PageRank (converged =", pr$converged, ", iter =", pr$iter, ") ===\n")
pr_df <- data.frame(cell = cell_names,
                    positive = round(pr$positive, 5),
                    negative = round(pr$negative, 5),
                    net = round(pr$net, 5),
                    total = round(pr$total, 5))
pr_df <- pr_df[order(-pr_df$net), ]
print(pr_df, row.names = FALSE)

# ---- Signed paths + LR annotation ----------------------------------------
paths <- annotate_paths(res$paths, res$edge_table, cell_names)

# Net-positive indirect paths (>=2 hops, even # of negatives).
ind_pos <- paths[paths$hop >= 2 & paths$net_sign == 1L, ]
ind_pos <- ind_pos[order(-ind_pos$contribution), ]
cat("\n=== Top net-POSITIVE indirect paths (hop>=2) ===\n")
print(head(ind_pos[, c("path_name","signs","hop","net_sign","contribution","lr_annotation")], 15),
      row.names = FALSE)

# Focused question: indirect positive routes that END at Tumor through an
# even number of inhibitory hops (the "enemy of my enemy" anti-tumor motif).
to_tumor <- ind_pos[ind_pos$target_name == "Tumor" &
                    grepl("-", ind_pos$signs, fixed = TRUE), ]
cat("\n=== Net-POSITIVE indirect paths ENDING at Tumor (with >=1 inhibitory hop) ===\n")
if (nrow(to_tumor) > 0) {
  print(head(to_tumor[, c("path_name","signs","hop","contribution","lr_annotation")], 15),
        row.names = FALSE)
} else {
  cat("(none found)\n")
}

# Specifically look for the Treg/CD8/Tumor motif seen in the synthetic gate.
motif <- paths[grepl("Treg", paths$path_name) & grepl("CD8", paths$path_name) &
               grepl("Tumor", paths$path_name) & paths$net_sign == 1L, ]
cat("\n=== Treg/CD8/Tumor net-positive paths (synthetic motif in real data) ===\n")
if (nrow(motif) > 0) {
  motif <- motif[order(-motif$contribution), ]
  print(head(motif[, c("path_name","signs","hop","contribution","lr_annotation")], 10),
        row.names = FALSE)
} else {
  cat("(none found)\n")
}

# ---- Persist -------------------------------------------------------------
# Note: the full path enumeration is large (~10^5 rows) and regenerable by
# re-running this script; we commit curated subsets instead.
write.csv(round(A_pos, 4), file.path(outdir, "A_pos.csv"))
write.csv(round(A_neg, 4), file.path(outdir, "A_neg.csv"))
write.csv(pr_df, file.path(outdir, "pagerank.csv"), row.names = FALSE)
write.csv(res$edge_table, file.path(outdir, "edge_table.csv"), row.names = FALSE)
sel_cols <- c("path_name","signs","hop","net_sign","contribution","lr_annotation")
write.csv(head(ind_pos[, sel_cols], 200),
          file.path(outdir, "indirect_positive_paths_top200.csv"), row.names = FALSE)
write.csv(to_tumor[, sel_cols],
          file.path(outdir, "indirect_positive_paths_to_tumor.csv"), row.names = FALSE)
if (nrow(motif) > 0)
  write.csv(motif[, sel_cols],
            file.path(outdir, "motif_treg_cd8_tumor.csv"), row.names = FALSE)
write.csv(data.frame(cell_type = names(tab), n = as.integer(tab)),
          file.path(outdir, "celltype_counts.csv"), row.names = FALSE)
write.csv(round(expr_ct[unique(c(lr_table$ligand, lr_table$receptor)), ], 4),
          file.path(outdir, "lr_gene_expression.csv"))

cat("\n=== DONE. Results in", outdir, "===\n")
