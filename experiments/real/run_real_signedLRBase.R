# run_real_signedLRBase.R
# Real-data signed CCI on GSE72056 using the curated signedLRBase
# (rule-based functional sign + OmniPath cross-check) instead of the 12-pair
# hand-list. Outputs to results/real_signedLRBase/ for comparison with run_real.R.
#
# Run: Rscript experiments/real/run_real_signedLRBase.R

suppressMessages({
  source("R/signed_lr_utils.R")
  source("R/signed_lrbase.R")
  source("R/gse72056.R")
})

set.seed(1)
infile <- "data/GSE72056_melanoma_single_cell_revised_v2.txt.gz"
outdir <- "results/real_signedLRBase"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- Cell-type expression (shared helper) --------------------------------
cat("Loading GSE72056 and assigning cell types ...\n")
ge <- load_gse72056_celltype_expr(infile, min_cells = 20L)
expr_ct <- ge$expr_ct
cell_names <- ge$cell_names
cat("Cell types:", paste(cell_names, collapse = ", "), "\n")
cat("Per-cell-type expression:", nrow(expr_ct), "x", ncol(expr_ct), "\n")

# ---- signedLRBase --------------------------------------------------------
lr_all <- load_signed_lrbase("data/signedLRBase/signedLRBase.csv")
present <- lr_all$ligand %in% rownames(expr_ct) & lr_all$receptor %in% rownames(expr_ct)
cat("\nsignedLRBase:", nrow(lr_all), "pairs;", sum(present),
    "usable in this dataset (both genes expressed-table present).\n")
lr_table <- lr_all[present, ]
cat("Usable sign distribution:\n"); print(table(sign = lr_table$sign))
cat("Usable category x sign:\n"); print(table(lr_table$category, lr_table$sign))

# ---- Signed CCI ----------------------------------------------------------
mats <- build_lr_matrices(lr_table, cell_names, expr = expr_ct)
res <- signed_cci(mats$lig_expr, mats$rec_expr, as.integer(lr_table$sign),
                  cell_names, max_hop = 3L, damping = 0.85)

A_pos <- res$cci$A_pos; A_neg <- res$cci$A_neg
cat("\n=== A_pos (sender -> receiver) ===\n"); print(round(A_pos, 2))
cat("\n=== A_neg (sender -> receiver) ===\n"); print(round(A_neg, 2))

pr <- res$pagerank
cat("\n=== Signed PageRank (converged =", pr$converged, ", iter =", pr$iter, ") ===\n")
pr_df <- data.frame(cell = cell_names,
                    positive = round(pr$positive, 5), negative = round(pr$negative, 5),
                    net = round(pr$net, 5), total = round(pr$total, 5))
pr_df <- pr_df[order(-pr_df$net), ]
print(pr_df, row.names = FALSE)

# ---- LR-name annotation (map pair index -> ligand_receptor) --------------
pair_names <- paste0(lr_table$ligand, "_", lr_table$receptor)
res$edge_table$pair_name <- pair_names[res$edge_table$lr_pair]
paths <- annotate_paths(res$paths, res$edge_table, cell_names)

ind_pos <- paths[paths$hop >= 2 & paths$net_sign == 1L, ]
ind_pos <- ind_pos[order(-ind_pos$contribution), ]
cat("\n=== Top net-POSITIVE indirect paths (hop>=2) ===\n")
print(head(ind_pos[, c("path_name","signs","hop","contribution","lr_annotation")], 15),
      row.names = FALSE)

motif <- paths[grepl("Treg", paths$path_name) & grepl("CD8", paths$path_name) &
               grepl("Tumor", paths$path_name) & paths$net_sign == 1L, ]
motif <- motif[order(-motif$contribution), ]
cat("\n=== Treg/CD8/Tumor net-positive paths ===\n")
print(head(motif[, c("path_name","signs","hop","contribution","lr_annotation")], 10),
      row.names = FALSE)

# ---- Persist -------------------------------------------------------------
write.csv(round(A_pos, 4), file.path(outdir, "A_pos.csv"))
write.csv(round(A_neg, 4), file.path(outdir, "A_neg.csv"))
write.csv(pr_df, file.path(outdir, "pagerank.csv"), row.names = FALSE)
write.csv(res$edge_table, file.path(outdir, "edge_table.csv"), row.names = FALSE)
sel <- c("path_name","signs","hop","net_sign","contribution","lr_annotation")
write.csv(head(ind_pos[, sel], 200),
          file.path(outdir, "indirect_positive_paths_top200.csv"), row.names = FALSE)
if (nrow(motif) > 0)
  write.csv(head(motif[, sel], 200), file.path(outdir, "motif_treg_cd8_tumor.csv"),
            row.names = FALSE)
write.csv(data.frame(cell_type = names(ge$counts), n = as.integer(ge$counts)),
          file.path(outdir, "celltype_counts.csv"), row.names = FALSE)

cat("\n=== DONE. Results in", outdir, "===\n")
