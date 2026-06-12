# lr_contribution.R
# Which ligand-receptor pairs actually DRIVE each signed edge?
#
# BuildSignedCCI's edge_table holds the per-LR score X_ijk = L_ik * R_jk for
# every (sender, receiver, LR pair). A_pos/A_neg are sums over same-sign LRs.
# Here we decompose each directed signed edge back into its LR contributions,
# rank them, and report the dominant LR per edge (the "contributing LR" detection).
#
# Run (after run_real.R):  Rscript experiments/real/lr_contribution.R

outdir <- "results/real"
et <- read.csv(file.path(outdir, "edge_table.csv"), stringsAsFactors = FALSE)

# LR-pair index -> name (must match the filtered lr_table order in run_real.R).
lr_names <- c("TGFB1_TGFBR1","IL10_IL10RA","CD274_PDCD1","FASLG_FAS",
              "LGALS9_HAVCR2","PVR_TIGIT","CD80_CD28","CD86_CD28",
              "IL12A_IL12RB1","IFNG_IFNGR1","CXCL9_CXCR3","CCL5_CCR5")
et$lr_name <- lr_names[et$lr_pair]
et$edge <- paste0(et$sender_name, "->", et$receiver_name)
et$polarity <- ifelse(et$sign > 0, "pos", "neg")

# ---- Per-edge LR ranking -------------------------------------------------
# Within each (edge, polarity), rank LR pairs by their score share.
split_key <- paste(et$edge, et$polarity, sep = "|")
et <- do.call(rbind, lapply(split(et, split_key), function(d) {
  d <- d[order(-d$score), ]
  d$share <- d$score / sum(d$score)
  d$rank <- seq_len(nrow(d))
  d
}))
rownames(et) <- NULL

write.csv(et[, c("edge","polarity","lr_name","score","share","rank",
                 "sender_name","receiver_name","sign")],
          file.path(outdir, "lr_contribution_by_edge.csv"), row.names = FALSE)

# Dominant (rank-1) LR per edge & polarity.
dom <- et[et$rank == 1, c("edge","polarity","lr_name","score","share")]
dom <- dom[order(-dom$score), ]
cat("=== Dominant LR pair per directed signed edge (top 20 by score) ===\n")
print(head(dom, 20), row.names = FALSE)
write.csv(dom, file.path(outdir, "dominant_lr_per_edge.csv"), row.names = FALSE)

# ---- LR ranking for the key (-)(-)=(+) motif edges -----------------------
motif_edges <- c("Tumor->Treg", "Treg->CD8", "CD8->Tumor",
                 "CAF->Treg", "Endo->CD8", "NK->Tumor")
cat("\n=== LR contributions on key motif edges (negative hops) ===\n")
for (e in motif_edges) {
  sub <- et[et$edge == e & et$polarity == "neg", ]
  if (nrow(sub) == 0) next
  sub <- sub[order(-sub$share), ]
  cat(sprintf("\n[%s]  (-)  A_neg = %.3f\n", e, sum(sub$score)))
  print(sub[, c("lr_name","score","share")], row.names = FALSE)
}

# ---- Global LR importance (total signed weight carried) ------------------
glob <- aggregate(score ~ lr_name + polarity, data = et, sum)
glob <- glob[order(-glob$score), ]
cat("\n=== Global LR-pair importance (total score across all edges) ===\n")
print(glob, row.names = FALSE)
write.csv(glob, file.path(outdir, "lr_global_importance.csv"), row.names = FALSE)

cat("\n=== DONE. LR-contribution tables in", outdir, "===\n")
