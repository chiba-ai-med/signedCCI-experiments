# sign_refinement_markers.R
# Simple marker-based "is the sign consistent in THIS dataset?" check for the
# context-dependent pairs of signedLRBase, on GSE72056 (melanoma TME).
#
# Idea (the simplest version of data-driven sign refinement): look at the
# RECEIVER cell type and ask whether it is in an ACTIVATED or a SUPPRESSED state,
# using small marker panels. If the receiver is activated -> incoming signals net
# to "+"; if suppressed/exhausted -> net "-". Compare this data-state sign to the
# prior (curated) sign.
#
# HONEST CAVEAT: this measures the receiver's OVERALL state, not the causal effect
# of each specific ligand (that needs CytoSig / DoRothEA-style footprints). It is a
# consistency check, not a per-pair causal correction.
#
# Run: Rscript experiments/real/sign_refinement_markers.R

suppressMessages({
  source("R/signed_lrbase.R")
  source("R/gse72056.R")
})

infile <- "data/GSE72056_melanoma_single_cell_revised_v2.txt.gz"
outdir <- "results/real_signedLRBase"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

ge <- load_gse72056_celltype_expr(infile, min_cells = 20L)
expr_ct <- ge$expr_ct                       # gene x cell type (mean log-expr)

# z-score each gene across cell types -> relative score per cell type.
z <- t(scale(t(expr_ct)))
z[is.na(z)] <- 0

# ---- Marker panels -------------------------------------------------------
panels <- list(
  Tnk_activation  = c("GZMA","GZMB","GZMK","PRF1","NKG7","GNLY","KLRD1","IFNG"),
  Tnk_exhaustion  = c("PDCD1","HAVCR2","LAG3","TIGIT","CTLA4","TOX","ENTPD1","BTLA"),
  Mac_m1          = c("TNF","IL1B","CXCL9","CXCL10","CXCL11","CD86","IL12B","NOS2"),
  Mac_m2          = c("MRC1","CD163","MSR1","IL10","ARG1","MARCO","CCL22")
)
score <- function(cell, genes) {
  g <- intersect(genes, rownames(z))
  if (length(g) == 0 || !cell %in% colnames(z)) return(NA_real_)
  mean(z[g, cell])
}

# ---- Receiver + panel assignment per context-dependent pair --------------
# lineage: "Tnk" (eval in CD8 / NK) or "Mac" (eval in Macrophage). "skip" = no clean panel.
recv_map <- list(
  IFNG_IFNGR1="Tnk", IFNG_IFNGR2="Tnk", IL1B_IL1R1="Tnk",
  TNF_TNFRSF1A="Tnk", TNF_TNFRSF1B="Tnk", TNFSF14_TNFRSF14="Tnk",
  TGFB1_TGFBR1="Tnk", TGFB1_TGFBR2="Tnk", TGFB2_TGFBR1="Tnk", TGFB3_TGFBR1="Tnk",
  PVR_CD96="Tnk", HMGB1_HAVCR2="Tnk", LGALS3_LAG3="Tnk", CD160_TNFRSF14="Tnk",
  `HLA-G_KIR2DL4`="NK",
  CSF1_CSF1R="Mac", IL4_IL4R="Mac", IL13_IL13RA1="Mac", VEGFA_FLT1="Mac",
  CXCL12_CXCR4="skip", TNFSF13B_TNFRSF13B="skip"
)

db <- load_signed_lrbase()
ctx <- db[as.logical(db$context_dependent), ]

rows <- lapply(seq_len(nrow(ctx)), function(i) {
  key <- paste0(ctx$ligand[i], "_", ctx$receptor[i])
  lineage <- recv_map[[key]]
  if (is.null(lineage) || lineage == "skip")
    return(data.frame(ligand=ctx$ligand[i], receptor=ctx$receptor[i],
                      prior_sign=ctx$sign[i], receiver=NA, act=NA, supp=NA,
                      data_sign=NA, agree=NA, stringsAsFactors=FALSE))
  if (lineage == "Mac") {
    recv <- "Macrophage"; a <- score(recv,panels$Mac_m1); s <- score(recv,panels$Mac_m2)
  } else {
    recv <- if (lineage=="NK") "NK" else "CD8"
    a <- score(recv,panels$Tnk_activation); s <- score(recv,panels$Tnk_exhaustion)
  }
  ds <- if (is.na(a)||is.na(s)) NA_integer_ else if (a>=s) 1L else -1L
  data.frame(ligand=ctx$ligand[i], receptor=ctx$receptor[i],
             prior_sign=ctx$sign[i], receiver=recv,
             act=round(a,3), supp=round(s,3), data_sign=ds,
             agree=ifelse(is.na(ds), NA, ds==ctx$sign[i]), stringsAsFactors=FALSE)
})
out <- do.call(rbind, rows)

cat("=== Receiver state scores (z-scored marker means) ===\n")
cat("CD8        : activation =", round(score("CD8",panels$Tnk_activation),3),
    " exhaustion =", round(score("CD8",panels$Tnk_exhaustion),3), "\n")
cat("NK         : activation =", round(score("NK",panels$Tnk_activation),3),
    " exhaustion =", round(score("NK",panels$Tnk_exhaustion),3), "\n")
cat("Macrophage : M1 =", round(score("Macrophage",panels$Mac_m1),3),
    " M2 =", round(score("Macrophage",panels$Mac_m2),3), "\n")

cat("\n=== Prior sign vs data-state sign (context-dependent pairs) ===\n")
print(out, row.names = FALSE)

cat("\nSummary: agree =", sum(out$agree, na.rm=TRUE),
    "/ disagree =", sum(out$agree==FALSE, na.rm=TRUE),
    "/ NA =", sum(is.na(out$agree)), "\n")
cat("\nDisagreements (prior sign not supported by receiver state in this dataset):\n")
print(out[which(out$agree==FALSE), c("ligand","receptor","prior_sign","receiver","act","supp","data_sign")],
      row.names = FALSE)

write.csv(out, file.path(outdir, "context_sign_marker_check.csv"), row.names = FALSE)
cat("\nWritten to", file.path(outdir, "context_sign_marker_check.csv"), "\n")
