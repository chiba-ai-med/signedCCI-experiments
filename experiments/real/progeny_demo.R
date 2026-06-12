# progeny_demo.R
# Demonstrate the GENERALISABLE (arbitrary-data) half of data-driven sign:
# a PRE-BUILT reference model (PROGENy) measures signed pathway activity on ANY
# dataset with no per-dataset marker curation. We then map pathway activity to
# the functional LR sign with a THIN prior (the residual knowledge step that
# cannot be removed: "pathway active" != "good for the cell").
#
# Pipeline step mapping:
#   (1) which genes to read   -> PROGENy's pre-built signed gene weights (not hand rules)
#   (2) measure activity       -> progeny() on this dataset (dataset-agnostic)
#   (3) activity -> +/- sign   -> thin prior table below (still needs knowledge)
#
# Run: Rscript experiments/real/progeny_demo.R

suppressMessages({
  library(progeny)
  source("R/gse72056.R")
})

infile <- "data/GSE72056_melanoma_single_cell_revised_v2.txt.gz"
outdir <- "results/real_signedLRBase"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

ge <- load_gse72056_celltype_expr(infile, min_cells = 20L)
expr_ct <- ge$expr_ct                       # genes x cell type (mean log-expr)

# ---- (1)+(2) Pre-built PROGENy model applied to THIS dataset --------------
# Pre-trained signed responsive-gene weights for 14 pathways; works on any human
# expression matrix. scale=TRUE -> relative activity across the cell types here.
cat("Running PROGENy (pre-built Human model) on per-cell-type expression ...\n")
act <- progeny(as.matrix(expr_ct), scale = TRUE, organism = "Human",
               top = 100, perm = 1)          # cell type x pathway
cat("\n=== PROGENy pathway activity (cell type x pathway, z-scaled) ===\n")
print(round(act, 2))

# ---- (3) THIN prior: pathway -> LR ligands + functional sign --------------
# This is the irreducible knowledge step. Each PROGENy pathway is tied to the
# signedLRBase ligands that drive it and to the functional sign of those pairs.
prior <- data.frame(
  pathway = c("TGFb","Trail","EGFR","VEGF","JAK-STAT","TNFa","NFkB"),
  lr_ligands = c("TGFB1/TGFB2/TGFB3","TNFSF10","EGF/TGFA/HBEGF/AREG","VEGFA",
                 "IFNG/IL2/IL12/IL15","TNF","TNF/IL1B"),
  functional_sign = c(-1,-1, 1, -1, 1, NA, 1),   # NA = genuinely context (TNFa: NFkB+ vs death-)
  note = c("suppressive cytokine","death ligand","growth/proliferation",
           "VEGFA/FLT1 immunosuppressive (myeloid)","activating cytokines",
           "context: NFkB-survival(+) vs apoptosis(-)","inflammatory activation"),
  stringsAsFactors = FALSE
)

# Pull each mapped pathway's activity across cell types.
get_path <- function(p) if (p %in% colnames(act)) round(act[, p], 2) else rep(NA, nrow(act))
cat("\n=== Mapping pre-built pathway activity -> LR functional sign (per cell type) ===\n")
for (i in seq_len(nrow(prior))) {
  p <- prior$pathway[i]
  vals <- get_path(p)
  cat(sprintf("\n[%s pathway]  ligands=%s  prior_sign=%s  (%s)\n",
              p, prior$lr_ligands[i],
              ifelse(is.na(prior$functional_sign[i]), "context", prior$functional_sign[i]),
              prior$note[i]))
  ord <- order(-vals)
  print(vals[ord])
}

# ---- Focused read-outs: where is each signed pathway most active? ---------
cat("\n=== Interpretation: top cell type per signed pathway (auto) ===\n")
for (i in seq_len(nrow(prior))) {
  p <- prior$pathway[i]
  if (!p %in% colnames(act)) next
  v <- act[, p]; top <- names(v)[which.max(v)]
  sgn <- ifelse(is.na(prior$functional_sign[i]), "context",
                ifelse(prior$functional_sign[i] > 0, "+1 (activating)", "-1 (suppressive/death)"))
  cat(sprintf("%-9s highest in %-11s (%.2f) -> sign %s  [%s]\n",
              p, top, max(v), sgn, prior$note[i]))
}

write.csv(round(act, 4), file.path(outdir, "progeny_pathway_activity.csv"))
write.csv(prior, file.path(outdir, "progeny_pathway_sign_prior.csv"), row.names = FALSE)
cat("\nWritten progeny_pathway_activity.csv + progeny_pathway_sign_prior.csv to", outdir, "\n")
cat("\nNote: PROGENy's model is PRE-BUILT, so steps (1)+(2) run on ANY human\n")
cat("dataset unchanged -- no per-dataset marker curation. Only step (3),\n")
cat("activity->sign, keeps a thin irreducible prior.\n")
