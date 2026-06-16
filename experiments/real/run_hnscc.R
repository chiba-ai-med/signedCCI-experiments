# run_hnscc.R
# Reproduce the signed-CCI context analysis on a SECOND dataset:
# GSE103322 (Puram et al. 2017, head & neck squamous cell carcinoma).
# Uses the identical shared analysis as melanoma (R/signed_context_analysis.R).
#
# Run: Rscript experiments/real/run_hnscc.R

source("R/gse103322.R")
source("R/signed_context_analysis.R")

ge <- load_gse103322_celltype_expr("data/GSE103322_HNSCC_all_data.txt.gz", min_cells = 20L)
cat("=== HNSCC (GSE103322) cell-type counts ===\n"); print(ge$counts)

run_signed_context_analysis(ge$expr_ct, ge$cell_names,
                            outdir = "results/real_hnscc", label = "HNSCC",
                            n_perm = 2000L)
