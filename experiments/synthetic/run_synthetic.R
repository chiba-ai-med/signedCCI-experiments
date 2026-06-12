# run_synthetic.R
# Toy immune x tumor signed cell-cell interaction (CCI) network.
# GATE TEST: verifies that the signed walkers reproduce (-)x(-)=(+)
# ("enemy of my enemy is my friend") on a known synthetic network.
#
# Run:
#   Rscript experiments/synthetic/run_synthetic.R

suppressMessages(source("R/signed_lr_utils.R"))

set.seed(1)
outdir <- "results/synthetic"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- Cell types ----------------------------------------------------------
cell_names <- c("Tumor", "Treg", "CD8", "NK", "DC")

# ---- Signed LR table (biological labels; sign = +1 activate / -1 inhibit) --
lr_table <- data.frame(
  sender   = c("Tumor","Treg","Treg","Tumor","Tumor","CD8","NK","DC","DC","CD8"),
  receiver = c("Treg","CD8","CD8","CD8","CD8","Tumor","Tumor","CD8","CD8","Tumor"),
  ligand   = c("TGFB1","IL10","TGFB1","CD274","FASLG","FASLG","FASLG","CD80","IL12A","IFNG"),
  receptor = c("TGFBR1","IL10RA","TGFBR1","PDCD1","FAS","FAS","FAS","CD28","IL12RB1","IFNGR1"),
  sign     = c( 1,      -1,    -1,     -1,     -1,    -1,    -1,    1,     1,       1),
  stringsAsFactors = FALSE
)
lr_table$pair_name <- paste0(lr_table$sender, "_", lr_table$ligand, "_",
                             lr_table$receptor, "_", lr_table$receiver)

cat("=== Signed LR table ===\n")
print(lr_table[, c("sender","receiver","ligand","receptor","sign")])

# ---- Build (cell type x LR pair) matrices (synthetic: isolate each edge) ---
m <- build_lr_matrices(lr_table, cell_names, expr = NULL)

# ---- Full signed-CCI pipeline --------------------------------------------
res <- signed_cci(m$lig_expr, m$rec_expr, m$lr_sign, cell_names,
                  max_hop = 3L, damping = 0.85)

A_pos <- res$cci$A_pos   # sender x receiver
A_neg <- res$cci$A_neg   # sender x receiver

cat("\n=== A_pos (sender row -> receiver col) ===\n"); print(A_pos)
cat("\n=== A_neg (sender row -> receiver col) ===\n"); print(A_neg)

# ---- Sanity on the constructed graph (BuildSignedCCI layout) -------------
stopifnot(A_pos["Tumor","Treg"]  > 0)   # Tumor -(+)-> Treg
stopifnot(A_neg["Treg","CD8"]    > 0)   # Treg  -(-)-> CD8
stopifnot(A_neg["Tumor","CD8"]   > 0)   # Tumor -(-)-> CD8
stopifnot(A_neg["CD8","Tumor"]   > 0)   # CD8   -(-)-> Tumor (FASLG)
stopifnot(A_pos["CD8","Tumor"]   > 0)   # CD8   -(+)-> Tumor (IFNG)
stopifnot(A_neg["NK","Tumor"]    > 0)   # NK    -(-)-> Tumor
stopifnot(A_pos["DC","CD8"]      > 0)   # DC    -(+)-> CD8
cat("\n[OK] A_pos/A_neg match the designed directed signed edges.\n")

# ---- Signed PageRank -----------------------------------------------------
pr <- res$pagerank
cat("\n=== Signed PageRank (converged =", pr$converged,
    ", iter =", pr$iter, ") ===\n")
pr_df <- data.frame(cell = cell_names,
                    positive = round(pr$positive, 5),
                    negative = round(pr$negative, 5),
                    net = round(pr$net, 5),
                    total = round(pr$total, 5))
pr_df <- pr_df[order(-pr_df$net), ]
print(pr_df, row.names = FALSE)

# ---- Signed paths --------------------------------------------------------
paths <- annotate_paths(res$paths, res$edge_table, cell_names)

get_path <- function(p, name, sign_str = NULL) {
  sel <- p[p$path_name == name, , drop = FALSE]
  if (!is.null(sign_str)) sel <- sel[sel$signs == sign_str, , drop = FALSE]
  sel
}

# (A) Required: 2-hop Treg -(-)-> CD8 -(-)-> Tumor has net_sign == +1.
p_tct <- get_path(paths, "Treg->CD8->Tumor", "-,-")
cat("\n=== [REQUIRED] 2-hop  Treg -(-)-> CD8 -(-)-> Tumor ===\n")
print(p_tct[, c("path_name","signs","net_sign","contribution","lr_annotation")],
      row.names = FALSE)
stopifnot(nrow(p_tct) == 1, p_tct$net_sign == 1L)
cat("[OK] (-)x(-) = (+):  enemy of my enemy is my friend.\n")

# (B) Required: 3-hop Tumor -(+)-> Treg -(-)-> CD8 -(-)-> Tumor net_sign == +1.
p_3 <- get_path(paths, "Tumor->Treg->CD8->Tumor", "+,-,-")
cat("\n=== [REQUIRED] 3-hop  Tumor -(+)-> Treg -(-)-> CD8 -(-)-> Tumor ===\n")
print(p_3[, c("path_name","signs","net_sign","contribution","lr_annotation")],
      row.names = FALSE)
stopifnot(nrow(p_3) == 1, p_3$net_sign == 1L)
cat("[OK] (+)x(-)x(-) = (+):  self-reinforcing anti-tumor loop.\n")

# (C) Sanity: a single (-) hop is net (-), a single (+) hop is net (+).
p1_neg <- get_path(paths, "Treg->CD8", "-")
p1_pos <- get_path(paths, "Tumor->Treg", "+")
stopifnot(nrow(p1_neg) >= 1, all(p1_neg$net_sign == -1L))
stopifnot(nrow(p1_pos) >= 1, all(p1_pos$net_sign == 1L))
cat("\n[OK] single-hop sign sanity:  (-) stays -,  (+) stays +.\n")

# ---- Show all even-negative (net positive) multi-hop indirect paths ------
ind <- paths[paths$hop >= 2 & paths$net_sign == 1L, ]
ind <- ind[order(-ind$contribution), ]
cat("\n=== Net-POSITIVE indirect paths (hop>=2), top by contribution ===\n")
print(head(ind[, c("path_name","signs","hop","net_sign","contribution")], 12),
      row.names = FALSE)

# ---- Persist results -----------------------------------------------------
write.csv(A_pos, file.path(outdir, "A_pos.csv"))
write.csv(A_neg, file.path(outdir, "A_neg.csv"))
write.csv(pr_df, file.path(outdir, "pagerank.csv"), row.names = FALSE)
write.csv(paths, file.path(outdir, "signed_paths.csv"), row.names = FALSE)
write.csv(res$edge_table, file.path(outdir, "edge_table.csv"), row.names = FALSE)

cat("\n=== ALL GATE ASSERTIONS PASSED ===\n")
cat("Results written to", outdir, "\n")
