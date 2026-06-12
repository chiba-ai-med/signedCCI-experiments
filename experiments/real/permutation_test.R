# permutation_test.R
# Significance of the signed-CCI structure on GSE72056 by a SIGN-PERMUTATION
# null model: keep expression, the LR scores X_ijk, and the weighted topology
# fixed; randomly reassign which LR pairs are +1 vs -1 (preserving their counts).
# Out-degree D_i = sum_j (A+ + A-)_ij is sign-invariant, so the unsigned graph is
# preserved and only the SIGN effect is tested.
#
# Path-contribution masses are computed via matrix powers of the doubled
# transition matrix M = [[P+,P-],[P-,P+]] (validated against the enumeration in
# symTensor::SignedPathContribution), which makes thousands of permutations cheap.
#
# Run: Rscript experiments/real/permutation_test.R

suppressMessages({
  source("R/signed_lr_utils.R")
  source("R/signed_lrbase.R")
  source("R/gse72056.R")
})

set.seed(1)
N_PERM <- 2000L
infile <- "data/GSE72056_melanoma_single_cell_revised_v2.txt.gz"
outdir <- "results/real_signedLRBase"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

ge <- load_gse72056_celltype_expr(infile, min_cells = 20L)
expr_ct <- ge$expr_ct
cell_names <- ge$cell_names
n <- length(cell_names)

lr_all <- load_signed_lrbase()
present <- lr_all$ligand %in% rownames(expr_ct) & lr_all$receptor %in% rownames(expr_ct)
lr <- lr_all[present, ]
K <- nrow(lr)
mats <- build_lr_matrices(lr, cell_names, expr = expr_ct)
lig <- mats$lig_expr; rec <- mats$rec_expr     # cell type x LR pair
sign0 <- as.integer(lr$sign)

# Per-LR score matrices Mk[i,j] = lig[i,k]*rec[j,k] (sender x receiver), flattened.
Xflat <- vapply(seq_len(K), function(k) as.vector(outer(lig[, k], rec[, k])),
                numeric(n * n))               # (n*n) x K

# Build doubled transition matrix M (2n x 2n) for a given sign vector.
# A_pos/A_neg are sender x receiver; walkers use A[to,from] = transpose.
build_M <- function(sgn) {
  Apos <- matrix(Xflat %*% (sgn == 1L),  n, n)
  Aneg <- matrix(Xflat %*% (sgn == -1L), n, n)
  AposT <- t(Apos); AnegT <- t(Aneg)           # to x from
  D <- colSums(AposT) + colSums(AnegT); D[D == 0] <- 1
  Pp <- sweep(AposT, 2, D, "/"); Pn <- sweep(AnegT, 2, D, "/")
  list(M = rbind(cbind(Pp, Pn), cbind(Pn, Pp)), Pp = Pp, Pn = Pn,
       Apos = Apos, Aneg = Aneg)
}

# Statistics from a sign vector.
stats <- function(sgn) {
  b <- build_M(sgn)
  # SignedPageRank net per cell type (use the package function for consistency).
  pr <- SignedPageRank(t(b$Apos), t(b$Aneg), damping = 0.85)
  net <- pr$net
  # Net-positive / negative path mass via matrix powers (start in + state).
  Mp <- b$M
  netpos <- netneg <- numeric(0)
  ML <- diag(2 * n)
  pos_pair_L2 <- NULL
  masses <- list()
  for (L in 1:3) {
    ML <- ML %*% Mp
    pos_block <- ML[1:n, 1:n]        # (+,t) <- (+,s): net-positive s->t mass
    neg_block <- ML[(n + 1):(2 * n), 1:n]
    masses[[L]] <- list(pos = sum(pos_block), neg = sum(neg_block))
    if (L == 2L) pos_pair_L2 <- pos_block   # source (col) -> target (row), net+
  }
  list(net = net, masses = masses, pos_pair_L2 = pos_pair_L2)
}

# ---- Observed -------------------------------------------------------------
obs <- stats(sign0)

# Cross-check M^L masses vs enumeration (observed only).
chk <- signed_cci(lig, rec, sign0, cell_names, max_hop = 3L)$paths
enum_pos2 <- sum(chk$contribution[chk$hop == 2 & chk$net_sign == 1L])
cat(sprintf("Validation: net-positive hop-2 mass  enumeration=%.5f  matrixpower=%.5f\n",
            enum_pos2, obs$masses[[2]]$pos))
enum_pos3 <- sum(chk$contribution[chk$hop == 3 & chk$net_sign == 1L])
cat(sprintf("Validation: net-positive hop-3 mass  enumeration=%.5f  matrixpower=%.5f\n",
            enum_pos3, obs$masses[[3]]$pos))

obs_fracpos <- sapply(2:3, function(L)
  obs$masses[[L]]$pos / (obs$masses[[L]]$pos + obs$masses[[L]]$neg))
names(obs_fracpos) <- c("hop2", "hop3")

# ---- Permutations ---------------------------------------------------------
cat(sprintf("\nRunning %d sign permutations ...\n", N_PERM))
null_net   <- matrix(NA_real_, N_PERM, n, dimnames = list(NULL, cell_names))
null_frac  <- matrix(NA_real_, N_PERM, 2, dimnames = list(NULL, c("hop2","hop3")))
null_pp2   <- array(NA_real_, c(N_PERM, n, n))
for (p in seq_len(N_PERM)) {
  s <- sample(sign0)
  st <- stats(s)
  null_net[p, ] <- st$net
  null_frac[p, ] <- sapply(2:3, function(L)
    st$masses[[L]]$pos / (st$masses[[L]]$pos + st$masses[[L]]$neg))
  null_pp2[p, , ] <- st$pos_pair_L2
}

# ---- p-values -------------------------------------------------------------
p_two <- function(obs_v, null_v) {
  pu <- (1 + sum(null_v >= obs_v)) / (length(null_v) + 1)
  pl <- (1 + sum(null_v <= obs_v)) / (length(null_v) + 1)
  min(1, 2 * min(pu, pl))
}
p_up  <- function(obs_v, null_v) (1 + sum(null_v >= obs_v)) / (length(null_v) + 1)

# (1) Per-cell-type SignedPageRank net (two-sided).
net_tab <- data.frame(
  cell = cell_names, net = round(obs$net, 5),
  null_mean = round(colMeans(null_net), 5),
  z = round((obs$net - colMeans(null_net)) / apply(null_net, 2, sd), 2),
  p = sapply(seq_len(n), function(i) p_two(obs$net[i], null_net[, i])))
net_tab$p_adj <- p.adjust(net_tab$p, "BH")
net_tab <- net_tab[order(net_tab$p_adj, -abs(net_tab$z)), ]

# (2) Net-positive indirect fraction (two-sided).
frac_tab <- data.frame(
  hop = c("hop2","hop3"),
  obs_frac_netpos = round(obs_fracpos, 4),
  null_mean = round(colMeans(null_frac), 4),
  p = c(p_two(obs_fracpos[1], null_frac[,1]), p_two(obs_fracpos[2], null_frac[,2])))

# (3) Directed net-positive 2-hop routes (source -> target), one-sided upper + BH.
pp <- expand.grid(target = seq_len(n), source = seq_len(n))
pp$obs <- mapply(function(t, s) obs$pos_pair_L2[t, s], pp$target, pp$source)
pp$null_mean <- mapply(function(t, s) mean(null_pp2[, t, s]), pp$target, pp$source)
pp$p <- mapply(function(t, s) p_up(obs$pos_pair_L2[t, s], null_pp2[, t, s]),
               pp$target, pp$source)
pp <- pp[pp$source != pp$target, ]
pp$p_adj <- p.adjust(pp$p, "BH")
pp$route <- paste0(cell_names[pp$source], " ~~> ", cell_names[pp$target])
pp <- pp[order(pp$p_adj, -pp$obs), ]

# ---- Report ---------------------------------------------------------------
cat("\n=== (1) SignedPageRank net vs sign-permutation null (two-sided) ===\n")
print(net_tab, row.names = FALSE)

cat("\n=== (2) Net-positive indirect connectivity fraction vs null ===\n")
print(frac_tab, row.names = FALSE)

cat("\n=== (3) Significant net-positive 2-hop routes (BH q<0.05), top 20 ===\n")
sig <- pp[pp$p_adj < 0.05, c("route","obs","null_mean","p","p_adj")]
sig$obs <- round(sig$obs, 5); sig$null_mean <- round(sig$null_mean, 5)
sig$p_adj <- signif(sig$p_adj, 3)
if (nrow(sig) > 0) print(head(sig, 20), row.names = FALSE) else cat("(none)\n")
cat(sprintf("\n%d of %d directed 2-hop routes significant (BH q<0.05).\n",
            sum(pp$p_adj < 0.05), nrow(pp)))

write.csv(net_tab, file.path(outdir, "perm_pagerank_net.csv"), row.names = FALSE)
write.csv(frac_tab, file.path(outdir, "perm_netpos_fraction.csv"), row.names = FALSE)
write.csv(pp[, c("route","obs","null_mean","p","p_adj")],
          file.path(outdir, "perm_netpos_routes.csv"), row.names = FALSE)
cat("\nWritten perm_*.csv to", outdir, "\n")
