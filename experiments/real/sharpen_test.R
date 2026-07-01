# sharpen_test.R
# Lever A: is the flat permutation-test result real, or an artifact of the dense
# (mean-expression) graph? Sharpen the graph by keeping each gene's expression
# ONLY in cell types where it is specifically elevated (z-score across cell types
# >= tau; else set to 0), then re-run the SIGN-PERMUTATION test (global signs).
# Compare p-values / #significant across tau (tau=NA = dense baseline).
#
# Run: Rscript experiments/real/sharpen_test.R

suppressMessages({
  source("R/signed_lr_utils.R")
  source("R/signed_lrbase.R")
  source("R/gse72056.R")
})
set.seed(1)
N_PERM <- 1000L
ge <- load_gse72056_celltype_expr("data/GSE72056_melanoma_single_cell_revised_v2.txt.gz", 20L)
expr0 <- ge$expr_ct; cell_names <- ge$cell_names; n <- length(cell_names)
lr_all <- load_signed_lrbase()

# Specificity gate: keep expr where z across cell types >= tau, else 0.
gate <- function(expr, tau) {
  if (is.na(tau)) return(expr)
  z <- t(scale(t(expr))); z[is.na(z)] <- -Inf
  expr * (z >= tau)
}

# Sign-permutation stats for a given expression matrix (global signs).
run_one <- function(expr, tau) {
  present <- lr_all$ligand %in% rownames(expr) & lr_all$receptor %in% rownames(expr)
  lr <- lr_all[present, ]; K <- nrow(lr)
  m <- build_lr_matrices(lr, cell_names, expr = expr)
  lig <- m$lig_expr; rec <- m$rec_expr; sign0 <- as.integer(lr$sign)
  Xf <- vapply(seq_len(K), function(k) as.vector(outer(lig[,k], rec[,k])), numeric(n*n))
  stat <- function(sgn) {
    Ap <- matrix(Xf %*% (sgn==1L), n, n); An <- matrix(Xf %*% (sgn==-1L), n, n)
    ApT <- t(Ap); AnT <- t(An); D <- colSums(ApT)+colSums(AnT); D[D==0] <- 1
    Pp <- sweep(ApT,2,D,"/"); Pn <- sweep(AnT,2,D,"/")
    M <- rbind(cbind(Pp,Pn),cbind(Pn,Pp)); fr <- numeric(2); ML <- diag(2*n)
    for (L in 1:3){ML <- ML%*%M; if(L>=2) fr[L-1] <- sum(ML[1:n,1:n])/(sum(ML[1:n,1:n])+sum(ML[(n+1):(2*n),1:n]))}
    net <- SignedPageRank(t(Ap), t(An), damping=0.85)$net
    list(fr=fr, net=net)
  }
  obs <- stat(sign0)
  nf <- matrix(NA,N_PERM,2); nnet <- matrix(NA,N_PERM,n)
  for (p in seq_len(N_PERM)){ s <- sample(sign0); st <- stat(s); nf[p,] <- st$fr; nnet[p,] <- st$net }
  p_two <- function(o,v){pu<-(1+sum(v>=o))/(length(v)+1);pl<-(1+sum(v<=o))/(length(v)+1);min(1,2*min(pu,pl))}
  q <- p.adjust(sapply(seq_len(n), function(i) p_two(obs$net[i], nnet[,i])), "BH")
  Adense <- matrix(Xf %*% rep(1,K), n, n)   # total magnitude (sign-free) for sparsity
  data.frame(
    tau = ifelse(is.na(tau),"dense",as.character(tau)),
    usable_pairs = K,
    nonzero_edge_frac = round(mean(Adense > 0), 3),
    fracpos_hop2 = round(obs$fr[1],4), p_hop2 = round(p_two(obs$fr[1], nf[,1]),4),
    fracpos_hop3 = round(obs$fr[2],4), p_hop3 = round(p_two(obs$fr[2], nf[,2]),4),
    n_sig_cells_q05 = sum(q < 0.05))
}

cat("Sharpening test on melanoma (global signs, N=", N_PERM, " perms)\n", sep="")
res <- do.call(rbind, lapply(c(NA, 0.5, 1.0, 1.5), function(tau) run_one(gate(expr0, tau), tau)))
cat("\n=== Sign-permutation significance vs specificity gate (tau) ===\n")
print(res, row.names = FALSE)
cat("\nRead: does a stricter gate (sparser graph) make the net-positive indirect\n")
cat("fraction and/or per-cell net significant, i.e. was the flat result dilution?\n")
dir.create("results/real_signedLRBase", showWarnings = FALSE, recursive = TRUE)
write.csv(res, "results/real_signedLRBase/sharpen_test.csv", row.names = FALSE)
cat("\nWritten results/real_signedLRBase/sharpen_test.csv\n")
