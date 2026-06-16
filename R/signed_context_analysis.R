# signed_context_analysis.R
# Shared analysis used for BOTH datasets so the methodology is identical:
#   global single-sign model vs per-receiver context-adjusted model, with a
#   sign-permutation null (N perms) for the net-positive indirect fraction and
#   the SignedPageRank net. Context rule: activating (+1) context-dependent pairs
#   flip to -1 when the receiver cell type is in a suppressed/exhausted state
#   (markers). See experiments/real/context_signs.R for the original (melanoma).

suppressMessages({
  source("R/signed_lr_utils.R")
  source("R/signed_lrbase.R")
})

#' @param expr_ct gene x cell type mean-expression matrix
#' @param cell_names cell-type names (columns of expr_ct)
#' @param outdir directory for CSV outputs
#' @param label short tag used in messages/filenames
#' @param n_perm number of sign permutations
#' @param seed RNG seed
run_signed_context_analysis <- function(expr_ct, cell_names, outdir,
                                        label = "data", n_perm = 2000L, seed = 1L) {
  set.seed(seed)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  n <- length(cell_names)

  lr_all <- load_signed_lrbase()
  present <- lr_all$ligand %in% rownames(expr_ct) & lr_all$receptor %in% rownames(expr_ct)
  lr <- lr_all[present, ]; K <- nrow(lr)
  mats <- build_lr_matrices(lr, cell_names, expr = expr_ct)
  lig <- mats$lig_expr; rec <- mats$rec_expr
  sign0 <- as.integer(lr$sign); ctx <- as.logical(lr$context_dependent)
  cat(sprintf("[%s] %d cell types; %d/%d signedLRBase pairs usable\n",
              label, n, K, nrow(lr_all)))

  # Receiver polarity from markers (activated +1 / suppressed -1; NA = unknown).
  z <- t(scale(t(expr_ct))); z[is.na(z)] <- 0
  panel <- list(
    act = c("GZMA","GZMB","GZMK","PRF1","NKG7","GNLY","KLRD1","IFNG"),
    exh = c("PDCD1","HAVCR2","LAG3","TIGIT","CTLA4","TOX","ENTPD1","BTLA"),
    m1  = c("TNF","IL1B","CXCL9","CXCL10","CXCL11","CD86","IL12B","NOS2"),
    m2  = c("MRC1","CD163","MSR1","IL10","ARG1","MARCO","CCL22"))
  sc <- function(ct, gg) { gg <- intersect(gg, rownames(z))
    if (!length(gg) || !ct %in% colnames(z)) NA else mean(z[gg, ct]) }
  polarity <- setNames(rep(NA_integer_, n), cell_names)
  for (ct in intersect(c("CD8","CD4conv","Treg","NK"), cell_names))
    polarity[ct] <- if (sc(ct, panel$act) >= sc(ct, panel$exh)) 1L else -1L
  if ("Macrophage" %in% cell_names)
    polarity["Macrophage"] <- if (sc("Macrophage", panel$m1) >= sc("Macrophage", panel$m2)) 1L else -1L

  ctx_sign <- function(perm_sign) {
    s <- matrix(perm_sign, nrow = K, ncol = n)
    for (j in seq_len(n)) {
      pj <- polarity[cell_names[j]]; if (is.na(pj)) next
      s[which(ctx & perm_sign == 1L & pj == -1L), j] <- -1L
    }
    s
  }
  build_AA <- function(sgn_mat) {
    Wpos <- t(rec) * (sgn_mat == 1L); Wneg <- t(rec) * (sgn_mat == -1L)
    list(Apos = lig %*% Wpos, Aneg = lig %*% Wneg)
  }
  stat <- function(A) {
    AposT <- t(A$Apos); AnegT <- t(A$Aneg)
    D <- colSums(AposT) + colSums(AnegT); D[D == 0] <- 1
    Pp <- sweep(AposT, 2, D, "/"); Pn <- sweep(AnegT, 2, D, "/")
    M <- rbind(cbind(Pp, Pn), cbind(Pn, Pp)); fr <- numeric(2); ML <- diag(2 * n)
    for (L in 1:3) { ML <- ML %*% M
      if (L >= 2) fr[L-1] <- sum(ML[1:n,1:n]) / (sum(ML[1:n,1:n]) + sum(ML[(n+1):(2*n),1:n])) }
    list(fr = fr, net = SignedPageRank(t(A$Apos), t(A$Aneg), damping = 0.85)$net)
  }

  glob <- stat(build_AA(matrix(sign0, K, n)))
  cont <- stat(build_AA(ctx_sign(sign0)))

  # NEED: flips
  sm <- ctx_sign(sign0)
  flips <- which(matrix(sign0, K, n) == 1L & sm == -1L, arr.ind = TRUE)
  pair_names <- paste0(lr$ligand, "_", lr$receptor)

  # Permutation null (both models; null re-applies the context rule).
  p_two <- function(o, nl) { pu <- (1+sum(nl>=o))/(length(nl)+1); pl <- (1+sum(nl<=o))/(length(nl)+1); min(1, 2*min(pu,pl)) }
  nullG <- matrix(NA, n_perm, 2); nullC <- matrix(NA, n_perm, 2)
  nGnet <- nCnet <- matrix(NA, n_perm, n, dimnames = list(NULL, cell_names))
  for (p in seq_len(n_perm)) {
    s <- sample(sign0)
    g <- stat(build_AA(matrix(s, K, n))); c2 <- stat(build_AA(ctx_sign(s)))
    nullG[p,] <- g$fr; nullC[p,] <- c2$fr; nGnet[p,] <- g$net; nCnet[p,] <- c2$net
  }
  frac <- data.frame(model = c("global","context"),
    fracpos_hop2 = round(c(glob$fr[1], cont$fr[1]), 4),
    null_hop2 = round(c(mean(nullG[,1]), mean(nullC[,1])), 4),
    p_hop2 = c(p_two(glob$fr[1], nullG[,1]), p_two(cont$fr[1], nullC[,1])),
    fracpos_hop3 = round(c(glob$fr[2], cont$fr[2]), 4),
    null_hop3 = round(c(mean(nullG[,2]), mean(nullC[,2])), 4),
    p_hop3 = c(p_two(glob$fr[2], nullG[,2]), p_two(cont$fr[2], nullC[,2])))
  net <- data.frame(cell = cell_names,
    net_global = round(glob$net, 4), q_global = p.adjust(sapply(seq_len(n), function(i) p_two(glob$net[i], nGnet[,i])), "BH"),
    net_context = round(cont$net, 4), q_context = p.adjust(sapply(seq_len(n), function(i) p_two(cont$net[i], nCnet[,i])), "BH"))
  net <- net[order(net$q_context), ]

  cat(sprintf("\n[%s] Receiver polarity (activated +1 / suppressed -1):\n", label)); print(polarity)
  cat(sprintf("\n[%s] NEED: %d activating context-edge signs flipped to - by suppressed receivers\n",
              label, nrow(flips)))
  if (nrow(flips) > 0) { print(table(receiver = cell_names[flips[,2]]))
    cat("flipped pairs:", paste(unique(pair_names[flips[,1]]), collapse = ", "), "\n") }
  cat(sprintf("\n[%s] IMPROVEMENT: net-positive indirect fraction vs sign-perm null\n", label))
  print(frac, row.names = FALSE)
  cat(sprintf("\n[%s] SignedPageRank net significance (BH q), global vs context\n", label))
  print(net, row.names = FALSE)

  write.csv(frac, file.path(outdir, "context_vs_global_fracpos.csv"), row.names = FALSE)
  write.csv(net,  file.path(outdir, "context_vs_global_net.csv"), row.names = FALSE)
  invisible(list(frac = frac, net = net, polarity = polarity, n_flips = nrow(flips)))
}
