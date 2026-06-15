# context_signs.R
# (Q1) Show hop-1/2/3 signed paths explicitly.
# (Q2) Make the sign context-dependent per receiver and test whether it is
#      (a) NEEDED and (b) an IMPROVEMENT, vs the global single-sign model.
#
# Context rule (the IFNG-on-exhausted-CD8 logic, made operational):
#   For context_dependent pairs whose prior is ACTIVATING (+1), the sign becomes
#   -1 when the RECEIVER cell type is in a suppressed/exhausted state (measured
#   from markers), because an activating signal does not land as "+" on an
#   exhausted cell. Suppressive (-1) priors are robust and kept. Non-context
#   pairs keep their prior sign. This turns the global sign s_k into s_{k,j}
#   (k = LR pair, j = receiver).
#
# Run: Rscript experiments/real/context_signs.R

suppressMessages({
  source("R/signed_lr_utils.R")
  source("R/signed_lrbase.R")
  source("R/gse72056.R")
})
set.seed(1)
N_PERM <- 2000L
outdir <- "results/real_signedLRBase"
ge <- load_gse72056_celltype_expr("data/GSE72056_melanoma_single_cell_revised_v2.txt.gz", 20L)
expr_ct <- ge$expr_ct; cell_names <- ge$cell_names; n <- length(cell_names)

lr_all <- load_signed_lrbase()
present <- lr_all$ligand %in% rownames(expr_ct) & lr_all$receptor %in% rownames(expr_ct)
lr <- lr_all[present, ]; K <- nrow(lr)
mats <- build_lr_matrices(lr, cell_names, expr = expr_ct)
lig <- mats$lig_expr; rec <- mats$rec_expr
sign0 <- as.integer(lr$sign)
ctx <- as.logical(lr$context_dependent)
pair_names <- paste0(lr$ligand, "_", lr$receptor)

# ---- (Q1) hop-1/2/3 paths (global signs) ---------------------------------
res0 <- signed_cci(lig, rec, sign0, cell_names, max_hop = 3L)
res0$edge_table$pair_name <- pair_names[res0$edge_table$lr_pair]
P <- annotate_paths(res0$paths, res0$edge_table, cell_names)
strip <- function(d) { d$lr_annotation <- gsub("\\[[^]]*\\]", "", d$lr_annotation); d }
cat("===== (Q1) Signed paths by hop (global signs) =====\n")
for (h in 1:3) {
  cat(sprintf("\n--- hop %d : top net-POSITIVE ---\n", h))
  d <- P[P$hop == h & P$net_sign == 1L, ]
  d <- strip(d[order(-d$contribution), ])
  print(head(d[, c("path_name","signs","contribution","lr_annotation")], 5), row.names = FALSE)
  cat(sprintf("--- hop %d : top net-NEGATIVE ---\n", h))
  d <- P[P$hop == h & P$net_sign == -1L, ]
  d <- strip(d[order(-d$contribution), ])
  print(head(d[, c("path_name","signs","contribution","lr_annotation")], 3), row.names = FALSE)
}

# ---- Receiver polarity from markers (activated +1 / suppressed -1) --------
z <- t(scale(t(expr_ct))); z[is.na(z)] <- 0
panel <- list(
  act = c("GZMA","GZMB","GZMK","PRF1","NKG7","GNLY","KLRD1","IFNG"),
  exh = c("PDCD1","HAVCR2","LAG3","TIGIT","CTLA4","TOX","ENTPD1","BTLA"),
  m1  = c("TNF","IL1B","CXCL9","CXCL10","CXCL11","CD86","IL12B","NOS2"),
  m2  = c("MRC1","CD163","MSR1","IL10","ARG1","MARCO","CCL22"))
sc <- function(ct, g) { g <- intersect(g, rownames(z)); if (!length(g)||!ct%in%colnames(z)) NA else mean(z[g,ct]) }
polarity <- setNames(rep(NA_integer_, n), cell_names)
for (ct in c("CD8","CD4conv","Treg","NK"))
  polarity[ct] <- if (sc(ct,panel$act) >= sc(ct,panel$exh)) 1L else -1L
polarity["Macrophage"] <- if (sc("Macrophage",panel$m1) >= sc("Macrophage",panel$m2)) 1L else -1L
cat("\n===== Receiver polarity (markers) =====\n")
print(polarity)

# ---- Build per-(pair,receiver) sign masks --------------------------------
# posmask[k,j] / negmask[k,j]: which (pair,receiver) contribute to A+ / A-.
ctx_sign <- function(perm_sign) {
  s <- matrix(perm_sign, nrow = K, ncol = n)            # default: global per pair
  for (j in seq_len(n)) {
    pj <- polarity[cell_names[j]]
    if (is.na(pj)) next
    flip <- which(ctx & perm_sign == 1L & pj == -1L)    # activating prior + suppressed receiver
    s[flip, j] <- -1L
  }
  s
}
build_AposAneg <- function(sgn_mat) {
  Wpos <- t(rec) * (sgn_mat == 1L)                      # K x n
  Wneg <- t(rec) * (sgn_mat == -1L)
  list(Apos = lig %*% Wpos, Aneg = lig %*% Wneg)        # n x n (sender x receiver)
}
# Doubled transition + net-positive path mass via matrix powers.
masses_and_net <- function(A) {
  AposT <- t(A$Apos); AnegT <- t(A$Aneg)
  D <- colSums(AposT) + colSums(AnegT); D[D==0] <- 1
  Pp <- sweep(AposT,2,D,"/"); Pn <- sweep(AnegT,2,D,"/")
  M <- rbind(cbind(Pp,Pn), cbind(Pn,Pp))
  fr <- numeric(2); ML <- diag(2*n)
  for (L in 1:3) { ML <- ML %*% M
    if (L>=2) fr[L-1] <- sum(ML[1:n,1:n]) / (sum(ML[1:n,1:n]) + sum(ML[(n+1):(2*n),1:n])) }
  pr <- SignedPageRank(t(A$Apos), t(A$Aneg), damping=0.85)
  list(fracpos = fr, net = pr$net)                      # fr = c(hop2,hop3)
}

glob <- masses_and_net(build_AposAneg(matrix(sign0, K, n)))
cont <- masses_and_net(build_AposAneg(ctx_sign(sign0)))

# ---- (Q2a) NEED: how many activating edges flip under context ------------
sm <- ctx_sign(sign0)
flips <- which(matrix(sign0,K,n)==1L & sm==-1L, arr.ind=TRUE)
cat("\n===== (Q2a) NEED: activating context-pair signs flipped to - by receiver state =====\n")
if (nrow(flips)>0) {
  ft <- data.frame(pair=pair_names[flips[,1]], receiver=cell_names[flips[,2]])
  cat(nrow(ft), "edge-signs flipped (activating signal onto a suppressed receiver).\n")
  print(table(receiver=ft$receiver))
  cat("flipped pairs (unique):", paste(unique(ft$pair), collapse=", "), "\n")
} else cat("none\n")

# ---- (Q2b) IMPROVEMENT: permutation test, global vs context --------------
p_up  <- function(o,nl) (1+sum(nl>=o))/(length(nl)+1)
p_two <- function(o,nl){pu<-(1+sum(nl>=o))/(length(nl)+1);pl<-(1+sum(nl<=o))/(length(nl)+1);min(1,2*min(pu,pl))}
nullG <- matrix(NA,N_PERM,2); nullC <- matrix(NA,N_PERM,2)
nullGnet <- nullCnet <- matrix(NA,N_PERM,n,dimnames=list(NULL,cell_names))
cat(sprintf("\nRunning %d permutations (both models) ...\n", N_PERM))
for (p in seq_len(N_PERM)) {
  s <- sample(sign0)
  g <- masses_and_net(build_AposAneg(matrix(s,K,n)))
  c2 <- masses_and_net(build_AposAneg(ctx_sign(s)))
  nullG[p,] <- g$fracpos; nullC[p,] <- c2$fracpos
  nullGnet[p,] <- g$net;  nullCnet[p,] <- c2$net
}
cmp <- data.frame(
  model = c("global","context"),
  fracpos_hop2 = round(c(glob$fracpos[1], cont$fracpos[1]),4),
  null_hop2    = round(c(mean(nullG[,1]), mean(nullC[,1])),4),
  p_hop2       = c(p_two(glob$fracpos[1],nullG[,1]), p_two(cont$fracpos[1],nullC[,1])),
  fracpos_hop3 = round(c(glob$fracpos[2], cont$fracpos[2]),4),
  null_hop3    = round(c(mean(nullG[,2]), mean(nullC[,2])),4),
  p_hop3       = c(p_two(glob$fracpos[2],nullG[,2]), p_two(cont$fracpos[2],nullC[,2])))
cat("\n===== (Q2b) IMPROVEMENT: net-positive indirect fraction vs sign-perm null =====\n")
print(cmp, row.names = FALSE)

netcmp <- data.frame(cell=cell_names,
  net_global=round(glob$net,4),  p_global=sapply(seq_len(n),function(i)p_two(glob$net[i],nullGnet[,i])),
  net_context=round(cont$net,4), p_context=sapply(seq_len(n),function(i)p_two(cont$net[i],nullCnet[,i])))
netcmp$q_global  <- p.adjust(netcmp$p_global,"BH")
netcmp$q_context <- p.adjust(netcmp$p_context,"BH")
netcmp <- netcmp[order(netcmp$q_context),]
cat("\nSignedPageRank net significance, global vs context (BH q):\n")
print(netcmp, row.names = FALSE)

write.csv(cmp, file.path(outdir,"context_vs_global_fracpos.csv"), row.names=FALSE)
write.csv(netcmp, file.path(outdir,"context_vs_global_net.csv"), row.names=FALSE)
cat("\nWritten context_vs_global_*.csv to", outdir, "\n")
