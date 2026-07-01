# run_immunotherapy.R
# Lever B: does the signed immune-immune communication structure differ between
# checkpoint-immunotherapy RESPONDERS and NON-RESPONDERS? (GSE120575, no tumor
# cells.) Headline statistic = "suppression index" = share of total signed edge
# weight that is negative. Significance = per-patient index, Mann-Whitney R vs NR.
#
# Run: Rscript experiments/real/run_immunotherapy.R

suppressMessages({
  source("R/gse120575.R"); source("R/signed_lrbase.R"); source("R/signed_lr_utils.R")
})
set.seed(1)
outdir <- "results/real_immunotherapy"; dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

lr_all <- load_signed_lrbase()
d <- if (file.exists("/tmp/gse120575_d.rds")) readRDS("/tmp/gse120575_d.rds") else
       load_gse120575(extra_genes = unique(c(lr_all$ligand, lr_all$receptor)))

# Signed statistics for a per-cell-type expression matrix.
signed_stats <- function(expr_ct, cell_names) {
  n <- length(cell_names)
  present <- lr_all$ligand %in% rownames(expr_ct) & lr_all$receptor %in% rownames(expr_ct)
  lr <- lr_all[present, ]; K <- nrow(lr)
  m <- build_lr_matrices(lr, cell_names, expr = expr_ct)
  lig <- m$lig_expr; rec <- m$rec_expr; s <- as.integer(lr$sign)
  Ap <- matrix(0, n, n); An <- matrix(0, n, n)
  for (k in seq_len(K)) { Mk <- outer(lig[,k], rec[,k]); if (s[k]==1) Ap <- Ap+Mk else An <- An+Mk }
  supp <- sum(An) / (sum(Ap) + sum(An) + 1e-12)          # negative share
  ApT <- t(Ap); AnT <- t(An); D <- colSums(ApT)+colSums(AnT); D[D==0] <- 1
  Pp <- sweep(ApT,2,D,"/"); Pn <- sweep(AnT,2,D,"/")
  M <- rbind(cbind(Pp,Pn),cbind(Pn,Pp)); ML <- M %*% M
  fr2 <- sum(ML[1:n,1:n]) / (sum(ML[1:n,1:n]) + sum(ML[(n+1):(2*n),1:n]))
  list(supp = supp, fracpos_hop2 = fr2, n_types = n)
}

group_stat <- function(mask, min_cells = 20L) {
  ce <- celltype_mean_expr(d, mask, min_cells)
  if (length(ce$cell_names) < 3) return(NULL)
  c(signed_stats(ce$expr_ct, ce$cell_names), list(n_cells = sum(mask, na.rm=TRUE)))
}

# ---- Group-level: Responder vs Non-responder (all, and Pre-only) ----------
cat("=== Group-level signed suppression (global signs) ===\n")
grp <- list(
  "Responder (all)"     = !is.na(d$response) & d$response=="Responder",
  "Non-responder (all)" = !is.na(d$response) & d$response=="Non-responder",
  "Responder (Pre)"     = !is.na(d$response) & d$response=="Responder"     & d$timepoint=="Pre",
  "Non-responder (Pre)" = !is.na(d$response) & d$response=="Non-responder" & d$timepoint=="Pre")
gt <- lapply(grp, group_stat)
gtab <- do.call(rbind, lapply(names(gt), function(nm) if (is.null(gt[[nm]])) NULL else
  data.frame(group=nm, n_cells=gt[[nm]]$n_cells, n_types=gt[[nm]]$n_types,
             supp_index=round(gt[[nm]]$supp,4), fracpos_hop2=round(gt[[nm]]$fracpos_hop2,4))))
print(gtab, row.names = FALSE)

# ---- Per-patient suppression index + Mann-Whitney R vs NR -----------------
pts <- unique(d$patient)
prow <- lapply(pts, function(p) {
  mask <- d$patient == p & !is.na(d$response)
  if (sum(mask) < 80) return(NULL)
  st <- group_stat(mask, min_cells = 15L); if (is.null(st)) return(NULL)
  data.frame(patient=p, response=d$response[which(mask)[1]], n_cells=sum(mask),
             supp_index=st$supp, fracpos_hop2=st$fracpos_hop2)
})
pt <- do.call(rbind, prow)
cat(sprintf("\n=== Per-patient (>=80 cells, >=3 types): %d patients (%d R, %d NR) ===\n",
            nrow(pt), sum(pt$response=="Responder"), sum(pt$response=="Non-responder")))
pt <- pt[order(pt$response, -pt$supp_index), ]
print(pt, row.names = FALSE)

wtest <- function(col) {
  R  <- pt[[col]][pt$response=="Responder"]; NR <- pt[[col]][pt$response=="Non-responder"]
  w <- suppressWarnings(wilcox.test(NR, R))            # NR vs R
  cat(sprintf("\n%s: median NR=%.4f vs R=%.4f, Wilcoxon p=%.4f\n",
              col, median(NR), median(R), w$p.value))
}
wtest("supp_index"); wtest("fracpos_hop2")

write.csv(gtab, file.path(outdir, "group_suppression.csv"), row.names = FALSE)
write.csv(pt,   file.path(outdir, "per_patient_suppression.csv"), row.names = FALSE)
cat("\nWritten to", outdir, "\n")
