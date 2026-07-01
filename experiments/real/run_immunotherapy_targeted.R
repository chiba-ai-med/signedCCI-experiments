# run_immunotherapy_targeted.R
# Targeted follow-up to the null global suppression index: test specific,
# biologically-motivated signed statistics between responders / non-responders
# (GSE120575), per patient, Mann-Whitney. Scale-robust (fractions / PageRank).
#
#   cd8_supp_frac : share of signals CD8 RECEIVES that are suppressive
#   treg_cd8_frac : Treg's suppressive contribution as a share of CD8's input
#   cd8_net       : CD8 SignedPageRank net (positive - negative), normalized
#   treg_supp_out : share of Treg's OUTPUT that is suppressive
#
# Run: Rscript experiments/real/run_immunotherapy_targeted.R

suppressMessages({
  source("R/gse120575.R"); source("R/signed_lrbase.R"); source("R/signed_lr_utils.R")
})
set.seed(1)
outdir <- "results/real_immunotherapy"; dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
lr_all <- load_signed_lrbase()
d <- if (file.exists("/tmp/gse120575_d.rds")) readRDS("/tmp/gse120575_d.rds") else
       load_gse120575(extra_genes = unique(c(lr_all$ligand, lr_all$receptor)))

targeted <- function(mask) {
  ct <- d$celltype[mask]; tab <- table(ct)
  keep <- names(tab)[tab >= 10]
  if (!all(c("Treg","CD8") %in% keep)) return(NULL)
  E <- sapply(keep, function(k) rowMeans(d$expr[, mask & d$celltype==k, drop=FALSE]))
  E <- as.matrix(E); cn <- keep; n <- length(cn)
  present <- lr_all$ligand %in% rownames(E) & lr_all$receptor %in% rownames(E)
  lr <- lr_all[present, ]; m <- build_lr_matrices(lr, cn, expr = E)
  lig <- m$lig_expr; rec <- m$rec_expr; s <- as.integer(lr$sign)
  Ap <- matrix(0,n,n,dimnames=list(cn,cn)); An <- Ap
  for (k in seq_len(nrow(lr))) { Mk <- outer(lig[,k], rec[,k]); if (s[k]==1) Ap<-Ap+Mk else An<-An+Mk }
  cd8_in_pos <- sum(Ap[,"CD8"]); cd8_in_neg <- sum(An[,"CD8"])
  pr <- SignedPageRank(t(Ap), t(An), damping=0.85); names(pr$net) <- cn
  data.frame(
    cd8_supp_frac = cd8_in_neg/(cd8_in_pos+cd8_in_neg+1e-12),
    treg_cd8_frac = An["Treg","CD8"]/(cd8_in_pos+cd8_in_neg+1e-12),
    cd8_net       = pr$net["CD8"],
    treg_supp_out = sum(An["Treg",])/(sum(Ap["Treg",])+sum(An["Treg",])+1e-12))
}

pts <- unique(d$patient)
rows <- lapply(pts, function(p) {
  mask <- d$patient==p & !is.na(d$response)
  if (sum(mask) < 80) return(NULL)
  st <- targeted(mask); if (is.null(st)) return(NULL)
  cbind(data.frame(patient=p, response=d$response[which(mask)[1]], n_cells=sum(mask)), st)
})
pt <- do.call(rbind, rows)
cat(sprintf("Patients with Treg & CD8 (>=10 each) and >=80 cells: %d (%d R, %d NR)\n",
            nrow(pt), sum(pt$response=="Responder"), sum(pt$response=="Non-responder")))

stats <- c("cd8_supp_frac","treg_cd8_frac","cd8_net","treg_supp_out")
res <- do.call(rbind, lapply(stats, function(s) {
  R <- pt[[s]][pt$response=="Responder"]; NR <- pt[[s]][pt$response=="Non-responder"]
  w <- suppressWarnings(wilcox.test(NR, R))
  data.frame(stat=s, median_NR=round(median(NR),4), median_R=round(median(R),4),
             direction=ifelse(median(NR)>median(R),"NR>R","NR<R"), p=w$p.value)
}))
res$p_adj <- p.adjust(res$p, "BH")
cat("\n=== Targeted per-patient comparison (Mann-Whitney NR vs R) ===\n")
print(transform(res, p=round(p,4), p_adj=round(p_adj,4)), row.names = FALSE)
cat("\n(Hypotheses: non-responders would have higher cd8_supp_frac / treg_cd8_frac /\n",
    "treg_supp_out and lower cd8_net.)\n")
write.csv(pt,  file.path(outdir,"per_patient_targeted.csv"), row.names=FALSE)
write.csv(res, file.path(outdir,"targeted_tests.csv"), row.names=FALSE)
cat("\nWritten to", outdir, "\n")
