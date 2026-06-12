# build_signedLRBase.R
# Build a SIGNED ligand-receptor base for cell-cell interaction:
#   functional-outcome sign (+1 = activation/survival/proliferation,
#   -1 = immunosuppression/exhaustion/death), assigned by curated functional
#   rules, then CROSS-CHECKED against OmniPath's molecular signalling sign
#   (is_stimulation/is_inhibition).
#
# Key point: the molecular sign (does L activate or inhibit R's signalling) is
# NOT the same as the functional outcome sign. E.g. FASLG->FAS is molecularly
# AGONIST (+) but the cellular outcome is death (-). We keep the functional sign
# as authoritative and record agreement/conflict with OmniPath for transparency.
#
# Run: Rscript experiments/signedLRBase/build_signedLRBase.R

suppressMessages(library(data.table))

indir  <- "experiments/signedLRBase"
outdir <- "data/signedLRBase"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cur <- fread(file.path(indir, "signed_lr_curated.csv"))
cat("Curated pairs:", nrow(cur), "\n")
stopifnot(all(cur$sign %in% c(-1, 1)))

# ---- Fetch OmniPath signed signalling network ----------------------------
op_url <- paste0("https://omnipathdb.org/interactions",
                 "?genesymbols=yes&fields=is_stimulation,is_inhibition,curation_effort",
                 "&datasets=omnipath,ligrecextra")
cat("Fetching OmniPath signed interactions ...\n")
op <- tryCatch(
  fread(op_url, showProgress = FALSE),
  error = function(e) {
    cat("fread(URL) failed (", conditionMessage(e), "); trying download.file\n")
    tmp <- tempfile(fileext = ".tsv")
    download.file(op_url, tmp, quiet = TRUE)
    fread(tmp, showProgress = FALSE)
  }
)
cat("OmniPath interactions:", nrow(op), "rows;",
    "cols:", paste(head(colnames(op), 12), collapse = ", "), "\n")

# Molecular sign per (source_genesymbol -> target_genesymbol).
op[, mol_sign := fifelse(is_stimulation == 1 & is_inhibition == 0,  1L,
                  fifelse(is_inhibition == 1 & is_stimulation == 0, -1L, 0L))]
op_key <- op[, .(db_stim = max(is_stimulation),
                 db_inhib = max(is_inhibition),
                 db_sign = {
                   s <- max(is_stimulation); i <- max(is_inhibition)
                   if (s == 1 && i == 0) 1L else if (i == 1 && s == 0) -1L else 0L
                 }),
             by = .(ligand = source_genesymbol, receptor = target_genesymbol)]

# ---- Join curated functional sign with OmniPath molecular sign -----------
m <- merge(cur, op_key, by = c("ligand", "receptor"), all.x = TRUE, sort = FALSE)

m[, agreement := fifelse(is.na(db_sign), "db_absent",
                  fifelse(db_sign == 0, "db_ambiguous",
                   fifelse(db_sign == sign, "match", "conflict")))]

# Confidence / interpretation of the rule-vs-DB comparison.
# Inhibitory pathways work BY molecular agonism of a suppressive receptor, so a
# (functional = -1) vs (molecular = +1) conflict is EXPECTED there, not an error
# -- it is exactly why the molecular sign cannot be used directly. We mark those
# as expected_override; only conflicts OUTSIDE the inhibitory categories are
# genuinely worth review.
expected_override_cats <- c("checkpoint_inhibitory", "nk_inhibitory",
                            "suppressive_cytokine", "death")
m[, confidence := fifelse(agreement == "match", "high",
                   fifelse(agreement == "conflict" & category %in% expected_override_cats,
                           "expected_override",
                    fifelse(agreement == "conflict", "review", "curated_only")))]
m[, source := "curated+omnipath"]

# Final sign = functional (curated) sign, authoritative by design.
setcolorder(m, c("ligand","receptor","sign","category","context_dependent",
                 "db_sign","agreement","confidence","source","notes"))

fwrite(m, file.path(outdir, "signedLRBase.csv"))

# ---- Report --------------------------------------------------------------
cat("\n=== signedLRBase built:", nrow(m), "pairs ===\n")
cat("\nSign distribution:\n"); print(table(sign = m$sign))
cat("\nCategory x sign:\n"); print(table(m$category, m$sign))
cat("\nAgreement with OmniPath molecular sign:\n"); print(table(m$agreement))
cat("\nConfidence:\n"); print(table(m$confidence))

cat("\n=== Conflicts (functional sign != OmniPath molecular sign) ===\n")
print(m[agreement == "conflict",
        .(ligand, receptor, sign, db_sign, category, confidence)])

# ---- Validation: known anchors -------------------------------------------
anchor <- function(lg, rc, expected) {
  s <- m[ligand == lg & receptor == rc, sign]
  stopifnot(length(s) == 1, s == expected)
}
anchor("CD274", "PDCD1", -1)   # PD-L1/PD-1 inhibitory
anchor("CD80",  "CD28",   1)   # costimulatory
anchor("CD80",  "CTLA4", -1)   # same ligand, inhibitory receptor
anchor("FASLG", "FAS",   -1)   # death (molecular agonist override)
anchor("IFNG",  "IFNGR1", 1)   # activating
anchor("PVR",   "TIGIT", -1)   # checkpoint
anchor("PVR",   "CD226",  1)   # same ligand, activating receptor
cat("\n[OK] anchor signs validated (PD-L1/PD-1=-1, CD80/CD28=+1, CD80/CTLA4=-1,",
    "FASLG/FAS=-1, IFNG/IFNGR1=+1, PVR/TIGIT=-1, PVR/CD226=+1).\n")

cat("\nWritten to", file.path(outdir, "signedLRBase.csv"), "\n")
