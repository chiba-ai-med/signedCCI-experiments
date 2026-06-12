# 01_download.R
# Download the GSE72056 melanoma TME single-cell dataset (Tirosh et al. 2016)
# into data/ (git-ignored). Idempotent: skips if the file already exists.
#
# Run: Rscript experiments/real/01_download.R

url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/",
              "GSE72056/suppl/GSE72056_melanoma_single_cell_revised_v2.txt.gz")
dir.create("data", showWarnings = FALSE)
dest <- "data/GSE72056_melanoma_single_cell_revised_v2.txt.gz"

if (file.exists(dest) && file.info(dest)$size > 1e7) {
  cat("Already present:", dest, "(",
      round(file.info(dest)$size / 1e6, 1), "MB )\n")
} else {
  cat("Downloading", url, "\n")
  options(timeout = 1200)
  download.file(url, dest, mode = "wb")
  cat("Saved", dest, "(", round(file.info(dest)$size / 1e6, 1), "MB )\n")
}
