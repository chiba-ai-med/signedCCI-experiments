# signed_lrbase.R
# Loader for the signed ligand-receptor base (data/signedLRBase/signedLRBase.csv).
# Returns a data.frame compatible with build_lr_matrices() (ligand, receptor, sign).

#' Load the signedLRBase.
#'
#' @param path CSV produced by experiments/signedLRBase/build_signedLRBase.R.
#' @param exclude_context if TRUE, drop context-dependent pairs.
#' @param categories optional character vector to keep only these categories.
#' @param drop_confidence optional character vector of confidence levels to drop
#'   (e.g. "review"). "expected_override" is kept by default -- those are the
#'   molecular!=functional cases the base exists to handle.
#' @return data.frame with at least ligand, receptor, sign (+ metadata columns).
load_signed_lrbase <- function(path = "data/signedLRBase/signedLRBase.csv",
                               exclude_context = FALSE,
                               categories = NULL,
                               drop_confidence = NULL) {
  stopifnot(file.exists(path))
  db <- read.csv(path, stringsAsFactors = FALSE)
  stopifnot(all(c("ligand", "receptor", "sign") %in% colnames(db)))
  stopifnot(all(db$sign %in% c(-1, 1)))
  if (exclude_context && "context_dependent" %in% colnames(db))
    db <- db[!as.logical(db$context_dependent), ]
  if (!is.null(categories))
    db <- db[db$category %in% categories, ]
  if (!is.null(drop_confidence) && "confidence" %in% colnames(db))
    db <- db[!db$confidence %in% drop_confidence, ]
  db$pair_name <- paste0(db$ligand, "_", db$receptor)
  db
}
