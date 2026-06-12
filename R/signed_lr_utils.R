# signed_lr_utils.R
# Helpers for signed cell-cell interaction (CCI) analysis with
# scTensor::BuildSignedCCI + symTensor::SignedPageRank / SignedPathContribution.
#
# Orientation note (verified against the rikenbit sources):
#   * BuildSignedCCI returns A[i, j] with i = sender (row), j = receiver (col),
#     i.e. A[from, to].
#   * SignedPageRank / SignedPathContribution normalise with colSums() and use
#     P[j, current] = A[j, current] / out_deg[current], i.e. they expect
#     A[to (row), from (col)].
#   => The BuildSignedCCI output must be TRANSPOSED before being handed to the
#      symTensor walkers. signed_cci() does this for you.

suppressMessages({
  library(scTensor)   # BuildSignedCCI
  library(symTensor)  # SignedPageRank, SignedPathContribution
})

#' Build (cell type x LR pair) ligand/receptor matrices and a sign vector.
#'
#' Two modes:
#'   * Synthetic (expr = NULL): each row of `lr_table` becomes one LR-pair
#'     column with the ligand expressed ONLY in its `sender` cell type and the
#'     receptor expressed ONLY in its `receiver` cell type. This encodes the
#'     intended directed signed edge exactly (no cross-talk), which is what the
#'     toy gate test needs.
#'   * Real (expr given): `expr` is a (gene x cell type) expression matrix.
#'     For each LR pair, lig_expr[, k] = expr[ligand, ] and
#'     rec_expr[, k] = expr[receptor, ]. BuildSignedCCI then forms every
#'     sender/receiver combination from the real expression.
#'
#' @param lr_table data.frame with columns: ligand, receptor, sign, and
#'   (synthetic mode) sender, receiver. Optional `pair_name`.
#' @param cell_names character vector of cell-type names (defines matrix rows).
#' @param expr optional (gene x cell type) numeric matrix; colnames must be a
#'   superset of `cell_names`, rownames must include all ligand/receptor genes.
#' @return list(lig_expr, rec_expr, lr_sign, pair_names) with the expression
#'   matrices laid out as (cell type x LR pair).
build_lr_matrices <- function(lr_table, cell_names, expr = NULL) {
  stopifnot(all(c("ligand", "receptor", "sign") %in% colnames(lr_table)))
  K <- nrow(lr_table)
  n <- length(cell_names)

  pair_names <- if (!is.null(lr_table$pair_name)) {
    as.character(lr_table$pair_name)
  } else {
    paste0(lr_table$ligand, "_", lr_table$receptor)
  }

  lig_expr <- matrix(0, nrow = n, ncol = K,
                     dimnames = list(cell_names, pair_names))
  rec_expr <- matrix(0, nrow = n, ncol = K,
                     dimnames = list(cell_names, pair_names))

  if (is.null(expr)) {
    # Synthetic: isolate each directed edge.
    stopifnot(all(c("sender", "receiver") %in% colnames(lr_table)))
    lig_level <- if (!is.null(lr_table$lig_level)) lr_table$lig_level else rep(1, K)
    rec_level <- if (!is.null(lr_table$rec_level)) lr_table$rec_level else rep(1, K)
    for (k in seq_len(K)) {
      s <- as.character(lr_table$sender[k])
      r <- as.character(lr_table$receiver[k])
      stopifnot(s %in% cell_names, r %in% cell_names)
      lig_expr[s, k] <- lig_level[k]
      rec_expr[r, k] <- rec_level[k]
    }
  } else {
    # Real: pull per-cell-type expression for each LR gene.
    expr <- expr[, cell_names, drop = FALSE]
    for (k in seq_len(K)) {
      lg <- as.character(lr_table$ligand[k])
      rg <- as.character(lr_table$receptor[k])
      if (!lg %in% rownames(expr) || !rg %in% rownames(expr)) {
        next  # gene missing -> leave column at 0 (no edge); caller can filter
      }
      lig_expr[, k] <- expr[lg, ]
      rec_expr[, k] <- expr[rg, ]
    }
  }

  list(lig_expr = lig_expr, rec_expr = rec_expr,
       lr_sign = as.integer(lr_table$sign), pair_names = pair_names)
}

#' Run the full signed-CCI pipeline: BuildSignedCCI -> (transpose) ->
#' SignedPageRank + SignedPathContribution.
#'
#' @param lig_expr,rec_expr (cell type x LR pair) matrices.
#' @param lr_sign integer vector (+1/-1), length ncol(lig_expr).
#' @param cell_names character vector naming the cell types (matrix rows).
#' @param max_hop max path length for SignedPathContribution.
#' @param damping PageRank damping factor.
#' @return list with cci (BuildSignedCCI output, sender x receiver),
#'   A_pos_t/A_neg_t (transposed = to x from), pagerank, paths (with named
#'   source/target), edge_table (with cell names), cell_names.
signed_cci <- function(lig_expr, rec_expr, lr_sign, cell_names,
                       max_hop = 3L, damping = 0.85) {
  cci <- BuildSignedCCI(lig_expr, rec_expr, lr_sign)

  # Name the sender x receiver matrices.
  dimnames(cci$A_pos) <- list(cell_names, cell_names)
  dimnames(cci$A_neg) <- list(cell_names, cell_names)

  # Transpose for the symTensor walkers (they expect A[to, from]).
  A_pos_t <- t(cci$A_pos)
  A_neg_t <- t(cci$A_neg)

  pagerank <- SignedPageRank(A_pos_t, A_neg_t, damping = damping)
  names(pagerank$positive) <- cell_names
  names(pagerank$negative) <- cell_names
  names(pagerank$net)      <- cell_names
  names(pagerank$total)    <- cell_names

  paths <- SignedPathContribution(A_pos_t, A_neg_t, max_hop = max_hop)
  if (nrow(paths) > 0) {
    paths$source_name <- cell_names[paths$source]
    paths$target_name <- cell_names[paths$target]
    paths$path_name <- vapply(strsplit(paths$path, "->"), function(idx) {
      paste(cell_names[as.integer(idx)], collapse = "->")
    }, character(1))
  }

  # Edge table with names.
  et <- cci$edge_table
  if (nrow(et) > 0) {
    et$sender_name <- cell_names[et$sender]
    et$receiver_name <- cell_names[et$receiver]
  }

  list(cci = cci, A_pos_t = A_pos_t, A_neg_t = A_neg_t,
       pagerank = pagerank, paths = paths, edge_table = et,
       cell_names = cell_names)
}

#' Annotate paths with the LR pairs that realise each signed hop.
#'
#' For every consecutive (from -> to) hop in a path, look up the LR pairs in
#' `edge_table` whose (sender, receiver) match the communication direction and
#' whose sign matches the hop's sign. Returns the input data.frame with an
#' added `lr_annotation` column (one entry per path, hops joined by " | ").
annotate_paths <- function(paths, edge_table, cell_names) {
  if (nrow(paths) == 0) return(paths)
  pair_lab <- if (!is.null(edge_table$pair_name)) {
    edge_table$pair_name
  } else if (!is.null(edge_table$lr_pair)) {
    paste0("LR", edge_table$lr_pair)
  } else {
    rep("?", nrow(edge_table))
  }

  paths$lr_annotation <- vapply(seq_len(nrow(paths)), function(r) {
    idx <- as.integer(strsplit(paths$path[r], "->")[[1]])      # node sequence
    sgn <- strsplit(paths$signs[r], ",")[[1]]                  # per-hop signs
    hops <- vapply(seq_along(sgn), function(h) {
      from <- idx[h]; to <- idx[h + 1L]
      want_sign <- if (sgn[h] == "+") 1L else -1L
      sel <- which(edge_table$sender == from &
                   edge_table$receiver == to &
                   edge_table$sign == want_sign)
      lrs <- if (length(sel)) paste(unique(pair_lab[sel]), collapse = "/") else "?"
      sprintf("%s-%s>%s[%s]", cell_names[from], sgn[h], cell_names[to], lrs)
    }, character(1))
    paste(hops, collapse = " | ")
  }, character(1))

  paths
}
