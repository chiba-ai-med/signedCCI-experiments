# signedCCI-experiments

Signed cell-cell interaction (CCI) experiments with the updated
[`symTensor`](https://github.com/rikenbit/symTensor) and
[`scTensor`](https://github.com/rikenbit/scTensor).

The central idea: assign a sign to each ligand-receptor (LR) pair
— activation `(+)`, inhibition / cell-death induction `(-)` — build a **signed**
CCI graph, and use Signed PageRank / Signed Path Contribution to recover
**indirect** interactions via the `(-)×(-)=(+)` property
("enemy of my enemy is my friend"), e.g. `Tumor → Treg → CD8 → Tumor`.

## Environment

conda env `signedcci` (R 4.5.3). Built **without apt**; uses the **rikenbit**
development versions (not the Bioconductor release):

```bash
# 1. Bootstrap scTensor's heavy Bioconductor dependency tree via bioconda
mamba create -n signedcci -c conda-forge -c bioconda \
  bioconductor-sctensor=2.20.0 r-remotes r-biocmanager git

# 2. rikenbit/symTensor (must precede scTensor; scTensor 2.23.1 depends on it)
Rscript -e 'remotes::install_github("rikenbit/symTensor", ref="main", upgrade="never")'

# 3. rikenbit/scTensor 2.23.1 — overwrites the conda 2.20.0 bootstrap
Rscript -e 'remotes::install_github("rikenbit/scTensor", ref="master", upgrade="never", build_vignettes=FALSE)'
```

Result: `symTensor` 0.1.0 + `scTensor` 2.23.1. The signed-CCI functions
(`BuildSignedCCI`, `SignedPageRank`, `SignedPathContribution`) are pure base R;
the heavy Bioconductor deps are only needed for a consistent install.

## API and a key orientation note

- `scTensor::BuildSignedCCI(lig_expr, rec_expr, lr_sign)` → `A_pos`/`A_neg`/`X_pos`/`X_neg`/`edge_table`,
  with `A[sender(row), receiver(col)]` (= `A[from, to]`).
- `symTensor::SignedPageRank(A_pos, A_neg, ...)` and
  `symTensor::SignedPathContribution(A_pos, A_neg, max_hop, ...)` expect
  `A[to(row), from(col)]`.

⚠️ **These are transposed relative to each other** — transpose `t(A_pos)`, `t(A_neg)`
before passing `BuildSignedCCI` output to the walkers.
`R/signed_lr_utils.R::signed_cci()` handles this.

## Layout

```
R/signed_lr_utils.R                     # build_lr_matrices / signed_cci / annotate_paths
experiments/synthetic/run_synthetic.R   # toy immune×tumor GATE test
results/synthetic/                      # A_pos/A_neg, pagerank, signed_paths, edge_table CSVs
```

## Synthetic gate test

```bash
/home/koki/anaconda3/envs/signedcci/bin/Rscript experiments/synthetic/run_synthetic.R
```

5 cell types (`Tumor, Treg, CD8, NK, DC`) with a biologically-labelled signed LR
table. All gate assertions pass:

- 2-hop `Treg ─(-)→ CD8 ─(-)→ Tumor` → `net_sign = +1`
- 3-hop `Tumor ─(+)→ Treg ─(-)→ CD8 ─(-)→ Tumor` → `net_sign = +1`
  (self-reinforcing anti-tumor loop)
- single `(-)` hop stays `-`, single `(+)` hop stays `+`
- `A_pos`/`A_neg` match the designed directed signed edges; Signed PageRank converges

This confirms `(-)×(-)=(+)` indirect propagation on a known network.

## Next

Real data (e.g. GSE72056, Tirosh melanoma TME): assign cell types
(`Tumor/Treg/CD8/NK/...`), aggregate per-cell-type expression, run
`build_lr_matrices(expr, lr_table)` → `signed_cci()`, and check whether
`Treg→CD8→Tumor`-type indirect positive paths appear in real data.
