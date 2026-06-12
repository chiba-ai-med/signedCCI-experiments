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

## Real data (GSE72056, Tirosh melanoma TME)

```bash
/home/koki/anaconda3/envs/signedcci/bin/Rscript experiments/real/01_download.R
/home/koki/anaconda3/envs/signedcci/bin/Rscript experiments/real/run_real.R
```

`experiments/real/run_real.R` reads the 23,684 gene × 4,645 cell matrix, assigns
cell types from the GEO annotation (malignant → `Tumor`; non-malignant codes →
`Bcell/Macrophage/Endo/CAF/NK`; T cells split into `CD8`/`Treg`/`CD4conv` by
`CD8A/CD8B` and `FOXP3`), aggregates per-cell-type mean log-expression, then runs
`build_lr_matrices(expr, lr_table)` → `signed_cci()` over 12 signed LR pairs
(checkpoints/death `(-)`: TGFB1, IL10, PD-L1/PD-1, FASLG, LGALS9, PVR;
costim/activation `(+)`: CD80/CD86, IL12, IFNG, CXCL9, CCL5).

**Result:** the `(-)×(-)=(+)` motif is recovered in real data — the top
net-positive indirect path is `Tumor ─(-)→ Treg ─(-)→ CD8` (contribution 0.019),
and several double-negative routes terminate at `Tumor`
(`CAF→Treg→Tumor`, `Endo→CD8→Tumor`, …). Outputs in `results/real/`.

**Caveat:** mean-expression aggregation yields a near-fully-connected signed
graph (each hop is realised by many LR pairs at once), so signed propagation
holds mechanically but cell-type separation in net PageRank is weak
(net ≈ 0.008–0.017). Sharper specificity would need LR-specificity weighting,
edge thresholding, or per-edge dominant-LR selection — left as future work.

## Future extensions

- LR-specificity weighting / thresholding for sharper real-data signals.
- End-to-end signed random walk at the LR level (Hypergraph / Tensor PageRank),
  per the symTensor/scTensor extension spec (not in this minimal implementation).
