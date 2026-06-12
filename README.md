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

### Contributing LR detection

```bash
/home/koki/anaconda3/envs/signedcci/bin/Rscript experiments/real/lr_contribution.R
```

`BuildSignedCCI`'s `edge_table` keeps the per-LR score `X_ijk = L_ik·R_jk`, so each
signed edge can be decomposed back into its LR contributions and ranked
(`results/real/lr_contribution_by_edge.csv`, `dominant_lr_per_edge.csv`,
`lr_global_importance.csv`). Examples:

- `CD8 ─(-)→ Tumor` cytotoxicity is driven by **FASLG/FAS** (40%) and LGALS9 (30%);
- `Treg ─(-)→ CD8` suppression by **IL10** (41%), LGALS9/TIM3 (38%), PD-L1/PD-1 (12%);
- globally, `CCL5/CCR5` carries the most positive weight and
  `LGALS9/HAVCR2 (TIM3)`, `IL10/IL10RA` the most negative.

This delivers both halves of the goal: detecting signed CCI **and** the LR pairs
responsible for each signed interaction.

## signedLRBase — a rule-based signed LR database

The 12-pair hand-list above was replaced by a curated, reproducible **signedLRBase**.
Existing resources (SIGNOR, OmniPath, CellPhoneDB v5) carry a *molecular signalling*
sign (does the ligand activate or inhibit the receptor's signalling), which is **not**
the *functional outcome* sign this project needs: `FASLG→FAS` is a molecular agonist
(+) but the cellular outcome is death (−). So we assign the **functional sign**
(+1 = activation/survival/proliferation, −1 = immunosuppression/exhaustion/death) by
curated rules, then cross-check against OmniPath.

```bash
/home/koki/anaconda3/envs/signedcci/bin/Rscript experiments/signedLRBase/build_signedLRBase.R
```

- `experiments/signedLRBase/signed_lr_curated.csv` — 91 immune/TME LR pairs with a
  functional sign + category (checkpoint_inhibitory / nk_inhibitory / death /
  suppressive_cytokine = −; costimulation / activating_cytokine / chemokine /
  growth_survival / nk_activating = +) and `context_dependent` flags.
- `experiments/signedLRBase/build_signedLRBase.R` — pulls OmniPath's molecular sign
  (`is_stimulation`/`is_inhibition`) and records agreement → `data/signedLRBase/signedLRBase.csv`.
- `R/signed_lrbase.R::load_signed_lrbase()` — loader (compatible with `build_lr_matrices`).

**Validation of the molecular≠functional gap:** of 91 pairs, 50 match OmniPath's
molecular sign, **26 are expected overrides** (all inhibitory pathways that work *by*
molecular agonism of a suppressive receptor — PD-L1/PD-1, CTLA4, TIGIT, HLA/NKG2A,
TRAIL, TGFβ, IL10, …), 15 are DB-ambiguous, and **0 are genuine review conflicts**.
i.e. using OmniPath's molecular sign directly would mis-sign ~29% of the negative pairs
— which is exactly why signedLRBase exists.

### Real data with signedLRBase

```bash
/home/koki/anaconda3/envs/signedcci/bin/Rscript experiments/real/run_real_signedLRBase.R
```

Re-runs GSE72056 with the 91-pair signedLRBase (88 usable; shared cell-typing helper
`R/gse72056.R`). The `(-)×(-)=(+)` motif persists (`Tumor→Treg→CD8`, `Treg→CD8→Tumor`
as net-positive double-negative paths); the richer death/NK-inhibitory coverage makes
NK a strong negative hub. Outputs in `results/real_signedLRBase/`.

## Context-dependence of the sign (scope & limits)

The functional sign of an LR pair is **not a global constant** — it is a function of
receiver cell type, cell state, co-signals, tissue, dose and time. No static database
can enumerate all contexts, so signedLRBase is a **prior** (default sign + `category` +
`context_dependent` + `confidence`), not ground truth, scoped to the immune×tumor TME.

Two experiments probe context-dependence on GSE72056:

**1. Marker check (rule-based, does NOT generalise).**
`experiments/real/sign_refinement_markers.R` scores each receiver's activation vs
exhaustion/M1-vs-M2 state with hand-picked panels. In this melanoma data CD8 is
**exhausted** (exhaustion z 1.57 > activation 1.23 — the classic Tirosh finding), so
nominally-activating signals onto CD8 (IFNG, TNF, IL1B) read as suppressed: 10/19
context pairs disagree with the prior. Caveat: this measures the receiver's *overall
state*, not per-ligand causation (all CD8-receiver pairs collapse to one score) — and
the hand panels must be re-curated for every new dataset, so it does not generalise.

**2. Pre-built reference (PROGENy, generalises to arbitrary data).**
`experiments/real/progeny_demo.R` applies PROGENy's pre-trained signed pathway model
— no per-dataset curation — and the signed activity localises sensibly:
TGFb (−) in CAF, EGFR (+) in Tumor, VEGF (−) & JAK-STAT/NFkB (+) in Macrophage,
Trail (−) in B/NK. The three steps: (1) *which genes* and (2) *measure activity* are
dataset-agnostic via the pre-built model; only (3) *activity → +/- sign* keeps a thin,
irreducible prior (`progeny_pathway_sign_prior.csv`) — because "pathway active" ≠ "good
for the cell" (same logic as FASLG being a molecular agonist but functionally death).

**Takeaway:** signs move in real data, but a single "true" sign is ill-posed; treat the
sign as prior + dataset-local evidence + uncertainty. For arbitrary data the
generalisable recipe is *pre-built reference (CytoSig / DoRothEA / PROGENy) + a thin
sign-mapping prior*. NicheNet cannot supply the sign — its ligand→target potential is
*unsigned* (activation and inhibition contribute additively).

## Future extensions

- LR-specificity weighting / thresholding for sharper real-data signals.
- End-to-end signed random walk at the LR level (Hypergraph / Tensor PageRank),
  per the symTensor/scTensor extension spec (not in this minimal implementation).
