# DM-CART: Dirichlet-Multinomial Classification and Regression Trees

Recursive partitioning driven by the **marginal Dirichlet-Multinomial
likelihood**, applied to U.S. natality data to model the full discretized
birth-weight distribution rather than a binary low-birth-weight indicator.

**Manuscript:** *Investigating Determinants of Birth Weight Using Bayesian
Tree-Based Nonparametric Modeling* — Adam M. Kurth (Brown University) and
P. Richard Hahn (Arizona State University).

This work extends the 2020–2021 analysis in
[Kurth (2025), MS thesis, Arizona State University](https://hdl.handle.net/2286/R.2.N.201810).
Three things changed: the cohort is now **2023–2024**, maternal race is carried
at its **full seven-level resolution** instead of a Black/other binary, and
marital status keeps an explicit **unknown** level instead of dropping the
11.1% of records where it is not reported. The class space grows from 128
profiles to **672**, which is affordable only because the Dirichlet prior
stabilizes the sparse cells that result.

---

## What the method does

Within a maternal-infant class $i$, category counts are multinomial given
$\boldsymbol{\theta}_i$, with a conjugate Dirichlet prior on
$\boldsymbol{\theta}_i$. Integrating $\boldsymbol{\theta}_i$ out gives a
closed-form marginal likelihood, and a candidate split of node $t$ is scored by

$$\Delta\ell = \ell(t_L) + \ell(t_R) - \ell(t)$$

Replacing an impurity measure with a likelihood is what carries the prior into
the splitting decision: a split that isolates a handful of births is not
rewarded for the extreme category proportions it produces, because those
proportions are shrunk toward $\boldsymbol{\alpha}$ before being scored.

Two consequences matter in practice:

- **Sparse cells are the design point, not a tolerated edge case.** A
  categorical predictor can be disaggregated to a resolution plain multinomial
  CART cannot sustain.
- **Categorical splits are searched exhaustively.** A $k$-level unordered
  predictor admits $2^{k-1}-1$ bipartitions; the conventional search evaluates
  only the $k-1$ that are contiguous in the declared level order. For race at
  seven levels that is 6 of 63, and *which* 6 depends on nothing but the order
  the levels were declared in. All 63 are evaluated here.

Uncertainty comes from a **two-tier parametric bootstrap** over $B = 10{,}000$
trees. Tier 1 redraws how the $T$ births distribute across the $I$ classes;
tier 2 redraws, within each class, how its births distribute across categories
using the posterior mean $\hat{\pi}_{i,k} = (n_{i,k}+\alpha_k)/(N_i+\alpha_0)$.

---

## Quick start

```bash
git clone https://github.com/adamkurth/dm-cart.git
cd dm-cart
```

```r
install.packages(c("rpart", "rpart.plot", "data.table", "ggplot2", "dplyr",
                   "tidyr", "patchwork", "reshape2", "xtable", "rprojroot"))
```

Place the two raw natality CSVs in `data/` (see [Data](#data) below), then:

```bash
./run-all.sh
```

That reproduces every number, figure, and table in the manuscript from the raw
files, start to finish, in roughly **10 minutes** on a current laptop — most of
it stage 2 reading the ~2 GB natality CSVs and stage 3 fitting the 10,000-tree
ensemble.

To smoke-test the wiring in well under a minute:

```bash
./run-all.sh --b 200
```

---

## `run-all.sh`

Four stages, each consuming the previous one's output. Finished stages are
skipped unless `--force` is given, so an interrupted run resumes.

| Stage | Script | Reads | Writes |
|---|---|---|---|
| 1 | `quantile.R` | prior-year CSV | decile cutpoints + Dirichlet prior |
| 2 | `data.R` | target-year CSV | class design matrix + counts |
| 3 | `dm-cart.R` | counts | fitted trees + bootstrap ensemble |
| 4 | `visual.R` | saved `.RData` | figures and tables in `results/<cohort>/article/` |
| 5 | `verify-numbers.R` | saved `.RData` | checks every number quoted in the manuscript |

| Flag | Effect |
|---|---|
| *(none)* | full pipeline, `B = 10000`, skipping finished stages |
| `--b N` | smaller ensemble — use for smoke runs |
| `--figures-only` | redraw article outputs from saved results; touches nothing else |
| `--force` | redo every stage, ignoring saved outputs |
| `--cohort YYYY-YYYY` | a different cohort (default `2023-2024`) |
| `-h`, `--help` | usage |

Only stage 3 is expensive. Stages 4 and 5 read nothing but saved `.RData` and
are safe to re-run at any time.

**Stage 5 is the one that catches drift.** The manuscript quotes about forty
numbers that the pipeline produces rather than the author typing. A re-fit moves
interval endpoints in the third decimal and selection percentages by a few
tenths of a point, and a stale number in the text is invisible unless something
checks it. `verify-numbers.R` compares each claim against the saved output and
fails loudly:

```
  lbw high stability %                       98.9             98.7             MISMATCH
```

Run it on its own at any time with `Rscript verify-numbers.R`.

**A guard worth knowing about.** Stage 3 refuses to overwrite a saved ensemble
with a *smaller* one:

```
REFUSING: saved results use B = 10000, larger than the requested B = 200.
Re-run with --force to overwrite, or --figures-only to keep them.
```

This exists because it is otherwise very easy to destroy a multi-hour
publication run while testing a figure change.

Scripts can also be run individually, in stage order:

```bash
Rscript quantile.R 2023-2024
Rscript data.R     2023-2024
DMCART_B=10000 Rscript dm-cart.R 2023-2024
Rscript visual.R   2023-2024        # full model
Rscript visual.R   2023-2024 lbw    # LBW-only model
Rscript verify-numbers.R
```

The bootstrap is seeded (`set.seed()` near the top of `dm-cart.R`), so a re-run
on the same input reproduces the published numbers exactly. Change the seed and
every interval endpoint moves; stage 5 will tell you which sentences need
updating.

---

## Files

| File | Role |
|---|---|
| `util.R` | Path resolution for every cohort and model arm, directory creation, and the tree-rendering helpers (`dm.tree.coords`, `draw.dm.tree`, `plot.dm.tree`). Sourced by everything else. |
| `quantile.R` | Prior-year deciles of the sub-2.5 kg distribution → cutpoints and the Dirichlet hyperparameter vector. |
| `data.R` | Raw natality CSV → the $I \times 7$ class design matrix and the $I \times K$ counts matrix. Handles the NCHS "not stated" sentinels. |
| `dm-cart.R` | The DM split rule, the fit, the bootstrap ensemble, variable importance, risk ranking, and the extreme-profile analysis. |
| `visual.R` | Draws every article figure and writes every article table fragment, reading only saved results. |
| `verify-numbers.R` | Checks every number quoted in the manuscript against the saved output. Run before sending the paper anywhere. |
| `run-all.sh` | Driver for the four stages. |
| `tex/article/` | The manuscript. |


### The custom `rpart` method

```r
dm.method <- list(init = myinit, eval = myeval, split = mysplit, method = "dm")
```

- `myinit(y, offset, parms, wt)` — sets up the DM tree; `parms` carries
  $\boldsymbol{\alpha}$.
- `myeval(y, wt, parms)` — the DM marginal log-likelihood at a node.
- `mysplit(y, wt, x, parms, continuous)` — $\Delta\ell$ over candidate splits;
  for categorical predictors, over all $2^{k-1}-1$ bipartitions.

**There is no `pred` element, deliberately.** `rpart`'s prediction path is not
usable with a user-written method: `predict.rpart()` returns the $n \times 1$
`yval` rather than the count vector, and calling `rpart:::pred.rpart()`
directly **segfaults**. Leaf probabilities are read straight off the fitted
object instead:

```r
dm.leaf.probs <- function(fit, alphavec, cat.names = colnames(Y.df)) {
  counts <- fit$frame$yval2[fit$where, , drop = FALSE]
  t(apply(counts, 1, function(z) (z + alphavec) / sum(z + alphavec)))
}
```

Note `fit$where` is a **row index into `fit$frame`**, not a node id. Several
other `rpart` internals are similarly counter-intuitive and are documented at
the point of use in `dm-cart.R` — in particular, `splits[, "index"]` for a
categorical split is a row number into `fit$csplit` (1 = left, 3 = right,
2 = absent), not a cutpoint.

---

## Data

The analysis uses the **U.S. Natality public-use micro-data** published by the
National Center for Health Statistics, available without restriction from the
[NCHS Vital Statistics Online Data Portal](https://www.cdc.gov/nchs/data_access/vitalstatsonline.htm)
(also mirrored by [NBER](https://www.nber.org/research/data/vital-statistics-natality-birth-data)).

The CSVs are ~2 GB each and are **not** in this repository. Place them in
`data/` named exactly:

```
data/natality2023us-original.csv     # prior year -- cutpoints and prior only
data/natality2024us-original.csv     # target year -- the counts being modeled
```

`util.R` checks for both and lists what is missing before anything runs.

Estimating the cutpoints and prior on 2023 while modeling counts on 2024 keeps
the two disjoint — the prior carries no information about the sample it is
applied to.

The files contain no direct identifiers, and NCHS suppresses every geographic
identifier before release. No ethics approval was required for this secondary
analysis.

### Fields read

Only these columns are read (`data.table::fread(select = ...)`, which is what
takes stage 2 from minutes to about 17 seconds):

| Field | Role |
|---|---|
| `dbwt` | birth weight in grams — the response, before discretization |
| `sex` | infant sex |
| `mager` | maternal age |
| `meduc` | maternal education |
| `precare5` | month prenatal care began |
| `cig_0` | cigarettes per day before pregnancy |
| `dmar` | marital status |
| `mracehisp` | maternal race and Hispanic origin |

**Sentinel codes are converted to `NA` before recoding, not after.** NCHS uses
in-range values for "not stated": `dbwt = 9999`, `meduc = 9`, `precare5 = 5`,
and `cig_0 >= 99` (98 is a legitimate maximum count, so the test is `>=`, not
`==`). Recoding first would silently fold "not stated" into a substantive
level — `cig_0 = 99` would become "smoker".

### Predictors

Seven predictors give $2^5 \cdot 3 \cdot 7 = 672$ classes.

| Predictor | Source | Levels | Encoding |
|---|---|---|---|
| `sex` | `sex` | 2 | 1 = male, 0 = female |
| `mager` | `mager` | 2 | 1 = over 33 years, 0 = 33 or younger |
| `meduc` | `meduc` | 2 | 1 = high school completed or beyond (`meduc >= 3`), 0 = otherwise |
| `precare1st` | `precare5` | 2 | 1 = care began in the first trimester, 0 = otherwise |
| `cig_0` | `cig_0` | 2 | 1 = any pre-pregnancy smoking, 0 = none |
| `dmar` | `dmar` | 3 | married; unmarried; **unknown** |
| `race` | `mracehisp` | 7 | White, Black, Hispanic, Asian, AIAN, NHOPI, multiracial |

Two of these differ from the thesis and are the point of the revision:

- **`race` is not binary.** Collapsing to Black/other, or to an "other"
  bucket, merges groups with different risk profiles. The fitted trees contain
  groupings such as {White, Hispanic, AIAN} against {Asian, NHOPI, multiracial}
  that a binary contrast cannot express and a contiguous search could not reach.
- **`dmar` keeps `unknown` as a level.** It is missing for 11.1% of records and
  the missingness is geographically clustered, but the public-use file
  suppresses every geographic identifier — so the mechanism is unobservable and
  no conditional imputation is testable. Carrying `unknown` explicitly assumes
  nothing and lets the fitted model report where the group belongs.

`meduc >= 3` is "high school completed **or beyond**". The thesis used
`meduc == 3`, which is high school and *nothing further* — it put college
graduates in the low-education group.

### Response

Birth weight is discretized into $K$ ordered categories: ten deciles of the
sub-2.5 kg distribution, plus a single normal-birth-weight category. The full
model keeps all births ($K = 11$); the LBW-only model restricts to births below
the threshold ($K = 10$), where the deciles are conditional on being LBW.

Cutpoints and prior mass below are the 2023 values actually used for the 2024
fit (`data/rebin/2023-2024/`). Deciles are not exactly equal because birth
weight is recorded in whole grams.

| Category | Range (g) | Prior mass, full model | Prior mass, LBW-only |
|---|---|---|---|
| Q1 | 227–1191 | 0.87% | 10% |
| Q2 | 1191–1660 | 0.88% | 10% |
| Q3 | 1660–1910 | 0.86% | 10% |
| Q4 | 1910–2070 | 0.89% | 10% |
| Q5 | 2070–2185 | 0.85% | 10% |
| Q6 | 2185–2272 | 0.86% | 10% |
| Q7 | 2272–2350 | 0.93% | 10% |
| Q8 | 2350–2410 | 0.96% | 10% |
| Q9 | 2410–2460 | 0.76% | 10% |
| Q10 | 2460–2500 | 0.83% | 10% |
| NBW | > 2500 | 91.31% | — |

The Q1 interval is closed on the left; every other interval is left-open,
right-closed. `data.R` reconciles the binned total against the record count and
stops if they disagree.

---

## Outputs

```
results/2023-2024/
├── article/          figures and table fragments the manuscript \input{}s
├── fullmodel/        full-model ensemble, risk tables, diagnostic plots
└── lbwmodel/         LBW-only equivalents
```

The four manuscript figures are drawn at **7.0 inches wide** — the one-column
measure of the journal's LaTeX class — so `\includegraphics[width=\textwidth]`
reproduces them at scale 1.0 and the point sizes set in `visual.R` are the ones
that reach the page. `ARTICLE.FIG.WIDTH` at the top of `visual.R` is the single
place to change this if the target layout changes.

### Headline results (2024 cohort, $B = 10{,}000$)

| | |
|---|---|
| Records | 3,638,436 → 3,476,907 retained (95.6%) |
| Classes | 672, of which 18 hold no births |
| Fitted trees | 23 terminal nodes (full), 8 (LBW-only) |
| Root split | maternal race, in 100% of replicates in **both** arms |
| Highest risk | class 181, $\Pr(\mathrm{LBW}) = 0.284\ [0.239, 0.330]$, $N = 1{,}710$ |
| Lowest risk | class 40, $\Pr(\mathrm{LBW}) = 0.054\ [0.051, 0.056]$, $N = 346{,}814$ |
| Spread | 5.3-fold between extremes |

The headline finding is that the determinants of **crossing** the 2.5 kg
threshold and the determinants of **severity given** low birth weight are
different sets of variables. Marital status, smoking, education and prenatal
care are selected in essentially every full-model replicate and in almost none
of the LBW-only replicates (3.9%, 1.9%, 1.0%, 6.9%). Race is the exception: it
survives the restriction and is the root split in both arms. A binary LBW
outcome cannot express this distinction, because it collapses the second
question entirely.

---

## Reproducibility notes

- **Extreme profiles are derived, not assumed.** The clinically intuitive
  "everything adverse" profile does not top the ranking in either sex — its
  female variant ranks 5th and its male variant 13th, and both fall below the
  support threshold. Fixing comparison groups in advance would have reported a
  rare profile as the high-risk group.
- **Ties are resolved at the resolution of the estimate.** Every class in a
  terminal node shares a posterior, so extremes are not unique. Classes within
  5% of the interval width of the extremum are treated as tied, and the tied
  class with the most observed births is reported.
- **The support threshold is a judgment, not a derived quantity.** Classes
  with fewer than 1,000 observed births are plotted but excluded from the
  ranking, because such a class contributes almost nothing to its terminal node
  and its estimate is dominated by its neighbors and the prior. Figure 3 shows
  the excluded classes rather than hiding them.
- **No fabricated fallbacks.** Every figure and table is drawn from saved
  bootstrap output. If the output is missing, the script stops rather than
  substituting placeholder values.

---

## Citation

```bibtex
@article{kurth-dmcart,
  author  = {Kurth, Adam M. and Hahn, P. Richard},
  title   = {Investigating Determinants of Birth Weight Using Bayesian
             Tree-Based Nonparametric Modeling},
  year    = {2026},
  note    = {Manuscript in preparation}
}
```

The archived code release that produced the published results:

```bibtex
@software{kurth-dmcart-code,
  author    = {Kurth, Adam M.},
  title     = {adamkurth/dm-cart: v1.0.0},
  year      = {2026},
  publisher = {Zenodo},
  version   = {v1.0.0},
  doi       = {10.5281/zenodo.22089989},
  url       = {https://doi.org/10.5281/zenodo.22089989}
}
```

The earlier thesis:

```bibtex
@mastersthesis{kurth-thesis2025,
  author  = {Kurth, Adam},
  title   = {Investigating Determinants of Birth Weight Using Bayesian
             Tree-Based Nonparametric Modeling},
  school  = {Arizona State University},
  address = {Tempe, AZ},
  type    = {master's thesis},
  year    = {2025},
  url     = {https://hdl.handle.net/2286/R.2.N.201810}
}
```

## Contact

Adam M. Kurth — adam_kurth@brown.edu — [github.com/adamkurth](https://github.com/adamkurth)
