# DM-CART Revision Plan — 2023/2024 cohort, expanded race

Status: proposal. Nothing in this document has been applied to the codebase yet.

Every claim below was verified by running code against this repo and against the raw
NCHS files in `data/`, not by reading alone. Verification commands are noted inline.
The professor's audit is folded in; where I disagree with it, I say so and why.

---

## 0. Executive summary

The audit's read of the **statistical core is correct**: `log.dm.likelihood()`,
`mysplit()`'s ΔL rule, the two-tier bootstrap, and the informed-prior construction in
`quantile.R` all match Ch. 2–3 as written. That part needs no repair.

But the audit stopped one layer short in two places, and both are fatal to the numbers
currently on disk:

1. **`mypred()` misroutes rows for *every* factor split, not just race.** The audit
   flagged the categorical *search*; the categorical *prediction* path is separately and
   more severely broken. Measured: 16 of 20 rows land in the wrong leaf, max probability
   error 0.61. Every element of `Yhat` — and therefore every high/low-risk number,
   Table 3.2, Table 3.3 — is computed from this.
2. **The prior year and the target year are the same year.** `get_paths()` computes
   `target.year` and then never returns it. The no-double-dipping design — the thing
   Section 2.x explicitly claims — is not in effect in any run that has happened.

There is also a data-coding tier the audit could not have caught without the raw files:
four NCHS sentinel "unknown" codes flow into the model as real values, and ~20% of
records are being dropped silently and non-randomly.

On race: **I recommend against the audit's item #1** (revert to binary). Details in §3.

Fix order is Tier 0 → 3 below. Tier 0 must land before any number is quoted anywhere.

---

## Tier 0 — invalidates every current result

### 0.1 `mypred()` sends rows to the wrong leaves  🔴

`dm-cart.R:154–239`. The manual tree walk does:

```r
split_val <- fit_copy$splits[split_rows[1], "index"]
go_left   <- as.numeric(cur_val) <= split_val
```

Two independent errors:

* For a **categorical** split, `splits[, "index"]` is *not* a cutpoint — it is the **row
  number of `fit$csplit`**, the table that actually encodes which levels go left
  (`1` = left, `3` = right, `2` = absent). Comparing a factor's integer code against a
  `csplit` row index is meaningless. Every predictor in this model is a factor, so this
  is wrong at every split, including the binary ones — it is not a race-only problem.
* `split_rows[1]` takes the *first* `splits` row bearing that variable name, so a
  variable reused at several nodes is always routed by its first split.

Measured on a controlled tree (probe script, 4 leaves / 20 classes):

```
rows where the current traversal != correct leaf: 16 of 20
max absolute probability error:                   0.612
row 1 correct: 0.2314 0.2675 0.2548 0.2463
row 1 current: 0.7400 0.0912 0.0874 0.0814
```

**Fix.** No traversal is needed at all. `X.matrix` *is* the complete enumeration of
classes, and each bootstrap tree is fit on exactly those rows, so rpart has already done
the routing: `fit$where` gives the frame **row index** for every training row.

```r
dm.leaf.probs <- function(fit, alphavec) {
  cnt <- fit$frame$yval2[fit$where, , drop = FALSE]   # leaf count vectors, in row order
  p   <- t(apply(cnt, 1, function(z) (z + alphavec) / sum(z + alphavec)))
  colnames(p) <- names(alphavec)
  p
}
```

Two traps to avoid:

* `predict.rpart()` on a custom method returns an `n x 1` yval, not the count matrix — useless here.
* `rpart:::pred.rpart()` **segfaults** on custom-method fits. Verified — R aborts. Do not use it.

If you ever need to score genuinely new rows, the only correct route is `csplit` +
`splits[, "index"]` per *node* (not per variable). For this design you never need to.

### 0.2 Prior year == target year → double-dipping  🔴

`util.R:38–58` computes `target.year` inside the branch and then returns `year = year`,
the *requested* year. `cohort.dir` likewise is computed and never returned.

```
$ Rscript -e 'source("util.R"); p <- get_paths(2023); cat(p$year, p$prior.year, is.null(p$cohort))'
2023 2023 TRUE
```

So `TARGET.YEAR <- init.paths$year` in all three scripts resolves to **2023**, the same
year as `PRIOR.YEAR`. The quantile cutpoints and the Dirichlet prior are built from the
same file the counts are built from. The on-disk artifacts confirm this is what actually
ran — `birthweight_data_rebin_2023.RData`, `bootstrap_tree_results_2023.RData`; there is
no 2024 output anywhere. The 2024 CSV has never been touched.

As a side effect `sprintf("... %s ...", init.paths$cohort)` returns `character(0)`, so the
cohort banner in all three scripts prints nothing — which is why this was never noticed.

**Fix.** Return the resolved cohort fields and consume them:

```r
list(
  ...,
  requested.year = requested.year,
  year           = target.year,   # was: year
  target.year    = target.year,
  prior.year     = prior.year,
  cohort         = cohort.dir,
  ...
)
```

and in each script use `TARGET.YEAR <- init.paths$target.year`. Then add a guard that
refuses to run when they collide:

```r
stopifnot(PRIOR.YEAR < TARGET.YEAR)
```

Note this renames output files (`*_2023.RData` → `*_2024.RData` for target-year
artifacts). Existing `results/2023-2024/` contents should be deleted, not merged — they
are double-dipped.

### 0.3 NCHS sentinel "unknown" codes enter the model as data  🔴

`na.omit()` catches `NA`, but NCHS encodes unknown as in-range integers. Counts below are
from a 2,000,000-record sample of `natality2024us-original.csv`:

| field | unknown code | count / 2M | currently becomes | should be |
|---|---|---|---|---|
| `dbwt` | `9999` | 1,648 | **counted as Normal (>2500 g)** | dropped |
| `cig_0` | `99` | 11,202 | smoker = 1 | dropped |
| `meduc` | `9` | 52,163 | "not HS" = 0 | dropped |
| `precare5` | `5` | 35,060 | "not adequate" = 0 | dropped |
| `mracehisp` | `8` | 22,804 | — already handled correctly | — |

`dbwt == 9999` is the worst of these: it inflates the NBW category, which is the
denominator for the entire LBW story.

**Fix.** Convert sentinels to `NA` *before* any recode, and print a missingness table:

```r
dat$dbwt[dat$dbwt == 9999]         <- NA
dat$cig_0[dat$cig_0 >= 99]         <- NA
dat$meduc[dat$meduc == 9]          <- NA
dat$precare5[dat$precare5 == 5]    <- NA
dat$mracehisp[dat$mracehisp == 8]  <- NA
```

### 0.4 `dmar` is blank on ~20% of records and they are dropped silently  🔴

| year | blank `dmar` |
|---|---|
| 2024 (2M sample) | 402,966 (20.1%) |
| 2023 (300k sample) | 121,159 (40.4%) |

These are **not** junk records. Of 220,254 blank-`dmar` rows checked, 220,194 have a
valid `dbwt`, 220,254 have `mracehisp`, 220,254 have `meduc`. Blank `dmar` co-occurs
exactly with blank `mar_p`, i.e. it is a reporting-area pattern, not a per-record lapse.
Because these CSVs are ordered by state, the loss is **geographically clustered** — the
class-level counts are silently reweighted toward reporting states.

`ifelse(dmar == 1, 1, 0)` maps `NA` → `NA` → `na.omit()` drops the row, with no message.

**Fix (needs your call — see §5, Q1).** At minimum: make it explicit and quantify it in
the paper. Options are (a) drop with a documented missingness table and a sensitivity
arm, or (b) carry "unknown" as a third `dmar` level, which the Tier-1 categorical fix
makes tractable. A referee will ask about a 20% non-random deletion; better to answer it
in the methods than in review.

---

## Tier 1 — methods

### 1.1 Categorical split search is under-searching  🟠

`dm-cart.R:130,148`. `direction <- ux` returns the sorted raw level codes, so rpart only
evaluates the k−1 splits that are *contiguous in arbitrary code order*. For 5-level race
that is 4 of the 15 possible bipartitions, and the 4 it checks are determined by the
order the levels happen to be declared in.

Quantified on a case whose true optimum is `{Black, Asian} | {White, Hispanic, Other}`
(non-contiguous by code):

```
ordinal (current code): best goodness =  35.86
exhaustive:             best goodness = 174.95   ← 4.9x
```

**Fix, and it is verified to work inside rpart's user-method API.** Search all
2^(k−1)−1 bipartitions internally, then hand rpart a `direction` permutation with the
winning left-set placed first, so the optimum is reachable as a contiguous split under
that ordering:

```r
mysplit <- function(y, wt, x, parms, continuous = FALSE) {
  ux <- sort(unique(x)); k <- length(ux)
  if (k < 2) return(list(goodness = numeric(0), direction = ux))   # guards the 1:0 bug

  parent <- -log.dm.likelihood(colSums(y), alpha = alphavec)
  gof <- function(L) parent +
    log.dm.likelihood(colSums(y[ L, , drop = FALSE]), alpha = alphavec) +
    log.dm.likelihood(colSums(y[!L, , drop = FALSE]), alpha = alphavec)

  if (continuous) {
    return(list(goodness = pmax(sapply(seq_len(k - 1), function(i) gof(x <= ux[i])), 0),
                direction = rep(1, k - 1)))
  }
  if (k <= MAX.EXHAUSTIVE.LEVELS) {                 # 12 -> 2047 partitions, ~ms
    bg <- -Inf; bset <- ux[1]
    for (m in seq_len(2^(k - 1) - 1)) {
      inL <- as.logical(bitwAnd(m, 2^(0:(k - 1))))
      g <- gof(x %in% ux[inL]); if (g > bg) { bg <- g; bset <- ux[inL] }
    }
    ord <- c(bset, setdiff(ux, bset))
  } else {
    ord <- ux[order(sapply(ux, function(u) {        # fallback: order by P(LBW) profile
      s <- colSums(y[x == u, , drop = FALSE]); sum(s[lbw.cols]) / max(sum(s), 1)
    }))]
  }
  list(goodness  = pmax(sapply(seq_len(k - 1), function(i) gof(x %in% ord[seq_len(i)])), 0),
       direction = ord)
}
```

Confirmed end-to-end: rpart accepts the permuted `direction` and writes the intended
partition into `csplit` (`improve = 174.95`, left = `{B,A}`, right = `{W,H,O}`).

Two properties worth stating in the paper:

* **Binary predictors are untouched.** k = 2 → exactly one partition → identical to the
  current code. So this changes nothing about the thesis's binary results; it only makes
  the multi-level race variable legitimate.
* Cost is negligible at the level counts in play (k = 7 → 63 partitions per node).

### 1.2 `meduc == 3` should be `meduc >= 3`  🟠

The audit's suspicion is confirmed against the actual distribution — the field follows
the standard 2003-revision coding:

```
1 =  64,593 (<=8th)      4 = 327,212 (some college)   7 = 207,581 (master's)
2 = 146,959 (9-12 no dip) 5 = 162,627 (associate)      8 =  66,054 (doctorate/prof)
3 = 530,542 (HS grad/GED) 6 = 442,269 (bachelor's)     9 =  52,163 (unknown)
```

`meduc == 3` labels a mother with a bachelor's degree as **not** having completed high
school. Table 2.1 says "HS completed / otherwise". Use `meduc >= 3` and drop code 9.

This changes the meaning of an existing predictor, so the binary-race reproduction arm
will not numerically match the thesis after this fix. That is correct — the thesis
number was wrong. Worth a footnote.

### 1.3 `precare5 == 1` is defensible, but rename it  🟡

Coding is `1` = 1st trimester, `2` = 2nd, `3` = 3rd, `4` = no prenatal care, `5` =
unknown. So `== 1` correctly means "prenatal care began in the first trimester" — but
that is *initiation timing*, not adequacy (Kotelchuck/APNCU is the adequacy index, and it
is not in this extract). Rename the predictor to `precare1st` and restate the Table 2.1
gloss as "care began in 1st trimester" so the paper does not overclaim. Drop code 5.

---

## 2. Tier 2 — reproducibility gaps (audit items 3–7, all confirmed)

| # | Item | Verified | Fix |
|---|---|---|---|
| 3 | Depth-controlled comparison (§3.2.4, Figs 3.4–3.5) never runs | `awk` over `dm-cart.R`: **0 active lines after 412** | rewrite inside the model loop |
| 4 | `tree.depth.results` has no writer | only writer is in the dead block | write it; see below |
| 5a | High/low risk conditions on 3 of 7 predictors | `dm-cart.R:325–326` | exact class lookup, below |
| 5b | `risk_summary` never saved | `bootstrap.tree.results` holds only 3 elements | persist it |
| 6 | `top.vars.used` computed, never saved | 2021 results *do* have `top.var.df` — it was lost in the rewrite | restore `top.var.df` |
| 7 | `B <- 1000` | `dm-cart.R:265` | see below |

### 2.1 Depth data must come from the real ensemble

The plot in `visual.R` wants a per-split depth distribution. Take it from the bootstrap
trees you already build, inside the existing loop — no refit needed:

```r
node.ids   <- as.numeric(rownames(dm.tree.b$frame))
node.depth <- floor(log2(node.ids))                       # heap indexing: correct
internal   <- dm.tree.b$frame$var != "<leaf>"
depth.rows[[b]] <- data.frame(b        = b,
                              Variable = as.character(dm.tree.b$frame$var[internal]),
                              Depth    = node.depth[internal])
```

**Do not** copy the dead code's `floor(log2(which(frame$var == "cig_0")))` — `which()`
returns *row positions* in the frame, not node ids. Those two coincide only for a
perfectly full tree, so that line was quietly wrong wherever it appeared.

The maxdepth ∈ {2,3,4,5} comparison is a separate, cheap refit of the *observed* counts;
keep it, but persist to `paths$tree.depth.results(TARGET.YEAR)` with both objects:

```r
save(depth_df, pred.stats.compare, var.across.depths, tree.sum, file = paths$tree.depth.results(TARGET.YEAR))
```

(`depth_df` with columns `Variable`, `Depth` is exactly what `visual.R` expects — keeping
that name means `visual.R` needs no adaptation.)

### 2.2 High/low risk must be a class, not an average over classes

Current code selects **32 of 320 rows** (5 race × 3 conditions fixed leaves 2^4 = 16 free
combinations × ...), then averages them — which is why the thesis's sharp contrast came
out blended. Specify the full 7-dimensional profile and look up the single row:

```r
find.class <- function(profile) {
  hit <- Reduce(`&`, Map(function(v, val) as.character(X.matrix[[v]]) == as.character(val),
                         names(profile), profile))
  stopifnot(sum(hit) == 1)      # fails loudly if the design changes under you
  which(hit)
}
high.profile <- list(sex="1", dmar="0", mager="0", meduc="0", precare5="0", cig_0="1", race="Black")
low.profile  <- list(sex="0", dmar="1", mager="0", meduc="1", precare5="1", cig_0="0", race="White")
```

Then compute percentile intervals from `Yhat` over that single index and save a tidy
frame under the name `visual.R` already looks for:

```r
risk_summary <- rbind(
  data.frame(Subgroup="High-Risk Subgroup", Category=cat.labels,
             Mean=hi$means, Lower=hi$lwr, Upper=hi$upr),
  data.frame(Subgroup="Low-Risk Subgroup",  Category=cat.labels,
             Mean=lo$means, Lower=lo$lwr, Upper=lo$upr))
```

Note `Lower`/`Upper` rather than `SE`: `dm-cart.R` uses percentile bootstrap intervals
while `visual.R` reconstructs `Mean ± 1.96·SE`. Two different intervals for the same
quantity, in the same paper. Publish the percentile interval and drop the normal
approximation.

The two profiles above need your sign-off — see §5, Q2.

### 2.3 `B`

Make it an env-overridable constant so the final run is a flag, not an edit:

```r
B <- as.integer(Sys.getenv("DMCART_B", unset = if (interactive()) 200 else 10000))
```

Runtime scales with I: at 448 classes and B = 10,000 this is the long pole. Budget for it
and run it once, at the end, after Tier 0–1 land.

---

## 3. Race — expansion design (I disagree with audit item #1)

The audit recommends reverting to binary `mrace15` to match the thesis. That is the right
advice for *reproducing* the thesis and the wrong advice for *superseding* it, which is
what you're doing. Your own text already flags this as the weak point:

> `chapter2.tex:45` — "…the primary dichotomy, *though this choice is entirely
> arbitrary*. (NEED TO REVISIT)"

So: keep the expansion, make it a switch, and keep binary as a reproduction/sensitivity arm.

```r
RACE.SCHEME <- "full7"   # "binary" | "collapsed5" | "full7"
```

| scheme | groups | I | purpose |
|---|---|---|---|
| `binary` | Black / non-Black | 2^7 = **128** | reproduces thesis; sensitivity arm |
| `collapsed5` | White, Black, Hispanic, Asian, Other | 2^6·5 = **320** | current `data.R` behaviour |
| `full7` | full `mracehisp` | 2^6·7 = **448** | **recommended default** |

`full7` = White NH, Black NH, AIAN NH, Asian NH, NHOPI NH, Multiracial NH, Hispanic.
Feasibility, from the 2M-record 2024 sample scaled to ~3.6M births/yr:

| group | /2M | ≈/yr | ≈ mean births per class (÷64) |
|---|---|---|---|
| White NH | 954,674 | 1.72M | 26,900 |
| Hispanic | 540,238 | 973k | 15,200 |
| Black NH | 285,363 | 514k | 8,030 |
| Asian NH | 128,725 | 232k | 3,620 |
| Multiracial NH | 50,688 | 91k | 1,430 |
| AIAN NH | 11,639 | 21k | 328 |
| NHOPI NH | 5,869 | 10.6k | **165** |

The thin tail (AIAN, NHOPI) is real, and it is also exactly the argument for this paper:
**the Dirichlet prior is what makes the disaggregation admissible.** A plain multinomial
CART would produce unstable per-class estimates at n≈165 spread over 10–11 categories;
the informed DM prior shrinks them toward the marginal LBW profile and the parametric
bootstrap reports honest uncertainty on them. Currently "we can handle sparsity" is
asserted; with `full7` it is demonstrated. That is a substantially stronger contribution
than the binary version, and it retires the "entirely arbitrary" caveat.

Two things this obligates:

1. **Tier 1.1 is a prerequisite, not an option.** With 7 levels the ordinal search sees
   6 of 63 partitions. Publishing `full7` without the exhaustive fix would be worse than
   publishing binary.
2. **Report the partitions, don't just report "race was used."** The scientifically
   interesting output is now *which grouping of races* the DM criterion selects, and how
   stable that grouping is across the 10,000 bootstrap trees. That comes straight out of
   `csplit` and is new content — see §4.

Collapsing AIAN/NHOPI/Multiracial into "Other" (`collapsed5`) is defensible but it is the
same arbitrariness the thesis apologises for, one level up. If a referee pushes back on
sparsity, `collapsed5` is the fallback, and having it behind a switch means producing that
arm is a one-line change.

Out of scope but worth knowing: `mracehisp` already folds Hispanic origin into race.
`mrace15` (15 bridged categories) × `mhisp_r` (6) is present in the extract if you ever
want a finer cross — I'd not go there for this paper.

---

## 4. Visualization refresh — bootstrap output only

Your constraint: **only outputted bootstrap data.** Concretely that means:

* **`visual.R` does not currently parse.** Line 18 is `} else {»` — a stray `»`.
  Verified: `parse("visual.R")` → `18:9: unexpected input`. Fix first; everything else in
  this section is moot until it does.
* **Delete every fallback branch.** The `risk_summary` fallback (`c(rep(0.016,10), 0.840)`,
  `visual.R:236–244`) is still live and still fabricates the headline result. The
  `depth_df` random-`sample()` fallback is already commented out — delete it rather than
  leave it as a comment someone re-enables. Replace both with:

  ```r
  if (!file.exists(f)) stop(sprintf("Missing %s — run: Rscript dm-cart.R %d-%d", f, PRIOR.YEAR, TARGET.YEAR))
  ```
* **Load into an environment, never the global.** `visual.R` does `load(...)` at top level
  after `rm(list=ls())`, so a stale object from a previous sourcing can satisfy an
  `exists()` check. Use `e <- new.env(); load(f, envir=e)` and read explicitly.
* **Percentile intervals everywhere**, from `Yhat` — not `Mean ± 1.96·SE` (§2.2).
* **Drop the hardcoded header.** `MODEL_YEAR <- 2024; WITH_2_5KG <- FALSE` at
  `visual.R:22–23` means the script only ever renders one of the four model/cohort
  combinations, and silently. Take `commandArgs()` like the other three scripts and loop
  over `PROCESS.WITH.2.5KG`.
* **Also fix**: the caption at `visual.R:319–322` describes the *old* binary risk profiles
  ("White", no race term in the high-risk line). Generate captions from the profile lists
  in §2.2 so they cannot drift from the model again.

New plots the expanded race variable earns, all from bootstrap output:

1. **Race-partition stability.** Across B trees, tabulate the `csplit` row wherever the
   split variable is `race`; plot frequency of each realized bipartition. Answers "does
   the DM criterion consistently separate the same race groups?" — this is the paper's
   new headline figure.
2. **Per-race LBW probability with percentile bands**, one panel per race group,
   marginalizing `Yhat` over the other 6 predictors. Shows the tail groups' wider
   intervals honestly rather than hiding them.
3. Keep the existing depth-distribution and variable-frequency figures; they become
   correct for free once §2.1 writes real data.

---

## 5. Decisions I need from you

**Q1 — `dmar` missingness (§0.4).** ~20% of 2024 and ~40% of the 2023 head sample,
state-clustered. (a) drop + documented missingness table + sensitivity arm, or (b) carry
`unknown` as a third level of `dmar` (now cheap, given Tier 1.1)? I lean (b) for the main
model with (a) as the sensitivity arm — dropping a fifth of the data non-randomly is the
kind of thing that sinks a submission.

**Q2 — risk profiles (§2.2).** The thesis's classes i=69 / i=28 are indices into the
128-class binary design and do not survive the re-encoding. Confirm the two 7-dimensional
profiles, or give me the ones you want. I've assumed high = male, unmarried, ≤33, no HS,
no 1st-trimester care, smoker, Black; low = female, married, ≤33, HS+, 1st-trimester care,
non-smoker, White.

**Q3 — `meduc >= 3` (§1.2).** Confirms the audit; changes an existing predictor's meaning,
so the binary arm will no longer reproduce the thesis's exact numbers. Confirm you want
the corrected definition (I recommend yes, with a footnote).

---

## 6. Sequencing

Tier 0 first, then Tier 1, then re-run end to end; Tier 2–4 are cosmetic until the numbers
are trustworthy.

1. `util.R` — return `target.year` / `cohort`; add the `PRIOR.YEAR < TARGET.YEAR` guard. *(§0.2)*
2. `data.R` — sentinel `NA`s, `dmar` policy, `meduc >= 3`, `RACE.SCHEME` switch, missingness table. *(§0.3, §0.4, §1.2, §3)*
3. `quantile.R` — no logic change; drop the duplicate 2.1 GB `read.csv` (lines 40 and 46 read the same file twice).
4. `dm-cart.R` — `dm.leaf.probs()`, exhaustive `mysplit()`, depth capture, exact-class risk, persist `top.var.df` + `risk_summary` + `tree.depth.results`, `B` via env. *(§0.1, §1.1, §2)*
5. Re-run: `quantile.R` → `data.R` → `dm-cart.R` with `DMCART_B=200` to smoke-test, then `DMCART_B=10000`.
6. `visual.R` — parse fix, remove fallbacks, percentile bands, new race figures. *(§4)*
7. `chapter2.tex` / `chapter3.tex` — I = 448 not 128, race table, `precare1st` gloss,
   `meduc` footnote, missingness paragraph, and retire the "entirely arbitrary" caveat.

Performance note, orthogonal to correctness: `read.csv()` on a 2.1 GB file is the slowest
step in the pipeline by a wide margin. `data.table::fread(path, select = raw.cols)` is
installed and reads only the 9 needed columns — worth doing while `data.R` is open anyway.

---

## Appendix — what the audit got right, unchanged

Recorded so it doesn't get re-litigated:

* `log.dm.likelihood()` = Eq. 2.3, exactly.
* `mysplit()`'s `parent.dev - child.dev` = the ΔL rule of §2.4.4.
* Two-tier bootstrap = Eqs. 3.1–3.2, including holding `multinomial.probs` fixed across
  replicates (no double Monte Carlo).
* `rowMeans` vs `rowSums` for `row.probs` — genuinely equivalent after normalization, as
  the audit noted. Not a bug.
* `quantile.R`'s informed prior = Table 2.2 / Fig 3.1.
* `sex`, `dmar`, `mager` (>33), `cig_0` (>0) encodings = Table 2.1.

Minor items, no action required beyond a cleanup pass:

* `final.data$x.y` (mean log birthweight) is never NA-filled while `x.x` is set to `0`;
  `0` is also wrong for an empty class — `NA` is the honest value. Neither feeds the model.
* `myinit()` declares `parms = list(alpha = 1)` while `myeval`/`mysplit` read the global
  `alphavec` by closure. Works, but pass `alphavec` through `parms` so the method is
  self-contained and testable.
* `mysplit`'s `1:(length(ux) - 1)` iterates `i = 1, 0` when a node has one unique x value.
  rpart doesn't appear to call it that way, but the `k < 2` guard in §1.1 removes the
  landmine.
