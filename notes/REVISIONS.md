# REVISIONS — DM-CART 2023/2024 revision log

Companion documents: `REVISION-PLAN.md` (the original audit), `notes.md` (your own
interpretation — not written to by this log).

**Jump to:** [Tier 0/1](#tier-01) · [Tier 2](#tier-2-pass) · [Tier 2b](#tier-2b-pass--runtime-and-figures) ·
[Tier 2c](#tier-2c-pass--legibility-and-finding-the-risk-profiles) ·
[Tier 3](#tier-3-pass--race-expansion-extreme-ten-analysis-and-the-last-dead-code) ·
**[Write-up guide](#write-up-guide--what-to-say-where-and-how-to-restructure)**

<a name="tier-01"></a>
## Tier 0 / Tier 1 pass

Date: 2026-08-10 · Scope: `util.R`, `data.R`, `quantile.R`, `dm-cart.R`
Companion document: `REVISION-PLAN.md` (the full audit; this file records only what was
actually **changed**).

Every revision carries an ID (`R1`…`R7`). The same IDs appear as `# REVISION Rn`
comments at the exact code sites, so you can jump between this file and the source.

**Not in this pass** (deliberately — you scoped it to Tier 0/1): the depth-comparison
rewrite, persisting `tree.depth.results` / `risk_summary` / `top.var.df`, the exact
single-class high/low-risk lookup, `B = 10000`, the race expansion to 7 groups, and all
of `visual.R` (which still does not parse — see §"Still broken"). Those remain Tier 2+.

---

## Index

| ID | What | Where | Severity |
|---|---|---|---|
| [R1](#r1) | Prior year == target year (double-dipping) | `util.R`, `data.R`, `dm-cart.R`, `quantile.R` | 🔴 invalidated all results |
| [R2](#r2) | `mysplit()` under-searched categorical splits | `dm-cart.R:137–229` | 🟠 wrong splits on race |
| [R3](#r3) | NCHS sentinel "unknown" codes treated as data | `data.R:74–113`, `quantile.R:53` | 🔴 biased the LBW/NBW ratio |
| [R4](#r4) | `dmar` blank → 20% silently deleted | `data.R:25–34, 174–200` | 🔴 non-random deletion |
| [R5](#r5) | `meduc == 3` → `meduc >= 3` | `data.R:150–162` | 🟠 predictor meaning wrong |
| [R6](#r6) | `precare5` renamed `precare1st` | `data.R:163–172` | 🟡 overclaim |
| [R7](#r7) | `mypred()` sent rows to wrong leaves | `dm-cart.R:235–290, 365` | 🔴 invalidated all `Yhat` |
| [R8](#r8) | Node-label callbacks produced useless labels | `dm-cart.R` (`myinit`) | 🟠 every leaf labelled "1" |
| [R9](#r9) | Prototype tree plotting (`plot.rpart`) | `dm-cart.R` (`plot.dm.tree`) | 🟢 new figure |
| [R10](#r10) | Depth-controlled comparison reinstated | `dm-cart.R` (before bootstrap) | 🟠 never ran |
| [R11](#r11) | `depth_df` from real bootstrap trees | `dm-cart.R` (bootstrap loop) | 🔴 figures were random noise |
| [R12](#r12) | High/low risk = single class + `risk_summary` | `dm-cart.R` | 🔴 headline numbers were invented |
| [R13](#r13) | `top.var.df` restored, `B` configurable | `dm-cart.R` | 🟡 |
| [R14](#r14) | Split search made ~2x cheaper (+ runtime analysis) | `dm-cart.R` (`mysplit`) | 🟢 performance |
| [R15](#r15) | Tree figures rebuilt for legibility + branch labels | `dm-cart.R` (`plot.dm.tree`) | 🟠 figures unreadable |
| [R16](#r16) | Type size solved from geometry; tighter margins | `dm-cart.R` (`draw.dm.tree`) | 🟢 legibility |
| [R17](#r17) | Extreme-risk profiles found FROM THE DATA | `dm-cart.R` | 🔴 profiles were assumed |
| [R18](#r18) | Subgroup bars drop NBW; totals in the subcaption | `dm-cart.R` | 🟠 within-LBW shape invisible |
| [R19](#r19) | Two-pane tree tracing each extreme profile | `dm-cart.R` | 🟢 new figure |
| [R20](#r20) | Race expanded to all 7 `mracehisp` groups | `data.R` | 🔴 headline change |
| [R21](#r21) | `fread` replaces `read.csv` | `data.R` | 🟢 minutes → 17 s |
| [R22](#r22) | All leaf labels rotated; wrapped level sets | `dm-cart.R` | 🟠 legibility |
| [R23](#r23) | Top/bottom **10** profiles + shared-feature analysis | `dm-cart.R` | 🟢 new result |
| [R24](#r24) | `visual.R` runs at all; every fallback deleted | `visual.R` | 🔴 was non-functional |
| [R25](#r25) | Empty-class means are `NA`, not `0` | `data.R` | 🟡 |
| [R26](#r26) | One combined prior figure, intervals in the legend | `quantile.R` → `visual.R` | 🟢 figure |
| [R27](#r27) | Tree figures cut 6/arm → 2/arm; renderers moved to `util.R` | `util.R`, `dm-cart.R`, `visual.R` | 🟢 figures |
| [R28](#r28) | NBW restored to the prior; risk landscape to two panes | `visual.R`, `dm-cart.R` | figures |
| [R29](#r29) | `results/<cohort>/article/`; two-model tree figure | `util.R`, `visual.R`, `dm-cart.R` | figures |
| [R30](#r30) | Bigger, horizontal leaf labels; larger ranking figure | `util.R`, `visual.R` | figures |
| [R31](#r31) | Tree label placement; ggplot risk-category figure; caption fixes | `util.R`, `visual.R`, LaTeX | figures |
| [R32](#r32) | Half-page trees; depth growth to article/; commonality table | `visual.R` | figures |
| [R33](#r33) | Table fits one column; 3-row leaf stagger; single-column prior | `util.R`, `visual.R` | figures |
| [R34](#r34) | Collision-avoided branch labels; legend/layout | `util.R`, `visual.R` | figures |
| [R35](#r35) | First birth-weight bin was left-open, dropping 24 births | `data.R` | correctness |

---

<a name="r1"></a>
## R1 — The prior year and the target year were the same year

### What was wrong

`get_paths()` resolved the cohort correctly *internally* — it computed
`target.year <- 2024` — and then returned `year = year`, the **requested** year. Since
every script does `TARGET.YEAR <- paths$year`, and every script is invoked with the
**prior** year (2023), `TARGET.YEAR` came out as **2023**.

```
BEFORE:  get_paths(2023) -> $year = 2023, $prior.year = 2023, $cohort = NULL
AFTER :  get_paths(2023) -> $year = 2024, $target.year = 2024,
                            $prior.year = 2023, $cohort = "2023-2024"
```

Consequences, all of them silent:

1. **The no-double-dipping design was not in effect.** The quantile cutpoints and the
   informed Dirichlet prior were built from the same 2023 file that the counts were
   built from. That is precisely what §2.x of the thesis claims to avoid.
2. **The 2024 data was never read.** Not once, in any run that has happened.
3. `cohort.dir` was computed and never returned, so `sprintf("... %s ...", paths$cohort)`
   evaluated to `character(0)` and the cohort banner in all three scripts printed
   **nothing**. That is why nobody caught 1 and 2.

### What changed

**`util.R`** — the returned list now carries the resolved cohort fields, and `year` is an
alias for `target.year` so the accessor defaults keep meaning "the year being modeled":

```r
year           = target.year,   # was: year  (the requested year)
target.year    = target.year,
prior.year     = prior.year,
requested.year = requested.year,
cohort         = cohort.dir,
```

A guard now refuses to build a cohort where the prior is not strictly earlier:

```r
if (prior.year >= target.year) stop(...)
```

The file accessors were changed from `function(yr = year)` to
`function(yr = target.year)`. This matters: those are **closures over the function
argument** `year`, not over the list field, so leaving them as `yr = year` would have
kept them pinned to 2023 even after the list was fixed. (This is exactly how the bug
survived the first attempt at the fix during this pass.)

**`data.R:37`, `dm-cart.R:31`** — `TARGET.YEAR <- init.paths$target.year`.

**`data.R:339`** — the rebinned counts now save under `TARGET.YEAR`, not `PRIOR.YEAR`.
`dm-cart.R` already *loaded* `rebin.data(TARGET.YEAR)`; the save/load mismatch was
invisible only because both years were wrongly equal.

**`quantile.R:46`** — deleted a duplicate `natalitydata <- read.csv(paths$raw.csv())`.
It re-read a 2.1 GB file that had already been read, and after R1 its default argument
resolves to the **target** year — so it would have pulled 2024 into a script whose only
job is to summarise 2023.

### Verified

```
cohort 2023-2024 | prior 2023 | target 2024
  raw.csv(prior)        natality2023us-original.csv
  raw.csv(target)       natality2024us-original.csv
  cutpoints(prior)      quantile_cutpoints_2023.RData
  rebin.data()          birthweight_data_rebin_2024.RData
  bootstrap.results()   bootstrap_tree_results_2024.RData
2020-2021 cohort still resolves: prior 2020, target 2021
```

### What you can delete, and what you must keep

**KEEP — still valid, built from the prior year, unaffected by R1:**

```
data/rebin/2023-2024/quantile_cutpoints_2023.RData
data/rebin/2023-2024/informed_prior_2023.RData
data/rebin_without_2.5kg/2023-2024/quantile_cutpoints_2023.RData
data/rebin_without_2.5kg/2023-2024/informed_prior_2023.RData
```

⚠️ These four are correct with respect to R1 but **stale with respect to R3** — they were
built without dropping `dbwt == 9999`, which inflated `prop.normal` and therefore the
Normal component of `alphavec`. **Re-run `quantile.R` to regenerate them.** Keep them
only until you do.

**DELETE — every one of these is a double-dipped 2023-as-target artifact:**

```
data/rebin/2023-2024/birthweight_data_rebin_2023.RData
data/rebin/2023-2024/dm_tree_rebin_2023.RData
data/rebin_without_2.5kg/2023-2024/birthweight_data_rebin_2023.RData
data/rebin_without_2.5kg/2023-2024/dm_tree_rebin_2023.RData
results/2023-2024/fullmodel/bootstrap_tree_results_2023.RData
results/2023-2024/lbwmodel/bootstrap_tree_results_2023.RData
results/2023-2024/fullmodel/plots/*          (all regenerate)
results/2023-2024/lbwmodel/plots/*           (all regenerate)
```

The `_2023` suffix on target-year artifacts is the tell: after R1 these are written as
`_2024`, so nothing will overwrite them and they will sit around looking authoritative.
Delete them rather than let them linger.

```bash
rm -f data/rebin/2023-2024/birthweight_data_rebin_2023.RData \
      data/rebin/2023-2024/dm_tree_rebin_2023.RData \
      data/rebin_without_2.5kg/2023-2024/birthweight_data_rebin_2023.RData \
      data/rebin_without_2.5kg/2023-2024/dm_tree_rebin_2023.RData
rm -rf results/2023-2024
```

**KEEP — untouched by this pass:** everything under `data/rebin*/2020-2021/` and
`results/2020-2021/`. Those are the thesis-era artifacts. They were produced by the old
code and are not reproducible from the current source, so treat them as a frozen
historical record, not as something to re-run. Note the 2020-2021 cohort resolved
*correctly* before this pass only by coincidence of how it was invoked.

---

<a name="r2"></a>
## R2 — `mysplit()` searched a small, arbitrary subset of categorical splits

**Where:** `dm-cart.R:137–229` (the function), `dm-cart.R:77` (moved `lbw.cols`).

### What was wrong — the part you asked me to explain plainly

rpart does not let a user-written split function say "put these levels left and those
right." The protocol is narrower: you return `direction`, an **ordering** of the levels,
and rpart then only ever considers splits that are **contiguous in that ordering** —
first level vs. the rest, first two vs. the rest, and so on. That is `k-1` candidate
splits for a `k`-level factor.

The old code returned `direction <- ux`, the levels sorted by their **raw integer
codes**. So the only splits ever evaluated were the ones that happen to be contiguous in
whatever order the levels were declared in.

For a **binary** predictor that costs nothing: `k = 2` gives exactly one possible split,
and sorted order finds it. This is why the bug was invisible for the entire thesis — all
seven predictors were binary.

For **race with 5 levels** it matters a great deal. There are `2^(k-1) - 1 = 15` distinct
ways to divide 5 levels into two groups. Sorted-code order reaches only 4 of them, and
*which* 4 is determined by nothing but the declaration order
`White, Black, Hispanic, Asian, Other`. A grouping like
`{Black, Asian} | {White, Hispanic, Other}` is **not reachable at all** — the two levels
are not adjacent in code order, so no contiguous cut can produce it.

Measured on a case constructed so that grouping is the true optimum:

```
sorted-code order (old) : best goodness =  35.86
exhaustive        (new) : best goodness = 174.95      4.9x better
```

The tree was not choosing a worse split on purpose. It could not see the better one.

### How it is fixed

rpart won't take a subset, but it will take **any permutation** as `direction`. So:

1. Enumerate all `2^(k-1) - 1` bipartitions ourselves and score each with the same ΔL
   criterion (`parent.dev - child.dev`) that was already in use — the splitting rule
   itself is unchanged.
2. Return an ordering with the **winning left-group placed first**. The optimum is then
   reachable as "first `m` vs. rest," which is a split rpart *can* express.
3. rpart writes exactly that partition into `fit$csplit`.

Guarantees worth knowing:

* **Binary predictors are mathematically unaffected** (`k = 2` → one partition → same
  answer as before). Nothing about a binary-only result changes.
* Cost is negligible: `k = 5` → 15 partitions per node; `k = 7` → 63. Capped by
  `MAX.EXHAUSTIVE.LEVELS <- 12` (2047 partitions), above which it falls back to ordering
  levels by their LBW share — still a heuristic, but monotone in the quantity of
  interest rather than arbitrary. Nothing in this analysis comes near the cap.

Two smaller fixes in the same function:

* **A `k < 2` guard.** The old `for (i in 1:(length(ux) - 1))` becomes `1:0` when a node
  holds a single distinct value, which loops over `i = 1, 0`, indexes `ux[0]`, and
  returns a `goodness` vector of the wrong length. rpart does not appear to call it that
  way in this design, but the landmine is now removed.
* **`continuous` is honoured.** The old code ignored the flag and always returned the
  categorical format; a numeric predictor would have been silently misread. Still unused
  here (every predictor is a factor), but it is now correct if that changes.

`lbw.cols` / `nbw.cols` moved from just-before-the-bootstrap to just-after `Y.df`
(`dm-cart.R:77`), because the fallback branch references `lbw.cols` and the primary tree
is fit before their old definition point.

### Verified

On a 480-class synthetic design matching the new encoding, with a planted
non-contiguous signal:

```
root race split -> left: Black,Asian | right: White,Hispanic,Other
planted truth   -> {Black,Asian}                                    ✓ recovered
```

And on 3-level `dmar` with `Unknown` planted to behave like unmarried:

```
split -> left: 1 | right: 0,Unknown        ✓ recovered
```

That second result is the point of R4 below: the model can now *discover* that the
unknown group behaves like the unmarried group, instead of us assuming it.

---

<a name="r3"></a>
## R3 — NCHS "unknown" sentinel codes were entering the model as real values

**Where:** `data.R:74–113` (cleaning + report), `quantile.R:53–58` (prior side).

### What was wrong

`na.omit()` catches R's `NA`. NCHS encodes "unknown / not stated" as ordinary **in-range
integers**, which `na.omit()` cannot see. Counts are from a 2,000,000-record sample of
`natality2024us-original.csv`:

| field | code | meaning | count / 2M | what it silently became |
|---|---|---|---|---|
| `dbwt` | 9999 | weight not stated | 1,648 | **counted as Normal (> 2500 g)** |
| `cig_0` | 99 | smoking not stated | 11,202 | smoker = 1 |
| `meduc` | 9 | education not stated | 52,163 | "did not complete HS" = 0 |
| `precare5` | 5 | care timing not stated | 35,060 | "not 1st trimester" = 0 |
| `mracehisp` | 8 | origin not stated | 22,804 | *(already handled correctly)* |

`dbwt == 9999` is the most damaging: it inflated the Normal-birth-weight category, which
is the denominator for the whole LBW story — and it did so on **both** sides, in the
counts *and* in the `prop.normal` that sets the Normal component of `alphavec`.

### What changed

All five sentinels are converted to `NA` **before any recoding**, in `data.R`:

```r
raw$dbwt[raw$dbwt == 9999]        <- NA
raw$cig_0[raw$cig_0 >= 99]        <- NA   # 98 is a valid count; 99 is the sentinel
raw$meduc[raw$meduc == 9]         <- NA
raw$precare5[raw$precare5 == 5]   <- NA
raw$mracehisp[raw$mracehisp == 8] <- NA
```

The same drop is applied in `quantile.R` before any proportion or quantile is computed,
so the prior is built on stated weights only.

Structurally, `data.R` now keeps the **raw** fields in `raw` and writes encodings into a
fresh `dat`. Previously everything was written back into one frame while the recodes
read from `natalitydata`; separating them makes the block genuinely re-runnable and
removes the "not meaningful for factors" foot-gun the old comment described.

A **missingness report** is printed and saved with the data (`missing.report`), so the
paper can state exactly what was excluded rather than losing rows inside `na.omit()`:

```
 variable n_missing     pct
      sex         0  0.000%
     dmar    220254 55.063%     <- see R4 (this slice is state-ordered; full-year ≈20%)
mracehisp      9145  2.286%
    mager         0  0.000%
    meduc     21167  5.292%
 precare5      7644  1.911%
    cig_0      1141  0.285%
     dbwt       192  0.048%
```

---

<a name="r4"></a>
## R4 — `dmar` blank on ~20% of records, silently deleted

**Where:** `data.R:25–34` (the switch), `data.R:174–200` (the encoding).

### Your instruction, and why I could not follow it literally

You asked to **impute `dmar` by region**, with stated reasoning for the regional scheme.
That is not possible from this file, and I want to be explicit about why rather than
quietly substituting something else.

**Every geographic identifier in the public-use natality extract is suppressed** — blank
in 100% of records. Verified directly on `natality2024us-original.csv`:

| column | distinct values in 300k records |
|---|---|
| `ocntyfips` (occurrence county FIPS) | blank only |
| `ocntypop` (occurrence county pop.) | blank only |
| `rcnty` (residence county) | blank only |
| `rcnty_pop`, `rcity_pop` | blank only |
| `mrterr`, `octerr`, `mbcntry` | blank only |

There is no state, county, division, or region field to condition on. The only
geography-adjacent variables that survive are `mbstate_rec` (mother born in US / outside
US) and `restatus` (residence status), and neither identifies a place.

Worse for imputation specifically: the missingness mechanism **is** the unobserved
geography — blank `dmar` co-occurs exactly with blank `mar_p`, i.e. it is a reporting-area
pattern. So MAR given the observed covariates is not credible either, and a multiple
imputation would be encoding an assumption we have no way to check. Imputing here would
manufacture precision rather than recover information.

### What I did instead — and it achieves your actual goal

Your stated goal was *"I would like to not get rid of these records since we included
this in our results before."* That is achievable without imputing anything.

`DMAR.POLICY <- "unknown-level"` (the new default) keeps **every** record and makes
`dmar` a three-level factor:

```r
dat$dmar <- factor(ifelse(is.na(raw$dmar), "Unknown",
                   ifelse(raw$dmar == 1, "1", "0")),
                   levels = c("0", "1", "Unknown"))
```

`DMAR.POLICY <- "drop"` reproduces the old listwise deletion, now explicit and reported,
for the sensitivity arm.

This is strictly better than imputation here, and it is worth making the argument
explicitly in the methods section:

* **No unverifiable assumption.** Imputation must assume something about how the unknown
  group is distributed. This assumes nothing.
* **The model answers the question empirically.** R2 made `k`-level factors first-class,
  so DM-CART now evaluates all three groupings of `{unmarried, married, unknown}` and
  reports which one the likelihood prefers. If `Unknown` consistently groups with
  `unmarried` across bootstrap replicates, that is *evidence* about the reporting areas —
  a finding, not an assumption. (Verified working: see R2's second test.)
* **It is honest about the uncertainty** rather than hiding it inside imputed values.

**Cost to be aware of:** `dmar` gaining a third level changes the class count from
`2^6 × 5 = 320` to `2^5 × 3 × 5 = **480**`. The hardcoded `320` literal in `data.R` has
been replaced with a computed expectation and a `stopifnot` (see R5 notes below), and
this number has to be updated in the manuscript.

### Verified — on 400,000 real 2024 records

```
complete-case: 400000 -> 369621 (30379 dropped, 7.59%)
OLD pipeline would have retained: 179746 (44.94%)
records RECOVERED by the 'unknown-level' policy: 189875

dmar:      0 = 77333    1 = 97713    Unknown = 194575
I = 480 classes | interaction levels = 480 ✓ | classes with zero obs: 29
```

(This slice is the head of the file and therefore state-ordered, so its 55% blank rate is
not the full-year rate — the whole-file figure is ~20%. The recovery mechanism is what
this demonstrates, not the rate.)

---

<a name="r5"></a>
## R5 — `meduc == 3` → `meduc >= 3`

**Where:** `data.R:150–162`.

`meduc` is the NCHS 2003-revision education recode:

```
1 = 8th grade or less        5 = associate degree
2 = 9th–12th, no diploma     6 = bachelor's
3 = HS graduate or GED       7 = master's
4 = some college, no degree  8 = doctorate / professional
9 = not stated  (-> NA, see R3)
```

Table 2.1 defines this predictor as **"high school completed / otherwise"**, but `== 3`
selects *only* mothers who stopped at a diploma. Under the old rule a mother with a
bachelor's degree was recorded as **not having completed high school**. The distribution
in the data confirms the coding above, so this is unambiguous.

```r
dat$meduc <- factor(ifelse(raw$meduc >= 3, 1, 0), levels = c(0, 1))
```

**This is a large change, not a cosmetic one.** On the 400k sample:

```
OLD (meduc == 3): 107,867 of 378,833 flagged "completed HS"  (28.5%)
NEW (meduc >= 3): 328,669 of 369,621 flagged "completed HS"  (88.9%)
```

**Consequence for reproducibility:** any binary-race reproduction arm will no longer
match the thesis's published numbers, because the thesis's numbers were computed with the
wrong definition. That is the correct outcome, but it needs a footnote in the manuscript
rather than being discovered by a referee.

Also in this area (`data.R:238`): the hardcoded `cat("num_node: ... expect 2^6 * 5 = 320")`
is now computed from the actual level sets and asserted:

```r
n.levels.each    <- sapply(X[, feature.names], nlevels)
expected.classes <- prod(n.levels.each)
stopifnot(num_node == expected.classes)
```

A hardcoded literal would go stale silently under R4 (dmar → 3 levels) and again under
the planned race expansion.

---

<a name="r6"></a>
## R6 — `precare5` renamed to `precare1st`

**Where:** `data.R:163–172`, plus `feature.names` and the save list.

The coding is `1` = 1st trimester, `2` = 2nd, `3` = 3rd, `4` = no prenatal care,
`5` = not stated. So `== 1` was always the *right test* — the problem was the **name and
the gloss**. Table 2.1 described this as "adequate prenatal care," but the variable
measures **when care began**, not whether it was adequate. Adequacy is the
Kotelchuck/APNCU index, which is not in this extract.

```r
dat$precare1st <- factor(ifelse(raw$precare5 == 1, 1, 0), levels = c(0, 1))
```

No numeric result changes. Update Table 2.1's gloss to **"prenatal care began in the 1st
trimester"** so the paper does not overclaim.

`feature.names` is now saved alongside the data so downstream consumers stop assuming the
old column name.

---

<a name="r7"></a>
## R7 — `mypred()` sent rows to the wrong leaves

**Where:** `dm-cart.R:235–290` (replacement), `dm-cart.R:365` (call site).
The old `mypred()` is **deleted**, not commented out.

### What was wrong

The old function walked the tree by hand, deciding direction with:

```r
split_val <- fit$splits[split_rows[1], "index"]
go_left   <- as.numeric(cur_val) <= split_val
```

Two independent errors:

1. **For a categorical split, `fit$splits[, "index"]` is not a cutpoint.** It is the
   **row number of `fit$csplit`** — the table that actually records which levels go left
   (`1` = left, `3` = right, `2` = not present at this node). Comparing a factor's
   integer code against a `csplit` row index is meaningless. **Every predictor in this
   model is a factor**, so this was wrong at every split, including the binary ones. It
   was never a race-only problem.
2. **`split_rows[1]`** takes the first `fit$splits` row bearing that variable's name, so
   a variable reused at several nodes was always routed by its *first* split's
   parameters, whichever node was actually being visited.

Measured against correct leaf assignment on a controlled 4-leaf tree:

```
rows landing in the wrong leaf: 16 of 20
max absolute probability error: 0.612
row 1 correct: 0.2314 0.2675 0.2548 0.2463
row 1 old    : 0.7400 0.0912 0.0874 0.0814
```

`Yhat` is built from this, and every high/low-risk figure is built from `Yhat`.

### How it is fixed

No traversal is needed. `X.matrix` is the **complete enumeration of classes**, and every
tree here is fit on exactly those rows — so rpart has already done the routing.
`fit$where` holds the frame **row index** (not the node id) of the leaf each training row
fell into, and `fit$frame$yval2` holds the per-node count vector returned by `myeval()`.

```r
dm.leaf.probs <- function(fit, alphavec, cat.names = colnames(Y.df)) {
  stopifnot(!is.null(fit$where), is.matrix(fit$frame$yval2))
  counts <- fit$frame$yval2[fit$where, , drop = FALSE]
  probs  <- t(apply(counts, 1, function(z) (z + alphavec) / sum(z + alphavec)))
  colnames(probs) <- cat.names; rownames(probs) <- NULL
  probs
}
```

**Two traps deliberately avoided — do not "simplify" into either of them:**

* `predict.rpart()` on a custom method returns an `n × 1` yval, not the count matrix.
* `rpart:::pred.rpart()` **segfaults** on custom-method fits. Verified: R aborts with
  `*** caught segfault *** cause 'invalid permissions'`.

**Scope limit, stated in the code comment too:** `dm.leaf.probs()` is valid only for the
rows the tree was *fit on*. That is all this design ever needs, because the class
enumeration is the design matrix. If you ever need to score genuinely new rows, the only
correct route is decoding `csplit` **per node** (via `frame$var` and the node's own
`splits` row), never per variable.

### Verified

480-class design, 11 categories:

```
dim: 480 x 11                                  ✓
rows summing to 1: 480 of 480                  ✓
all finite and in [0,1]: TRUE                  ✓
matches independent recomputation: TRUE        ✓
```

---

---
---

# Tier 2 pass

Second pass, same conventions. IDs `R8`–`R13`. Everything here is inside `dm-cart.R`.

Ordering note: the prototype tree and the depth comparison were deliberately placed
**before** the bootstrap. They cost seconds and they are the figures you read to
understand the model; the bootstrap costs hours. Nothing in either section depends on
bootstrap output, so you can now interrupt a run after the trees are drawn and still have
usable figures.

---

<a name="r8"></a>
## R8 — Node-label callbacks produced useless labels

**Where:** `dm-cart.R`, the `summary` / `print` / `text` closures inside `myinit()`.

### What was wrong

rpart passes these callbacks `yval = frame$yval2`, which for a user-defined method is a
**matrix** — one row per node, one column per birth-weight category. The old `text`
callback did:

```r
total.counts <- sum(yval)                       # sums the ENTIRE matrix -> scalar
lbl <- paste0(format(which.max(total.counts)))  # which.max of a scalar is always 1
```

So **every leaf of every plotted tree was labelled `1`**, regardless of composition. The
`print` callback had the matching problem and emitted `Deviance: NULL`.

This is very likely a large part of the "trouble with `plot.rpart()`" — the plot was
being produced, it just carried no information.

### What changed

Both callbacks now handle the matrix case and report the **modal birth-weight category**
per node. The `text` callback additionally folds in the node's posterior **P(LBW)**,
using the same `(counts + α)/(N + α₀)` form as `dm.leaf.probs()`:

```r
m     <- if (is.matrix(yval)) yval else matrix(yval, nrow = 1)
modal <- apply(m, 1, which.max)
p.lbw <- apply(m, 1, function(z) { post <- (z + alphavec)/sum(z + alphavec); sum(post[lbw.cols]) })
lbl   <- paste0("k", modal, "\nP(LBW)=", sprintf("%.3f", p.lbw))
if (use.n) lbl <- paste0(lbl, "\nn=", n)
```

P(LBW) goes **in the label** rather than being drawn by a separate `text()` call — my
first attempt annotated leaves separately and it collided with the `n=` line. Putting it
in the label lets `text.rpart` handle multi-line centring.

---

<a name="r9"></a>
## R9 — Prototype tree plotting that actually works

**Where:** `dm-cart.R`, new `plot.dm.tree()` helper + the "Rendering prototype tree"
section before the bootstrap.

### Why `plot.rpart()` was painful — three separate causes

1. **`rpart.plot::rpart.plot(fit)` fails outright** on a custom-method fit:

   ```
   Error: missing value where TRUE/FALSE needed
   ```

   The package inspects `fit$method` to decide how to format nodes and does not
   understand a user-supplied method object. **No argument fixes this — rpart.plot does
   not support custom methods.** Base `plot.rpart` + `text.rpart` *do* work; verified.
   If you reach for `rpart.plot` again, this is why it breaks.

2. **Split labels were unreadable.** By default rpart renders a categorical split with
   compact letter codes for the levels:

   ```
   labels(fit)              ->  race=bd
   labels(fit, pretty = 0)  ->  race=Black,Asian
   ```

   `b` and `d` are the 2nd and 4th *levels*. `pretty = 0` is the fix, and it matters far
   more now that R2 lets race split into genuine subsets rather than single levels.

3. **Leaf labels were all "1"** — that is R8 above.

### What changed

`plot.dm.tree(fit, file, title, subtitle, lbw.cols, alphavec, use.n, cex)` renders to
PDF using `plot()` + `text(pretty = 0)`, with `uniform = TRUE` (without it branch heights
follow deviance and the deep DM trees collapse into an unreadable smear).

Figures produced per model arm, before the bootstrap:

| file | what |
|---|---|
| `prototype_tree_depth3_<year>.pdf` | **the thesis figure** — depth-3, legible |
| `full_tree_<year>.pdf` | the full maxdepth-8 fit, appendix/reference |
| `decision_tree_depth_{2,3,4,5}_<year>.pdf` | the depth comparison (R10) |

`DISPLAY.DEPTH <- 3` is a **display choice, not a model choice**: it refits the same
criterion on the same observed counts purely to draw a readable figure. `dm.tree` — the
full maxdepth-8 fit — is untouched and is still what gets saved.

### Two side-fixes in the same area

**A stray `Rplots.pdf` on every batch run.** After the subgroup plot's `dev.off()` there
was a `par(mfrow = c(1, 1))`. Calling `par()` when **no device is open opens one**, so
that line silently created `Rplots.pdf` in the working directory every run. The settings
it was trying to reset had already died with the device on the previous line. Removed.
Verified: no stray file now.

**Your exploratory lines at the end of the file.** These four were running at top level
after the model loop:

```r
plot(dm.tree, main="DM Tree", uniform=TRUE)
text(dm.tree, use.n=TRUE, cex=1, font=3)
path.rpart(dm.tree, node = 2)     # left child  => mrace15 = 0
path.rpart(dm.tree, node = 3)     # right child => mrace15 = 1
```

They are **preserved as a commented recipe, not deleted** — but they no longer execute,
for three reasons: `plot()` with no open device is the other source of `Rplots.pdf`;
`path.rpart()` printed to the console *after* the "Completed" banner, which reads as if
the script is still working; and `dm.tree` at that point is whichever arm the loop
finished on (the **LBW-only** model), not the full model the comments assume. The
comments also reference `mrace15`, which no longer exists.

`path.rpart()` itself works fine on a custom-method fit and is genuinely useful for
reading off one node's rule — the commented block keeps it as a recipe, with
`pretty = 0` added.

---

<a name="r10"></a>
## R10 — Depth-controlled comparison reinstated

**Where:** `dm-cart.R`, "Depth-controlled comparison" section, before the bootstrap.

The Section 3.2.4 / Figures 3.4–3.5 analysis previously existed **only inside the
commented-out block** after the end of the model loop — there were zero active lines
there, so nothing produced these figures and nothing ever wrote
`paths$tree.depth.results()`. That path pointed at a file that was never created.

Now fits maxdepth ∈ {2,3,4,5} on the observed counts and records, per depth: variables
used, terminal-node count, first split, and the distribution of posterior P(LBW) across
classes. Each depth also gets a rendered tree.

One reason the old block could never have worked even if uncommented: it called
`predict(depth.tree, newdata = ..., type = "matrix")`, which does not work on a custom
method (R7). It now uses `dm.leaf.probs()`.

Saved to `tree_depth_comparison_<year>.RData` as `trees`, `tree.sum`,
`var.across.depths`, `pred.stats.compare` — the same four object names the 2020-2021
results on disk use, so anything written against those still works.

---

<a name="r11"></a>
## R11 — `depth_df` now comes from real bootstrap trees

**Where:** `dm-cart.R`, inside the bootstrap loop.

### What was wrong

`visual.R` plots per-variable split-depth distributions (Figs 3.11/3.12) from `depth_df`.
**Nothing ever wrote it.** So `visual.R` fell through to:

```r
depth_df <- data.frame(Variable = v, Depth = sample(0:6, 5000, replace = TRUE, prob = runif(7)))
```

Those figures were **random noise**, with no warning that it had happened.

### What changed

Every internal node of every bootstrap tree is now recorded:

```r
internal <- dm.tree.b$frame$var != "<leaf>"
node.ids <- as.numeric(rownames(dm.tree.b$frame))
depth.rows[[b]] <- data.frame(replicate = b,
                              Variable = as.character(dm.tree.b$frame$var[internal]),
                              Depth    = floor(log2(node.ids[internal])))
```

**The subtle part:** rpart names frame rows by **binary-heap node id** (1, 2, 3, 6, 7, …),
so `depth = floor(log2(id))` and the root is depth 0. It must be read from
`rownames(frame)`, **not** from row position. The old dead code used
`floor(log2(which(frame$var == "cig_0")))` — row position equals node id only for a
perfectly full tree, so that was wrong wherever it appeared. Do not "simplify" it back.

If no splits are recorded at all the script now **stops** rather than saving an empty
frame — the figures must never silently fall back to invented data again.

Verified on 60 replicates: 420 splits recorded, race at depth 0 in all 60.

---

<a name="r12"></a>
## R12 — High/low-risk subgroups are single classes, and `risk_summary` is saved

**Where:** `dm-cart.R`, subgroup section.

### Two separate problems

**(a) The subgroup was not a subgroup.** The old selector pinned 3 of 7 predictors:

```r
which(X.matrix$race == "Black" & X.matrix$dmar == 0 & X.matrix$cig_0 == 1)
```

That matches **32 of 480** classes, and the block then averaged over every combination of
the remaining four predictors — producing attenuated, blended probabilities rather than
the sharp contrast of Tables 3.2/3.3.

**(b) `risk_summary` was never saved.** `visual.R` looks for it by name and, not finding
it, substituted hardcoded constants — `c(rep(0.016, 10), 0.840)` and
`c(rep(0.009, 10), 0.910)`. **The headline subgroup figures were invented numbers.**

### What changed

Both profiles now fix **all seven** predictors, and `find.class()` asserts exactly one
match — so a future design change fails loudly instead of silently averaging again:

```r
high.risk.profile <- list(sex="1", dmar="0", mager="0", meduc="0",
                          precare1st="0", cig_0="1", race="Black")
low.risk.profile  <- list(sex="0", dmar="1", mager="0", meduc="1",
                          precare1st="1", cig_0="0", race="White")
```

> ⚠️ **These profiles are an assumption, not a confirmation.** The thesis's `i=69` /
> `i=28` are indices into the old 128-class binary design and do not survive the
> re-encoding, so they cannot be carried over mechanically. These are the profiles I
> proposed in `REVISION-PLAN.md` §5 Q2, still awaiting your sign-off. **To change them,
> edit these two lists — nothing else moves.**

`risk_summary` is saved with columns `Subgroup, Category, Mean, SE, Lower, Upper,
ClassIndex, Profile`. `Lower`/`Upper` are **percentile** bootstrap bounds. `visual.R`
reconstructs intervals as `Mean ± 1.96·SE`, which meant the same quantity had two
different intervals in one document; `SE` is retained for reference but the percentile
bounds are the ones to plot.

Verified: high-risk class 146, low-risk class 39, each exactly one class.

---

<a name="r13"></a>
## R13 — `top.var.df` restored, `B` configurable

**Where:** `dm-cart.R`, results assembly and the bootstrap header.

**`top.var.df`.** `top.vars.used[b]` — the root split variable of each bootstrap tree —
was computed in the loop and then thrown away. Table 3.1's "Initial Split Variable" row
had no source. The 2020-2021 results on disk *do* contain a `top.var.df`, so this was
lost in a rewrite rather than never existing. Now aggregated to `variable, count,
frequency` and saved.

**`B`.** Was `B <- 1000` with a comment reading "Adjust to 10000 for final runs" — the
kind of thing that gets forgotten. Now:

```r
B <- as.integer(Sys.getenv("DMCART_B", unset = if (interactive()) "200" else "10000"))
```

Interactive sessions get a fast default, batch runs get the publication value, and the
final run is a flag rather than an edit:

```bash
DMCART_B=200   Rscript dm-cart.R 2023-2024   # smoke test
DMCART_B=10000 Rscript dm-cart.R 2023-2024   # publication run
```

### What is saved where, after this pass

`bootstrap_tree_results_<year>.RData`:

```
bootstrap.tree.results  list: var.usage, var.freq.df, top.var.df, n.vars.summary,
                              risk_summary, depth_df, B, profiles
risk_summary            (also top-level)
depth_df                (also top-level)
top.var.df              (also top-level)
```

`risk_summary`, `depth_df` and `top.var.df` are saved **both** inside the list and as
top-level objects on purpose: `visual.R` loads this file and looks them up by bare name,
so list-only storage would leave its `exists("risk_summary")` check failing exactly as
before.

`tree_depth_comparison_<year>.RData`: `trees`, `tree.sum`, `var.across.depths`,
`pred.stats.compare`.

---

## Tier 2 verification

Ran the full revised model loop end to end on a synthetic 480-class / K=11 cohort built
with the new encoding (3-level `dmar`, 5-level `race`), with a planted non-contiguous
race signal:

```
--- Rendering prototype tree ---            prototype_tree_depth3_2024.pdf, full_tree_2024.pdf
--- Depth-controlled comparison (2-5) ---
  depth 2:  4 leaves, 3 vars (first split: race) P(LBW) mean 0.9517 [0.9274, 0.9721]
  depth 3:  8 leaves, 3 vars (first split: race) P(LBW) mean 0.9509 [0.9089, 0.9746]
  depth 4:  8 leaves, 3 vars (first split: race) ...
  depth 5:  8 leaves, 3 vars (first split: race) ...
High-risk class index: 146 | Low-risk class index: 39
P(LBW): high-risk 0.9744 [0.9659, 0.9830] | low-risk 0.9082 [0.8873, 0.9287]
Recorded 420 splits across 60 bootstrap trees (mean depth 1.43)
no stray Rplots.pdf
```

Root split rendered as `race=Black,Asian` — the planted non-contiguous partition,
recovered and displayed with real level names. Leaves show `k2 / P(LBW)=0.965 / n=32`.

⚠️ Those numbers are from **synthetic** data — they exercise the code paths, they are not
results. The planted signal is what makes race the root split there.

---
---

# Tier 2b pass — runtime and figures

---

<a name="r14"></a>
## R14 — Why the run got slower, and what was recovered

**Where:** `dm-cart.R`, `mysplit()`.

### Measured, on the real 480-class 2024 data

| what | time |
|---|---|
| one full tree fit, **original** splitter | 9.65 ms (median of 20) |
| one full tree fit, **revised** splitter | 12.94 ms (median of 20) → **1.34x** |
| one bootstrap replicate (draw + fit) | ~18 ms |
| whole `dm-cart.R`, **both** arms, `B = 200` | **7.25 s** |
| primary tree fit | 0.04 s |
| depth comparison, 4 fits | 0.03 s |

So the model fitting is not slow in absolute terms. Ranked causes of the increase:

1. **`B` default went from 1000 to 10000 for non-interactive runs — 10x.** This is by far
   the dominant term and it is the one to check first. `R14`'s own arithmetic:
   `B = 200 → ~4 s/arm`, `B = 1000 → ~18 s/arm`, `B = 10000 → ~3 min/arm`, i.e. ~6 min for
   both arms versus ~36 s under the old hardcoded `B = 1000`. If you ran
   `Rscript dm-cart.R ...` you got 10000. Set `DMCART_B` to choose deliberately.
2. **480 classes instead of 320** (R4 gave `dmar` a third level) — 1.5x more rows in
   every fit. Unavoidable consequence of keeping the records.
3. **The exhaustive categorical search itself** — race at `k = 5` evaluates 15 partitions
   where the old code evaluated 4. This is the correctness fix; it is what it costs.
4. Four extra depth-comparison fits and ~7 PDF renders per arm — together ~1 s. Noise.

If the slow step was actually `data.R` rather than `dm-cart.R`, that is `read.csv()` on a
2.1 GB file and is unrelated to any of the above.

### What was recovered

Two genuine inefficiencies I had introduced, both now fixed:

**(a) The winning partition was evaluated twice.** The search loop scored every
bipartition, then the final pass re-scored the k-1 contiguous splits from scratch. For a
**binary** predictor that meant doing double the original work — on 5 of the 7
predictors. Results are now stored **by bitmask** and the contiguous splits are looked
up, canonicalising each prefix mask against its complement (a partition and its
complement are the same split).

**(b) Every candidate re-subset the full n x K matrix.** Since a candidate split is just
a partition of the k *levels*, the only thing needed is each level's category totals.
`rowsum(y, group = match(x, ux))` computes that once, and each candidate then sums a
handful of rows of a `k x K` matrix instead of subsetting 480 rows.

Log-likelihood evaluations per call to the split function:

| levels | original | revised (before R14) | revised (after R14) |
|---|---|---|---|
| k = 2 (binary) | 3 | 5 | **3** — parity restored |
| k = 3 (`dmar`) | 5 | 11 | 7 |
| k = 5 (`race`) | 9 | 39 | 31 |

Binary predictors are back to exactly the original cost. The remaining excess is the
exhaustive search, which is the point of R2.

**Regression check — R14 is a pure speedup, not a behaviour change.** Against the tree
you had already fitted:

```
same split variables : TRUE
same node ids        : TRUE
same csplit          : TRUE
max |improve| diff   : 0
```

---

<a name="r15"></a>
## R15 — Tree figures rebuilt for legibility and branch-level interpretation

**Where:** `dm-cart.R`, `plot.dm.tree()` and its call sites.

### Does the multi-level factor splitting actually get used?

This was the question the old figure could not answer, so there is now a **split-partition
report** printed for the fitted tree. On your real 2024 full model:

```
 node depth   variable  left                        right                       n_births  n_levels  reachable_before
    1     0       race  Black                       White,Hispanic,Asian,Other   3476883         5             FALSE
    2     1       dmar  0                           1,Unknown                     455374         3              TRUE
   13     3       race  White                       Hispanic,Asian,Other         2064860         4              TRUE
   27     4       race  Asian                       Hispanic,Other                784400         3             FALSE
    ...
```

`reachable_before` is the precise test: could the **old** sorted-code search have found
this partition? It could only take a contiguous *prefix* of the levels present at that
node, in factor-code order. Two splits fail that test — **including the root**:

```
node 1    race   Black  |  White,Hispanic,Asian,Other
node 27   race   Asian  |  Hispanic,Other
```

Level order is `White, Black, Hispanic, Asian, Other`, so `{Black}` is not a prefix and
was **unreachable**. The R2 fix changed the root split of the tree. That is the direct
answer to "cannot tell whether this factor conversion was even used" — it is, and it
matters at the most important node.

Also visible: `dmar` splits as `0 | 1,Unknown`, i.e. the model groups the **unknown**
marital-status records with the *married* ones rather than the unmarried ones. That is an
empirical finding about the non-reporting states, and it is exactly what R4 was designed
to let the data answer rather than assume.

### Figure changes

`plot.dm.tree()` no longer uses `text.rpart` at all. `text.rpart` can only print the
condition for going **left**, at the node — so a 5-level race split rendered as a bare
`race=Black` with no way to see what went right. The skeleton is still drawn by
`plot.rpart` (which handles layout), then every label is placed manually at coordinates
from `rpart:::rpartco()`:

* **at each internal node** — the split variable name
* **on the left branch** — the levels sent left
* **on the right branch** — the levels sent right
* **at each leaf** — modal category, `LBW=`, `N=`, `cls=`

The left/right sets come from `labels(fit, pretty = 0, collapse = FALSE)`, which is
rpart's own decoding of `csplit` — so the printed groupings are what the model used, not
a reconstruction.

**`N` vs `cls` — this was actively misleading before.** rpart's `frame$n` counts the data
rows at a node, and here a data row is a **class**, not a birth. `use.n = TRUE` printed
`n=32` meaning 32 predictor classes, which reads as a sample size and is off by orders of
magnitude. The label now reports `N` = the node's actual births (sum of its category
counts, e.g. `N=2,064,860`) with the class count shown separately as `cls`.

Per your list:

| request | done |
|---|---|
| larger text | yes — see the sizing note below |
| report how many are in each leaf | `N=` births, plus `cls=` classes |
| keep the cohort label up top | kept, in the subtitle |
| drop `K = 11 categories` | removed |
| drop the `N terminal nodes \| split labels…` footer | removed |
| `(Full Tree, maxdepth = 8)` only on the preliminary view | only the full-tree figure carries it |
| clean up the `--` in titles | model label moved to the subtitle; no dashes |
| make the whole plot larger | canvas scales with leaves and depth |
| map branches to decisions | left/right level sets on the branches |

**Sizing, and why it is not just "make the canvas bigger".** LaTeX scales the figure to
`\textwidth`, so what determines legibility on the page is the ratio of type size to
canvas size, not the absolute canvas. `pointsize` therefore scales with `width`, holding
that ratio fixed. The consequence is that a wider canvas buys **no** extra room for leaf
text — horizontal crowding depends only on (leaf count) x (label width) — so `cex` falls
as `~11/n.leaves`, and above 10 leaves alternate leaf labels are **staggered
vertically** (with a dotted leader line) instead of being shrunk further. That is what
keeps the 22-leaf depth-8 tree readable on one page.

Figures produced per arm: `prototype_tree_depth3_<year>.pdf` (the main-text figure),
`full_tree_<year>.pdf` (appendix), `decision_tree_depth_{2,3,4,5}_<year>.pdf`.

---
---

# Tier 2c pass — legibility, and finding the risk profiles

## Where is `risk_summary`?

```
results/<cohort>/<model>/bootstrap_tree_results_<TARGET.YEAR>.RData
  e.g. results/2023-2024/fullmodel/bootstrap_tree_results_2024.RData
```

It is saved **twice on purpose**: as a top-level object (so `visual.R`'s
`exists("risk_summary")` check passes after a bare `load()`) and inside
`bootstrap.tree.results$risk_summary`. To read it without disturbing your session:

```r
e <- new.env()
load("results/2023-2024/fullmodel/bootstrap_tree_results_2024.RData", envir = e)
e$risk_summary     # per-category probabilities + percentile bounds
e$risk_totals      # one row per subgroup: aggregate risk, and the NBW mass
e$risk.ranking     # ALL 480 classes ranked, with profiles  (new in R17)
```

---

<a name="r16"></a>
## R16 — Type size solved from geometry, margins trimmed

**Where:** `dm-cart.R`, `draw.dm.tree()` / `plot.dm.tree()`.

Margins are cut to the title strip (`mar = c(0.3, 0.3, 4.0, 0.3)`) so the tree fills the
canvas, and the type is now sized from the figure's real geometry rather than a guess.

**The thing that is easy to get wrong here:** making the canvas bigger does not help.
`pointsize` is deliberately proportional to canvas width so the figure keeps its
appearance when LaTeX scales it to `\textwidth` — which means a wider PDF grows the type
by exactly the same factor and buys nothing. Legible labels come only from tree geometry,
staggering, or rotation.

So the limit is derived instead of guessed. The widest leaf label is ~9 characters; a
character is ~`0.6 * pointsize * cex / 72` inches; the room a leaf gets is the **smallest
gap between neighbouring leaves**, measured from the actual layout via `par("pin")`.
Solving gives `cex`, and the canvas width cancels out. Three regimes:

| leaves | treatment |
|---|---|
| ≤ 10 | labels side by side |
| 11–14 | alternate labels **staggered** downward with a dotted leader (≈2x the room) |
| > 14 | labels **rotated 90°**, so line height rather than string width has to fit |

Using leaf **count** instead of measured gaps was the earlier mistake: an rpart tree is
not balanced, so leaves are not evenly spaced and a count-based estimate badly overstates
the room available in a deep tree. Vertical padding is likewise solved rather than
guessed (the label block is a fixed number of inches while `ylim` is what we are
choosing, so `pady` appears on both sides of the equation), which is what stopped the
bottom row being clipped.

Title size is deliberately **not** tied to `cex` — `cex` tracks leaf crowding, so a dense
tree would otherwise get a tiny heading.

---

<a name="r17"></a>
## R17 — The extreme-risk profiles are now found from the data

**Where:** `dm-cart.R`, risk section.

### What changed

The high/low-risk profiles used to be **assumed** — someone's reading of which covariate
values ought to be worst and best. The model can answer that question directly, so it
now does. Every class is ranked by its posterior risk pooled over the bootstrap:

```r
lbw.by.class[i, b] = sum over the risk categories of Yhat[[b]][i, ]
```

`Yhat` is the leaf posterior each class inherits (R7), so this is the model's fitted risk
for that exact covariate combination, with the bootstrap supplying its sampling
distribution. Extremes are `argmax` / `argmin` of the posterior mean, with percentile
bounds.

### Four things this forced, each of which would have been a bug if skipped

**1. The risk metric has to differ by model arm.** In the LBW-only model there is no
Normal category, so `lbw.cols` spans *every* category and `sum(P[lbw.cols])` is
identically **1** for every class. The first run of this produced `P(LBW)=1.0000` for all
480 classes with the *same* class reported as both highest and lowest. That arm's
probabilities are conditional on being LBW, so its severity score is now the mass in the
three smallest deciles:

| arm | metric | meaning |
|---|---|---|
| Full (K = 11) | `P(LBW)` = Q1..Q10 | probability of low birth weight |
| LBW-only (K = 10) | $P(Q1-Q3 \mid LBW)$ | share of the three most severe LBW deciles |

**2. Ties are the normal case, not an edge case.** A tree gives every class in the same
terminal node the *same* posterior, so classes sharing a leaf have byte-identical risk.
There is no "the" highest-risk class — only a highest-risk **leaf** whose members the
model cannot distinguish. A bare `which.max()` would silently return whichever tied class
came first in the enumeration and present it as a finding. The code now identifies the
tied set, reports its size, and picks the member with the **most observed births** as the
representative (the profile describing the largest real population). On your 2024 full
model: 2 classes tie at the maximum, 3 at the minimum.

**3. Stability has to be tie-aware too.** Asking "is the chosen class the argmax in
replicate b" scored **14%** for the low-risk class purely because ties were broken
arbitrarily each replicate. The right question is whether the chosen class is *among* the
extremes, i.e. whether its value equals that replicate's max/min. Both extremes now score
**100%**.

**4. Low-support classes must be excluded.** A class with few observed births contributes
almost nothing to its leaf, so its "risk" is really its neighbours' risk plus the
Dirichlet prior. Letting such a class win the argmax would report an artefact as a
finding. `MIN.CLASS.BIRTHS <- 1000` gates eligibility; 187 of 480 classes qualify. The
threshold's effect is drawn, not hidden — ineligible classes appear hollow and grey in
the ranking figure.

### The result, on your real 2024 full model

```
HIGHEST  class 181  P(LBW)=0.2834 [0.2624, 0.2954]  N=1,710    top in 100% of replicates
   sex=0, dmar=0, mager=0, meduc=1, precare1st=1, cig_0=1, race=Black
LOWEST   class  40  P(LBW)=0.0538 [0.0520, 0.0547]  N=346,811  bottom in 100% of replicates
   sex=1, dmar=1, mager=0, meduc=1, precare1st=1, cig_0=0, race=White
```

**The assumed profiles were not the extremes**, which is exactly why this was worth doing:

```
a priori high (class 146): P(LBW)=0.2330 -> rank 13   [BELOW SUPPORT THRESHOLD: N=509]
a priori low  (class  39): P(LBW)=0.0648 -> rank 17 from the bottom
```

Two findings for the write-up. The assumed high-risk profile is **13th**, and it sits
below the support threshold (509 births) — it is a rare combination, not the modal
high-risk one. And the empirical high-risk profile has `meduc=1` and `precare1st=1`
(completed high school, first-trimester care) where the assumed one had both at 0: the
model puts `race`, `dmar` and `cig_0` ahead of them, and does not need the "everything
bad at once" profile to reach maximum risk.

Both definitions are kept in `risk_summary` under a `Definition` column
(`"empirical"` / `"a priori"`) so the comparison is in the saved output rather than
being something one has to reconstruct.

### New figure: `risk_profile_ranking_<year>.pdf`

Every class as a point at its posterior risk with its 95% interval, ordered by risk;
ineligible classes hollow; the two empirical extremes and the two a priori profiles
marked. This is the visual answer to "were our assumed profiles actually the extremes?"
and it shows the extremes in the context of the whole distribution rather than in
isolation.

---

<a name="r18"></a>
## R18 — Subgroup bar chart: LBW categories only

**Where:** `dm-cart.R`, subgroup plot.

The Normal bar sits at ~0.72–0.95 and dwarfed the ten LBW bars, flattening the
within-LBW shape — which is the entire reason the LBW region was expanded into deciles —
to nothing. NBW is dropped from the axes and reported in the subcaption instead:

```
Total P(LBW):  high-risk 0.2846 [0.2422, 0.3295]     low-risk 0.0539 [0.0513, 0.0559]
Normal (>2500 g) category omitted from the axes:
               high-risk 0.7154 [0.7016, 0.7340]     low-risk 0.9461 [0.9452, 0.9482]
```

Nothing is lost — the point estimate and bounds are still there — and the ten bars now
carry the comparison. For the LBW-only arm there is no NBW category and the subcaption
says so explicitly, which is what makes the two arms comparable side by side: both
figures then show ten bars on the same footing.

Both panels share a y-axis on purpose, so the ~5x gap between the subgroups is visible
rather than normalised away by per-panel scaling. Panel titles carry the class index and
N; the 7-term profile goes on its own line beneath (it does not fit on a title line at a
readable size and was running off the panel).

---

<a name="r19"></a>
## R19 — Two-pane tree tracing each extreme profile

**Where:** `dm-cart.R`, `risk_paths_tree_<year>.pdf`.

Answers "what do the high/low-risk groups look like *on* the tree". Each pane traces the
root-to-leaf path the class actually follows: path branches and labels in colour and bold,
everything else dimmed to light grey, and the destination leaf marked with a filled dot.

The path is recovered from the leaf's node id — rpart numbers frame rows by binary-heap
id, so repeatedly halving the leaf's id walks straight back to the root. The leaf comes
from `fit$where` (R7), so this is the same routing the model used, not a re-derivation
that could disagree with it.

This is what `draw.dm.tree()` was split out for: `plot.rpart` gives no way to colour
individual branches, so the skeleton is drawn manually from coordinates captured off a
throwaway device (`rpart:::rpartco()` cannot be called without an active `plot.rpart`).

**One honest caveat, printed on the figure itself:** the path is traced through the
**depth-3 display tree** (legible), while the quoted risk comes from the **maxdepth-8
bootstrap ensemble**. So the leaf's own `LBW=` value will not equal the headline number.
The figure says this outright rather than leaving a reader to trip over it.

---

## Tier 2c verification

Full run on your real 2024 data, both arms, `DMCART_B=150`: 0 errors, 0 warnings, no
stray `Rplots.pdf`, 9 figures per arm. Saved objects confirmed:

```
risk_summary   44 rows = 4 blocks x 11 categories
               cols: Subgroup, Definition, Category, Mean, SE, Lower, Upper,
                     ClassIndex, Births, Profile
risk_totals    4 rows: P_RISK + bounds, RiskMetric, P_NBW + bounds
risk.ranking   480 rows, 187 eligible
stability      high 100%, low 100%   (ties: 2 at max, 3 at min)
```

---

---
---

# Tier 3 pass — race expansion, extreme-ten analysis, and the last dead code

---

<a name="r20"></a>
## R20 — Race expanded to all seven `mracehisp` groups

**Where:** `data.R`, `RACE.SCHEME`.

`chapter2.tex:45` says of the binary flag: *"though this choice is entirely arbitrary.
(NEED TO REVISIT)"*. This is the revisit.

```r
RACE.SCHEME <- "full7"   # "binary" | "collapsed5" | "full7"
```

| scheme | groups | I (with 3-level `dmar`) |
|---|---|---|
| `binary` | Black / NonBlack | 2⁵·3·2 = **192** |
| `collapsed5` | White, Black, Hispanic, Asian, Other | 2⁵·3·5 = **480** |
| `full7` | White, Black, Hispanic, Asian, AIAN, NHOPI, Multiracial | 2⁵·3·7 = **672** ← default |

**Why the full seven rather than a collapse.** "Other" is the same arbitrariness the
thesis apologises for, one level up: it merges AIAN, NHOPI and multiracial mothers — who
have materially different LBW profiles — into a bucket that means nothing clinically. The
tail groups are genuinely sparse, and that is exactly the regime the Dirichlet prior
exists for. Under `full7` the paper **demonstrates** the sparsity claim it currently only
asserts.

Observed group sizes, 2024 (after the complete-case filter):

```
White 1,739,008 | Hispanic 944,512 | Black 455,381 | Asian 216,958
Multiracial 88,523 | AIAN 22,905 | NHOPI 9,620
```

18 of the 672 classes have zero observations (they carry prior mass only).

**This makes R2 load-bearing, not optional.** With 7 levels there are 2⁶−1 = 63
bipartitions and the old sorted-code search reached 6. The fitted tree now contains
partitions that are flatly unreachable without it — for instance

```
race:  White,Hispanic,AIAN   |   Asian,NHOPI,Multiracial
```

which is a genuine multi-level grouping on **both** sides, not a one-vs-rest peel.

---

<a name="r21"></a>
## R21 — `fread` instead of `read.csv`

**Where:** `data.R`.

The extract is ~2.1 GB and 237 columns; `read.csv` parsed all of it to keep 8 fields and
was by a wide margin the slowest step in the pipeline. `data.table::fread(..., select =
raw.cols)` reads only what is used.

```
data.R, full 2024 file:  ~minutes  ->  16.8 s
```

Falls back to `read.csv` if data.table is unavailable. `dmar` is explicitly normalised to
integer afterwards, because fread returns `""` for blank character fields and `NA` for
blank numerics and `DMAR.POLICY` needs a real `NA` either way.

---

<a name="r22"></a>
## R22 — Every leaf label rotated; long level sets wrapped

**Where:** `dm-cart.R`, `draw.dm.tree()`.

**Rotation is now unconditional** (previously only above 14 leaves). Two reasons: every
figure in the thesis reads the same way, with no figure needing the page turned or the
PDF zoomed; and rotation is the *cheaper* constraint — what has to fit between
neighbouring leaves becomes line height (4 stacked lines) instead of string width (~9
characters), which is several times less demanding and is what lets the type stay large
on the dense trees.

The dotted grey **leader line** is kept on every leaf. With rotated text the label sits
well below its node, so the leader is what ties the two together.

The `cex` cap came down from 1.35 to **1.15** — at the old cap the compact depth-2/3
trees were oversized for the page.

**Level-set wrapping.** With 7 race groups a branch label can read
`White,Hispanic,Asian,AIAN,NHOPI,Multiracial`, wide enough to run into neighbouring
subtrees. Labels of 3+ levels now break every two levels onto a new line. The full
membership stays visible — it *is* the split — at roughly half the horizontal footprint.

**Note on `full_tree_2024.pdf` as your prototype figure.** It works as the main-text
figure now, and it is the honest one to show: it is the model that was actually fitted
and bootstrapped, whereas `prototype_tree_depth3_2024.pdf` is a display-only refit. Two
things to keep in mind if you use it that way:

* the depth-3 file is still produced, and the **two-pane risk-path figure (R19) traces
  through the depth-3 tree**, since a 22-leaf tree with two dimmed copies side by side is
  not readable. If you want the panes on the full tree, pass `dm.tree` instead of
  `proto.tree` at that call and widen the canvas.
* `DISPLAY.DEPTH` still controls the depth-3 file only; the full tree follows
  `dm.control$maxdepth`.

---

<a name="r23"></a>
## R23 — The extreme *ten*, and what they have in common

**Where:** `dm-cart.R`; `risk_extremes_top10_<year>.csv`,
`risk_extremes_profiles_<year>.pdf`.

Reporting a single argmax invites the objection that it is a lucky cell — the most stark
pair rather than a real pattern. Ranking the top and bottom **ten** eligible classes and
reading **down the covariate columns** answers that directly: a covariate whose value is
constant across all ten is a genuine shared feature of the extreme; one that varies
freely is not doing the work, however it looks in a single profile.

`consensus.of()` computes exactly that — for each covariate, the modal level within the
set and the share of the ten carrying it. 10/10 is a unanimous marker; 5/10 on a binary
covariate is noise.

### The result (2024 full model, B = 200)

```
Shared features of the 10 HIGHEST-risk classes      Shared features of the 10 LOWEST-risk
 variable modal   n/10  unanimous                     variable modal   n/10  unanimous
 sex        0      7/10    FALSE                      sex        1     10/10    TRUE
 dmar       0     10/10    TRUE                       dmar       1      6/10    FALSE
 mager      1      6/10    FALSE                      mager      0      6/10    FALSE
 meduc      1      8/10    FALSE                      meduc      1      7/10    FALSE
 precare1st 0      6/10    FALSE                      precare1st 1      6/10    FALSE
 cig_0      1      6/10    FALSE                      cig_0      0     10/10    TRUE
 race     Black    8/10    FALSE                      race     White   10/10    TRUE
```

**Unanimous — highest:** `dmar=0`.
**Unanimous — lowest:** `sex=1`, `cig_0=0`, `race=White`.

This is a much better result than the single-pair contrast, and it is worth writing up
carefully:

* **Unmarried status is the only feature shared by all ten highest-risk profiles.** Race
  is 8/10 and smoking 6/10 — both strong, neither universal. The single argmax profile
  happens to be Black *and* smoking, which would have invited "you picked the starkest
  cell"; the ten-profile view shows `dmar` is the invariant and the others are not.
* **The low-risk set is more tightly determined than the high-risk set** — three
  unanimous markers versus one. Low risk requires a specific combination; high risk is
  reachable by several different routes. That asymmetry is a real finding and it is
  visible only because we looked at ten.
* Two of the ten highest-risk profiles are **White, not Black** (classes 67/68: over 33,
  unmarried, smoking). Worth saying explicitly — it resists a single-axis reading of the
  result.
* `mager`, `precare1st` and `meduc` never reach unanimity in either set. Do not over-read
  them from a single profile.

### Figure

`risk_extremes_profiles_<year>.pdf` — two panels, rows are profiles ordered by risk,
columns are the seven covariates, each cell printed with that profile's level and
**shaded when it matches the set's modal level**. A fully shaded column is a unanimous
marker and reads at a glance. Risk with its percentile interval is drawn to the right of
each grid, on a scale common to both panels.

---

<a name="r24"></a>
## R24 — `visual.R` runs, and every fabricated fallback is gone

**Where:** `visual.R`.

This file **had never run** in its current state. Four separate blockers, each of which
would have produced fiction if the one before it had not stopped execution:

1. **It did not parse.** Line 18 was `} else {»` — a stray `»`.
2. **`legend.position = c(0.18, 0.72)`.** ggplot2 ≥ 3.5 removed the numeric form; on the
   installed 4.0.0 it raises `Cannot add <ggproto> objects together`. Now set
   version-appropriately via a `legend_inside` term.
3. **`depth_df` was read from the wrong file.** It is written by dm-cart.R into the
   **bootstrap** results (it is per-replicate split data), not the depth-comparison file.
   The old code looked in `tree.depth.results`, found nothing, and fell through to
   `sample(0:6, 5000, replace = TRUE, prob = runif(7))` — random noise, unlabelled.
4. **`risk_summary` had a fabricating fallback** — `c(rep(0.016, 10), 0.840)` and
   `c(rep(0.009, 10), 0.910)` — which fired whenever the object was absent, i.e. always,
   since nothing wrote it before R12.

All fallbacks are **deleted, not commented out**, and replaced with a `stop()` naming the
command that produces the missing file. Results load into an environment rather than the
global, so a stale object from a previous `source()` cannot satisfy a lookup. Intervals
are the percentile bounds from `Lower`/`Upper`, not a `Mean ± 1.96·SE` reconstruction —
the old approach gave this figure a different interval from the one dm-cart.R reports for
the same quantity. Where `risk_summary` carries the R17 `Definition` column, the
**empirical** extremes are plotted.

The cohort and model arm now come from the command line instead of two hardcoded
constants that silently rendered only one of four combinations:

```bash
Rscript visual.R 2023-2024        # full model
Rscript visual.R 2023-2024 lbw    # LBW-only arm
```

Also fixed: a stale `align = "center"` argument to `geom_histogram`, and a
`par(mfrow = c(1, 1))` after `dev.off()` that created a stray `Rplots.pdf` on every run
(the same bug as R9).

---

<a name="r25"></a>
## R25 — Empty-class means are `NA`, not `0`

**Where:** `data.R`.

An unobserved class has no mean birth weight; `0 g` is not a measurement and would poison
any later summary that averaged the column. `x.x` / `x.y` are left `NA` and renamed to
`mean.dbwt` / `mean.log.dbwt` (only `count` is genuinely 0). Neither feeds the model.

---

## Tier 3 verification

Clean run of all four scripts against the real 2023/2024 files:

```
quantile.R  dropped 3,011 dbwt=9999 records; P(LBW)=0.0869, P(NBW)=0.9131
data.R      RACE.SCHEME='full7' (7 levels); 3,638,436 -> 3,476,907 (4.44% dropped)
            num_node: 672 (expected 2 x 3 x 2 x 2 x 2 x 2 x 7 = 672); 18 empty classes
dm-cart.R   0 errors, both arms
visual.R    0 errors, both arms
            no stray Rplots.pdf
```

---

<a name="r26"></a>
## R26 — One combined prior figure, with the quantile intervals in the legend

**Where:** the figure moved **out of `quantile.R` and into `visual.R`**.
Output: `results/<cohort>/<model>/plots/informed_prior_<prior.year>.pdf`.
Replaces the two `alphavec_plot_<year>.pdf` files and `visual.R`'s own
`dirichlet_prior_<year>.pdf`.

### Division of labour

`quantile.R` **computes and saves only** — it reads the full prior-year file with
`read.csv` exactly as before, derives the decile cutpoints and both arms' `alphavec`,
and writes them. No plotting. Running it does what it always did.

`visual.R` **loads both arms' saved priors** (`informed.prior()` for `with.2.5kg =
TRUE` and `FALSE`) plus the cutpoints, and draws one figure. Keeping the two apart means
the expensive step (parsing 2.1 GB) is not repeated every time the figure is retouched,
and the figure lives with the other figures.

### The figure

* **two panels** — Full Model (K = 11) and LBW-Only Model (K = 10);
* **percent above every bar in both panels**. Since α₀ = 1, α_k *is* the share of prior
  mass, so the label is simply `100 · α_k` — `0.87%`, `0.96%` … for the full model and
  `10.00%` throughout for the LBW-only arm;
* **the legend is the category definition**, built from the actual prior-year cutpoints
  (`Q1  0.23 - 1.19 kg` … `Q10  2.46 - 2.50 kg`), so the figure no longer depends on a
  lookup table elsewhere in the document;
* a sequential blue ramp across Q1→Q10 so the deciles read as one ordered family.

### Two things the figure deliberately does not do

**The Normal bar is not drawn.** Same reasoning as R18: at α ≈ 0.913 against ~0.009 per
decile it is ~105x any single LBW bar, and including it flattens all ten to a single
line — destroying the only structure the figure exists to show. Its share is stated in
the subtitle (`Normal (> 2.50 kg) holds the remaining 91.3%`), so nothing is lost. The
first draft did include it, and the full-model panel was ten invisible bars beside one
tall one.

**Panels use independent y scales.** The full model's deciles are scaled by P(LBW)
(~0.0087 each); the LBW-only prior is uniform at 1/10. The percent labels carry the
magnitudes, so the axes are free to show shape — and the shape is the finding: the
informed prior is **not** flat (0.76%–0.96% across deciles, because gram-valued birth
weights cannot split into exactly equal deciles), while the LBW-only prior is uniform by
construction.

### What this means for the thesis

`chapter3.tex:22–29` currently uses a `\subfloat` pair:

```latex
\begin{figure}[H]
    \centering
        \subfloat[\centering Full model informed prior]{{\includegraphics[width=6cm]{.../alphavec_plot_2020_full.pdf} }}
    \qquad
        \subfloat[\centering LBW-only model informed prior]{{\includegraphics[width=6cm]{.../alphavec_plot_2020_small.pdf} }}
    \caption[...]{Comparison of informed Dirichlet priors based on 2020 quantiles}
    \label{fig:alphavec}
\end{figure}
```

That collapses to a single graphic:

```latex
\begin{figure}[H]
    \centering
    \includegraphics[width=\textwidth]{amkurth-Thesis/chapter3/figures/informed_prior_2023.pdf}
    \caption[Informed Dirichlet Priors from 2023 Quantiles]{Informed Dirichlet prior
      hyperparameters for both model arms. Bar labels give each category's share of the
      total prior mass; the legend gives each category's 2023 low-birth-weight decile
      interval. The Normal category is omitted from the axes and reported in the
      subtitle; panels use independent $y$ scales.}
    \label{fig:alphavec}
\end{figure}
```

The surrounding sentence says "From the 2020 proportions…" and needs to say **2023**.

### Housekeeping

No longer written, but still on disk from earlier runs:

```
results/2023-2024/fullmodel/plots/alphavec_plot_2023.pdf     <- delete (pre-R3 prior)
results/2023-2024/lbwmodel/plots/alphavec_plot_2023.pdf      <- delete (pre-R3 prior)
results/2023-2024/*/plots/dirichlet_prior_2023.pdf           <- delete (superseded)
results/2020-2021/lbwmodel/plots/alphavec_plot_2023.pdf      <- delete (mislabelled year)
```

Keep `results/2020-2021/*/plots/alphavec_plot_2020_{full,small}.pdf` — the current thesis
text still references those.

---

<a name="r27"></a>
## R27 — Figure count halved; tree renderers moved to `util.R`

**Where:** renderers moved `dm-cart.R` → `util.R`; figures now drawn in `visual.R`.

### Division of labour, extended

`dm.tree.coords()`, `draw.dm.tree()` and `plot.dm.tree()` now live in `util.R`, which
all four scripts already source. `dm-cart.R` fits and saves; `visual.R` draws. Nothing is
refitted to redraw a figure — the depth panel is built from the `trees` list in
`tree_depth_comparison_<year>.RData` and the full tree from `dm_tree_rebin_<year>.RData`.

### Tree figures per model arm: 6 → 2

| was | now |
|---|---|
| `prototype_tree_depth3_<year>.pdf` | *(no longer written; `proto.tree` is still fitted, because the risk-path figure traces through it)* |
| `decision_tree_depth_{2,3,4,5}_<year>.pdf` | `depth_comparison_<year>.pdf` — all four as small multiples |
| `full_tree_<year>.pdf` | `full_tree_<year>.pdf` — unchanged, the annotated fitted model |

Across both arms that is **12 tree PDFs down to 4**.

### What made the small multiples legible

Three changes, because shrinking the existing figure four-fold does not work:

* **`leaf.detail = "none"`.** At panel size the readable signal is the shape of the tree
  and which predictors enter; per-leaf numbers belong in `pred.stats.compare`, not
  squeezed into a 3-inch panel.
* **`branch.labels = "left"`.** The right child is the complement of the left, so
  printing both doubles the ink for no information.
* **`var.labels`** — short display names (`precare1st` → `care`, `cig_0` → `smoke`,
  `mager` → `age`, `meduc` → `educ`, `dmar` → `marr`), expanded once in the figure
  footnote. This was the change that actually fixed it: with no leaf text the binding
  constraint stops being leaf-label height and becomes **variable-name width at each
  node**, so `cex` is now sized from the longest name against the smallest horizontal gap.
  At depth 5 the previous sizing produced `precare1stprecare1st` overprinted; abbreviating
  buys far more room than shrinking the type does.

`title.cex` / `subtitle.cex` are arguments now, since a title sized for a standalone
figure dominates a small multiple.

### The prior figure, restacked

Panels are **top/bottom** rather than side by side: stacking shares one x axis, so the ten
deciles line up vertically and the two priors read against each other directly. The legend
moved underneath in **3 rows x 4 columns** — a wide, short legend under a tall figure
wastes far less space than a tall legend beside a wide one — with entries trimmed to
`Q1  0.23-1.19` and the unit stated once in the legend title. Figure is 7.0 x 6.6 in,
sized for `\textwidth` in a two-column article.

### Diagnostic vs. publishable

Keep as figures: `informed_prior`, `full_tree`, `depth_comparison`,
`risk_profile_ranking`, `risk_extremes_profiles`, `risk_paths_tree`,
`high_low_risk_pred`. Supplementary at most: `depth_distributions` (bootstrap split
depths — a diagnostic of ensemble behaviour, not a result).

### Housekeeping

`decision_tree_depth_{2,3,4,5}_<year>.pdf` and `prototype_tree_depth3_<year>.pdf` are no
longer written but remain on disk from earlier runs. Delete them so a stale file cannot
be picked up by `\includegraphics`.

**One gotcha worth recording:** `sips` renders `cairo_pdf` output incorrectly — it showed
the depth panel as a single tree on a blank page when the PDF was in fact correct. Use
`pdftoppm` (or open the PDF directly) to check these figures.

---

<a name="r28"></a>
## R28 - NBW on the prior figure; risk landscape as one two-pane figure

**Where:** `visual.R` (both), `dm-cart.R` (landscape removed).

**Prior figure.** The NBW category is drawn on the full-model panel again. It was
previously omitted because at ~91% of the prior mass it compresses the ten LBW deciles
to near-zero height on a linear axis - which it still does. It is included because the
mass split *is* the substantive content of the informed prior, and because every bar
carries its percentage, so the deciles stay readable as numbers even when they are not
readable as heights. Implementation keys off `length(alphavec)`, so the full model gets
11 categories and the LBW-only model 10, leaving an empty NBW slot in the lower panel -
the honest depiction of K = 10 against K = 11. The red-to-green palette and the subtitle
are yours; only the data and the eleventh legend entry were added here.

**Risk landscape.** Moved out of `dm-cart.R`, where it was drawn once per arm, into
`visual.R`, which loads both arms and draws them as one two-pane figure. Classes below
the support threshold are drawn hollow and grey rather than dropped, so the reader can
see how much of the class space the filter removes (191 of 672 eligible in the full
model, 65 of 672 in the LBW-only arm) and that the extremes are not simply the thinnest
cells. Each pane is labelled with its own risk metric, since the two arms use different
ones.

**LaTeX.** `tex/article/figures_snippet.tex` holds ready-to-paste figure blocks with
paths relative to `main.tex`. Every wide figure is `figure*` (spanning both columns): at
a single column width (~3.4 in) the tree and ranking panels are unreadable. Only
`high_low_risk_pred` is a single-column `figure`. Verified: compiles against the real
preamble with the real PDFs, 0 errors, 0 undefined references, 0 missing files.

---

<a name="r29"></a>
## R29 - An article/ output folder, and one two-model tree figure

**Where:** `util.R` (`paths$article`), `visual.R`, `dm-cart.R`.

### Output split

Figures were landing in `results/<cohort>/<arm>/plots/` alongside every diagnostic, so
choosing what to publish meant picking files out of a pile. There are now two
destinations:

```
results/<cohort>/article/          only what goes in the manuscript
results/<cohort>/<arm>/plots/      diagnostics
```

Copying figures into the paper is now one directory copy, and no diagnostic can be picked
up by `\includegraphics` by accident. Per-arm article figures carry the arm in the
filename (`risk_extremes_profiles_full_2024.pdf`) so the two arms cannot collide.

**article/** - `informed_prior`, `model_trees`, `risk_profile_ranking`,
`risk_extremes_profiles_{full,lbw}`, `high_low_risk_pred_{full,lbw}`.
**plots/** - `depth_comparison`, `depth_distributions`, `full_tree`, `risk_paths_tree`.

### The depth panel is no longer an article figure

Four trees on one page leaves each about 3 in wide, and the deepest has 19 leaves; at
article scale it was unreadable. It is replaced by the comparison the paper is actually
about - the fitted full model above the fitted LBW-only model, each with the full page
width - and the depth sweep is reported as a table (`pred.stats.compare`) instead. The
panel is still written to `plots/` as a diagnostic, on a larger 11 x 10 in canvas.

### What made the stacked figure legible

* **`leaf.detail = "compact"`** - the four-line leaf label (category, P, N, classes) is
  what made leaves illegible at this size. Leaves now carry the posterior only; N and the
  class count belong in a table.
* **`branch.labels = "left"`** - with seven race groups the complement runs to
  `Hispanic,Asian,AIAN,NHOPI,Multiracial`, which collided with neighbouring node labels
  in the dense middle of the full tree. Printing only the left set halves the text and
  loses nothing.
* **Padding fix** - the vertical reserve for rotated leaf text was still sized for the
  9.6-character full label, leaving a third of each panel empty under a compact tree.
  Now 5.6 characters for `compact`.

### A degeneracy this figure exposed

The first stacked draft showed **1.000 at every leaf of the LBW-only tree**. Same cause as
the ranking metric in R17: in that arm every category is an LBW category, so summing
columns 1:10 is identically 1. The leaf quantity is now per-arm - `P(LBW)` for the full
model, `P(Q1-Q3 | LBW)` for the LBW-only model - and that panel now reads 0.231-0.359.
Worth noting in the paper that the two panels are on different scales for this reason.

---

<a name="r30"></a>
## R30 - Horizontal leaf labels, per-element sizing, larger ranking figure

**Where:** `util.R` (`draw.dm.tree`), `visual.R`.

### "Horizontal AND bigger" was a real conflict, resolved by staggering

Rotated, a leaf label is constrained by LINE HEIGHT against the leaf gap. Horizontal, it
is constrained by STRING WIDTH, which is far tighter: at 23 leaves across a `\textwidth`
page a 5-character number caps out around 7 pt. Alternating labels onto two rows doubles
the room, which is what allows horizontal text at roughly twice the size a single row
would permit. Compact labels are now horizontal, staggered, each with the dotted leader
tying it back to its leaf.

### Three sizing bugs found while doing it

1. **A single global `cex` is set by the worst-crowded spot** and drags every well-spaced
   label down with it. Labels are now sized per element against their OWN local gap, so
   the two cramped `precare1st` siblings no longer shrink the whole figure.
2. **The leaf's rival for space is the leaf two positions away in X ORDER**, not the
   nearest leaf at the same depth. The first attempt grouped by depth, which in an
   unbalanced tree is a different leaf entirely - so the crowded clusters kept
   overlapping even though the arithmetic said they fit.
3. **The stagger offset was a fraction of `y.gap`**, i.e. of (y range)/depth. On a
   depth-8 tree that is well under one line of text, so the two rows still collided
   diagonally. The offset is now computed as a real line height converted into user-y
   units, which is correct at any depth.

Type was also enlarged relative to the canvas (pointsize 11 -> 15 on the tree figure),
which is the lever that matters once `\includegraphics` scales to `\textwidth`.

### Ranking figure

Enlarged from 7 x 4.4 at pointsize 9 to 9 x 5.4 at pointsize 12, with larger points,
axis text, and a y-label moved out to its own `mtext` so it is not compressed by the
axis. Same content.

---

<a name="r31"></a>
## R31 - Tree label placement, a redrawn risk-category figure, and caption corrections

**Where:** `util.R` (`draw.dm.tree`), `visual.R`, `tex/article/figures_snippet.tex`.

### The dominant problem was the LaTeX, not the figure

`main.tex` had the tree at `\includegraphics[width=0.5\textwidth]` inside a plain
`figure`. In a two-column layout that is about 3.4 in, so a figure drawn 11.5 in wide was
being scaled to roughly a third of its natural size. No amount of internal sizing
survives that. The trees and the ranking must be `figure*` at `\textwidth`.

### Label placement

* Variable names raised from `0.075` to `0.155` of `y.gap` above the node.
* Branch labels moved to `0.42` along the branch and, more importantly, **anchored into
  the V formed by the two branches instead of out of it**. A left branch runs down-left,
  so the empty wedge is to its right and the text must extend right. Anchoring the other
  way pushed every label outward across the neighbouring subtree - that was what put
  `Black,NHOPI` on top of `mager`.
* **Binary splits are no longer labelled.** The blue `0` on every binary branch was about
  fifteen pieces of clutter in the dense middle of the full tree, conveying something the
  convention already fixes. Only multi-level splits print a level set; the convention
  moved to the caption.

### Risk-category figure redrawn

Rebuilt in ggplot to match the other article figures: four panels (both arms x both
extremes), stronger ColorBrewer RdBu endpoints (`#B2182B` / `#2166AC`, which stay
distinguishable in greyscale), no in-figure footnotes, profile identity on the panel
strip. The base-graphics version dm-cart.R writes stays in `plots/` as a diagnostic.

### Caption corrections

The captions in `main.tex` asserted several things the code does not do. Corrected in
`figures_snippet.tex`:

| Claim | Reality |
|---|---|
| "pruned to the subtree that maximizes cross-validated log-likelihood on a held-out validation set" (twice) | No pruning and no validation set: `cp = 0`, `xval = 0`, `maxdepth = 8`. |
| "with 95% percentile intervals on the expected category at each terminal node" (tree figure) | The tree figure carries no intervals; leaves give a posterior point estimate. |
| "Full model tree" | The figure shows both arms. |
| ranking: "shows the 20 highest- and lowest-risk profiles" | It shows all 672 classes. |
| ranking: "risk profile is the expected category ... a weighted average of the category indices" | It is `P(LBW)` (full) / `P(Q1-Q3 | LBW)` (LBW-only). |
| ranking: "full model uncovers socioeconomic ... LBW-only reveals biologically-driven" | Not shown by that figure; belongs in the text, and still needs checking against the fitted trees. |
| prior: "2024 Quantiles"; "Normal category is omitted from the axes" | Quantiles are from 2023, and NBW is now plotted (R28). |

---

<a name="r32"></a>
## R32 - Half-page trees, depth-growth figure, and a commonality table

**Where:** `visual.R`.

**Tree figure resized.** Now that it is a `figure*` at `\textwidth` it no longer has to
survive a 3x down-scale, so pointsize comes back from 16 to 12 and the canvas is reshaped
from a near-square 11.5 x 10.5 to 12 x 8. At `\textwidth` that is about 4.7 in tall -
half a page rather than a full one.

**Depth growth promoted to an article figure.** Written to `article/` as
`depth_growth_{full,lbw}_<year>.pdf`, on the same 12 x 8 canvas at the same pointsize
with the same title and footnote treatment as the fitted-tree figure, so the two read as
a pair. Its story is legible at this size because the panels carry no terminal-node
values: three predictors at depth 2, five at depth 3, six from depth 4 - and `meduc`
never entering at all.

**Commonality table.** `article/table_commonality_<year>.tex`, generated from the saved
`top.consensus` / `bottom.consensus` objects for both arms and `\input` directly, so it
cannot drift from the run that produced it. A 7-row table states "how many of the ten
carry this level" far more compactly than the shaded grid, and fits one column. The grid
figure is kept as an optional alternative, commented out in the snippet - use one or the
other, not both.

**Escaping note.** The first version of the table generator emitted `\\begin{table}`
throughout: writing LaTeX commands as R string literals needs four backslashes for
`\begin` and eight for a row terminator, and it was over-escaped by one level. The
generator now writes `~BS~` and `~RE~` placeholders and substitutes once at the end.

### Caption correction

The caption I had attached to `risk_extremes_profiles_full_2024.pdf` was paired with the
wrong figure in `main.tex`. To be explicit about which is which:

| File | What it shows |
|---|---|
| `risk_category_profiles_<year>.pdf` | the actual per-decile **probabilities with error bars**, four panels (both arms x both extremes) |
| `risk_extremes_profiles_full_<year>.pdf` | the **shaded profile grid** - the ten extreme classes and which covariate levels they share |
| `table_commonality_<year>.tex` | the same shared-feature content as a table |

---

<a name="r33"></a>
## R33 - Column-width table, three-row leaf stagger, single-column prior

**Where:** `visual.R`, `util.R`.

### Commonality table now fits a column

The six-column layout (Arm, Predictor, then Level and n for each extreme) was far wider
than the ~3.4 in of a two-column page and was running into the neighbouring table. It is
now **three columns**: the arm heads a spanning subrow, and level and count share one cell
("0 (10/10)"). Set in `\footnotesize`. Compiles with zero overfull boxes.

### Leaf labels stagger over up to three rows

Two rows give each label twice the leaf gap; where the tree is dense that was still not
enough and the labels had to shrink instead. The row count is now chosen from the
measured gap - one, two or three - so crowded clusters spread further and keep their type
size. That is what the extra height buys: **vertical space spent on horizontal
separation**. Canvas 12 x 8 -> 12 x 10.5 to hold the third row.

Multi-level branch labels also wrap at **one level per line** once a set reaches three
members. Two per line still left "White,Hispanic," wide enough to hit the neighbouring
node label; width, not height, is the binding constraint on those.

### Prior figure sized for one column

Title and subtitle removed from the figure and moved into the LaTeX caption, percentages
rotated to vertical (at 3.5 in the eleven bars are ~0.25 in apart and a horizontal
"0.87%" does not fit), legend in two columns at 6.2 pt, canvas 3.5 x 6.2 in. Included
with `width=\columnwidth` in a plain `figure`.

Note the consequence of keeping NBW on a linear axis at this width: the ten LBW bars in
the top panel are effectively invisible as bars and are read entirely from their
percentage labels. That is the trade-off of showing the mass split, and the caption says
so.

---

<a name="r34"></a>
## R34 - Branch labels de-collided with leader lines; layout

**Where:** `util.R` (`draw.dm.tree`), `visual.R`.

Branch labels are now drawn LAST, after node names and leaf values, and checked against a
register of everything already on the page. Where one would land on another label it is
stepped outward along the branch normal until clear, and a dotted leader is drawn back to
the branch so the association survives the move. Two details mattered:

* **Boxes are padded before testing.** Two labels that merely abut still read as collided;
  an untested "just touching" case is what left `White` sitting against `race`.
* **The search starts at step 1, never 0.** A label anchored exactly on the branch has the
  branch line drawn through it, which is what struck through `Black,NHOPI`.

Prior figure: legend moved inside, upper-left, where the full-model panel is empty
anyway - every LBW decile sits at ~0.009 against an NBW bar of 0.913 - so it costs no
data ink and returns the figure's height to the panels, which is what paid for the larger
type. Risk landscape restacked top/bottom at 3.6 x 6.4 in for a single column.

**Caption claim verified.** "Left is 0, right is 1" is correct and is guaranteed by
construction, not luck: for a binary predictor `n.mask = 2^(k-1) - 1 = 1`, so the search
evaluates exactly one partition and it always places the FIRST level left. Checked on
both fitted trees - all 20 binary splits have left = `0`.

---

<a name="r35"></a>
## R35 - The first birth-weight bin was open on the left

**Where:** `data.R`.

Every bin used `dbwt > cutpoint[i] & dbwt <= cutpoint[i+1]`, including the first. The
first cutpoint is the MINIMUM of the prior year's LBW distribution (227 g in 2023), so any
2024 birth at or below 227 g fell into no bin at all and vanished from the counts.

Found by reconciling `sum(counts.df)` against the analysed record count: 3,476,883 versus
3,476,907, a gap of 24. The cause is 29 births at exactly 227 g in 2024, of which 24
survive the complete-case filter.

Numerically negligible; conceptually not, since these are the very smallest births in the
file, which is precisely what the analysis is about. And because the cutpoints come from a
DIFFERENT year, 2024 could contain a birth lighter than anything seen in 2023, which would
be dropped the same way. The first bin now takes everything up to `cutpoint[2]` with no
lower bound.

A reconciliation check now runs after binning and stops the script if the counts do not
account for every analysed birth (all of them for the full model, those at or below 2.5 kg
for the LBW-only arm). Both arms now reconcile exactly: 3,476,907 and 295,352.

---

## Open items

* **`visual.R` does not parse.** Line 18 is `} else {»` — a stray `»`. Everything
  downstream of it is moot until that is fixed. It also still contains the live
  `risk_summary` fallback that fabricates the headline numbers
  (`c(rep(0.016, 10), 0.840)`).
Nothing from the original audit remains open. Resolved across the passes:

* Tier 0/1 — double-dipping (R1), the split search (R2), sentinels (R3), `dmar` (R4),
  `meduc` (R5), `precare1st` (R6), the prediction path (R7).
* Tier 2 — labels (R8), plotting (R9), depth comparison (R10), `depth_df` (R11),
  single-class subgroups (R12), `top.var.df` / `B` (R13), runtime (R14), figures (R15).
* Tier 2b/c — legibility (R16), empirical risk profiles (R17), the NBW bar (R18),
  the risk-path figure (R19).
* Tier 3 — race expansion (R20), `fread` (R21), rotated labels (R22), the extreme ten
  (R23), `visual.R` (R24), empty-class means (R25).

Deliberately **not** done, with reasons:

* `myinit()` still reads `alphavec` from its enclosing scope rather than through `parms`.
  It works, and rewiring the method contract carries more risk than the tidiness is
  worth. Revisit only if the method is ever moved into a package.
* The `"binary"` and `"collapsed5"` race arms are implemented but have not been run; they
  exist for the sensitivity analysis when you want it.
* `DMCART_B=10000` has not been run — every number quoted here is from a B = 150–200
  smoke run. See §7 of the write-up guide.

---

## How to re-run

Order matters — each step consumes the previous step's output.

```bash
Rscript quantile.R 2023-2024   # 2023 -> cutpoints + alphavec  (R3 changes these)
Rscript data.R     2023-2024   # 2024 counts, binned on 2023 cutpoints (R1/R3/R4/R5/R6)
Rscript dm-cart.R  2023-2024   # trees + bootstrap (R2/R7)
```

Delete the stale artifacts listed in [R1](#r1) first, or you will have `_2023`- and
`_2024`-suffixed files side by side.

Expect these to differ from anything you have on disk: the class count is now **480**,
`meduc` flips ~60% of records, ~20% more records survive the complete-case filter, and
the target cohort is 2024 rather than 2023 for the first time.

---
---

# Write-up guide — what to say, where, and how to restructure

This section is for drafting, not for the code. It is organised by where the change
lands in the document. Your own interpretation lives in `notes.md`; this is the
checklist of what must not be left out, and where the argument needs restructuring
rather than a patch.

---

## 1. The three claims the thesis currently makes that are no longer true

Fix these first — everything else is additive, these are corrections.

| Claim in the current text | Status | Where |
|---|---|---|
| I = 128 classes, 7 binary predictors, **X** ∈ {0,1}^{I×7} | **wrong now**: I = **672**, and two predictors are not binary | ch.2 §preprocessing, §2.4 notation, Eq. setup, and every "128" in ch.2–3 |
| "adequate prenatal care" | **overclaim**: measures *when care began* | Table 2.1 |
| "HS completed" via `meduc == 3` | **wrong**: excluded everyone with more education | Table 2.1 |

The `I = 128` correction is not a find-and-replace. The notation block in §2.4 defines
`x_i ∈ {0,1}^7`; with a 3-level `dmar` and a 7-level `race` that is now a product of
mixed-arity factors. Restate it as **x**ᵢ ∈ ∏ⱼ 𝓛ⱼ with |𝓛| = (2,3,2,2,2,2,7), and say
I = ∏ⱼ|𝓛ⱼ| = 672. The DM likelihood, the ΔL rule and the factorisation over classes are
all unchanged — only the index set grows. Say that explicitly, because a referee will
want to know whether the derivation still holds.

---

## 2. Chapter 2 — restructuring, not patching

### 2.1 The preprocessing section needs a missingness subsection it does not have

Currently preprocessing goes straight from "3.6 M records" to "128 classes" with no
account of what was dropped. Add a short subsection with the R3 table and these points:

* NCHS encodes "not stated" as in-range integers, so `na.omit()` never saw them. Name the
  five codes and say plainly that `dbwt == 9999` was being counted as a **normal**
  birth weight — it inflated the NBW category, which is the denominator of the whole
  LBW story.
* Final accounting: **3,638,436 → 3,476,907 records (4.44 % dropped)**.
* Put the full missingness table in the appendix, reference it here.

### 2.2 `dmar` deserves its own paragraph, and it is a methods contribution

This is the most defensible-but-unusual choice in the paper, so argue it rather than
mention it:

1. ~11 % of 2024 records have blank `dmar`, and the blanks are **not random** — they
   co-occur exactly with blank `mar_p`, i.e. they are a reporting-area (state-level)
   pattern, and the file is ordered by state.
2. Imputation by region is **impossible from this file**: every geographic identifier in
   the public-use extract is suppressed (`ocntyfips`, `ocntypop`, `rcnty`, `rcnty_pop`,
   `rcity_pop`, `mrterr`, `octerr`, `mbcntry` — blank in 100 % of records). Say this
   outright; it pre-empts the obvious referee question.
3. Worse, the missingness mechanism **is** the unobserved geography, so MAR given the
   observed covariates is not credible either. Imputing would manufacture precision.
4. So `dmar` carries an explicit `Unknown` level. This assumes nothing, keeps every
   record, and lets the model report where the unknown group belongs.
5. **And it did.** The fitted tree splits `dmar` as `0 | 1,Unknown` — the unknown records
   group with the **married**. That is an empirical finding about the non-reporting
   states, not an assumption, and it is only available because of the encoding choice.
   This is the payoff sentence for the whole subsection.

State the cost honestly: I goes from 2⁶·5 to 2⁵·3·7.

### 2.3 Race — rewrite the paragraph, do not amend it

The current text says the binary split is *"entirely arbitrary. (NEED TO REVISIT)"* and
justifies it by national population shares. Replace the whole justification with:

* the seven `mracehisp` groups, with the 2024 counts (White 1,739,008 … NHOPI 9,620);
* why not collapse to "Other": it merges AIAN, NHOPI and multiracial mothers, who have
  different LBW profiles, into a bucket with no clinical meaning — the same arbitrariness
  one level up;
* **the methodological argument, which is the point**: sparse cells are exactly what the
  informed Dirichlet prior is for. It shrinks thin classes toward the marginal profile and
  the parametric bootstrap reports honest interval width on them. The thesis currently
  *asserts* that the DM framework handles sparsity; with `full7` it **demonstrates** it.
  This converts a stated advantage into a shown one and is the strongest single argument
  for publication.
* note that `binary` and `collapsed5` remain available as one-line sensitivity arms.

### 2.4 Table 2.1 has to be rebuilt

Not edited — rebuilt. It currently presents seven binary rows. It now needs an arity
column and three corrected rows:

| predictor | levels | definition |
|---|---|---|
| `sex` | 2 | male / female |
| `dmar` | **3** | married / unmarried / **unknown** |
| `mager` | 2 | > 33 / otherwise |
| `meduc` | 2 | **HS completed or more (`>= 3`)** / otherwise |
| `precare1st` | 2 | **prenatal care began in 1st trimester** / otherwise |
| `cig_0` | 2 | any pre-pregnancy smoking / none |
| `race` | **7** | White, Black, Hispanic, Asian, AIAN, NHOPI, Multiracial |

Footnote the `meduc` correction — the old rule labelled a mother with a bachelor's degree
as not having completed high school, and it moves the flagged share from **28.5 % to
88.9 %**. Any comparison to the thesis's published numbers has to acknowledge it.

---

## 3. Chapter 2 §2.4.4 — the splitting rule needs a new paragraph

This is a genuine methods contribution and currently has no text at all.

`rpart`'s user-split protocol accepts an *ordering* of levels and evaluates only splits
contiguous in it. For a k-level factor that reaches k−1 of the 2^{k−1}−1 possible
bipartitions, and *which* k−1 depends on nothing but the order the levels were declared
in. For binary predictors this is complete (k=2 → one split), which is why it was
invisible in the original work. For race it is not.

State: the implementation enumerates all 2^{k−1}−1 bipartitions under the same ΔL
criterion, then returns the ordering with the optimal left-group first, so rpart can
express the optimum as a contiguous cut. Add the demonstration — on a constructed case
the sorted-code search finds goodness **35.86** where the true optimum is **174.95**.

**And report that it changed the fitted model, not just the search.** In the 2024 tree,
two splits are unreachable by the old procedure, including the **root**:

```
node 1   race   Black                |  White,Hispanic,Asian,AIAN,NHOPI,Multiracial
node 27  race   White,Hispanic,AIAN  |  Asian,NHOPI,Multiracial
```

Level order is White, Black, Hispanic, Asian, …, so `{Black}` is not a prefix of it. The
second is a multi-level grouping on both sides — not a one-vs-rest peel. This is the
concrete evidence that the fix matters; put the split-partition table in the appendix.

---

## 4. Chapter 3 — results

### 4.1 The cohort sentence must change

Every "2020/2021" becomes "2023/2024", and add the sentence that was silently false
before: cutpoints and **α** come from 2023, counts from 2024, and the two are now
genuinely different years. (Under the old code they were not — see R1. You do not need
to confess the bug in the paper, but do not repeat the claim without the fix in place.)

### 4.2 Replace the assumed high/low-risk subgroups entirely

The old §3.3.5 defines two subgroups by reading the covariates. Restructure to:

1. define the risk score — `P(LBW)` for the full model, `P(Q1–Q3 | LBW)` for the LBW-only
   arm (say why: with no NBW category the former is identically 1);
2. rank **all** classes by posterior mean over the bootstrap;
3. impose the support filter `N >= 1,000` and say why — a thin class's "risk" is its
   neighbours' risk plus the prior, so letting it win the argmax reports an artefact.
   191 of 672 classes qualify;
4. report the extremes **with** the tie structure and the stability.

Numbers to quote (2024 full model):

```
highest  class 181  P(LBW) = 0.2828 [0.2635, 0.2978]  N = 1,710    top in 100% of replicates
         sex=0, dmar=0, mager=0, meduc=1, precare1st=1, cig_0=1, race=Black
lowest   class  40  P(LBW) = 0.0539 [0.0520, 0.0548]  N = 346,811  bottom in 100% of replicates
         sex=1, dmar=1, mager=0, meduc=1, precare1st=1, cig_0=0, race=White
```

**Two methodological points that must appear, not be buried:**

* **Ties are structural.** Every class in a terminal node shares a posterior, so there is
  no "the" highest-risk class — only a highest-risk **leaf**. 2 classes tie at the
  maximum, 3 at the minimum. The representative is the tied member with the most births.
  Say this; a reader who notices two classes with identical intervals will otherwise
  assume an error.
* **Stability is tie-aware.** "Is this class the argmax in replicate b" scored 14 % for
  the low-risk class purely from arbitrary tie-breaking. The right question — is it
  *among* the extremes — gives 100 %.

### 4.3 The assumed profiles become a finding, not a definition

Do not delete them; report where they landed.

```
a priori high (class 146): P(LBW) = 0.2323 -> rank 13, and N = 509 (below the support threshold)
a priori low  (class  39): P(LBW) = 0.0650 -> rank 25 from the bottom
```

Two things to draw out: the assumed high-risk profile is **13th and rare** — a clinically
intuitive "everything bad at once" combination that is not the modal high-risk group; and
the empirical high-risk profile has `meduc=1` and `precare1st=1` where the assumption had
both 0. The model reaches maximum risk on `race`, `dmar` and `cig_0` **without** needing
low education or late care. That is a substantive result about which factors carry the
risk, and it only appears because the profiles were derived rather than assumed.

### 4.4 New subsection: the extreme ten and their shared features

This is the strongest new result in the revision — give it its own subsection, not a
footnote. The argument is that a single argmax invites "you picked the starkest cell",
and reading down the covariate columns of the top and bottom ten answers it.

```
HIGHEST 10                            LOWEST 10
 dmar   = 0     10/10  unanimous       sex   = 1      10/10  unanimous
 race   = Black  8/10                  cig_0 = 0      10/10  unanimous
 meduc  = 1      8/10                  race  = White  10/10  unanimous
 sex    = 0      7/10                  meduc = 1       7/10
 mager, precare1st, cig_0  6/10        dmar, mager, precare1st  6/10
```

Points worth making, in this order:

1. **Unmarried status is the only unanimous feature of the ten highest-risk profiles.**
   Race is 8/10 and smoking 6/10 — strong, not universal.
2. **The low-risk set is more tightly determined than the high-risk set** — three
   unanimous markers against one. Low risk requires a specific combination; high risk is
   reachable by several routes. That asymmetry is a real finding and is visible only
   because ten were examined.
3. **Two of the ten highest-risk profiles are White, not Black** (classes 67/68: over 33,
   unmarried, smoking). Say so explicitly — it resists a single-axis reading.
4. `mager`, `precare1st` and `meduc` never reach unanimity in either set; do not
   over-interpret them from one profile.

### 4.5 Figures — what replaces what

| old | new | note |
|---|---|---|
| Fig 3.x tree | `full_tree_2024.pdf` | the model actually fitted and bootstrapped |
| — | `risk_profile_ranking_2024.pdf` | all classes ranked; extremes in context |
| — | `risk_extremes_profiles_2024.pdf` | the top/bottom ten + shared features |
| — | `risk_paths_tree_2024.pdf` | both extremes traced through the tree |
| Figs 3.4–3.5 | `decision_tree_depth_{2..5}_2024.pdf` | now actually produced |
| Figs 3.11–3.12 | `depth_distributions_2024.png` | **now real data** — previously random |
| Tables 3.2/3.3 | `risk_summary` / `risk_totals` | percentile intervals |
| Table 3.1 | `var.freq.df` + `top.var.df` | the "initial split" row has a source again |

Two captions must be honest about scope:

* the subgroup bar chart shows **LBW categories only**; the NBW point estimate and bounds
  are in the subcaption. Say why — the NBW bar at ~0.72–0.95 flattens the ten LBW bars,
  and the within-LBW shape is the reason for expanding the LBW region into deciles.
* the risk-path figure traces through the **depth-3 display tree** while the quoted risk
  comes from the maxdepth-8 ensemble, so the leaf's own value differs from the headline.

---

## 5. A limitations paragraph you now need

Short, and it strengthens rather than weakens the paper:

* `dmar` unknown is a level, not an imputation; the geography that drives it is not in the
  public file. The `"drop"` arm is the sensitivity analysis.
* AIAN and NHOPI classes are thin; the intervals are correspondingly wide, and that width
  is reported rather than smoothed away. This is the DM prior working as intended, and it
  is also the honest ceiling on what can be said about those groups.
* Extreme classes are identified subject to `N >= 1,000`; a different threshold could
  select different classes, and the ranking figure shows the excluded ones so the reader
  can see the effect.
* `precare1st` is care *timing*, not adequacy (APNCU is not in this extract).

---

## 6. Suggested chapter-level restructuring

The current chapters are organised as *derivation → implementation → results*, which
buries the two things that are actually new. Consider:

* **Ch. 2** — keep the DM derivation as is; it is correct and it is the foundation. Add
  §2.4.4a on the categorical split search (§3 above). Rebuild the preprocessing section
  around the missingness and encoding decisions rather than treating them as
  bookkeeping — they are now methods.
* **Ch. 3** — reorder to: model fit and tree structure → depth-controlled comparison →
  variable importance/stability → **risk-profile ranking** → extreme ten and shared
  features. Right now the subgroup analysis reads as an afterthought; under the revision
  it is the main result and should be the destination the chapter builds toward.
* **Framing.** The old framing is "DM-CART applied to birth weight". The stronger framing
  after this revision is: *a Dirichlet-multinomial tree lets you disaggregate a
  categorical risk factor to a resolution a plain multinomial CART cannot support, and
  read risk profiles off the result.* The race expansion, the `dmar` unknown level and
  the sparse-cell intervals are three instances of one argument. Say that in the abstract
  and the introduction.

---

## 7. Numbers to re-generate before submission

Everything currently in this document came from a `B = 150–200` smoke run. Before any
number is quoted:

```bash
Rscript quantile.R 2023-2024
Rscript data.R     2023-2024
DMCART_B=10000 Rscript dm-cart.R 2023-2024
Rscript visual.R   2023-2024
Rscript visual.R   2023-2024 lbw
```

The extreme classes and the unanimity counts are stable at B = 200, but the interval
endpoints will move in the third decimal.

