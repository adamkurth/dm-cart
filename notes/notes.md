## Revisions Notes
1. `mysplit()` was changed to handle multi-level categorical predictors correctly. The old version searched a small subset of possible splits, which was fine for binary predictors but not for multi-level ones. Examples of the multi-level groupings (splits) would not have been considered, and the resulting tree would have been suboptimal. The new version exhaustively searches all possible partitions of the levels by enumerating over $2^{k-1} - 1$ possible splits for a $k$-level categorical predictor. This is now correct using $\Delta\ell$ criterion as before, but the search space is now complete.
2. Missingness: before using `na.omit()` instantly dropped any rows with missing values —— proceeding with complete cases. From a sample of 2M records, we found a damaging error in how these codes were considered. `dbwt` = 9999 (birthweight not stated in record), was counted as >2500 g into the NBW category! This was the most damaging, however, smoking `cig_0` (99) became `smoker==1` (not true), `meduc` (9) became `did not complete HS` (not true), `precare5` (5) became `not 1st trimester` (not true), and `mracehisp` (8) were also misclassified. The fix was to recode these sentinel values to NA before any further processing.

Now: we convert these all to NAs before any recoding is done in the data.R script (and likewise done in quantile.R), and we create a missingness report for 2024 dataset (which should be included in the thesis appendix). The report is also saved to a file for later reference. 
```
--- Missingness in raw target-year fields (year 2024 ) ---
  variable n_missing     pct
       sex         0  0.000%
      dmar    402966 11.075%
 mracehisp     35805  0.984%
     mager         0  0.000%
     meduc     69489  1.910%
  precare5     70222  1.930%
     cig_0     18039  0.496%
      dbwt      3106  0.085%
Total records read: 3638436
```
2. `dmar` (marital status) was missing in about ~20% of records. The missingness in this variable was not random, it varied by geography. By first thought was to impute the missing values, but since geographic region (state/county/region) is not available in this public dataset (what so ever). So I could not impute conditioning on region, nor would a MAR assumption be testable. Instead, I encoded `dmar` as a 3-level factor: 0 = married, 1 = unmarried, 2 = unknown. This retained all of the records, did not impose any unverified assumptions of MAR (imputing must assume something about how the group is distribute, this approach assumes nothing), and the model answers the question empirically (i.e. likelihood reports which group is prefered. If in the bootstrap, the unknown group is consistently grouped with the married or unmarried, that is evidence about the reporting areas— this is a finding not an assumption), and lastly it's honest about the uncertainty of missingness rather than hiding it behind an imputation. The cost, however, is that we gain a 3rd level and changes this from $2^6 * 5 = 320$ to $2^5 * 5 * 3 = 480$ 
3. `meduc` recorded to be >= 3 (meaning high school or above) this changes from:
```
OLD (meduc == 3): 107,867 of 378,833 flagged "completed HS"  (28.5%)
NEW (meduc >= 3): 328,669 of 369,621 flagged "completed HS"  (88.9%)
```
A binary-race reproduction arm will no longer be able to reproduce the old results, beacuse the thesis's numbers were computed using the wrong definition of "completed HS". The new definition is correct, and the old one was wrong. 

4. `precare5` renamed to `precare1st`. Coding is `1` = 1st trimester, `2` = 2nd, `3` = 3rd, `4` = no care, `5` = not stated. So `==1` was always the right test, but the name was misleading. Table 2.1 describes the coding as "adequate prenatal care", but the variable measures when the care started, not whether it was adequate. Adequacy is the Kotelchuck/APNCU index, which is not in this extract. (no numerical results change, just make sure this is changed in the thesis text.)
5. `mypred()` sent the rows (classes) to the wrong leaves. The old function manually traversed by hand: 
```r
split_val <- fit$splits[split_rows[1], "index"]
go_left   <- as.numeric(cur_val) <= split_val
```
where (1) for a categorical split, `fit$splits[,"index"]` is not a cutpoint. Its the row number of `fit$csplit`, the table that actually records which levels go left (`=1`) and which go right (`=3`), `=2` not present at this node. Comparing the factors against a `csplit` row index is meaningless since every predictor is in this model is a factor, so this was wrong at every split, including binary ones. It was never a race variable problem. (2) `split_rows[1]` takees the first from `fit$splits` row bearing the variable name, so a variable reused at several nodes was always routed by its first split's parameters, whichever node was actually being visited. 


This is fixed: No traversal is needed. `X.matrix` is the **complete enumeration of classes**, and every tree here is fit on exactly those rows — so rpart has already done the routing. `fit$where` holds the frame **row index** (not the node id) of the leaf each training row fell into, and `fit$frame$yval2` holds the per-node count vector returned by `myeval()`.

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

**Scope limit, stated in the code comment too:** `dm.leaf.probs()` is valid only for the rows the tree was *fit on*. That is all this design ever needs, because the class enumeration is the design matrix. If you ever need to score genuinely new rows, the only correct route is decoding `csplit` **per node** (via `frame$var` and the node's own `splits` row), never per variable.

480-class design, 11 categories:

```
dim: 480 x 11                                  ✓
rows summing to 1: 480 of 480                  ✓
all finite and in [0,1]: TRUE                  ✓
matches independent recomputation: TRUE        ✓
```


6. Fixed `rpart.plot()`: previously this handed a `total.counts <- sum(yval); lbl <- paste0(format(which.max(total.counts))` which total.counts was a single number, so `which.max()` always returned 1 and no information was being passed to the label. Now it handles a matrix, where we can report modal birth-weight category per node. the `text` callback additionally folds the node's posterior P(LBW) using the same $(\mathbf{y} + \alpha) / (N+ \alpha_0)$ form as `dm.leaf.probs()`, so the node label is now informative. The `rpart.plot()` call is now:
```r
m     <- if (is.matrix(yval)) yval else matrix(yval, nrow = 1)
modal <- apply(m, 1, which.max)
p.lbw <- apply(m, 1, function(z) { post <- (z + alphavec)/sum(z + alphavec); sum(post[lbw.cols]) })
lbl   <- paste0("k", modal, "\nP(LBW)=", sprintf("%.3f", p.lbw))
if (use.n) lbl <- paste0(lbl, "\nn=", n)
```
now P(LBW) is in the label rather than drawn seperately as `text()` call. 

Additionally, `rpart.plot::rpart.plot(fit)` fails because of the custom method fit. We made our own `dm.rpart.plot()` function 
...

7. Depth-controlled comparison: Before the bootstrap stage, Section 3.2.4/Figures 3.4-3.5 analysis we made a `dm.leaf.probs()` function which outputs the posterior probabilities of each class. At each depth = {2,3,4,5} records, per depth: variables used, terminal-node count, first split, and the distribution of posterior P(LBW) across classes. 

8. The depth dataframe (depth_df) comes now from real bootstrap results. Previously, if performed, only the first split or a subset of nodes might have been recorded. Now, every internal node of each bootstrap is recorded:
```r 
internal <- dm.tree.b$frame$var != "<leaf>"
node.ids <- as.numeric(rownames(dm.tree.b$frame))
depth.rows[[b]] <- data.frame(replicate = b,
                              Variable = as.character(dm.tree.b$frame$var[internal]),
                              Depth    = floor(log2(node.ids[internal])))
```
`rpart` names frame rows by binary-heap node id (root = 1, left child = 2*parent, right child = 2*parent + 1) etc, so `depth=floor(log2(node.ids)` and the root depth 0. It must be read from `rownames(frame)`, not row position 


9. High/Low-risk subgroups are single classes, `risk_summary` is saved: (1) previously the call to grab the subgroups was `which(X.matrix$race == "Black" & X.matrix$dmar == 0 & X.matrix$cig_0 == 1)` matched 32 / 480 classes, then averaged over all combinations of the remaining 4 predictors which produced blended probabilities rather than sharp contrasts of Tables 3.2/3.3. UPDATE!! (2) `risk_summary` was not saved. (3) the risk profiles were assumed to be high/low risk based on interpretation of the profiles, but this is not actually the case and must confirm that these are in fact the highest/lowest risk subgroups

10. the high/low risk subgroups analysis: this used to show the *assumed* (reading the covariates as what ought to be best/worst) highest/lowest risk subgroups but now this is computed empirically from the bootstrap results. Every class is ranked by the posterior risk pooled over all the bootstrap replicates `lbw.by.class[i,b] = sum over all risk categories of Yhat[[b]][i,]` where `Yhat` is the leaf posterior each class inherits, so this is the model's fitted risk for the exact combination of covariates, with the bootstrap supplying the sampling distribution. Extremes are `argmax`/`argmin` of the posterior mean with percentile bounds. 
    
This forced 4 things: (1) the risk metric has to differ by model arm: In the LBW-only model, there is no normal category, so the `lbw.cols` span every category and `sum(P(LBW))=1` is 1 for all classes in the LBW-only model. The first run of produced this for all classes with same classes reported as highest/lowest. That arm's conditional on being LBW, so its severity score is now the mass in the smallest 3 deciles: 

 
| arm | metric | meaning |
|---|---|---|
| Full (K = 11) | `P(LBW)` = Q1..Q10 | probability of low birth weight |
| LBW-only (K = 10) | $P(Q1-Q3 \mid LBW)$ | share of the three most severe LBW deciles |


(2) A tree gives every class in the same terminal node the *same* posterior, so classes sharing a leaf have the same risk technically. However, these classes share the same posterior probabilies *by leaf* so many classes are tied for highest/lowest risk. The code now identifies the tied set, reports its size, and picks the member with the most observed births as the representative (the profile describing the largest real population). On your 2024 full model: 2 classes tie at the maximum, 3 at the minimum.

(3) Stability has to be tie-aware too: Asking "is the chosen class the argmax in replicate b?" is not enough, scored 14% for the low-risk class purely because the ties were broken differently in the bootstrap. The right question is whether the chosen class is *among* the extreme, i.e. whether its value equals the min/max of the bootstrap replicate. Both extremes are now 100%. 

(4) low support subclasses are excluded from the high/low risk analysis: A class with few observed births contributes very little to its leaf, so its "risk" is really its neighbor's plus the dirichlet prior. LEtting such a class win the argmax would report artefacts rather than a real finding. `MIN.CLASS.BIRTHS <- 1000` now gates the high/low risk analysis for classes to qualify. 

```
HIGHEST  class 181  P(LBW)=0.2834 [0.2624, 0.2954]  N=1,710    top in 100% of replicates
   sex=0, dmar=0, mager=0, meduc=1, precare1st=1, cig_0=1, race=Black
LOWEST   class  40  P(LBW)=0.0538 [0.0520, 0.0547]  N=346,811  bottom in 100% of replicates
   sex=1, dmar=1, mager=0, meduc=1, precare1st=1, cig_0=0, race=White
```
two findings for the write-up:. The assumed high-risk profile is 13th, and sits below the threshold of 1000 births (at 509), so it is not a real finding. The empirical high-risk profile has `meduc=1` and `precare1st=1` (completed HS, first-trimester care), ...


Figure ranking: every class has a posterior risk with 95% interval, ordered by risk; ineligible classes are greyed out in the plot; the two empirical extremes and the two a priori profiles marked. This is the visual answer to "were our assumed profiles actually the extremes?" and it shows the extremes in the context of the whole distribution rather than in isolation.


11.  Race expanded to full 7-level factor: `White`, `Black`, `Hispanic`, `Asian`, `AIAN`, `NHOPI`, `Multiracial` (2^5 * 3 * 7 = 672 classes). If we used `other` instead of the 5-factor method, this would merge AIAN (American Indian/Alaska Native) and NHOPI (Native Hawaiian/Other Pacific Islander) into a single "other" category, which would be misleading and resulting in meaning nothing clinically. These *are* distinct populations with different risk profiles. The tail groups are genuinely sparse, and this is exactly what the Dirichlet prior is designed to handle. 

The counts: are 
```
White 1,739,008 | Hispanic 944,512 | Black 455,381 | Asian 216,958
Multiracial 88,523 | AIAN 22,905 | NHOPI 9,620
```
with only 18 profiles yielding 0 observations (carrying prior mass only). 

This makes the revision of the `mysplit()` function even more important, because the old version would have missed many of the possible splits for these multi-level race predictors. The new version exhaustively searches all possible partitions of the levels by enumerating over $2^{k-1} - 1 = 2^{7-1} - 1 = 63$ bipartitions and the old code reached 6. This represents a genuine multi-level grouping on both sides, not one-vs-rest. 


12.  Extreme 10 risk subgroups visual: Among eligible ($<1000$ births) highest/lowest 10 are ranked by posterior mean risk, with 95% intervals also showing their profile composition and commonality. A covariate that is constant across all 10 profiles is noted. Among the highest risk profiles, `dmar=0` (unmarried) is the only constant, where in the lowest risk profiles, `sex=1` (male), `race=White`, `cig_0=0` (non-smoker), were all constant. This is a visual answer to "what are the commonalities among the extreme risk profiles?" and it shows the extremes in the context of the whole distribution rather than in isolation. `consensus.of()` function calculates the n/10, and marks this as unanimous if n=10. 

- Unmarried status is the only unanimous covariate among the highest risk profiles. Race 8/10, and smoking 6/10. Both are strong, but not universal from our findings. 
- Lowest risk was more tightly determined than the highest risk set, 3 unanimous markers: male, white, non-smoker. The lowest risk profiles are more homogeneous than the highest risk profiles, and this means that there are many more ways to be high-risk than low-risk. This is a finding, not an assumption. (MAKE SURE THIS IS CORRECT)
- 2/10 of the highest risk profiles are "White", not any other race category. These are classes 67/68: over 33, unmarried, smoking.

13. limitations paragraph: 
    
- `dmar` unknown is now a level, not imputed. the geographi that drives it is not in the public file. The "drop" arms is in the sensitivity analysis (??) 
- AIAN/NHOPI classes are thinl intervals are correspondingly wide, and that width is reported rather than smoothed away. the DM prior is working as intended, and handling also the honest ceiling on what cna be said about those groups. 
- Extreme classes are identified to `N>=1000`; a different threshold would yield different select different extreme classes, and the ranking figure shows the excluded ones so the reader can see what is excluded.
- `precare1st` is care *timing* not adequacy of (APNCU is not this extract)

14. Thesis changes: current organization is *derivations -> implementation -> results* which buries two things that are now NEW. consider: 
   - Ch 2 keep DM derivation as is, this is the foundation. Add 2.4.4a on the category split search (`mysplit()`) (S3 method in REVISIONS.md). Rebuild the preprocessing section around missingness and encoding decisions rather than treating them as bookkeeping, they are now methods of modeling.
   - Ch 3 reorder to: model fit and tree structure -> depth controlled comparison -> variable importance / stability -> **risk profile ranking** -> extreme 10 and shared features. Right now the subgroup analysis reads as an afterthought, but under the revision its the main result and should be the destination that the chapter builds toward. The depth-controlled comparison is a sanity check, and the variable importance/stability is a check on the model's robustness, but the risk profile ranking is the main result.
   - Framing: old: DM-CART applied to LBW but the stronger framing after this revision is: *DM Trees let you disaggregate a categorical risk factor to a resolution a plain multinomial CART cannot support, and read risk profiles off the result.* The race expansion, the `dmar` unknown level, sparse cell intervalsare three instances of one argument. Say this in the abstact and the introduction.












