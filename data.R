# data.R
# 1) Collapse predictors into a single node variable (xnode) for each observation
# 2) Aggregate node statistics and construct design/counts matrices
# 3) Rebin target year data against prior year cutpoints for both Full and LBW-only models
#
# Usage via Command Line:
#   Rscript data.R 2023-2024
#   Rscript data.R 2023
#
# Default (Interactive / RStudio):
#   Defaults to cohort prior year 2023 (target year 2024)
library(rprojroot)
source(file.path(rprojroot::find_root(rprojroot::is_git_root), "util.R"))

# ----- Step 1: Parse Command Line Arguments -----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  arg_str <- args[1]
  PRIOR.YEAR <- as.integer(sub("-.*", "", arg_str))
} else {
  PRIOR.YEAR <- 2023
}

# ----- Analysis switches (see REVISIONS.md) -----
# REVISION R4: how to handle records where `dmar` (marital status) is blank.
#   "unknown-level" (default) -- keep the records; `dmar` becomes a THREE-level
#                                factor Married / Unmarried / Unknown.
#   "drop"                    -- listwise deletion (the old behaviour, but now
#                                explicit and reported). Use for the sensitivity arm.
# Imputation by state/region is NOT offered because the public-use natality
# file suppresses every geographic identifier -- see REVISIONS.md R4 for the
# verification of that claim.
DMAR.POLICY <- "unknown-level"
stopifnot(DMAR.POLICY %in% c("unknown-level", "drop"))

# Resolve target year and cohort paths
# REVISION R1: `$target.year`, not `$year`. `get_paths()` used to return the
# REQUESTED year here, which made TARGET.YEAR == PRIOR.YEAR == 2023 and meant
# the prior was built from the same data it was applied to (double-dipping),
# with the 2024 file never being read. See REVISIONS.md R1.
init.paths  <- get_paths(year = PRIOR.YEAR, with.2.5kg = TRUE)
TARGET.YEAR <- init.paths$target.year   # 2024 -- the cohort being modeled
PRIOR.YEAR  <- init.paths$prior.year    # 2023 -- cutpoints + alphavec ONLY
raw.csv.path <- init.paths$raw.csv(TARGET.YEAR)

cat(sprintf("\n==================================================\n"))
cat(sprintf(" Processing Data for Cohort: %s (Prior: %d, Target: %d)\n", 
            init.paths$cohort, PRIOR.YEAR, TARGET.YEAR))
cat(sprintf(" Target Raw Data Source: %s\n", raw.csv.path))
cat(sprintf("==================================================\n\n"))

# ----- Step 2: Read Raw Target Data (Once for efficiency) -----
if (!file.exists(raw.csv.path)) {
  stop(sprintf("Raw target data file not found: %s", raw.csv.path))
}

# REVISION R21: data.table::fread instead of read.csv.
# The natality extract is ~2.1 GB and 237 columns wide; read.csv parsed all of
# it to keep 8 fields and was by a wide margin the slowest step in the whole
# pipeline. fread with `select` reads only the columns we use. Same values,
# same types (kept as a plain data.frame so nothing downstream changes).
raw.cols <- c('sex', 'dmar', 'mracehisp', 'mager', 'meduc', 'precare5', 'cig_0', 'dbwt')
if (requireNamespace("data.table", quietly = TRUE)) {
  natalitydata <- as.data.frame(
    data.table::fread(raw.csv.path, select = raw.cols, showProgress = FALSE))
} else {
  natalitydata <- read.csv(raw.csv.path)[, raw.cols]
}

# BOY: sex
# MARRIED: dmar
# OVER33: mager
# HIGH SCHOOL: meduc
# FULL PRENATAL: precare5
# SMOKER: cig_0
# BIRTH WEIGHT: dbwt (response variable)

# Select raw features (raw.cols defined above, at the read)
raw <- natalitydata[, raw.cols]

# ==========================================================================
# REVISION R3: NCHS sentinel "unknown" codes -> NA, BEFORE any recoding.
# ==========================================================================
# `na.omit()` only catches R's NA. NCHS encodes "unknown / not stated" as
# ordinary IN-RANGE integers, so every one of these was previously flowing
# into the model as if it were a real observed value:
#
#   field       code   what it really means    what it USED to become
#   ---------   ----   ---------------------   -------------------------------
#   dbwt        9999   weight not stated       counted as Normal (> 2500 g)  <-- worst
#   cig_0       99     smoking not stated      smoker = 1
#   meduc       9      education not stated    "did not complete HS" = 0
#   precare5    5      care timing not stated  "not 1st trimester" = 0
#   mracehisp   8      origin not stated       (already handled correctly)
#
# `dbwt == 9999` is the most damaging: it inflated the Normal-birth-weight
# category, which is the denominator for the entire LBW story.
# Counts are reported in the missingness table below.
# --------------------------------------------------------------------------
raw$dbwt[raw$dbwt == 9999]        <- NA   # weight not stated
raw$cig_0[raw$cig_0 >= 99]        <- NA   # 99 = not stated (98 is a valid max count)
raw$meduc[raw$meduc == 9]         <- NA   # education not stated
raw$precare5[raw$precare5 == 5]   <- NA   # month care began not stated
raw$mracehisp[raw$mracehisp == 8] <- NA   # origin unknown or not stated
raw$sex[raw$sex == ""]            <- NA   # defensive; not observed in 2023/2024
# fread returns "" for blank character fields and NA for blank numerics; dmar is
# numeric, but normalise defensively so DMAR.POLICY sees a real NA either way.
if (is.character(raw$dmar)) raw$dmar[raw$dmar == ""] <- NA
raw$dmar <- suppressWarnings(as.integer(raw$dmar))

# `dmar` is blank (read as NA by read.csv) on a large, geographically
# clustered share of records -- handled by DMAR.POLICY below, not here.

# -------- Missingness report (REVISION R3) --------
# Printed and saved so the paper can state exactly what was excluded and why,
# rather than silently losing rows inside na.omit().
missing.report <- data.frame(
  variable  = raw.cols,
  n_missing = sapply(raw.cols, function(v) sum(is.na(raw[[v]]))),
  pct       = sapply(raw.cols, function(v) 100 * mean(is.na(raw[[v]]))),
  row.names = NULL
)
cat("\n--- Missingness in raw target-year fields (year", TARGET.YEAR, ") ---\n")
print(transform(missing.report, pct = sprintf("%.3f%%", pct)), row.names = FALSE)
cat(sprintf("Total records read: %d\n", nrow(raw)))
save(missing.report, file = sprintf("missing_report_%d.RData", TARGET.YEAR))

# ==========================================================================
# Covariate construction
# ==========================================================================
# All encodings below read from `raw` (sentinel-cleaned, never re-encoded) and
# write into a fresh `dat`. Keeping the source and the encoded output in
# separate objects makes this block safe to re-run in a live session --
# comparing an already-factored column with `>`/`==` throws "not meaningful
# for factors" and silently yields all-NA, which is what caused the earlier
# warning when everything was written back into one data frame.
dat <- data.frame(dbwt = raw$dbwt)

# ==========================================================================
# REVISION R20: race is now a configurable scheme, defaulting to all 7 groups
# ==========================================================================
# The thesis used a binary Black / non-Black flag and says of it, in
# chapter2.tex:45, "*though this choice is entirely arbitrary*. (NEED TO
# REVISIT)". This is that revisit.
#
#   "binary"      Black vs all other            I = 2^5 * 3 * 2 =  192
#   "collapsed5"  White/Black/Hispanic/Asian/Other  I = 2^5 * 3 * 5 =  480
#   "full7"       every mracehisp group         I = 2^5 * 3 * 7 =  672   <- default
#
# (Class counts include the 3-level `dmar` from R4.)
#
# Why full7 is the default rather than a collapse: "Other" as a bucket is the
# same arbitrariness the thesis apologises for, one level up -- it merges AIAN,
# NHOPI and multiracial mothers, who have materially different LBW profiles,
# into a category that means nothing clinically. The tail groups ARE sparse
# (NHOPI averages a few hundred births per class), but that is precisely the
# regime the Dirichlet prior exists for: it shrinks thin cells toward the
# marginal profile and the parametric bootstrap reports honest width on them.
# Under this scheme the paper DEMONSTRATES the sparsity claim it currently only
# asserts. `collapsed5` remains one line away if a referee wants it.
#
# NOTE: this makes the exhaustive categorical search (R2) load-bearing, not
# optional. With 7 levels there are 2^6 - 1 = 63 bipartitions and the old
# sorted-code search reached only 6 of them.
RACE.SCHEME <- "full7"
stopifnot(RACE.SCHEME %in% c("binary", "collapsed5", "full7"))

race.map <- list(
  # mracehisp: 1 White NH, 2 Black NH, 3 AIAN NH, 4 Asian NH, 5 NHOPI NH,
  #            6 multiracial NH, 7 Hispanic, 8 not stated (-> NA above)
  binary = c(`1` = "NonBlack", `2` = "Black", `3` = "NonBlack", `4` = "NonBlack",
             `5` = "NonBlack", `6` = "NonBlack", `7` = "NonBlack"),
  collapsed5 = c(`1` = "White", `2` = "Black", `3` = "Other", `4` = "Asian",
                 `5` = "Other", `6` = "Other", `7` = "Hispanic"),
  full7 = c(`1` = "White", `2` = "Black", `3` = "AIAN", `4` = "Asian",
            `5` = "NHOPI", `6` = "Multiracial", `7` = "Hispanic")
)
race.levels <- list(
  binary     = c("NonBlack", "Black"),
  collapsed5 = c("White", "Black", "Hispanic", "Asian", "Other"),
  full7      = c("White", "Black", "Hispanic", "Asian", "AIAN", "NHOPI", "Multiracial")
)
race.labels <- race.map[[RACE.SCHEME]]
dat$race <- factor(race.labels[as.character(raw$mracehisp)],
                   levels = race.levels[[RACE.SCHEME]])
cat(sprintf("\nRACE.SCHEME = '%s' (%d levels)\n", RACE.SCHEME, nlevels(dat$race)))

# -------- Binary predictors --------
dat$sex   <- factor(ifelse(raw$sex == "M", 1, 0), levels = c(0, 1))  # 1 = male
dat$mager <- factor(ifelse(raw$mager > 33, 1, 0), levels = c(0, 1))  # 1 = over 33

# `meduc == 3` -> `meduc >= 3`.
# meduc is the NCHS 2003-revision education recode:
#   1 = 8th grade or less        5 = associate degree
#   2 = 9th-12th, no diploma     6 = bachelor's
#   3 = HS graduate or GED       7 = master's
#   4 = some college, no degree  8 = doctorate / professional
#   9 = not stated (-> NA above)
# Table 2.1 defines this predictor as "high school completed / otherwise",
# but `== 3` selects ONLY mothers who stopped at a diploma -- it labelled a
# mother with a bachelor's degree as not having completed high school. The
# correct encoding of "completed HS" is `>= 3`. See REVISIONS.md R5.
dat$meduc <- factor(ifelse(raw$meduc >= 3, 1, 0), levels = c(0, 1))  # 1 = completed HS or more

# REVISION R6: `precare5` renamed to `precare1st`.
# precare5 codes WHEN prenatal care began, not whether it was adequate:
#   1 = 1st trimester   2 = 2nd   3 = 3rd   4 = no prenatal care
#   5 = not stated (-> NA above)
# `== 1` is the correct test, but the old name/gloss ("adequate prenatal
# care") overclaims: adequacy is the Kotelchuck/APNCU index, which is not in
# this extract. The predictor now says what it measures. See REVISIONS.md R6.
dat$precare1st <- factor(ifelse(raw$precare5 == 1, 1, 0), levels = c(0, 1))  # 1 = care began in 1st trimester

dat$cig_0 <- factor(ifelse(raw$cig_0 > 0, 1, 0), levels = c(0, 1))   # 1 = any pre-pregnancy smoking

# -------- REVISION R4: marital status (dmar) --------
# `dmar` is blank on ~20% of 2024 records. Those records are NOT junk: they
# carry valid dbwt, race and education. Blank `dmar` co-occurs exactly with
# blank `mar_p`, i.e. it is a reporting-area (state-level) pattern, and
# because the file is ordered by state the loss is geographically clustered.
# Under the old `ifelse(dmar == 1, 1, 0)` these became NA and were deleted by
# na.omit() with no message -- a non-random 20% deletion, unreported.
#
# Imputation by region was considered and is NOT possible from this file:
# every geographic identifier in the public-use extract is suppressed (blank
# in 100% of records) -- ocntyfips, ocntypop, rcnty, rcnty_pop, rcity_pop,
# mrterr, octerr, mbcntry. There is no state, county, division or region
# field to condition an imputation on, and the missingness mechanism is
# precisely the unobserved geography, so MAR given the observed covariates is
# not credible either. See REVISIONS.md R4 for the verification.
#
# Default instead keeps every record and lets the model speak: `dmar` carries
# an explicit "Unknown" level. DM-CART then reports empirically which side of
# a split the unknown group falls on, instead of assuming it.
if (DMAR.POLICY == "unknown-level") {
  dat$dmar <- factor(ifelse(is.na(raw$dmar), "Unknown",
                     ifelse(raw$dmar == 1, "1", "0")),
                     levels = c("0", "1", "Unknown"))
} else {  # "drop": listwise deletion, for the sensitivity arm
  dat$dmar <- factor(ifelse(raw$dmar == 1, 1, 0), levels = c(0, 1))
}

# Order columns so the response is last and predictors read in Table 2.1 order
feature.names <- c("sex", "dmar", "mager", "meduc", "precare1st", "cig_0", "race")
dat <- dat[, c(feature.names, "dbwt")]
p <- length(feature.names)

# -------- Complete-case filter (now explicit and reported) --------
n.before <- nrow(dat)
dat <- na.omit(dat)
n.after <- nrow(dat)
cat(sprintf("\nComplete-case filter: %d -> %d records (%d dropped, %.2f%%)\n",
            n.before, n.after, n.before - n.after,
            100 * (n.before - n.after) / n.before))
cat(sprintf("DMAR.POLICY = '%s'\n", DMAR.POLICY))
missing.report$policy_dmar <- DMAR.POLICY
missing.report$n_analysed  <- n.after

X <- dat[, feature.names]
y <- log(dat$dbwt)

# generalized class index: works for a mix of binary + multi-level factors
# drop=FALSE is required -- it fixes the level set (and therefore the
# integer coding) to the full theoretical set of combinations, regardless
# of whether every combination is actually observed in X. Using drop=TRUE
# here would silently renumber classes if any cell were empty, misaligning
# xnode between X and X.all with no warning.

xnode <- as.integer(interaction(X[, feature.names], drop = FALSE))
num_node <- length(levels(interaction(X[, feature.names], drop = FALSE)))

# generate all possible combinations across the actual levels of each predictor
X.all <- expand.grid(lapply(X[, feature.names], levels))
colnames(X.all) <- feature.names
for (col in feature.names) {
  X.all[[col]] <- factor(X.all[[col]], levels = levels(X[[col]]))
}
X.all$xnode <- as.integer(interaction(X.all[, feature.names], drop = FALSE))

# REVISION R4/R5: the expected class count is now COMPUTED from the actual
# level sets rather than hardcoded as "2^6 * 5 = 320". It changes with
# DMAR.POLICY (dmar gains an "Unknown" level -> 3 levels) and would change
# again under the planned race expansion, so a hardcoded literal here would
# go stale silently.
n.levels.each <- sapply(X[, feature.names], nlevels)
expected.classes <- prod(n.levels.each)
cat("\nPredictor level counts:\n")
print(n.levels.each)
cat(sprintf("num_node: %d (expected %s = %d)\n",
            num_node, paste(n.levels.each, collapse = " x "), expected.classes))
stopifnot(num_node == expected.classes)

# -------------------
# Sanity check: confirm X and X.all agree on what each xnode integer means.
# Cheap insurance against Bugs 4/5 recurring silently after future edits.
check.classes <- sample(unique(xnode), min(5, length(unique(xnode))))
for (cls in check.classes) {
  from.X    <- unique(X[xnode == cls, feature.names])
  from.Xall <- X.all[X.all$xnode == cls, feature.names]
  stopifnot(nrow(from.X) == 1, all(as.character(from.X) == as.character(from.Xall)))
}
cat(sprintf("xnode alignment check passed for %d sampled classes\n", length(check.classes)))
# -------------------

cat("classes with zero observations:", sum(table(factor(xnode, levels=1:num_node)) == 0), "\n")
cat("race distribution:\n"); print(table(dat$race))
cat("dmar distribution:\n"); print(table(dat$dmar))



# Aggregations
agg.dbwt <- aggregate(dat$dbwt, list(xnode = xnode), mean)
agg.log.dbwt <- aggregate(y, list(xnode = xnode), mean)
count.observations <- aggregate(list(count = rep(1, nrow(dat))), list(xnode = xnode), sum)

final.data <- merge(X.all, agg.dbwt, by = "xnode", all.x = TRUE)
final.data <- merge(final.data, agg.log.dbwt, by = "xnode", all.x = TRUE)
final.data <- merge(final.data, count.observations, by = "xnode", all.x = TRUE)

# REVISION R25: an unobserved class has NO mean birth weight -- 0 g is not a
# measurement, it is a missing value, and writing 0 there would poison any
# later summary that averaged this column. `x.x` (mean dbwt) and `x.y` (mean
# log dbwt) are therefore left as NA; only `count` is genuinely 0.
# (Neither mean feeds the model; the counts matrix is what DM-CART consumes.)
names(final.data)[names(final.data) == "x.x"] <- "mean.dbwt"
names(final.data)[names(final.data) == "x.y"] <- "mean.log.dbwt"
final.data$count[is.na(final.data$count)] <- 0

X.matrix <- final.data[, feature.names]
Y.matrix <- data.frame(count = final.data$count)



# ----- Step 3: Loop Across Both Rebin Types (With & Without 2.5kg) -----
for (PROCESS.WITH.2.5KG in c(TRUE, FALSE)) {

  model_label <- if (PROCESS.WITH.2.5KG) "Full Model (rebin)" else "LBW-Only Model (rebin_without_2.5kg)"
  cat(sprintf("\n--------------------------------------------------\n"))
  cat(sprintf(" Rebinning Data for: %s\n", model_label))
  cat(sprintf("--------------------------------------------------\n"))

  paths <- get_paths(year = PRIOR.YEAR, with.2.5kg = PROCESS.WITH.2.5KG)
  ensure_dirs(paths)

  cutpoints_file <- paths$quantile.cutpoints(PRIOR.YEAR)
  if (!file.exists(cutpoints_file)) {
    stop(sprintf("Cutpoints file not found: %s. Please run quantile.R first!", cutpoints_file))
  }
  
  load(cutpoints_file)
  prior.cutpoints <- cut.points
  probs <- seq(0, 1, by = 0.1)

  # Bin counts for LBW categories (Q1..Q10)
  counts <- list()
  for (i in 1:(length(prior.cutpoints) - 1)) {
    # REVISION R35: the first bin is closed on the LEFT.
    # The rule was `dbwt > cutpoint[1]` for every bin, so any 2024 birth at or
    # below the 2023 minimum landed in no bin at all and vanished from the
    # counts -- 29 births at exactly 227 g, of which 24 survived the
    # complete-case filter, which is why sum(counts.df) came out 24 short of
    # the analysed N. Small, but they are the very smallest births in the file:
    # precisely the ones this analysis is about. Because the cutpoints come
    # from a DIFFERENT year, 2024 can also contain a birth lighter than
    # anything seen in 2023, so the first bin takes everything up to
    # cutpoint[2] rather than being bounded below at all.
    in.bin <- if (i == 1) {
      dat$dbwt <= prior.cutpoints[i + 1]
    } else {
      dat$dbwt > prior.cutpoints[i] & dat$dbwt <= prior.cutpoints[i + 1]
    }
    bin.count <- aggregate(list(n = in.bin), list(xnode = xnode), sum)
    merged <- merge(data.frame(xnode = X.all$xnode), bin.count, by = "xnode", all.x = TRUE)
    merged$n[is.na(merged$n)] <- 0

    name <- paste0(
      "counts_",
      formatC(probs[i],   format="f", digits=2), "_",
      formatC(probs[i+1], format="f", digits=2),
      "_quantile_",
      round(prior.cutpoints[i]   / 1000, 2), "_",
      round(prior.cutpoints[i+1] / 1000, 2),
      "kg"
    )
    counts[[name]] <- merged$n
  }

  # Add Normal category (> 2.5kg) only for Full Model
  if (PROCESS.WITH.2.5KG) {
    above.count <- aggregate(list(n = (dat$dbwt > 2500)), list(xnode = xnode), sum)
    merged.above <- merge(data.frame(xnode = X.all$xnode), above.count, by = "xnode", all.x = TRUE)
    merged.above$n[is.na(merged.above$n)] <- 0
    counts[["counts_above_2.5kg"]] <- merged.above$n
  }

  counts.df <- do.call(cbind, counts)
  stopifnot(nrow(X.matrix) == num_node, nrow(counts.df) == num_node)

  # REVISION R35: the counts must account for every analysed birth. For the
  # full model that is all of them; for the LBW-only model it is exactly the
  # births at or below 2.5 kg. Checked rather than assumed.
  expect.total <- if (PROCESS.WITH.2.5KG) nrow(dat) else sum(dat$dbwt <= 2500)
  if (sum(counts.df) != expect.total)
    stop(sprintf("Counts do not reconcile: sum(counts.df) = %s but %s births were analysed.",
                 format(sum(counts.df), big.mark = ","),
                 format(expect.total, big.mark = ",")))
  cat(sprintf("Counts reconcile: %s births binned.\n",
              format(sum(counts.df), big.mark = ",")))

  # REVISION R1: saved under TARGET.YEAR, not PRIOR.YEAR. These are the
  # target-year (2024) counts; dm-cart.R already loads rebin.data(TARGET.YEAR),
  # so the old PRIOR.YEAR save was a name mismatch that only went unnoticed
  # because TARGET.YEAR was itself wrongly resolving to PRIOR.YEAR.
  output_rdata <- paths$rebin.data(TARGET.YEAR)
  save(
    X.matrix,
    Y.matrix,
    counts.df,
    final.data,
    prior.cutpoints,
    feature.names,    # REVISION R6: consumers must not assume "precare5"
    missing.report,   # REVISION R3: exclusions travel with the data
    DMAR.POLICY,      # REVISION R4: which arm produced this file
    file = output_rdata
  )

  cat("Saved rebinned data to:", output_rdata, "\n")
}

cat("\nCompleted data processing and rebinning for all model variants.\n\n")