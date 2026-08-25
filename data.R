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

DMAR.POLICY <- "unknown-level"
stopifnot(DMAR.POLICY %in% c("unknown-level", "drop"))
init.paths  <- get_paths(year = PRIOR.YEAR, with.2.5kg = TRUE)
TARGET.YEAR <- init.paths$target.year   # 2024 -- the cohort being modeled
PRIOR.YEAR  <- init.paths$prior.year    # 2023 -- cutpoints + alphavec ONLY
raw.csv.path <- init.paths$raw.csv(TARGET.YEAR)

# Both raw CSVs must be present before anything downstream runs: this script
# reads the target year, and quantile.R reads the prior year for the cutpoints.
# Reports rather than stops, so a missing file surfaces as a clear message
# instead of a read error several lines later.
check_raw_data <- function(paths) {
  needed  <- c(paths$raw.csv(paths$year), paths$raw.csv(paths$prior.year))
  missing <- needed[!file.exists(needed)]
  if (length(missing) > 0) {
    cat("Missing raw natality files -- download these before proceeding:\n")
    cat(paste(" -", missing), sep = "\n")
  } else {
    cat("All required raw natality files found for year", paths$year,
        "(and prior year", paths$prior.year, ").\n")
  }
  invisible(missing)
}
check_raw_data(init.paths)

cat(sprintf("\n==================================================\n"))
cat(sprintf(" Processing Data for Cohort: %s (Prior: %d, Target: %d)\n", 
            init.paths$cohort, PRIOR.YEAR, TARGET.YEAR))
cat(sprintf(" Target Raw Data Source: %s\n", raw.csv.path))
cat(sprintf("==================================================\n\n"))

# ----- Step 2: Read Raw Target Data (Once for efficiency) -----
if (!file.exists(raw.csv.path)) {
  stop(sprintf("Raw target data file not found: %s", raw.csv.path))
}

# data.table::fread instead of read.csv.
# The natality extract is ~2.1 GB and 237 columns wide
raw.cols <- c('sex', 'dmar', 'mracehisp', 'mager', 'meduc', 'precare5', 'cig_0', 'dbwt')
if (requireNamespace("data.table", quietly = TRUE)) {
  natalitydata <- as.data.frame(
    data.table::fread(raw.csv.path, select = raw.cols, showProgress = FALSE))
} else {
  natalitydata <- read.csv(raw.csv.path)[, raw.cols]
}
raw <- natalitydata[, raw.cols]

# NCHS sentinel "unknown" codes -> NA, BEFORE any recoding
raw$dbwt[raw$dbwt == 9999]        <- NA   # weight not stated
raw$cig_0[raw$cig_0 >= 99]        <- NA   # 99 = not stated (98 is a valid max count)
raw$meduc[raw$meduc == 9]         <- NA   # education not stated
raw$precare5[raw$precare5 == 5]   <- NA   # month care began not stated
raw$mracehisp[raw$mracehisp == 8] <- NA   # origin unknown or not stated
raw$sex[raw$sex == ""]            <- NA   # defensive; not observed in 2023/2024
if (is.character(raw$dmar)) raw$dmar[raw$dmar == ""] <- NA
raw$dmar <- suppressWarnings(as.integer(raw$dmar))

# `dmar` is blank (read as NA by read.csv) on a large, geographically
# clustered share of records -- handled by DMAR.POLICY below, not here.

# -------- Missingness report --------
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
missing.report$n_read <- nrow(raw)   # the flow table needs the denominator
# save(missing.report, file = sprintf("missing_report_%d.RData", TARGET.YEAR))

dat <- data.frame(dbwt = raw$dbwt)

# ==========================================================================
# race is now a configurable scheme, defaulting to all 7 groups
# ==========================================================================
# previously used binary Black / non-Black flag, this is revisited due to 
# shortcomings in the original draft 
#
#   "binary"      Black vs all other            I = 2^5 * 3 * 2 =  192
#   "collapsed5"  White/Black/Hispanic/Asian/Other  I = 2^5 * 3 * 5 =  480
#   "full7"       every mracehisp group         I = 2^5 * 3 * 7 =  672   <- default
# 
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

# `meduc >= 3` (previously `meduc == 3`) is the correct test for "completed high school or more"
# meduc is the NCHS 2003-revision education recode:
#   1 = 8th grade or less        5 = associate degree
#   2 = 9th-12th, no diploma     6 = bachelor's
#   3 = HS graduate or GED       7 = master's
#   4 = some college, no degree  8 = doctorate / professional
#   9 = not stated (-> NA above)
dat$meduc <- factor(ifelse(raw$meduc >= 3, 1, 0), levels = c(0, 1))  # 1 = completed HS or more

# `precare5` renamed to `precare1st`.
# precare5 codes when prenatal care began, not whether it was adequate:
#   1 = 1st trimester   2 = 2nd   3 = 3rd   4 = no prenatal care
#   5 = not stated (-> NA above)
# `== 1` is the correct test, not `<= 1` which would include the "no prenatal care" category.
dat$precare1st <- factor(ifelse(raw$precare5 == 1, 1, 0), levels = c(0, 1))  # 1 = care began in 1st trimester

dat$cig_0 <- factor(ifelse(raw$cig_0 > 0, 1, 0), levels = c(0, 1))   # 1 = any pre-pregnancy smoking

# -------- marital status (dmar) --------
# `dmar` is blank on ~11% of 2024 records. imputation  not possible since 
# every geographic identifier in the public-use extract is suppressed.
# default instead keeps every record and lets the model speak: `dmar` carries
# an explicit "Unknown" level.
if (DMAR.POLICY == "unknown-level") {
  dat$dmar <- factor(ifelse(is.na(raw$dmar), "Unknown",
                     ifelse(raw$dmar == 1, "1", "0")),
                     levels = c("0", "1", "Unknown"))
} else {  # "drop": listwise deletion, for the sensitivity arm
  dat$dmar <- factor(ifelse(raw$dmar == 1, 1, 0), levels = c(0, 1))
}
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


xnode <- as.integer(interaction(X[, feature.names], drop = FALSE))
num_node <- length(levels(interaction(X[, feature.names], drop = FALSE)))

# generate all possible combinations across the actual levels of each predictor
X.all <- expand.grid(lapply(X[, feature.names], levels))
colnames(X.all) <- feature.names
for (col in feature.names) {
  X.all[[col]] <- factor(X.all[[col]], levels = levels(X[[col]]))
}
X.all$xnode <- as.integer(interaction(X.all[, feature.names], drop = FALSE))

n.levels.each <- sapply(X[, feature.names], nlevels)
expected.classes <- prod(n.levels.each)
cat("\nPredictor level counts:\n")
print(n.levels.each)
cat(sprintf("num_node: %d (expected %s = %d)\n",
            num_node, paste(n.levels.each, collapse = " x "), expected.classes))
stopifnot(num_node == expected.classes)

# -------------------
# confirm X and X.all agree on what each xnode integer means.
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
names(final.data)[names(final.data) == "x.x"] <- "mean.dbwt"
names(final.data)[names(final.data) == "x.y"] <- "mean.log.dbwt"
final.data$count[is.na(final.data$count)] <- 0 # fill in zero counts for nodes with no observations

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

  # add NBW only for Full Model
  if (PROCESS.WITH.2.5KG) {
    above.count <- aggregate(list(n = (dat$dbwt > 2500)), list(xnode = xnode), sum)
    merged.above <- merge(data.frame(xnode = X.all$xnode), above.count, by = "xnode", all.x = TRUE)
    merged.above$n[is.na(merged.above$n)] <- 0
    counts[["counts_above_2.5kg"]] <- merged.above$n
  }

  counts.df <- do.call(cbind, counts)
  stopifnot(nrow(X.matrix) == num_node, nrow(counts.df) == num_node)

  # the counts must account for every analyzed birth, otherwise something went wrong in the aggregation
  expect.total <- if (PROCESS.WITH.2.5KG) nrow(dat) else sum(dat$dbwt <= 2500)
  if (sum(counts.df) != expect.total)
    stop(sprintf("Counts do not reconcile: sum(counts.df) = %s but %s births were analyzed.",
                 format(sum(counts.df), big.mark = ","),
                 format(expect.total, big.mark = ",")))
  cat(sprintf("Counts reconcile: %s births binned.\n",
              format(sum(counts.df), big.mark = ",")))

  # saved under TARGET.YEAR, not PRIOR.YEAR.
  output_rdata <- paths$rebin.data(TARGET.YEAR)
  save(
    X.matrix,
    Y.matrix,
    counts.df,
    final.data,
    prior.cutpoints,
    feature.names,
    missing.report,
    DMAR.POLICY,
    file = output_rdata
  )

  cat("Saved rebinned data to:", output_rdata, "\n")
}

cat("\nCompleted data processing and rebinning for all model variants.\n\n")