# quantile.R
# Computes quantile cutpoints and Dirichlet prior hyperparameters from raw natality data.
#
# Usage via Command Line:
#   Rscript quantile.R 2023-2024
#   Rscript quantile.R 2023
#
# Default (Interactive / RStudio):
#   Defaults to prior year 2023

library(rprojroot)
source(file.path(rprojroot::find_root(rprojroot::is_git_root), "util.R"))

PROCESS.WITH.2.5KG <- TRUE  # Set explicitly for downstream scripts

# ----- Step 1: Parse Command Line Arguments -----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  arg_str <- args[1]
  # Extracts starting year from "2023-2024" or plain "2023"
  PRIOR.YEAR <- as.integer(sub("-.*", "", arg_str))
} else {
  PRIOR.YEAR <- 2023
}

# Resolve raw CSV path once
init.paths <- get_paths(year = PRIOR.YEAR, with.2.5kg = TRUE)
raw.csv.path <- init.paths$raw.csv(PRIOR.YEAR)

cat(sprintf("\n==================================================\n"))
cat(sprintf(" Processing Prior Data for Year: %d\n", PRIOR.YEAR))
cat(sprintf(" Raw Data Source: %s\n", raw.csv.path))
cat(sprintf("==================================================\n\n"))

# ----- Step 2: Read Raw Data (Once for efficiency) -----
if (!file.exists(raw.csv.path)) {
  stop(sprintf("Raw data file not found: %s", raw.csv.path))
}

natalitydata <- read.csv(raw.csv.path)
all.weights <- natalitydata$dbwt

paths <- get_paths(year = PRIOR.YEAR, with.2.5kg = PROCESS.WITH.2.5KG)
ensure_dirs(paths)

# REVISION R1/R2: the duplicate `natalitydata <- read.csv(paths$raw.csv())`
# that used to sit here has been REMOVED. It re-read a 2.1 GB file that had
# already been read into `natalitydata` above, and `raw.csv()` now defaults to
# the TARGET year -- so after the R1 fix it would have silently pulled the
# 2024 file into a script whose entire job is to summarise the PRIOR year.
# Everything below uses `all.weights`, taken from the prior-year read above.

# REVISION R3: drop the NCHS "not stated" weight sentinel before computing
# any proportion or quantile. dbwt == 9999 was previously counted as a
# Normal (> 2500 g) birth, which biased prop.normal upward and therefore
# biased the Normal component of the informed Dirichlet prior.
n.dbwt.sentinel <- sum(all.weights == 9999, na.rm = TRUE)
all.weights <- all.weights[!is.na(all.weights) & all.weights != 9999]
cat(sprintf("Dropped %d records with dbwt = 9999 (weight not stated) from the prior.\n",
            n.dbwt.sentinel))

# Proportions & filter LBW weights
count.lbw <- sum(all.weights <= 2500, na.rm = TRUE)
count.normal <- sum(all.weights > 2500, na.rm = TRUE)
total.count <- length(all.weights)

prop.lbw <- count.lbw / total.count
prop.normal <- count.normal / total.count

cat(sprintf("Proportion of LBW (<= 2500g): %.4f\n", prop.lbw))
cat(sprintf("Proportion of Normal (> 2500g): %.4f\n", prop.normal))

dat.lte2.5 <- all.weights[all.weights <= 2500]
num.quantiles <- 10

# ----- Step 3: Loop Across Both Rebin Types (With & Without 2.5kg) -----
# REVISION R26: this script computes and SAVES only -- it no longer draws the
# prior figure. Both arms' priors are now shown together in one panel in
# visual.R, which reads the files written here. Keeping the plotting out of
# here means quantile.R does one job: process the full prior-year file and
# produce the cutpoints and alphavec.
for (PROCESS.WITH.2.5KG in c(TRUE, FALSE)) {

  model_label <- if (PROCESS.WITH.2.5KG) "Full Model (rebin)" else "LBW-Only Model (rebin_without_2.5kg)"
  cat(sprintf("\n--------------------------------------------------\n"))
  cat(sprintf(" Generating Outputs for: %s\n", model_label))
  cat(sprintf("--------------------------------------------------\n"))

  paths <- get_paths(year = PRIOR.YEAR, with.2.5kg = PROCESS.WITH.2.5KG)
  ensure_dirs(paths)

  # Calculate quantiles and bin boundaries
  cut.points <- c(
    min(dat.lte2.5, na.rm = TRUE),
    quantile(dat.lte2.5, probs = seq(1/num.quantiles, (num.quantiles-1)/num.quantiles, 1/num.quantiles), na.rm = TRUE),
    max(dat.lte2.5, na.rm = TRUE)
  )

  birth.bins <- cut(
    dat.lte2.5,
    breaks = cut.points,
    include.lowest = TRUE,
    right = TRUE,
    labels = FALSE
  )

  bin.counts <- table(factor(birth.bins, levels = 1:(length(cut.points) - 1)))

  # Build Prior Vector
  if (PROCESS.WITH.2.5KG) {
    lbw.props <- as.numeric(bin.counts) / sum(bin.counts)
    alphavec.lbw <- lbw.props * prop.lbw
    alphavec <- c(alphavec.lbw, prop.normal)
    category.labels <- c(paste0("Q", 1:10), "Normal")
  } else {
    alphavec <- rep(1/num.quantiles, num.quantiles)
    category.labels <- paste0("Q", 1:10)
  }

  alpha0 <- 1
  alphavec <- alpha0 * alphavec

  # Save Cutpoints & Priors to target directory
  cutpoints_file <- paths$quantile.cutpoints(PRIOR.YEAR)
  prior_file     <- paths$informed.prior(PRIOR.YEAR)

  save(cut.points, file = cutpoints_file)
  save(alphavec, file = prior_file)

  cat("Saved cutpoints to :", cutpoints_file, "\n")
  cat("Saved prior vector to:", prior_file, "\n")

}

cat("\nCompleted quantile generation for all model variants.\n\n")