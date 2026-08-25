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

# parse command line arguments for prior year
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  arg_str <- args[1]
  # Extracts starting year from "2023-2024" or plain "2023"
  PRIOR.YEAR <- as.integer(sub("-.*", "", arg_str))
} else {
  PRIOR.YEAR <- 2023
}

# resolve raw CSV path once
init.paths <- get_paths(year = PRIOR.YEAR, with.2.5kg = TRUE)
raw.csv.path <- init.paths$raw.csv(PRIOR.YEAR)

cat(sprintf("\n==================================================\n"))
cat(sprintf(" Processing Prior Data for Year: %d\n", PRIOR.YEAR))
cat(sprintf(" Raw Data Source: %s\n", raw.csv.path))
cat(sprintf("==================================================\n\n"))

# read raw natality data
if (!file.exists(raw.csv.path)) {
  stop(sprintf("Raw data file not found: %s", raw.csv.path))
}

natalitydata <- read.csv(raw.csv.path)
all.weights <- natalitydata$dbwt

paths <- get_paths(year = PRIOR.YEAR, with.2.5kg = PROCESS.WITH.2.5KG)
ensure_dirs(paths)

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

for (PROCESS.WITH.2.5KG in c(TRUE, FALSE)) {

  model_label <- if (PROCESS.WITH.2.5KG) "Full Model (rebin)" else "LBW-Only Model (rebin_without_2.5kg)"
  cat(sprintf("\n--------------------------------------------------\n"))
  cat(sprintf(" Generating Outputs for: %s\n", model_label))
  cat(sprintf("--------------------------------------------------\n"))

  paths <- get_paths(year = PRIOR.YEAR, with.2.5kg = PROCESS.WITH.2.5KG)
  ensure_dirs(paths)

  # calculate quantiles and bin boundaries
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