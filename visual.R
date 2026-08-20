# ==============================================================================
# visual.R
# Visualization pipeline for Dirichlet-Multinomial CART robustness outputs
# ==============================================================================
rm(list = ls())
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ------------------------------------------------------------------------------
# 1. INITIALIZE CONFIG & PATHS VIA UTIL.R
# ------------------------------------------------------------------------------
# Source util.R from root directory
util_path <- file.path(rprojroot::find_root(rprojroot::is_git_root), "util.R")
if (file.exists(util_path)) {
  source(util_path)
} else {
  stop("Could not locate `util.R` in project root.")
}

# REVISION R24: cohort and model arm come from the command line like the other
# three scripts, instead of two hardcoded constants that silently rendered only
# one of the four model/cohort combinations.
#   Rscript visual.R 2023-2024            # full model (default)
#   Rscript visual.R 2023-2024 lbw        # LBW-only arm
.args      <- commandArgs(trailingOnly = TRUE)
PRIOR.YEAR <- if (length(.args) > 0) as.integer(sub("-.*", "", .args[1])) else 2023
WITH_2_5KG <- !(length(.args) > 1 && tolower(.args[2]) %in% c("lbw", "lbwmodel", "false"))
.p0        <- get_paths(year = PRIOR.YEAR, with.2.5kg = WITH_2_5KG)
MODEL_YEAR <- .p0$target.year
cat(sprintf("visual.R: cohort %s | target %d | %s\n", .p0$cohort, MODEL_YEAR,
            if (WITH_2_5KG) "full model" else "LBW-only model"))


# MODEL_YEAR <- if (exists("YEAR", envir = .GlobalEnv)) get("YEAR", envir = .GlobalEnv) else 2024
# WITH_2_5KG <- if (exists("PROCESS.WITH.2.5KG", envir = .GlobalEnv)) get("PROCESS.WITH.2.5KG", envir = .GlobalEnv) else TRUE

paths      <- get_paths(year = MODEL_YEAR, with.2.5kg = WITH_2_5KG)
ensure_dirs(paths)

# REVISION R29: two destinations. `article.dir` holds only the figures meant
# for the manuscript; `output.dir` stays the diagnostic dump.
output.dir  <- paths$plots
article.dir <- paths$article
dir.create(article.dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 2. LOAD CUTPOINTS DATA
# ------------------------------------------------------------------------------
# Attempt to load cutpoints using util.R accessor (check prior year first, then current year)
cutpoints.file <- paths$quantile.cutpoints(paths$prior.year)

if (!file.exists(cutpoints.file)) {
  cutpoints.file <- paths$quantile.cutpoints(paths$year)
}

if (file.exists(cutpoints.file)) {
  load(cutpoints.file)
} else {
  stop(sprintf("Cutpoints file not found at: %s", cutpoints.file))
}

# Standardize variable name (handles 'cut.points' vs 'cutpoints')
if (!exists("cutpoints") && exists("cut.points")) {
  cutpoints <- cut.points
}

if (!exists("cutpoints") || is.null(cutpoints)) {
  stop("Neither `cutpoints` nor `cut.points` vector was found in the loaded .RData file.")
}
# ------------------------------------------------------------------------------
# 3. DYNAMIC BW CATEGORY LABELS BUILDER (FIXED)
# ------------------------------------------------------------------------------
format_bin_labels <- function(cuts, include_nbw = TRUE) {
  n.cuts <- length(cuts)
  n.bins <- n.cuts - 1
  labels <- character(n.bins)

  for (i in seq_len(n.bins)) {
    lower <- round(cuts[i])
    upper <- cuts[i + 1]

    if (is.infinite(upper)) {
      labels[i] <- sprintf("K%d: > %d g", i, lower)
    } else {
      labels[i] <- sprintf("K%d: %d–%d g", i, lower, round(upper))
    }
  }

  # Append Normal Birth Weight (NBW) category if cutpoints end at 2500g
  if (include_nbw && !is.infinite(cuts[n.cuts])) {
    nbw_label <- sprintf("K%d: > %d g", n.bins + 1, round(cuts[n.cuts]))
    labels <- c(labels, nbw_label)
  }

  return(labels)
}

# Derive bin labels matching WITH_2_5KG flag
cat.labels      <- format_bin_labels(cutpoints, include_nbw = WITH_2_5KG)
cat.labels        <- sub(":.*", "", cat.labels)
legend.txt.full <- paste(cat.labels, collapse = "\n")
legend.txt.lbw  <- paste(cat.labels[1:10], collapse = "\n")


# ------------------------------------------------------------------------------
# 1. INFORMED PRIOR PLOT (ALPHAVEC) -- BOTH MODEL ARMS IN ONE FIGURE
# ------------------------------------------------------------------------------
# REVISION R26: quantile.R used to draw one base-graphics PDF per model arm, so
# the two priors that are meant to be compared lived in separate files, and the
# axis read "Q1 ... Q10" with the birth-weight ranges defined only in a table
# elsewhere. quantile.R now computes and saves only; this section loads BOTH
# arms' saved priors and draws them together, with the legend giving each
# quantile category's definition so the figure stands on its own.
prior.arms <- list(
  "Full Model (K = 11)"     = get_paths(year = PRIOR.YEAR, with.2.5kg = TRUE),
  "LBW-Only Model (K = 10)" = get_paths(year = PRIOR.YEAR, with.2.5kg = FALSE)
)

all.cats <- c(paste0("Q", 1:10), "NBW")
prior.df <- do.call(rbind, lapply(names(prior.arms), function(arm) {
  pp <- prior.arms[[arm]]
  f  <- pp$informed.prior(pp$prior.year)
  if (!file.exists(f))
    stop(sprintf("Missing %s\n  Run:  Rscript quantile.R %d-%d", f, pp$prior.year, pp$target.year))
  e <- new.env(); load(f, envir = e)
  a <- e$alphavec
  # REVISION R28: the NBW category is now drawn on the full-model panel. Its
  # length drives this: the full model returns 11 alphas, the LBW-only model
  # 10, so `all.cats[seq_along(a)]` gives each arm exactly its own categories
  # and the LBW-only panel simply has no bar in the NBW slot -- the honest
  # depiction of K = 10 against K = 11. NBW holds ~91% of the full model's
  # prior mass, so on a linear axis the ten deciles are compressed; every bar
  # carries its percentage, so they stay readable as numbers.
  data.frame(model    = arm,
             Category = factor(all.cats[seq_along(a)], levels = all.cats),
             alpha    = a,
             stringsAsFactors = FALSE)
}))
prior.df$model <- factor(prior.df$model, levels = names(prior.arms))
# alpha0 = 1, so alpha_k IS the share of prior mass; percent is just x100.
prior.df$pct <- 100 * prior.df$alpha

# Legend entries ARE the category definitions, from the prior year's cutpoints.
cp.kg <- cutpoints / 1000
# R28: the 11th entry is the NBW label the palette above already allows for.
interval.labels <- c(sprintf("Q%-2d %.2f-%.2f", 1:10, cp.kg[1:10], cp.kg[2:11]),
                     sprintf("NBW > %.2f", cp.kg[length(cp.kg)]))
prior.df$Interval <- factor(interval.labels[as.integer(prior.df$Category)],
                            levels = interval.labels)

nbw.alpha <- local({
  pp <- prior.arms[["Full Model (K = 11)"]]
  e <- new.env(); load(pp$informed.prior(pp$prior.year), envir = e); e$alphavec[11]
})

# REVISION R27: panels stacked TOP/BOTTOM rather than side by side, and the
# legend moved beneath them in two rows. Stacking shares one x axis between the
# arms, so the ten deciles line up vertically and the two priors are read
# against each other directly; and a wide, short legend under a tall figure
# wastes far less space than a tall legend beside a wide one. Legend entries are
# trimmed to "Q1  0.23-1.19" with the unit stated once in the legend title --
# the figure needs to say what a category IS, not repeat "kg" ten times.
# Define the colors: 10 colors for LBW (Red to Yellow/Light Green) + 1 distinct color for NBW (Deep Green)
# Adjust the number '10' below if your number of LBW intervals changes.
# Define a blue gradient based on the exact number of labels in your intervals
# (This ensures it smoothly covers all LBW deciles + the NBW category)
prior_colors <- colorRampPalette(c("#DBEAFE", "#1E3A8A"))(length(interval.labels))

p_prior <- ggplot(prior.df, aes(x = Category, y = alpha, fill = Interval)) +
  geom_col(colour = "white", linewidth = 0.35, width = 0.82) +
  # Slightly increased the percentage label size on the bars
  geom_text(aes(label = sprintf("%.2f%%", pct)),
            angle = 90, hjust = -0.12, size = 3.2, colour = "#374151") + 
  facet_wrap(~ model, ncol = 1, scales = "free_y") +
  
  scale_fill_manual(values = setNames(prior_colors, interval.labels),
                    name = "Birth Weight Category") +
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.30))) +
  guides(fill = guide_legend(ncol = 2, byrow = TRUE, title.position = "top",
                             keyheight = unit(0.45, "cm"))) + # Bumped keyheight
  
  labs(x = NULL, y = expression(paste("Prior hyperparameter  ", alpha[k]))) +
  theme_minimal(base_size = 11) + # Increased base size
  theme(
    # INCREASED TEXT SIZES THROUGHOUT
    strip.text          = element_text(face = "bold", size = 11, hjust = 0,
                                       margin = margin(b = 4)),
    axis.title.y        = element_text(size = 10.5, colour = "#374151"),
    axis.text           = element_text(size = 9.5),
    
    panel.grid.major.x  = element_blank(),
    panel.grid.minor    = element_blank(),
    panel.grid.major.y  = element_line(colour = "#E5E7EB", linewidth = 0.35),
    panel.spacing       = unit(1.1, "lines"),
    
    legend.position       = "inside",
    legend.position.inside = c(0.01, 0.99),
    legend.justification  = c(0, 1),
    legend.background     = element_rect(fill = "white", colour = "#D1D5DB",
                                         linewidth = 0.25),
    
    # LEGEND SPECIFIC INCREASES
    legend.title          = element_text(face = "bold", size = 10.5), # Bumped to 10.5
    legend.text           = element_text(size = 9.5, family = "mono"), # Bumped to 9.5
    legend.key.size       = unit(0.45, "cm"), # Bumped to match larger text
    legend.margin         = margin(5, 6, 5, 6) # Expanded margins for breathing room
  )

ggsave(
  filename = file.path(article.dir, sprintf("informed_prior_%d.pdf", paths$prior.year)),
  plot     = p_prior,
  width    = 7.0,
  height   = 6.6,
  device   = grDevices::cairo_pdf
)

cat("Saved combined prior figure to:",
    file.path(article.dir, sprintf("informed_prior_%d.pdf", paths$prior.year)), "\n")
# ------------------------------------------------------------------------------
# 1b. TREE FIGURES -- ONE DEPTH-COMPARISON PANEL, ONE ANNOTATED FULL TREE
# ------------------------------------------------------------------------------
# REVISION R27: the thesis carried depth-2/3/4/5 as four separate full-page
# figures per model arm (eight in total), plus a full tree and a depth-3
# "display" tree. For an article that is six tree figures per arm to make two
# points. This block replaces them with:
#
#   depth_comparison_<year>.pdf  -- the four depths as small multiples, one
#                                   figure, showing how structure and variable
#                                   entry change with depth
#   full_tree_<year>.pdf         -- the fitted model, fully annotated
#
# Both are drawn from objects dm-cart.R already saved, so nothing is refitted.
library(rpart)

depth.file <- paths$tree.depth.results(paths$year)
tree.file  <- paths$dm.tree.rebin(paths$year)

if (file.exists(depth.file) && file.exists(tree.file)) {

  d.env <- new.env(); load(depth.file, envir = d.env)   # trees, tree.sum, ...
  t.env <- new.env(); load(tree.file,  envir = t.env)   # dm.tree
  depth.trees <- d.env$trees
  lbw.cols    <- 1:10
  # This arm's prior; the old global `alphavec` disappeared when the prior
  # section was rewritten to load both arms into local environments (R26).
  a.env <- new.env(); load(paths$informed.prior(paths$prior.year), envir = a.env)
  alphavec <- a.env$alphavec
  # Short display names for the panel figure only; the standalone full tree
  # keeps the real column names.
  short.vars <- c(sex = "sex", dmar = "marr", mager = "age", meduc = "educ",
                  precare1st = "care", cig_0 = "smoke", race = "race")

  # ---- (a) THE article tree figure: both arms, stacked ------------------
  # REVISION R29: the depth-2..5 small multiples were unreadable at article
  # size -- four trees on one page leaves each about 3 in wide, and the
  # deepest has 19 leaves. Replaced by the comparison the paper is actually
  # about: the fitted full model over the fitted LBW-only model, one above the
  # other, each getting the full page width.
  #
  # Leaf labels are `leaf.detail = "compact"` (posterior P(LBW) only). The
  # four-line label -- category, P(LBW), N, classes -- is what made the leaves
  # illegible at this size; N and the class count belong in a table, and the
  # posterior is the quantity the figure exists to show.
  # The leaf quantity has to differ by arm. In the LBW-only model every
  # category IS an LBW category, so summing columns 1:10 gives 1.000 at every
  # leaf -- the same degeneracy R17 fixed for the risk ranking, which would
  # have shown a tree of identical 1.000s. That arm reports the severity
  # score instead: the mass in the three smallest deciles.
  arms.for.fig <- list(
    "Full model (K = 11)"     = list(paths = get_paths(year = PRIOR.YEAR, with.2.5kg = TRUE),
                                     cols = 1:10, lab = "P(LBW)"),
    "LBW-only model (K = 10)" = list(paths = get_paths(year = PRIOR.YEAR, with.2.5kg = FALSE),
                                     cols = 1:3,  lab = "P(Q1-Q3 | LBW)")
  )
  tree.files <- vapply(arms.for.fig, function(z) z$paths$dm.tree.rebin(z$paths$target.year),
                       character(1))

  if (all(file.exists(tree.files))) {
    models.file <- file.path(article.dir, sprintf("model_trees_%d.pdf", paths$year))
    # REVISION R32: now that this is a figure* at \textwidth, the type no
    # longer has to fight a 3x down-scale, so pointsize comes back down to 12.
    # The canvas is also reshaped to a 12 x 8 (3:2) aspect: at \textwidth that
    # is about 4.7 in tall, i.e. half a page, rather than the near-square
    # 11.5 x 10.5 which filled one.
    grDevices::cairo_pdf(models.file, width = 12, height = 10.5, pointsize = 12)
    par(mfrow = c(2, 1), oma = c(2.4, 0, 2.6, 0))
    for (arm in names(arms.for.fig)) {
      z  <- arms.for.fig[[arm]]
      pp <- z$paths
      te <- new.env(); load(pp$dm.tree.rebin(pp$target.year), envir = te)
      ae <- new.env(); load(pp$informed.prior(pp$prior.year),  envir = ae)
      ft <- te$dm.tree
      draw.dm.tree(
        ft,
        title    = arm,
        subtitle = sprintf("%d terminal nodes | leaves give posterior %s",
                           sum(ft$frame$var == "<leaf>"), z$lab),
        lbw.cols = z$cols, alphavec = ae$alphavec,
        # branch.labels = "left": with seven race groups the complement set
        # runs to "Hispanic,Asian,AIAN,NHOPI,Multiracial", which collided with
        # neighbouring node labels in the dense middle of the full tree.
        # Printing only the left set halves the text and loses nothing, since
        # the right child is by definition everything else.
        leaf.detail = "compact", branch.labels = "left",
        title.cex = 1.20, subtitle.cex = 0.78, cex.max = 1.05
      )
    }
    mtext(sprintf("Fitted DM-CART trees, %d cohort", paths$year),
          side = 3, line = 0.8, outer = TRUE, cex = 1.20, font = 2)
    mtext(paste("Node labels give the split variable. Blue labels give the levels sent LEFT at a multi-level split;",
                "the right child is the complement."),
          side = 1, line = 1.0, outer = TRUE, cex = 0.72, col = "#4B5563")
    mtext("Binary splits are unlabelled: the left branch is always the reference level of Table 1.",
          side = 1, line = 0.1, outer = TRUE, cex = 0.72, col = "#4B5563")
    dev.off()
    cat("Saved two-model tree figure to:", models.file, "\n")
  }

  # ---- (b) diagnostics, NOT article figures -----------------------------
  # Kept because they are useful while working, written to plots/ so they
  # cannot be mistaken for manuscript figures. The depth sweep is reported in
  # the paper as a table (pred.stats.compare), not a figure.
  # REVISION R32: the depth sweep is an article figure again, written to
  # article/ and formatted to match the fitted-tree figure -- same 12 x 8
  # canvas, same pointsize, same title/footnote treatment -- so the two read
  # as a pair rather than as figures from different papers.
  panel.file <- file.path(article.dir,
                          sprintf("depth_growth_%s_%d.pdf",
                                  if (WITH_2_5KG) "full" else "lbw", paths$year))
  grDevices::cairo_pdf(panel.file, width = 12, height = 8, pointsize = 12)
  par(mfrow = c(2, 2), oma = c(2.4, 0, 3.0, 0))
  for (d in names(depth.trees)) {
    ft <- depth.trees[[d]]
    draw.dm.tree(
      ft,
      title    = sprintf("maxdepth = %s", d),
      subtitle = sprintf("%d terminal nodes | %d predictors used",
                         sum(ft$frame$var == "<leaf>"),
                         length(unique(as.character(ft$frame$var[ft$frame$var != "<leaf>"])))),
      lbw.cols = lbw.cols, alphavec = alphavec,
      leaf.detail = "none", branch.labels = "left",
      var.labels = short.vars, cex.max = 0.95,
      title.cex = 1.05, subtitle.cex = 0.70
    )
  }
  mtext(sprintf("Tree growth by maximum depth -- %s",
                if (WITH_2_5KG) "Full model" else "LBW-only model"),
        side = 3, line = 0.9, outer = TRUE, cex = 1.20, font = 2)
  # mtext(paste("Node labels abbreviated: marr = marital status, age = maternal age > 33, educ = high school completed,",
  #             "care = 1st-trimester care, smoke = any pre-pregnancy smoking."),
  #       side = 1, line = 1.0, outer = TRUE, cex = 0.72, col = "#4B5563")
  mtext("Blue labels give the levels sent LEFT; binary splits are unlabelled. Terminal-node values omitted at this scale.",
        side = 1, line = 0.1, outer = TRUE, cex = 0.72, col = "#4B5563")
  dev.off()
  cat("Saved depth-growth figure to:", panel.file, "\n")

  plot.dm.tree(
    t.env$dm.tree,
    file     = file.path(output.dir, sprintf("full_tree_%d.pdf", paths$year)),
    title    = sprintf("DM-CART fitted tree -- %s",
                       if (WITH_2_5KG) "full model" else "LBW-only model"),
    subtitle = sprintf("cohort %s | leaves report modal category, posterior P(LBW), births (N) and classes (cls)",
                       paths$cohort),
    lbw.cols = lbw.cols, alphavec = alphavec
  )

} else {
  message("Tree objects not found; skipping tree figures. Run dm-cart.R first.")
}


# ------------------------------------------------------------------------------
# 1c. RISK-PROFILE LANDSCAPE -- BOTH MODEL ARMS, TWO PANES
# ------------------------------------------------------------------------------
# REVISION R28: dm-cart.R used to draw this per arm, giving two separate
# single-pane figures. Both arms are loaded here and drawn as one two-pane
# figure, so the article carries one figure instead of two and the two risk
# landscapes can be compared directly.
#
# Every class is plotted at its posterior risk with its 95% percentile interval,
# ordered by risk. Classes below the support threshold are drawn hollow and grey
# rather than dropped: the reader can see how much of the class space the filter
# removes, and that the extremes were not simply the thinnest cells.
risk.arms <- list(
  "Full model"     = get_paths(year = PRIOR.YEAR, with.2.5kg = TRUE),
  "LBW-only model" = get_paths(year = PRIOR.YEAR, with.2.5kg = FALSE)
)

risk.files <- vapply(risk.arms, function(pp) pp$bootstrap.results(pp$target.year), character(1))
if (all(file.exists(risk.files))) {

  rank_pdf <- file.path(article.dir, sprintf("risk_profile_ranking_%d.pdf", paths$year))
  # REVISION R34: stacked top/bottom and sized for a single column, so the
  # figure fits the half-margin layout rather than needing the full width.
  grDevices::cairo_pdf(rank_pdf, width = 3.6, height = 6.4, pointsize = 9)
  par(mfrow = c(2, 1), mar = c(3.8, 4.2, 2.8, 0.8), oma = c(2.6, 0, 1.8, 0))

  for (arm in names(risk.arms)) {
    pp <- risk.arms[[arm]]
    e  <- new.env(); load(pp$bootstrap.results(pp$target.year), envir = e)
    rk <- e$risk.ranking[order(e$risk.ranking$p_lbw), ]
    btr <- e$bootstrap.tree.results
    metric   <- unique(e$risk_totals$RiskMetric)
    min.brth <- btr$min.class.births
    hi <- btr$profiles$empirical$high.index
    lo <- btr$profiles$empirical$low.index
    pos <- seq_len(nrow(rk))

    plot(pos, rk$p_lbw, type = "n", ylim = range(rk$lower, rk$upper),
         xlab = "", ylab = "", cex.axis = 1.05, las = 1)
    mtext(metric, side = 2, line = 2.8, cex = 0.95)
    grid(nx = NA, ny = NULL, col = "#E5E7EB", lty = "solid")
    # ineligible first, so the eligible points sit on top
    ie <- !rk$eligible
    segments(pos[ie], rk$lower[ie], pos[ie], rk$upper[ie], col = "#F3F4F6", lwd = 1)
    points(pos[ie], rk$p_lbw[ie], pch = 1, cex = 0.42, col = "#D1D5DB")
    segments(pos[!ie], rk$lower[!ie], pos[!ie], rk$upper[!ie], col = "#93C5FD", lwd = 1)
    points(pos[!ie], rk$p_lbw[!ie], pch = 19, cex = 0.62, col = "#1D4ED8")

    for (z in list(list(i = lo, col = "#3B82F6", side = 4),
                   list(i = hi, col = "#EF4444", side = 2))) {
      k <- match(z$i, rk$class)
      points(k, rk$p_lbw[k], pch = 21, bg = z$col, col = "black", cex = 1.9, lwd = 1.2)
      text(k, rk$p_lbw[k], sprintf("i=%d", z$i), pos = z$side, offset = 0.8,
           cex = 1.0, font = 2, col = z$col)
    }
    mtext("Class, ordered by posterior risk", side = 1, line = 2.3, cex = 0.9)
    title(arm, cex.main = 1.15, font.main = 2, line = 1.4)
    mtext(sprintf("%d of %d classes eligible", sum(rk$eligible), nrow(rk)),
          side = 3, line = 0.35, cex = 0.82, col = "#4B5563")
  }

  mtext("Posterior risk across all classes",
        side = 3, line = 0.2, outer = TRUE, cex = 1.05, font = 2)
  # Two lines: as one it overran the 7 in width at both edges.
  mtext(sprintf("Filled = at least %s observed births; hollow grey = below that threshold, excluded from the ranking.",
                format(min.brth, big.mark = ",")),
        side = 1, line = 0.9, outer = TRUE, cex = 0.68, col = "#4B5563")
  mtext("Bars are 95% percentile bootstrap intervals; red and blue mark the highest- and lowest-risk classes.",
        side = 1, line = 0.05, outer = TRUE, cex = 0.68, col = "#4B5563")
  dev.off()
  cat("Saved two-pane risk landscape to:", rank_pdf, "\n")

} else {
  message("Bootstrap results missing for one or both arms; skipping risk landscape.")
}


# ------------------------------------------------------------------------------
# 1d. RISK-CATEGORY PROFILE -- the two extreme classes, both arms
# ------------------------------------------------------------------------------
# REVISION R31: redrawn in ggplot to match the other article figures. The
# base-graphics version dm-cart.R writes is kept in plots/ as a diagnostic.
# Deliberately spare: no in-figure footnotes, one subtitle line, the profile
# on the panel strip. Everything else belongs in the caption.
if (all(file.exists(risk.files))) {

  rs <- do.call(rbind, lapply(names(risk.arms), function(arm) {
    pp <- risk.arms[[arm]]
    e  <- new.env(); load(pp$bootstrap.results(pp$target.year), envir = e)
    d  <- e$risk_summary
    if ("Definition" %in% names(d)) d <- d[d$Definition == "empirical", ]
    d  <- d[d$Category %in% paste0("Q", 1:10), ]
    d$arm <- arm
    d
  }))
  rs$Category <- factor(rs$Category, levels = paste0("Q", 1:10))
  rs$panel <- factor(
    sprintf("%s  --  %s (class %d)",
            rs$arm, sub("-Risk Subgroup", "-risk", rs$Subgroup), rs$ClassIndex),
    levels = unique(sprintf("%s  --  %s (class %d)",
                            rs$arm, sub("-Risk Subgroup", "-risk", rs$Subgroup), rs$ClassIndex)))
  rs$grp <- ifelse(grepl("High", rs$Subgroup), "Highest-risk class", "Lowest-risk class")

  # Stronger, print-safe contrast than the previous red/blue: the ColorBrewer
  # RdBu endpoints, which stay distinguishable in greyscale as well.
  risk.cols <- c("Highest-risk class" = "#B2182B", "Lowest-risk class" = "#2166AC")

  p_riskcat <- ggplot(rs, aes(x = Category, y = Mean, fill = grp)) +
    geom_col(width = 0.78, colour = "white", linewidth = 0.3) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.28,
                  linewidth = 0.4, colour = "#1F2937") +
    facet_wrap(~ panel, ncol = 2, scales = "free_y") +
    scale_fill_manual(values = risk.cols, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
    labs(
      title    = "Birth-weight category probabilities for the extreme-risk classes",
      subtitle = "Posterior means with 95% percentile bootstrap intervals; the normal-birth-weight category is excluded.",
      x = "Birth-weight decile", y = "Probability"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title          = element_text(face = "bold", size = 12, colour = "#111827"),
      plot.subtitle       = element_text(size = 8.6, colour = "#4B5563", margin = margin(b = 8)),
      plot.title.position = "plot",
      strip.text          = element_text(face = "bold", size = 9, hjust = 0,
                                         margin = margin(b = 4)),
      axis.title          = element_text(size = 9, colour = "#374151"),
      axis.text           = element_text(size = 8.2),
      panel.grid.major.x  = element_blank(),
      panel.grid.minor    = element_blank(),
      panel.grid.major.y  = element_line(colour = "#E5E7EB", linewidth = 0.35),
      panel.spacing       = unit(1.2, "lines"),
      legend.position     = "none"
    )

  riskcat.file <- file.path(article.dir, sprintf("risk_category_profiles_%d.pdf", paths$year))
  ggsave(riskcat.file, p_riskcat, width = 9.0, height = 5.6, device = grDevices::cairo_pdf)
  cat("Saved risk-category profile figure to:", riskcat.file, "\n")
}


# ------------------------------------------------------------------------------
# 1e. COMMONALITY TABLE -- what the extreme profiles share
# ------------------------------------------------------------------------------
# REVISION R32: the shared-feature analysis is emitted as a LaTeX fragment as
# well as a figure. In an article a 7-row table states "how many of the ten
# carry this level" far more compactly than a shaded grid, and it fits a single
# column. \input{} the fragment; it needs only booktabs, already in main.tex.
#
# NOTE ON ESCAPING: the template uses ~BS~ for a backslash and ~RE~ for a LaTeX
# row terminator, substituted at the end. Writing the commands directly means
# four backslashes for "\begin" and eight for a row end, which is unreadable
# and came out doubled on the first attempt.
if (all(file.exists(risk.files))) {

  rows <- list()
  for (arm in names(risk.arms)) {
    pp <- risk.arms[[arm]]
    e  <- new.env(); load(pp$bootstrap.results(pp$target.year), envir = e)
    hi <- e$top.consensus; lo <- e$bottom.consensus
    for (k in seq_len(nrow(hi))) {
      rows[[length(rows) + 1]] <- data.frame(
        arm      = arm,
        variable = gsub("_", "~BS~_", hi$variable[k], fixed = TRUE),
        hi.level = hi$modal[k], hi.n = hi$n[k], hi.of = hi$of[k],
        lo.level = lo$modal[k], lo.n = lo$n[k], lo.of = lo$of[k],
        stringsAsFactors = FALSE)
    }
  }
  comm <- do.call(rbind, rows)
  bold.if <- function(txt, n, of) ifelse(n == of, sprintf("~BS~textbf{%s}", txt), txt)

  # REVISION R33: three columns, not six. The previous layout carried an "Arm"
  # column plus separate Level and n columns for each extreme, which is far
  # wider than the ~3.4 in of a two-column page -- it was running into the
  # neighbouring table. The arm now heads a spanning subrow, and level and count
  # are combined into one cell ("0 (10/10)"), which fits a single column.
  tmpl <- c(
    "% Auto-generated by visual.R -- do not edit by hand.",
    "~BS~begin{table}[tb]",
    "~BS~centering~BS~footnotesize",
    paste0("~BS~caption[Shared features of the extreme-risk classes]{Shared features of ",
           "the ten highest- and ten lowest-risk classes. Each cell gives the level ",
           "carried by most of the ten and how many carry it; ~BS~textbf{bold} marks a ",
           "level shared by all ten. Classes with fewer than 1,000 observed births are ",
           "excluded.}"),
    "~BS~label{tab:commonality}",
    "~BS~begin{tabular}{@{}lcc@{}}",
    "~BS~toprule",
    "~BS~textbf{Predictor} & ~BS~textbf{Highest 10} & ~BS~textbf{Lowest 10} ~RE~",
    "~BS~midrule")
  for (a in unique(comm$arm)) {
    d <- comm[comm$arm == a, ]
    tmpl <- c(tmpl,
      sprintf("~BS~multicolumn{3}{@{}l@{}}{~BS~textit{%s}} ~RE~", a))
    for (k in seq_len(nrow(d))) {
      tmpl <- c(tmpl, sprintf("~BS~quad ~BS~texttt{%s} & %s (%d/%d) & %s (%d/%d) ~RE~",
        d$variable[k],
        bold.if(d$hi.level[k], d$hi.n[k], d$hi.of[k]), d$hi.n[k], d$hi.of[k],
        bold.if(d$lo.level[k], d$lo.n[k], d$lo.of[k]), d$lo.n[k], d$lo.of[k]))
    }
    if (a != tail(unique(comm$arm), 1)) tmpl <- c(tmpl, "~BS~addlinespace")
  }
  tmpl <- c(tmpl, "~BS~bottomrule", "~BS~end{tabular}", "~BS~end{table}")

  tex <- gsub("~RE~", "SLASHSLASH", tmpl, fixed = TRUE)
  tex <- gsub("~BS~", "SLASH", tex, fixed = TRUE)
  tex <- gsub("SLASHSLASH", "BACKBACK", tex, fixed = TRUE)
  tex <- gsub("SLASH", "BACK", tex, fixed = TRUE)
  tex <- gsub("BACKBACK", strrep(intToUtf8(92), 2), tex, fixed = TRUE)
  tex <- gsub("BACK", intToUtf8(92), tex, fixed = TRUE)

  comm.file <- file.path(article.dir, sprintf("table_commonality_%d.tex", paths$year))
  writeLines(tex, comm.file)
  cat("Saved commonality table to:", comm.file, "\n")
}

# ------------------------------------------------------------------------------
# 2: PREDICTOR VARIABLE SPLIT DEPTH DISTRIBUTIONS
# ------------------------------------------------------------------------------
# REVISION R24: depth_df is written by dm-cart.R into the BOOTSTRAP results
# file (it is per-replicate split data), not the depth-comparison file. The old
# code looked in tree.depth.results, found nothing, and fell through to
# `sample(0:6, 5000, replace = TRUE, prob = runif(7))` -- i.e. this figure was
# random noise with no warning. The fallback is deleted, not commented: if the
# data is not there the script stops and says how to produce it.
boot.file <- paths$bootstrap.results(paths$year)
if (!file.exists(boot.file))
  stop(sprintf("Missing %s\n  Run:  DMCART_B=10000 Rscript dm-cart.R %d-%d",
               boot.file, paths$prior.year, paths$target.year))
boot_env <- new.env()
load(boot.file, envir = boot_env)          # into an env, never the global
if (!exists("depth_df", envir = boot_env))
  stop(sprintf("`depth_df` not found in %s -- re-run dm-cart.R (REVISION R11).", boot.file))
depth_df <- get("depth_df", envir = boot_env)

mean_depths <- depth_df %>%
  group_by(Variable) %>%
  summarize(mean_depth = mean(Depth, na.rm = TRUE), .groups = "drop")

p_depths <- ggplot(depth_df, aes(x = Depth, fill = Variable)) +
  # R24: dropped the stale `align = "center"` argument (unknown to geom_histogram)
  geom_histogram(binwidth = 0.8, color = "white") +
  geom_vline(
    data        = mean_depths,
    aes(xintercept = mean_depth),
    linetype    = "dashed",
    color       = "black",
    linewidth   = 0.6
  ) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 3) +
  labs(
    title = "Distribution of Depths for Each Predictor Variable",
    x     = "Depth",
    y     = "Number of Splits"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(size = 14, face = "plain", hjust = 0),
    strip.background = element_blank(),
    strip.text       = element_text(face = "plain"),
    panel.grid.minor = element_blank(),
    legend.position  = "right"
  )

ggsave(
  filename = file.path(output.dir, sprintf("depth_distributions_%d.png", paths$year)),
  plot     = p_depths,
  width    = 10,
  height   = 8,
  dpi      = 300
)


# ------------------------------------------------------------------------------
# 6. PLOT 3 & 4: HIGH VS LOW RISK SUBGROUPS
# ------------------------------------------------------------------------------
# REVISION R24: the fallback that used to sit here fabricated the headline
# numbers -- c(rep(0.016, 10), 0.840) and c(rep(0.009, 10), 0.910) -- whenever
# `risk_summary` was absent, which was ALWAYS, because nothing wrote it before
# REVISION R12. It is deleted. Missing input is now an error.
#
# Intervals are the PERCENTILE bootstrap bounds carried in Lower/Upper, not a
# Mean +/- 1.96*SE reconstruction; the old approach gave this figure a different
# interval from the one dm-cart.R reports for the same quantity.
if (!exists("risk_summary", envir = boot_env))
  stop(sprintf("`risk_summary` not found in %s -- re-run dm-cart.R (REVISION R12/R17).", boot.file))
risk_summary <- get("risk_summary", envir = boot_env)

# R17 added a Definition column ("empirical" / "a priori"); plot the empirical
# extremes, falling back to whatever is present for older result files.
if ("Definition" %in% names(risk_summary) && any(risk_summary$Definition == "empirical"))
  risk_summary <- risk_summary[risk_summary$Definition == "empirical", ]

high_sub <- risk_summary %>% filter(Subgroup == "High-Risk Subgroup")
low_sub  <- risk_summary %>% filter(Subgroup == "Low-Risk Subgroup")

high.risk.probs <- high_sub$Mean
high.risk.lwr   <- high_sub$Lower
high.risk.upr   <- high_sub$Upper

low.risk.probs  <- low_sub$Mean
low.risk.lwr    <- low_sub$Lower
low.risk.upr    <- low_sub$Upper

# Explicitly omit NBW category (trims 11 categories to C1–C10)
n_cats <- min(10, length(high.risk.probs))

high.risk.probs <- high.risk.probs[1:n_cats]
high.risk.lwr   <- high.risk.lwr[1:n_cats]
high.risk.upr   <- high.risk.upr[1:n_cats]

low.risk.probs  <- low.risk.probs[1:n_cats]
low.risk.lwr    <- low.risk.lwr[1:n_cats]
low.risk.upr    <- low.risk.upr[1:n_cats]

cat.labels.C    <- paste0("C", 1:n_cats)

# Target output PDF path
pdf_out_path <- file.path(output.dir, sprintf("high_low_risk_pred_%d.pdf", paths$year))

pdf(file = pdf_out_path, width = 8.5, height = 8.5)

# Expanded bottom outer margin (oma = c(4.5, 0, 4, 0)) to hold two-line caption
par(mfrow = c(2, 1), 
    mar = c(3, 4.5, 3, 2),    
    oma = c(4.5, 0, 4, 0))      

# Dynamic Y-axis limits for LBW subset
y_max_high <- max(high.risk.upr, na.rm = TRUE) * 1.15
y_max_low  <- max(low.risk.upr, na.rm = TRUE) * 1.15

# TOP PANEL: High Risk Subgroup
high.bars <- barplot(high.risk.probs, 
                     main = "",
                     xlab = "", 
                     ylab = "Probability",
                     ylim = c(0, y_max_high),
                     names.arg = cat.labels.C,  
                     col = "red",
                     cex.names = 1.0,
                     cex.axis = 1.0,
                     border = "darkred")

title("High Risk Subgroup", line = 1, cex.main = 1.2)

arrows(high.bars, high.risk.lwr, 
       high.bars, high.risk.upr, 
       angle = 90, code = 3, length = 0.05, 
       col = "black", lwd = 1.5)

# BOTTOM PANEL: Low Risk Subgroup
low.bars <- barplot(low.risk.probs,
                    main = "", 
                    xlab = "", 
                    ylab = "Probability",
                    ylim = c(0, y_max_low),
                    names.arg = cat.labels.C,
                    col = "blue",
                    cex.names = 1.0,
                    cex.axis = 1.0,
                    border = "darkblue")

title("Low Risk Subgroup", line = 1, cex.main = 1.2)

arrows(low.bars, low.risk.lwr, 
       low.bars, low.risk.upr, 
       angle = 90, code = 3, length = 0.05, 
       col = "black", lwd = 1.5)

# MAIN TITLES
mtext(sprintf("Birth Weight Category Probabilities (%d)", paths$year), 
      side = 3, line = 2, outer = TRUE, cex = 1.4, font = 2)

mtext("Predicted probabilities with 95% confidence intervals", 
      side = 3, line = 0.5, outer = TRUE, cex = 1.1)

# TWO-LINE CAPTION AT BOTTOM
mtext("Low Risk Factors: Married, Non-smoker, White, 1st Trimester Prenatal Care", 
      side = 1, line = 1.2, outer = TRUE, cex = 0.85, col = "black")

mtext("High Risk Factors: Unmarried, Smoker, Delayed/No Prenatal Care", 
      side = 1, line = 2.6, outer = TRUE, cex = 0.85, col = "black")

dev.off()
# REVISION R24: `par(mfrow = c(1, 1))` removed from here. par() with no open
# device OPENS one, so this line created a stray Rplots.pdf on every batch run
# (same bug as REVISION R9 in dm-cart.R); the settings it reset died with the
# device closed on the line above.
