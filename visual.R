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

# cohort and model arm come from the command line like the other
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
ARTICLE.FIG.WIDTH <- 7.0

# two destinations. `article.dir` holds only the figures meant
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
# quantile.R used to draw one base-graphics PDF per model arm, so
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
  data.frame(model    = arm,
             Category = factor(all.cats[seq_along(a)], levels = all.cats),
             alpha    = a,
             stringsAsFactors = FALSE)
}))
prior.df$model <- factor(prior.df$model, levels = names(prior.arms))
# alpha0 = 1, so alpha_k IS the share of prior mass; percent is just x100.
prior.df$pct <- 100 * prior.df$alpha
cp.kg <- cutpoints / 1000
# the 11th entry is the NBW label the palette above already allows for.
interval.labels <- c(sprintf("Q%-2d %.2f-%.2f", 1:10, cp.kg[1:10], cp.kg[2:11]),
                     sprintf("NBW > %.2f", cp.kg[length(cp.kg)]))
prior.df$Interval <- factor(interval.labels[as.integer(prior.df$Category)],
                            levels = interval.labels)

nbw.alpha <- local({
  pp <- prior.arms[["Full Model (K = 11)"]]
  e <- new.env(); load(pp$informed.prior(pp$prior.year), envir = e); e$alphavec[11]
})

prior_colors <- colorRampPalette(c("#DBEAFE", "#1E3A8A"))(length(interval.labels))

p_prior <- ggplot(prior.df, aes(x = Category, y = alpha, fill = Interval)) +
  geom_col(color = "white", linewidth = 0.35, width = 0.82) +
  # Slightly increased the percentage label size on the bars
  geom_text(aes(label = sprintf("%.2f%%", pct)),
            angle = 90, hjust = -0.12, size = 2.3, color = "#374151") + 
  facet_wrap(~ model, ncol = 1, scales = "free_y") +
  
  scale_fill_manual(values = setNames(prior_colors, interval.labels),
                    name = "Birth Weight Category") +
  
  # Headroom for the rotated bar labels. The NBW bar reaches 0.913 of the
  # axis, and its label is the longest, so this is what stops ggplot
  # clipping the trailing '%' at the panel edge.
  scale_y_continuous(expand = expansion(mult = c(0, 0.58))) +
  guides(fill = guide_legend(ncol = 2, byrow = TRUE, title.position = "top",
                             keyheight = unit(0.26, "cm"),
                             keywidth  = unit(0.26, "cm"))) +
  
  labs(x = NULL, y = expression(paste("Prior hyperparameter  ", alpha[k]))) +
  theme_minimal(base_size = 8.5) +
  theme(
    # Sized against the 8 pt caption, not the 10 pt body: a figure whose
    # labels match the running text reads as oversized.
    strip.text          = element_text(face = "bold", size = 8.5, hjust = 0,
                                       margin = margin(b = 3)),
    axis.title.y        = element_text(size = 8, color = "#374151"),
    axis.text           = element_text(size = 7.2),
    
    panel.grid.major.x  = element_blank(),
    panel.grid.minor    = element_blank(),
    panel.grid.major.y  = element_line(color = "#E5E7EB", linewidth = 0.35),
    panel.spacing       = unit(1.1, "lines"),
    
    legend.position       = "inside",
    legend.position.inside = c(0.01, 0.99),
    legend.justification  = c(0, 1),
    legend.background     = element_rect(fill = "white", color = "#D1D5DB",
                                         linewidth = 0.25),
  
    # LEGEND SPECIFIC INCREASES
    legend.title          = element_text(face = "bold", size = 7.6),
    legend.text           = element_text(size = 6.6, family = "mono"),
    legend.key.size       = unit(0.26, "cm"),
    legend.margin         = margin(3, 4, 3, 4),
    legend.spacing.y      = unit(0.02, "cm")
  )

ggsave(
  filename = file.path(article.dir, sprintf("informed_prior_%d.pdf", paths$prior.year)),
  plot     = p_prior,
  width    = ARTICLE.FIG.WIDTH,
  height   = 3.6,
  device   = grDevices::cairo_pdf
)

cat("Saved combined prior figure to:",
    file.path(article.dir, sprintf("informed_prior_%d.pdf", paths$prior.year)), "\n")




# ------------------------------------------------------------------------------
# 1b. TREE FIGURES -- ONE DEPTH-COMPARISON PANEL, ONE ANNOTATED FULL TREE
# ------------------------------------------------------------------------------
# output here:
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
  a.env <- new.env(); load(paths$informed.prior(paths$prior.year), envir = a.env)
  alphavec <- a.env$alphavec
  # Short display names for the panel figure only
  short.vars <- c(sex = "sex", dmar = "marital", mager = "age", meduc = "educ",
                  precare1st = "care", cig_0 = "smoking", race = "race")
  # Single words instead for readability in full-tree figure
  plain.vars <- c(sex = "Sex", dmar = "Marital", mager = "Age",
                  meduc = "Education", precare1st = "Prenatal",
                  cig_0 = "Smoking", race = "Race")

  # ---- (a) THE article tree figure: both arms, stacked ------------------
  # Leaf labels are `leaf.detail = "compact"`
  arms.for.fig <- list(
    "Full model (K = 11)"     = list(paths = get_paths(year = PRIOR.YEAR, with.2.5kg = TRUE),
                                     cols = 1:10, lab = "P(LBW)"),
    "LBW-only model (K = 10)" = list(paths = get_paths(year = PRIOR.YEAR, with.2.5kg = FALSE),
                                     cols = 1:3,  lab = "P(Q1, Q_2, or Q3 | LBW)")
  )
  tree.files <- vapply(arms.for.fig, function(z) z$paths$dm.tree.rebin(z$paths$target.year),
                       character(1))

  if (all(file.exists(tree.files))) {
    models.file <- file.path(article.dir, sprintf("model_trees_%d.pdf", paths$year))
    grDevices::cairo_pdf(models.file, width = 14, height = 11, pointsize = 12)
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
        subtitle = sprintf("%d terminal nodes | leaves give posterior %s and the number of births",
                           sum(ft$frame$var == "<leaf>"), z$lab),
        lbw.cols = z$cols, alphavec = ae$alphavec,
        # branch.labels = "left": with seven race groups the complement set
        # runs to "Hispanic,Asian,AIAN,NHOPI,Multiracial", which is too long to fit in the figure. 
        leaf.detail = "compact", branch.labels = "left",
        var.labels = plain.vars,
        title.cex = 1.20, subtitle.cex = 0.95, cex.max = 1.05
      )
    }
    mtext(sprintf("Fitted DM-CART trees, %d cohort", paths$year),
          side = 3, line = 0.8, outer = TRUE, cex = 1.20, font = 2)
    mtext(paste("Node labels give the split variable. Blue labels give the levels sent LEFT at a multi-level split;",
                "the right child is the complement."),
          side = 1, line = 1.0, outer = TRUE, cex = 0.85, col = "#4B5563")
    mtext("Binary splits are unlabelled: the left branch is always the reference level of Table 1.",
          side = 1, line = 0.1, outer = TRUE, cex = 0.85, col = "#4B5563")
    dev.off()
    cat("Saved two-model tree figure to:", models.file, "\n")
  }

  # ---- (b) diagnostics -----------------------------
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
# dm-cart.R used to draw this per arm, giving two separate
# single-pane figures. Both arms are loaded here and drawn as one two-pane
# figure, so the article carries one figure instead of two and the two risk
# landscapes can be compared directly.
#
# Every class is plotted at its posterior risk with its 95% percentile interval,
# ordered by risk. Classes below the support threshold are drawn hollow and gray
# rather than dropped, so the reader can see how many classes were excluded from the ranking.
risk.arms <- list(
  "Full model"     = get_paths(year = PRIOR.YEAR, with.2.5kg = TRUE),
  "LBW-only model" = get_paths(year = PRIOR.YEAR, with.2.5kg = FALSE)
)

risk.files <- vapply(risk.arms, function(pp) pp$bootstrap.results(pp$target.year), character(1))
if (all(file.exists(risk.files))) {

  rank_pdf <- file.path(article.dir, sprintf("risk_profile_ranking_%d.pdf", paths$year))
  grDevices::cairo_pdf(rank_pdf, width = ARTICLE.FIG.WIDTH, height = 7.2, pointsize = 9)
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

    for (z in list(list(i = lo, col = "#3B82F6", side = 1),
                   list(i = hi, col = "#EF4444", side = 2))) {
      k <- match(z$i, rk$class)
      points(k, rk$p_lbw[k], pch = 21, bg = z$col, col = "black", cex = 1.9, lwd = 1.2)
      text(k, rk$p_lbw[k], sprintf("i=%d", z$i), pos = z$side, offset = 1.2,
           cex = 1.0, font = 2, col = z$col)
    }
    mtext("Class, ordered by posterior risk", side = 1, line = 2.3, cex = 0.9)
    title(arm, cex.main = 1.15, font.main = 2, line = 1.4)
    mtext(sprintf("%d of %d classes eligible", sum(rk$eligible), nrow(rk)),
          side = 3, line = 0.35, cex = 0.82, col = "#4B5563")
  }
  mtext("Posterior risk across all classes",
        side = 3, line = 0.2, outer = TRUE, cex = 1.05, font = 2)
  # mtext(sprintf("Filled = at least %s observed births; hollow gray = below that threshold, excluded from the ranking.",
                # format(min.brth, big.mark = ",")),
        # side = 1, line = 0.9, outer = TRUE, cex = 0.85, col = "#4B5563")
  # mtext("Bars are 95% percentile bootstrap intervals; red and blue mark the highest- and lowest-risk classes.",
  #       side = 1, line = 0.05, outer = TRUE, cex = 0.85, col = "#4B5563")
  dev.off()
  cat("Saved two-pane risk landscape to:", rank_pdf, "\n")

} else {
  message("Bootstrap results missing for one or both arms; skipping risk landscape.")
}


# ------------------------------------------------------------------------------
# 1d. RISK-CATEGORY PROFILE -- the two extreme classes, both arms
# ------------------------------------------------------------------------------
# Four panels: highest- and lowest-risk class in each arm. Only the ten LBW
# deciles are plotted. In the full model NBW holds ~91% of the mass, which on a
# shared axis would flatten the deciles to nothing.
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
  # Panel order is fixed by taking levels in the order the rows were built,
  # so the two arms stay in the same left-to-right order as everywhere else.
  rs$panel <- factor(
    sprintf("%s  --  %s (class %d)",
            rs$arm, sub("-Risk Subgroup", "-risk", rs$Subgroup), rs$ClassIndex),
    levels = unique(sprintf("%s  --  %s (class %d)",
                            rs$arm, sub("-Risk Subgroup", "-risk", rs$Subgroup), rs$ClassIndex)))
  rs$grp <- ifelse(grepl("High", rs$Subgroup), "Highest-risk class", "Lowest-risk class")

  risk.cols <- c("Highest-risk class" = "#B2182B", "Lowest-risk class" = "#2166AC")

  p_riskcat <- ggplot(rs, aes(x = Category, y = Mean, fill = grp)) +
    geom_col(width = 0.78, color = "white", linewidth = 0.3) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.28,
                  linewidth = 0.4, color = "#1F2937") +
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
      plot.title          = element_text(face = "bold", size = 12, color = "#111827"),
      plot.subtitle       = element_text(size = 8.6, color = "#4B5563", margin = margin(b = 8)),
      plot.title.position = "plot",
      strip.text          = element_text(face = "bold", size = 9, hjust = 0,
                                         margin = margin(b = 4)),
      axis.title          = element_text(size = 9, color = "#374151"),
      axis.text           = element_text(size = 8.2),
      panel.grid.major.x  = element_blank(),
      panel.grid.minor    = element_blank(),
      panel.grid.major.y  = element_line(color = "#E5E7EB", linewidth = 0.35),
      panel.spacing       = unit(1.2, "lines"),
      legend.position     = "none"
    )

  riskcat.file <- file.path(article.dir, sprintf("risk_category_profiles_%d.pdf", paths$year))
  ggsave(riskcat.file, p_riskcat, width = ARTICLE.FIG.WIDTH, height = 4.35,
         device = grDevices::cairo_pdf)
  cat("Saved risk-category profile figure to:", riskcat.file, "\n")
}

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

  # ~BS~ stands for a backslash, ~RE~ for a LaTeX row terminator; both are
  # substituted at the end. Writing the backslashes literally would need four
  # of them for \begin and eight for a row end, which is unreadable and was
  # the source of a doubling bug.
  tmpl <- c(
    "% Auto-generated by visual.R -- do not edit by hand.",
    "~BS~begin{table}[!htbp]",
    "~BS~centering~BS~small",
    paste0("~BS~caption[Shared features of the extreme-risk classes]{Shared features of ",
           "the ten highest- and ten lowest-risk classes. Each cell gives the level ",
           "carried by most of the ten and how many carry it; ~BS~textbf{bold} marks a ",
           "level shared by all ten. Classes with fewer than 1,000 observed births are ",
           "excluded.~BS~label{tab:commonality}}"),
        "~BS~begin{tabular*}{~BS~textwidth}{@{~BS~extracolsep~BS~fill}lcc@{}}",
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
  tmpl <- c(tmpl, "~BS~bottomrule", "~BS~end{tabular*}", "~BS~end{table}")

  # Two-stage substitution: the row terminator has to become two backslashes
  # and a lone ~BS~ one, so both are routed through distinct placeholders
  # before either is turned into a real backslash.
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
# 1f. RESULTS TABLES -- variable selection and the depth sweep
# ------------------------------------------------------------------------------
# Written from the same saved ensemble that produced the figures, so a table
# can never quote a different run from the text around it. Same ~BS~ / ~RE~
# escaping as above.
if (all(file.exists(risk.files))) {

  unesc <- function(v) {
    v <- gsub("~RE~", "SLASHSLASH", v, fixed = TRUE)
    v <- gsub("~BS~", "SLASH", v, fixed = TRUE)
    v <- gsub("SLASHSLASH", "BACKBACK", v, fixed = TRUE)
    v <- gsub("SLASH", "BACK", v, fixed = TRUE)
    v <- gsub("BACKBACK", strrep(intToUtf8(92), 2), v, fixed = TRUE)
    gsub("BACK", intToUtf8(92), v, fixed = TRUE)
  }
  texvar <- function(v) gsub("_", "~BS~_", v, fixed = TRUE)

  ## ---- (a) variable selection frequency + mean split depth ----------------
  freq <- lapply(names(risk.arms), function(arm) {
    pp <- risk.arms[[arm]]
    e  <- new.env(); load(pp$bootstrap.results(pp$target.year), envir = e)
    b  <- e$bootstrap.tree.results
    dep <- aggregate(Depth ~ Variable, e$depth_df, mean)
    m <- merge(b$var.freq.df, dep, by.x = "variable", by.y = "Variable", all.x = TRUE)
    m$root <- b$top.var.df$frequency[match(m$variable, b$top.var.df$variable)]
    m$root[is.na(m$root)] <- 0
    m <- m[order(-m$frequency, m$Depth), ]
    list(arm = arm, m = m, B = b$B, nvar = b$n.vars.summary[["Median"]],
         spt = nrow(e$depth_df) / b$B)
  })
  vars.order <- c("sex", "dmar", "race", "mager", "precare1st", "cig_0", "meduc")
  A <- freq[[1]]; Bm <- freq[[2]]
  # A predictor absent from an arm's summary yields no row at all, not an NA
  # row, so guard on length as well as on NA.
  pull <- function(z, v, col) {
    x <- z$m[[col]][match(v, z$m$variable)]
    if (length(x) == 0 || is.na(x)) NA_real_ else x
  }
  # Proportions render as percentages; depths stay on their natural scale.
  fpct <- function(x) if (is.na(x)) "---" else sprintf("%.1f~BS~%%", 100 * x)
  fdep <- function(x) if (is.na(x)) "---" else sprintf("%.2f", x)

  # Three stacked blocks share one set of columns: root-split share, overall
  # selection frequency, then mean split depth.
  t1 <- c("% Auto-generated by visual.R -- do not edit by hand.",
    "~BS~begin{table}[!htbp]", "~BS~centering~BS~small",
    paste0("~BS~caption[Variable importance in the full and LBW-only models]{Variable ",
           "importance across the ", format(A$B, big.mark = ","), "-tree bootstrap ensemble. ",
           "~BS~emph{Initial split} is the percentage of replicates in which the predictor is ",
           "chosen at the root; ~BS~emph{frequency} the percentage in which it appears anywhere ",
           "in the tree; ~BS~emph{mean depth} its average split depth, with the root at zero. ",
           "A dash marks a predictor that is never selected.~BS~label{tab:varfreq}}"),
        "~BS~begin{tabular*}{~BS~textwidth}{@{~BS~extracolsep~BS~fill}lcc@{}}", "~BS~toprule",
    sprintf(" & %s & %s ~RE~", A$arm, Bm$arm),
    "~BS~cmidrule(lr){2-2} ~BS~cmidrule(lr){3-3}",
    "~BS~multicolumn{3}{@{}l}{~BS~textit{Initial split variable}} ~RE~", "~BS~midrule")
  for (v in vars.order) {
    ra <- pull(A, v, "root"); rb <- pull(Bm, v, "root")
    # Skip predictors that never reach the root in either arm; a block of
    # zeroes says nothing the surrounding rows do not.
    if ((is.na(ra) || ra == 0) && (is.na(rb) || rb == 0)) next
    t1 <- c(t1, sprintf("~BS~quad ~BS~texttt{%s} & %s & %s ~RE~", texvar(v), fpct(ra), fpct(rb)))
  }
  t1 <- c(t1, "~BS~midrule",
    "~BS~multicolumn{3}{@{}l}{~BS~textit{Variable frequency}} ~RE~", "~BS~midrule")
  for (v in vars.order)
    t1 <- c(t1, sprintf("~BS~quad ~BS~texttt{%s} & %s & %s ~RE~", texvar(v),
                        fpct(pull(A, v, "frequency")), fpct(pull(Bm, v, "frequency"))))
  t1 <- c(t1, "~BS~midrule",
    "~BS~multicolumn{3}{@{}l}{~BS~textit{Mean split depth}} ~RE~", "~BS~midrule")
  for (v in vars.order)
    t1 <- c(t1, sprintf("~BS~quad ~BS~texttt{%s} & %s & %s ~RE~", texvar(v),
                        fdep(pull(A, v, "Depth")), fdep(pull(Bm, v, "Depth"))))
  t1 <- c(t1, "~BS~midrule",
    sprintf("~BS~quad Median predictors per tree & %d & %d ~RE~", A$nvar, Bm$nvar),
    sprintf("~BS~quad Mean splits per tree & %.1f & %.1f ~RE~", A$spt, Bm$spt),
    "~BS~bottomrule", "~BS~end{tabular*}", "~BS~end{table}")
  writeLines(unesc(t1), file.path(article.dir, sprintf("table_varfreq_%d.tex", paths$year)))

  # ---- (b) depth sweep ----------------------------------------------------
  dfile <- paths$tree.depth.results(paths$year)
  if (file.exists(dfile)) {
    de <- new.env(); load(dfile, envir = de)
    ps <- de$pred.stats.compare
    t2 <- c("% Auto-generated by visual.R -- do not edit by hand.",
      "~BS~begin{table}[!htbp]", "~BS~centering~BS~small",
      paste0("~BS~caption[Depth-controlled refits]{Depth-controlled refits of the ",
             if (WITH_2_5KG) "full" else "LBW-only", " model. Posterior risk is ",
             "summarized across all classes at each maximum depth.}"),
      sprintf("~BS~label{tab:depth-%s}", if (WITH_2_5KG) "full" else "lbw"),
      "~BS~begin{tabular*}{~BS~textwidth}{@{~BS~extracolsep~BS~fill}rrrrrr@{}}", "~BS~toprule",
      "~BS~textbf{Depth} & ~BS~textbf{Leaves} & ~BS~textbf{Vars} & ~BS~textbf{Mean} & ~BS~textbf{Min} & ~BS~textbf{Max} ~RE~",
      "~BS~midrule")
    for (k in seq_len(nrow(ps)))
      t2 <- c(t2, sprintf("%d & %d & %d & %.4f & %.4f & %.4f ~RE~", ps$depth[k],
                          ps$n.terminal.nodes[k], ps$n.variables[k],
                          ps$mean.lbw.prob[k], ps$min.lbw.prob[k], ps$max.lbw.prob[k]))
    t2 <- c(t2, "~BS~bottomrule", "~BS~end{tabular*}", "~BS~end{table}")
    writeLines(unesc(t2), file.path(article.dir,
                 sprintf("table_depth_%s_%d.tex",
                         if (WITH_2_5KG) "full" else "lbw", paths$year)))
  }

  # ---- (c) predictors, definitions and missingness, in one table ----------
  # This used to be two floats: a hand-written predictor table in main.tex and a
  # generated record-flow table. They shared a key column and were read
  # together, so they are merged here. Generated rather than hand-written so the
  # missingness counts cannot drift from the run that produced them.
  #
  # `variable` in missing.report holds the RAW NCHS field names, which differ
  # from the analysis names for two predictors, so the mapping is explicit.
  defs <- data.frame(
    field = c("sex", "mager", "meduc", "precare1st", "cig_0", "dmar", "race"),
    nchs  = c("sex", "mager", "meduc", "precare5",   "cig_0", "dmar", "mracehisp"),
    lev   = c("2", "2", "2", "2", "2", "3", "7"),
    def   = c("Male; female",
              "Maternal age $>33$~yr; $~BS~le 33$~yr",
              "High school completed or beyond; otherwise",
              "Prenatal care began in the first trimester; otherwise",
              "Any pre-pregnancy smoking; none",
              "Married; unmarried; unknown",
              "White, Black, Hispanic, Asian, AIAN, NHOPI, multiracial"),
    stringsAsFactors = FALSE)

  xe <- new.env(); load(paths$rebin.data(paths$year), envir = xe)
  mr <- xe$missing.report
  # n_read is preferred when data.R recorded it; the fallback reconstructs the
  # denominator from a field's count and percentage, for older saved runs.
  n.read <- if (!is.null(mr$n_read)) mr$n_read[1] else round(mr$n_missing[2] / (mr$pct[2] / 100))
  miss <- function(nchs) {
    k <- match(nchs, mr$variable)
    if (is.na(k) || mr$n_missing[k] == 0) return("0")
    sprintf("%s (%.2f~BS~%%)", format(mr$n_missing[k], big.mark = ","), mr$pct[k])
  }

  t3 <- c("% Auto-generated by visual.R -- do not edit by hand.",
    "~BS~begin{table}[!htbp]", "~BS~centering~BS~small",
    paste0("~BS~caption[Predictors, definitions and missingness]{Predictors, their ",
           "NCHS source fields, and per-field missingness for the ", paths$year,
           " cohort. Five predictors are binary; marital status and maternal race ",
           "are multi-level, giving $2^{5}~BS~cdot3~BS~cdot7=672$ classes. Counts are ",
           "records with the field missing after NCHS ``not stated'' sentinel codes ",
           "are converted to missing. Marital status is retained as an explicit ",
           "~BS~emph{unknown} level rather than excluded; the remaining fields are ",
           "handled by complete-case deletion.~BS~label{tab:vars}}"),
    "~BS~begin{tabular*}{~BS~textwidth}{@{~BS~extracolsep~BS~fill}llcp{5.0cm}r@{}}",
    "~BS~toprule",
    paste("~BS~textbf{Field} & ~BS~textbf{NCHS field} & ~BS~textbf{Levels} &",
          "~BS~textbf{Definition} & ~BS~textbf{Missing} ~RE~"),
    "~BS~midrule")
  for (k in seq_len(nrow(defs)))
    t3 <- c(t3, sprintf("~BS~texttt{%s} & ~BS~texttt{%s} & %s & %s & %s ~RE~",
                        texvar(defs$field[k]), texvar(defs$nchs[k]),
                        defs$lev[k], defs$def[k], miss(defs$nchs[k])))
  t3 <- c(t3, "~BS~midrule",
    "~BS~multicolumn{5}{@{}l}{~BS~textit{Outcome}} ~RE~",
    sprintf(paste("~BS~texttt{dbwt} & ~BS~texttt{dbwt} & 11 / 10 & Birth weight in grams,",
                  "discretized into ten LBW deciles plus NBW & %s ~RE~"), miss("dbwt")),
    "~BS~midrule",
    sprintf("~BS~multicolumn{4}{@{}l}{Records read} & %s ~RE~", format(n.read, big.mark = ",")),
    sprintf("~BS~multicolumn{4}{@{}l}{Analyzed} & %s (%.1f~BS~%%) ~RE~",
            format(mr$n_analysed[1], big.mark = ","), 100 * mr$n_analysed[1] / n.read),
    "~BS~bottomrule", "~BS~end{tabular*}", "~BS~end{table}")
  writeLines(unesc(t3), file.path(article.dir, sprintf("table_predictors_%d.tex", paths$year)))

  cat("Saved results tables (varfreq, depth, flow) to:", article.dir, "\n")
}


# ------------------------------------------------------------------------------
# 2: PREDICTOR VARIABLE SPLIT DEPTH DISTRIBUTIONS
# ------------------------------------------------------------------------------
boot.file <- paths$bootstrap.results(paths$year)
if (!file.exists(boot.file))
  stop(sprintf("Missing %s\n  Run:  DMCART_B=10000 Rscript dm-cart.R %d-%d",
               boot.file, paths$prior.year, paths$target.year))
boot_env <- new.env()
load(boot.file, envir = boot_env)          # into an env, never the global
if (!exists("depth_df", envir = boot_env))
  stop(sprintf("`depth_df` not found in %s -- re-run dm-cart.R.", boot.file))
depth_df <- get("depth_df", envir = boot_env)

mean_depths <- depth_df %>%
  group_by(Variable) %>%
  summarize(mean_depth = mean(Depth, na.rm = TRUE), .groups = "drop")

p_depths <- ggplot(depth_df, aes(x = Depth, fill = Variable)) +
  # dropped the stale `align = "center"` argument (unknown to geom_histogram)
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


if (!exists("risk_summary", envir = boot_env))
  stop(sprintf("`risk_summary` not found in %s -- re-run dm-cart.R.", boot.file))
risk_summary <- get("risk_summary", envir = boot_env)

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

# output
pdf_out_path <- file.path(output.dir, sprintf("high_low_risk_pred_%d.pdf", paths$year))

pdf(file = pdf_out_path, width = 8.5, height = 8.5)

# Expanded bottom outer margin (oma = c(4.5, 0, 4, 0)) to hold two-line caption
par(mfrow = c(2, 1), 
    mar = c(3, 4.5, 3, 2),    
    oma = c(4.5, 0, 4, 0))      

# dynamic Y-axis limits for LBW subset
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
