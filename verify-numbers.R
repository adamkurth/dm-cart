# verify-numbers.R
# ==============================================================================
# Checks every number quoted in tex/article/main.tex against the saved analysis
# output. Run it after any re-fit, before sending the manuscript anywhere:
#
#   Rscript verify-numbers.R
#
# Exits 0 if everything matches, 1 if anything does not.
#
# The manuscript quotes about three dozen figures that are
# produced by the pipeline, not typed by hand. This catches any accidental changes 
# to the analysis that would make the manuscript wrongly claim a number.
#
# The expected values below are what the manuscript currently claims. When a
# check fails, decide which is wrong: usually it is the text, and the fix is to
# update the manuscript and then update the value here to match.
# ==============================================================================

suppressMessages(library(rpart))

fails <- 0L
check <- function(label, computed, claimed) {
  agree <- identical(as.character(computed), as.character(claimed))
  if (!agree) fails <<- fails + 1L
  cat(sprintf("  %-42s %-16s %-16s %s\n", label, computed, claimed,
              if (agree) "ok" else "MISMATCH"))
}

root <- normalizePath(".")
p <- function(...) file.path(root, ...)

fe <- new.env(); load(p("results/2023-2024/fullmodel/bootstrap_tree_results_2024.RData"), envir = fe)
le <- new.env(); load(p("results/2023-2024/lbwmodel/bootstrap_tree_results_2024.RData"),  envir = le)
al <- new.env(); load(p("data/rebin/2023-2024/informed_prior_2023.RData"),                envir = al)
dt <- new.env(); load(p("data/rebin/2023-2024/birthweight_data_rebin_2024.RData"),        envir = dt)
tf <- new.env(); load(p("data/rebin/2023-2024/dm_tree_rebin_2024.RData"),                 envir = tf)
tl <- new.env(); load(p("data/rebin_without_2.5kg/2023-2024/dm_tree_rebin_2024.RData"),   envir = tl)

cat("\ncheck                                      computed         claimed          status\n")
cat(strrep("-", 92), "\n")

# ---- Section 3, Data ---------------------------------------------------------
mr <- dt$missing.report
check("records read",      format(mr$n_read[1], big.mark = ","),     "3,638,436")
check("records analyzed",  format(mr$n_analysed[1], big.mark = ","), "3,476,907")
check("retained %",        sprintf("%.1f", 100 * mr$n_analysed[1] / mr$n_read[1]), "95.6")
check("predictor classes", nrow(dt$Y.matrix), "672")
check("empty classes",     sum(rowSums(dt$Y.matrix) == 0), "18")

# ---- Section 5.1, informed prior --------------------------------------------
av <- al$alphavec
check("alpha_11",          sprintf("%.3f", av[11]),                  "0.913")
check("smallest decile %", sprintf("%.2f", 100 * min(av[1:10])),     "0.76")
check("largest decile %",  sprintf("%.2f", 100 * max(av[1:10])),     "0.96")
check("2023 LBW share %",  sprintf("%.1f", 100 * sum(av[1:10])),     "8.7")

# ---- Section 5.2, fitted trees ----------------------------------------------
check("terminal nodes, full", sum(tf$dm.tree$frame$var == "<leaf>"), "23")
check("terminal nodes, LBW",  sum(tl$dm.tree$frame$var == "<leaf>"), "8")

# Figure 2 labels come from the single fitted tree, not the ensemble; the
# caption quotes this value explicitly, so it is checked too.
leaf.p <- function(fit, alpha, cols) {
  cnt <- fit$frame$yval2[fit$where, , drop = FALSE]
  pr  <- t(apply(cnt, 1, function(z) (z + alpha) / sum(z + alpha)))
  rowSums(pr[, cols, drop = FALSE])
}
check("fitted-tree P(LBW), class 181",
      sprintf("%.3f", leaf.p(tf$dm.tree, av, 1:10)[181]), "0.286")

# ---- Sections 5.4 and 6, risk stratification --------------------------------
expected <- list(
  full = list(eligible = "191",
              high = list(class = "181", mean = "0.284", ci = "[0.262,0.297]",
                          n = "1,710",   stab = "100.0"),
              low  = list(class = "40",  mean = "0.054", ci = "[0.052,0.055]",
                          n = "346,814", stab = "100.0")),
  lbw  = list(eligible = "65",
              high = list(class = "139", mean = "0.358", ci = "[0.350,0.365]",
                          n = "2,783",   stab = "98.9"),
              low  = list(class = "327", mean = "0.230", ci = "[0.206,0.254]",
                          n = "3,011",   stab = "98.8")))

for (arm in c("full", "lbw")) {
  e  <- if (arm == "full") fe else le
  rk <- e$risk.ranking
  pr <- e$bootstrap.tree.results$profiles$empirical
  check(paste(arm, "eligible classes"), sum(rk$eligible), expected[[arm]]$eligible)
  for (k in c("high", "low")) {
    want <- expected[[arm]][[k]]
    idx  <- pr[[paste0(k, ".index")]]
    r    <- rk[rk$class == idx, ]
    check(sprintf("%s %s class",    arm, k), idx, want$class)
    check(sprintf("%s %s mean",     arm, k), sprintf("%.3f", r$p_lbw), want$mean)
    check(sprintf("%s %s interval", arm, k), sprintf("[%.3f,%.3f]", r$lower, r$upper), want$ci)
    check(sprintf("%s %s births",   arm, k), format(r$births, big.mark = ","), want$n)
    check(sprintf("%s %s stability %%", arm, k),
          sprintf("%.1f", 100 * pr[[paste0(k, ".stability")]]), want$stab)
  }
}

fold <- function(e, hi, lo) {
  rk <- e$risk.ranking
  sprintf("%.1f", rk$p_lbw[rk$class == hi] / rk$p_lbw[rk$class == lo])
}
check("full model fold change", fold(fe, 181, 40),  "5.3")
check("LBW-only fold change",   fold(le, 139, 327), "1.6")

# The a priori "everything adverse" profile is four classes, not one; the
# manuscript states their ranks and support, so both are checked.
X  <- dt$final.data
rk <- fe$risk.ranking
rk <- rk[order(-rk$p_lbw), ]; rk$rank <- seq_len(nrow(rk))
adv <- which(X$dmar == 0 & X$meduc == 0 & X$precare1st == 0 &
             X$cig_0 == 1 & X$race == "Black")
adv <- adv[order(rk$rank[match(adv, rk$class)])]
check("everything-adverse ranks",
      paste(rk$rank[match(adv, rk$class)], collapse = ","), "1,5,9,13")
check("everything-adverse births",
      paste(rk$births[match(adv, rk$class)], collapse = ","), "186,481,161,509")
check("everything-adverse all ineligible",
      all(!rk$eligible[match(adv, rk$class)]), "TRUE")

# Figure 4's decile profiles are quoted in the text, so they are checked too.
prof <- function(env, subgroup) {
  d <- env$risk_summary
  if ("Definition" %in% names(d)) d <- d[d$Definition == "empirical", ]
  d <- d[d$Subgroup == subgroup & d$Category %in% paste0("Q", 1:10), ]
  d$Mean[order(match(d$Category, paste0("Q", 1:10)))]
}
fh <- prof(fe, "High-Risk Subgroup"); fl <- prof(fe, "Low-Risk Subgroup")
lh <- prof(le, "High-Risk Subgroup"); ll <- prof(le, "Low-Risk Subgroup")
check("class 181 decile range", sprintf("%.3f-%.3f", min(fh), max(fh)), "0.019-0.033")
check("class 40 decile range",  sprintf("%.3f-%.3f", min(fl), max(fl)), "0.005-0.006")
check("class 139 Q1 / Q10",     sprintf("%.3f/%.3f", lh[1], lh[10]),    "0.134/0.079")
check("class 327 Q1 / peak",    sprintf("%.3f/%.3f", ll[1], max(ll)),   "0.071/0.128")
check("class 327 peak decile",  paste0("Q", which.max(ll)),             "Q8")

# Classes 67 and 68 are named in the shared-feature paragraph.
check("class 67/68 profile",
      paste(as.character(X$dmar[67]), as.character(X$mager[67]),
            as.character(X$cig_0[67]), as.character(X$race[67]), sep = "/"),
      "0/1/1/White")

cat(strrep("-", 92), "\n")
if (fails == 0L) {
  cat("All checks passed.\n\n")
} else {
  cat(sprintf("%d MISMATCH(ES). Update the manuscript, then update the expected\n", fails))
  cat("values in this script to match.\n\n")
  quit(status = 1L)
}
