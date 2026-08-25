# dm-cart.R
# Fits Dirichlet-Multinomial (DM) Decision Trees using custom rpart methods
# and runs multinomial parametric bootstrapping across model variations.
#
# Usage via Command Line:
#   Rscript dm-cart.R 2023-2024
#   Rscript dm-cart.R 2023
#
# Default (Interactive / RStudio):
#   Defaults to prior year 2023 (target year 2024)

library(rpart)
library(ggplot2)
library(xtable)
library(reshape2)
library(dplyr)
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
init.paths  <- get_paths(year = PRIOR.YEAR, with.2.5kg = TRUE)
TARGET.YEAR <- init.paths$target.year   # 2024 -- the cohort being modeled
PRIOR.YEAR  <- init.paths$prior.year    # 2023 -- cutpoints + alphavec ONLY

set.seed(08242026) # for reproducibility of the bootstrap

cat(sprintf("\n==================================================\n"))
cat(sprintf(" Running DM-CART Analysis for Cohort: %s\n", init.paths$cohort))
cat(sprintf(" Prior Year: %d  |  Target Year: %d\n", PRIOR.YEAR, TARGET.YEAR))
cat(sprintf("==================================================\n\n"))


# -------- Step 2: Loop across both model types (w/ 2.5kg & w/out 2.5kg) --------
for (PROCESS.WITH.2.5KG in c(TRUE, FALSE)){

    model.label <- if (PROCESS.WITH.2.5KG) "Full Model (rebin)" else "LBW-Only Model (rebin_without_2.5kg)"
    cat(sprintf("\n--------------------------------------------------\n"))
    cat(sprintf(" Processing Model: %s\n", model.label))
    cat(sprintf("--------------------------------------------------\n"))


    paths <- get_paths(year = PRIOR.YEAR, with.2.5kg = PROCESS.WITH.2.5KG)
    ensure_dirs(paths)

    # Load Prior Hyperparameters (alphavec)
    prior.file <- paths$informed.prior(PRIOR.YEAR)
    if (!file.exists(prior.file)) {
      stop(sprintf("Prior file not found: %s. Please run quantile.R first!", prior.file))
    }
    load(prior.file)

    # Load Rebinned Birth Weight Data
    rebin.file <- paths$rebin.data(TARGET.YEAR)
    if (!file.exists(rebin.file)) {
      stop(sprintf("Rebinned data file not found: %s. Please run data.R first!", rebin.file))
    }
    load(rebin.file)

    # Verify Loaded Data
    Y.df <- counts.df
    response.cols <- colnames(Y.df)

    # categorical index sets, for computing posterior P(LBW) in leaves
    lbw.cols <- 1:10                                                   # Q1..Q10, the LBW region
    nbw.cols <- if (PROCESS.WITH.2.5KG) 11:ncol(Y.df) else integer(0)  # Normal BW, full model only

    cat(sprintf("Loaded Dirichlet Hyperparameter vector (alphavec) length: %d\n", length(alphavec)))
    cat(sprintf("Counts DataFrame dimension: %d rows x %d categories\n", nrow(Y.df), ncol(Y.df)))

    #-----------------------------------------------------------
    # B. Simplified Dirichlet-Multinomial Functions
    #-----------------------------------------------------------
    log.dm.likelihood <- function(counts, alpha = 1) {
      N <- sum(counts)
      alpha_0 <- sum(alpha)

      term1 <- lgamma(alpha_0) - lgamma(N + alpha_0)
      term2 <- sum(lgamma(counts + alpha) - lgamma(alpha))

      ll <- term1 + term2
      if (is.na(ll) || is.infinite(ll)) {
        ll <- -Inf
      }
      return(ll)
    }


    #-----------------------------------------------------------
    # C. Custom rpart Method (DM with log-likelihood deviance)
    #-----------------------------------------------------------
    myinit <- function(y, offset, parms=NULL, wt=NULL) {
      list(
        method="dm",
        y = y,
        parms = list(alpha=1),
        numresp = ncol(y),
        numy = ncol(y),
        # rpart producing useless labels, handles matrix frame$yval2, so we override the default summary/print/text callbacks to produce more useful output 
        summary = function(yval, dev, wt, ylevel, digits) {
          cat("DM Model\n")
          cat("Deviance: ", format(dev, digits=digits), "\n")
          cat("Counts: ", paste(format(yval, digits=digits), collapse=", "), "\n")
        },
        print = function(yval, dev, wt, ylevel, digits) {
          # yval may arrive as a matrix (many nodes) or a single row.
          m <- if (is.matrix(yval)) yval else matrix(yval, nrow = 1)
          modal <- apply(m, 1, which.max)
          tot   <- rowSums(m)
          cat(paste0("modal category k=", modal, ", n=", format(tot, digits = digits)),
              sep = "\n")
        },
        # text() is called at each node to produce the label. It returns a character vector, one element per node. 
        # this is where we compute the posterior P(LBW) for the node and format it into the label.
        text = function(yval, dev, wt, ylevel, digits, n, use.n){
          m     <- if (is.matrix(yval)) yval else matrix(yval, nrow = 1)
          modal <- apply(m, 1, which.max)          # per-node modal category
          # compute posterior P(LBW) for each node, using the Dirichlet-Multinomial posterior mean:
          #   (counts + alpha) / sum(counts + alpha)
          p.lbw <- apply(m, 1, function(z) {
            post <- (z + alphavec) / sum(z + alphavec)
            sum(post[lbw.cols])
          })
          lbl <- paste0("k", modal, "\nP(LBW)=", sprintf("%.3f", p.lbw))
          if (use.n) lbl <- paste0(lbl, "\nn=", n)
          lbl
        }
      )
    }

    # myeval() is called at each node to compute the deviance of the node's counts under the DM model. It returns a list with the counts and the deviance.
    myeval <- function(y, wt = NULL, parms = 1) {
      counts <- colSums(y)
      dev <- -log.dm.likelihood(counts, alpha = alphavec)
      return(list(label = counts, deviance = dev))
    }

    MAX.EXHAUSTIVE.LEVELS <- 12   # 2^k -1 = 2^11-1 = 2047 partitions per node: still trivial

    mysplit <- function(y, wt, x, parms, continuous = FALSE) {
      ux <- sort(unique(x))
      k  <- length(ux)
      if (k < 2) return(list(goodness = numeric(0), direction = ux))
      # the candidate splits are the k-1 contiguous splits of the ordered levels, or all 2^(k-1)-1 non-empty subsets of the unordered levels. The goodness of each candidate split is computed as the log-likelihood of the left child plus the log-likelihood of the right child minus the log-likelihood of the parent. The direction is the order of the levels in the split.
      lev.sums   <- rowsum(y, group = match(x, ux), reorder = TRUE)   # k x K
      total      <- colSums(lev.sums)
      parent.ll  <- log.dm.likelihood(total, alpha = alphavec)
      
      # as before: L_left + L_right - L_parent, which equals parent.dev - child.dev.
      goodness.of <- function(in.left) {
        L <- colSums(lev.sums[in.left, , drop = FALSE])
        log.dm.likelihood(L, alpha = alphavec) +
          log.dm.likelihood(total - L, alpha = alphavec) - parent.ll
      }

      if (continuous) {
        # Ordered predictor: the k-1 cutpoint splits ARE the whole search space.
        goodness <- vapply(seq_len(k - 1),
                           function(i) goodness.of(seq_len(k) <= i), numeric(1))
        return(list(goodness = pmax(goodness, 0), direction = rep(1, k - 1)))
      }

      # ---- Unordered factor: choose the ordering to search under ----
      # The exhaustive search is only feasible for k <= 12, so for larger k we fall back to a heuristic ordering of the levels by their LBW share. The exhaustive search is guaranteed to find the best split.
      bit <- 2^(0:(k - 1)) # bitmask for the k levels, so bit[i] is the bit for level i
      if (k <= MAX.EXHAUSTIVE.LEVELS) {
        # results are stored in a table indexed by the bitmask of the levels sent left. The table is symmetric: the complement of a mask is the same split, so we only need to compute half of them. The best split is the one with the highest goodness (gvals).
        n.mask <- 2^(k - 1) - 1
        gvals  <- numeric(n.mask)
        for (m in seq_len(n.mask)) gvals[m] <- goodness.of(as.logical(bitwAnd(m, bit)))

        best.m   <- which.max(gvals)
        best.idx <- which(as.logical(bitwAnd(best.m, bit)))
        ord.idx  <- c(best.idx, setdiff(seq_len(k), best.idx))   # winner first
        # a prefix of `ord.idx` is itself a mask; canonicalize it against its complement (same split) so it indexes into the table we just built.
        full <- 2^k - 1
        goodness <- vapply(seq_len(k - 1), function(i) {
          mask <- sum(bit[ord.idx[seq_len(i)]])
          gvals[min(mask, full - mask)]
        }, numeric(1))

      } else {
        # Too many levels to enumerate. Fall back to ordering the levels by their LBW share, so that the contiguous splits rpart examines are at least monotone in the quantity of interest rather than arbitrary. This is a heuristic that may not find the best split, but it is fast and reasonable.
        row.tot   <- rowSums(lev.sums)
        lbw.share <- ifelse(row.tot == 0, 0,
                            rowSums(lev.sums[, lbw.cols, drop = FALSE]) / pmax(row.tot, 1))
        ord.idx   <- order(lbw.share)
        goodness  <- vapply(seq_len(k - 1),
                            function(i) goodness.of(seq_len(k) %in% ord.idx[seq_len(i)]),
                            numeric(1))
      }

      list(goodness = pmax(goodness, 0), direction = ux[ord.idx])
    }

    # replaces predict(..., type = "matrix") for a custom method, which does not work. This is the same computation as in mytext() above, but returns a matrix of posterior probabilities rather than a formatted label.
    dm.leaf.probs <- function(fit, alphavec, cat.names = colnames(Y.df)) {
      stopifnot(!is.null(fit$where), is.matrix(fit$frame$yval2))

      counts <- fit$frame$yval2[fit$where, , drop = FALSE]                      # leaf counts, per class row
      probs  <- t(apply(counts, 1, function(z) (z + alphavec) / sum(z + alphavec)))

      colnames(probs) <- cat.names
      rownames(probs) <- NULL
      probs
    }


    dm.method <- list(init=myinit, eval=myeval, split=mysplit, method="dm")

    # =========================================================================
    # plot.dm.tree()
    # =========================================================================
    # the default rpart plotting functions are not very useful for multi-level factor splits, because they only print the first level of the split and then "..." for the rest. We do not want to let this mislead the reader into thinking that the model is only peeling off one level at a time, when in fact it may be grouping multiple levels together. We also want to print the posterior P(LBW) and the number of births at each leaf, which are important for interpreting the model.
    #
    # we draw the skeleton of the plot using plot.rpart (handle layout) then place labels ourselves at coordinates from rpart:::rpartco(). The left/right level sets come from labels(fit, pretty = 0, collapse = FALSE), which is rpart's own decoding of csplit -- so the printed groupings are exactly what the model used, not a reconstruction.

    #     at each internal node : the SPLIT VARIABLE name
    #     on the left branch    : the levels sent left
    #     on the right branch   : the levels sent right
    #     at each leaf          : modal category, P(LBW), and N births
    #
    # - IMPORTANT: in rpart's frame$n counts the data rows at a node, where instead here the data rows is now a CLASS, not a birth. The label reports N = sum of the node's category counts (actual births), with the class count given separately as "cls". This is a more accurate representation of the sample size at each node, and avoids confusion with the total number of classes in the data.
    # - Tree-drawing helpers (dm.tree.coords / draw.dm.tree / plot.dm.tree) are in util.R
    # =========================================================================

    #-----------------------------------------------------------
    # D. Fit Primary DM-based rpart Tree
    #-----------------------------------------------------------
    cat("\n--- Fitting Primary DM Decision Tree ---\n")
    dm.control <- rpart.control(minsplit=2, cp=0, maxdepth=8, xval=0, usesurrogate = 0)
    dm.tree <- rpart(
      formula = Y.df ~ .,
      data = data.frame(X.matrix),
      method = dm.method,
      control = dm.control
    )

      # save main tree
      dm.tree.file <- paths$dm.tree.rebin(TARGET.YEAR)
      save(dm.tree, file = dm.tree.file)
      cat("Saved decision tree to:", dm.tree.file, "\n")

      cat("\n--- Rendering prototype tree ---\n")

      DISPLAY.DEPTH <- 3    # for prototype tree figure only, not the main fit 
      proto.tree <- rpart(
        formula = Y.df ~ .,
        data    = data.frame(X.matrix),
        method  = dm.method,
        control = rpart.control(minsplit = 2, cp = 0, maxdepth = DISPLAY.DEPTH,
                                xval = 0, usesurrogate = 0)
      )
      plot.dm.tree(
        dm.tree,
        file     = file.path(paths$plots, sprintf("full_tree_%d.pdf", TARGET.YEAR)),
        title    = sprintf("DM-CART (Full Tree, maxdepth = %d)", dm.control$maxdepth),
        subtitle = sprintf("%s | cohort %s | %d classes",
                           model.label, init.paths$cohort, nrow(X.matrix)),
        lbw.cols = lbw.cols, alphavec = alphavec
      )

      # ---- split-partition report --------------------------------
      report.splits <- function(fit, label) {
        fr <- fit$frame; leaf <- fr$var == "<leaf>"
        if (all(leaf)) { cat("  (no splits)\n"); return(invisible(NULL)) }
        lb  <- labels(fit, pretty = 0, collapse = FALSE)
        ids <- as.numeric(rownames(fr))
        out <- data.frame(
          node     = ids[!leaf],
          depth    = floor(log2(ids[!leaf])),
          variable = as.character(fr$var[!leaf]),
          left     = lb[!leaf, 1],
          right    = lb[!leaf, 2],
          n_births = round(rowSums(fr$yval2)[!leaf]),
          stringsAsFactors = FALSE
        )
        # Was this partition reachable by the OLD sorted-code search? That search could only take a contiguous PREFIX of the levels present at the node, ordered by their factor codes. Anything else is a split that only the exhaustive search can find. This is the concrete test of whether the multi-level factor handling is doing any work.
        xlev <- attr(fit, "xlevels")
        out$n_levels <- NA_integer_
        out$reachable_before <- NA
        for (r in seq_len(nrow(out))) {
          v <- out$variable[r]
          if (is.null(xlev[[v]])) next                       # numeric predictor
          Ls <- trimws(strsplit(out$left[r],  ",")[[1]])
          Rs <- trimws(strsplit(out$right[r], ",")[[1]])
          present <- c(Ls, Rs)
          ordered.present <- present[order(match(present, xlev[[v]]))]
          out$n_levels[r] <- length(present)
          m <- length(Ls)
          out$reachable_before[r] <-
            setequal(Ls, head(ordered.present, m)) || setequal(Rs, head(ordered.present, length(Rs)))
        }
        out <- out[order(out$depth, out$node), ]
        cat(sprintf("\nSplit partitions -- %s\n", label))
        print(out, row.names = FALSE)

        multi <- out[!is.na(out$n_levels) & out$n_levels > 2, ]
        new   <- multi[!multi$reachable_before, ]
        cat(sprintf("  %d splits on a >2-level factor; %d of them are NOT reachable by the\n",
                    nrow(multi), nrow(new)))
        cat("  old sorted-code search, i.e. they exist only because of \n")
        cat("  the new exhaustive search. These are:\n")
        if (nrow(new)) {
          for (r in seq_len(nrow(new)))
            cat(sprintf("    node %-4d %-10s  %s  |  %s\n",
                        new$node[r], new$variable[r], new$left[r], new$right[r]))
        } else cat("    (none at this fit)\n")
        invisible(out)
      }
      split.report <- report.splits(dm.tree, sprintf("full tree, %s", model.label))

      # ---- Depth-controlled comparison ----
      cat("\n--- Depth-controlled comparison (maxdepth 2-5) ---\n")
      depths    <- c(2, 3, 4, 5)
      trees     <- list()
      tree.sum  <- list()

      for (d in depths) {
        d.tree <- rpart(
          formula = Y.df ~ .,
          data    = data.frame(X.matrix),
          method  = dm.method,
          control = rpart.control(minsplit = 2, cp = 0, maxdepth = d,
                                  xval = 0, usesurrogate = 0)
        )
        trees[[as.character(d)]] <- d.tree

        vars.used <- unique(as.character(d.tree$frame$var[d.tree$frame$var != "<leaf>"]))

        probs   <- dm.leaf.probs(d.tree, alphavec)
        lbw.prob <- rowSums(probs[, lbw.cols, drop = FALSE])

        tree.sum[[as.character(d)]] <- list(
          depth            = d,
          variables        = vars.used,
          n.variables      = length(vars.used),
          n.terminal.nodes = sum(d.tree$frame$var == "<leaf>"),
          first.split      = if (length(vars.used)) vars.used[1] else "none",
          mean.lbw.prob    = mean(lbw.prob),
          median.lbw.prob  = median(lbw.prob),
          min.lbw.prob     = min(lbw.prob),
          max.lbw.prob     = max(lbw.prob),
          sd.lbw.prob      = sd(lbw.prob)
        )

        cat(sprintf("  depth %d: %2d leaves, %d vars (first split: %-10s) P(LBW) mean %.4f [%.4f, %.4f]\n",
                    d, tree.sum[[as.character(d)]]$n.terminal.nodes, length(vars.used),
                    tree.sum[[as.character(d)]]$first.split,
                    mean(lbw.prob), min(lbw.prob), max(lbw.prob)))
      }

      # Which variables survive at each depth
      var.across.depths <- data.frame(variable = colnames(X.matrix),
                                      stringsAsFactors = FALSE)
      for (d in depths) {
        var.across.depths[[paste0("depth_", d)]] <-
          as.integer(var.across.depths$variable %in% tree.sum[[as.character(d)]]$variables)
      }
      var.across.depths$total.usage <- rowSums(var.across.depths[, paste0("depth_", depths)])
      var.across.depths <- var.across.depths[order(-var.across.depths$total.usage), ]

      pred.stats.compare <- do.call(rbind, lapply(as.character(depths), function(d) {
        s <- tree.sum[[d]]
        data.frame(depth = s$depth, mean.lbw.prob = s$mean.lbw.prob,
                   min.lbw.prob = s$min.lbw.prob, max.lbw.prob = s$max.lbw.prob,
                   sd.lbw.prob = s$sd.lbw.prob,
                   n.terminal.nodes = s$n.terminal.nodes, n.variables = s$n.variables)
      }))

      cat("\nVariable usage across depths:\n"); print(var.across.depths, row.names = FALSE)

      #-----------------------------------------------------------
      # F. Multinomial Parametric Bootstrap Sampling
      #-----------------------------------------------------------
      # B is now controlled by the environment variable DMCART_B, which is set in the shell before running the script. This allows for different bootstrap sizes for interactive vs batch runs.
      # Use:
      #   DMCART_B=200 Rscript dm-cart.R 2023-2024      # smoke test
      #   DMCART_B=10000 Rscript dm-cart.R 2023-2024    # publication run
      B <- as.integer(Sys.getenv("DMCART_B",
                                 unset = if (interactive()) "200" else "10000"))
      stopifnot(is.finite(B), B >= 2)
      n.rows <- nrow(counts.df)
      # NBW / LBW split is now a single class, so we can compute the multinomial probabilities for each row of Y.df using the Dirichlet-Multinomial posterior mean. This is used to generate bootstrap samples of the counts, which are then used to fit new trees and record their structure and predictions.
      Yhat <- vector("list", length = B)

      multinomial.probs <- matrix(0, nrow = n.rows, ncol = ncol(Y.df))
      for (i in 1:n.rows) {
        cell.counts <- as.numeric(Y.df[i, ])
        multinomial.probs[i, ] <- (cell.counts + alphavec) / (sum(cell.counts) + sum(alphavec))
      }

      row.probs <- rowMeans(counts.df)
      row.probs <- row.probs / sum(row.probs)

      tree.structure <- list()
      top.vars.used <- character(B)
      n.vars.used <- integer(B)
      depth.rows <- vector("list", B)   # real split depths

      cat(sprintf("Executing Parametric Bootstrap (B = %d)...\n", B))
      pb <- txtProgressBar(0, B, style = 3)

      for (b in seq_len(B)) {
        n.star <- as.vector(rmultinom(1, size = sum(Y.df[]), prob = row.probs))
        boot.counts <- t(sapply(1:n.rows, function(j) rmultinom(1, size = n.star[j], prob = multinomial.probs[j, ])))
        colnames(boot.counts) <- colnames(Y.df)
        boot.df <- as.data.frame(boot.counts)

        model.data <- cbind(boot.df, X.matrix)
        r.formula <- as.formula(paste0("cbind(", paste(colnames(boot.df), collapse = ","), ") ~ ."))

        dm.tree.b <- rpart(
          formula = r.formula,
          data = model.data,
          method  = dm.method,
          control = dm.control
        )

        vars.used <- unique(dm.tree.b$frame$var[dm.tree.b$frame$var != "<leaf>"])
        n.vars.used[b] <- length(vars.used)
        top.vars.used[b] <- if (length(vars.used)) vars.used[1] else "none"

        tree.structure[[b]] <- list(
          variables = vars.used,
          frame = dm.tree.b$frame,
          splits = dm.tree.b$splits
        )

        # record the depth of every split, per replicate
        internal <- dm.tree.b$frame$var != "<leaf>"
        if (any(internal)) {
          node.ids <- as.numeric(rownames(dm.tree.b$frame))
          depth.rows[[b]] <- data.frame(
            replicate = b,
            Variable  = as.character(dm.tree.b$frame$var[internal]),
            Depth     = floor(log2(node.ids[internal])),
            stringsAsFactors = FALSE
          )
        }
        Yhat[[b]] <- dm.leaf.probs(dm.tree.b, alphavec)
        setTxtProgressBar(pb, b)
      }
      close(pb)

      # high/low-risk subgroups are now single CLASSES
      find.class <- function(profile, X = X.matrix) {
        hit <- Reduce(`&`, Map(function(v, val) as.character(X[[v]]) == as.character(val),
                               names(profile), profile))
        if (sum(hit) != 1) {
          stop(sprintf(paste0("Risk profile matched %d classes, expected exactly 1.\n",
                              "  profile: %s\n",
                              "  Does it name every predictor in X.matrix (%s)?"),
                       sum(hit),
                       paste(names(profile), unlist(profile), sep = "=", collapse = ", "),
                       paste(colnames(X), collapse = ", ")))
        }
        which(hit)
      }

      # dmar: "0" = unmarried, "1" = married ("Unknown" exists under DMAR.POLICY
      # = "unknown-level" but is deliberately not used for either risk profile).
      high.risk.profile <- list(sex = "1", dmar = "0", mager = "0", meduc = "0",
                                precare1st = "0", cig_0 = "1", race = "Black")
      low.risk.profile  <- list(sex = "0", dmar = "1", mager = "0", meduc = "1",
                                precare1st = "1", cig_0 = "0", race = "White")

      apriori.high.index <- find.class(high.risk.profile)
      apriori.low.index  <- find.class(low.risk.profile)
      describe.profile <- function(p) paste(names(p), unlist(p), sep = "=", collapse = ", ")
      


      # here we find the most extreme classes in the fitted model, which may or may not be the same as the apriori profiles. The apriori profiles are just a reference point, but the model may find a different class to be the highest or lowest risk. We report both the apriori and the fitted extremes.
      # 
      # pooled over the bootstrap replicates, we compute the posterior P(LBW) for each class, and then find the class with the highest and lowest mean P(LBW). We also compute the 95% credible interval for each class's P(LBW) across the bootstrap replicates. The classes with the highest and lowest mean P(LBW) are reported as the fitted high-risk and low-risk classes, respectively.
      # 
      #   lbw.by.class[i, b] = sum_{k in LBW} Yhat[[b]][i, k]
      #
      # Yhat is the leaf posterior probability vector for each class in each bootstrap replicate. The sum over the LBW columns gives the posterior P(LBW) for that class in that replicate. We then take the mean and quantiles across replicates to get the point estimate and credible interval for each class's P(LBW). 
      # We then take the argmax and argmin of the mean P(LBW) across classes to find the fitted high-risk and low-risk classes. We also report the stability of these extremes across replicates, which is the proportion of replicates in which the chosen class is among the extreme classes within a tolerance defined by TIE.REL. 
      # The tolerance is defined as a fraction of the width of the 95% credible interval for that class's P(LBW), so that we consider classes to be tied if their means are within that range.
      # Lastly, the minimum number of observed births per class is enforced to avoid reporting artifacts from classes with very few observations. Classes with fewer than MIN.CLASS.BIRTHS observed births are excluded from the ranking of extremes, but their stability is still reported. 
      MIN.CLASS.BIRTHS <- 1000

      # ---- the ranking metric has to differ by model arm ----------------------
      if (PROCESS.WITH.2.5KG) {
        risk.cols  <- lbw.cols          # Q1..Q10 out of 11
        risk.label <- "P(LBW)"
        risk.desc  <- "probability of low birth weight"
      } else {
        risk.cols  <- 1:3               # three smallest deciles, conditional on LBW
        risk.label <- "P(Q1, Q2 or Q3 | LBW)"
        risk.desc  <- "share of the three most severe LBW deciles"
      }

      class.births <- rowSums(Y.df)
      # I x B matrix of the posterior risk score per class per replicate
      lbw.by.class <- vapply(Yhat, function(P) rowSums(P[, risk.cols, drop = FALSE]),
                             numeric(n.rows))
      class.mean <- rowMeans(lbw.by.class)
      class.lwr  <- apply(lbw.by.class, 1, quantile, 0.025)
      class.upr  <- apply(lbw.by.class, 1, quantile, 0.975)

      eligible <- which(class.births >= MIN.CLASS.BIRTHS)
      if (length(eligible) < 2)
        stop(sprintf("Only %d classes clear MIN.CLASS.BIRTHS = %d; cannot rank risk profiles.",
                     length(eligible), MIN.CLASS.BIRTHS))

      # ties are defined at the resolution of the estimate 
      # i.e., if the 95% CI is [0.10, 0.12] then any class with mean P(LBW) within 5% of that range is considered tied.
      TIE.REL <- 0.05
      pick.extreme <- function(want.max) {
        v      <- class.mean[eligible]
        target <- if (want.max) max(v) else min(v)
        k      <- eligible[which(abs(v - target) == min(abs(v - target)))[1]]
        tol    <- max(TIE.REL * (class.upr[k] - class.lwr[k]), 1e-10)
        tied   <- eligible[abs(v - target) <= tol]
        list(index = tied[which.max(class.births[tied])], tied = tied, tol = tol)
      }
      hi.pick <- pick.extreme(TRUE);  lo.pick <- pick.extreme(FALSE)
      is.high.risk.indices <- hi.pick$index
      is.low.risk.indices  <- lo.pick$index
      # stability is measured against same tolerance, so it reports whether the chosen class is among the statistically indistinguishable extremes in a replicate -- not whether it wins a tie-break decided by noise.
      hi.col.max <- apply(lbw.by.class[eligible, , drop = FALSE], 2, max)
      lo.col.min <- apply(lbw.by.class[eligible, , drop = FALSE], 2, min)
      high.stability <- mean(abs(lbw.by.class[is.high.risk.indices, ] - hi.col.max) <= hi.pick$tol)
      low.stability  <- mean(abs(lbw.by.class[is.low.risk.indices,  ] - lo.col.min) <= lo.pick$tol)

      profile.of <- function(i) as.list(setNames(as.character(unlist(X.matrix[i, ])),
                                                 colnames(X.matrix)))
      empirical.high.profile <- profile.of(is.high.risk.indices)
      empirical.low.profile  <- profile.of(is.low.risk.indices)

      # ranking table
      risk.ranking <- data.frame(
        class    = seq_len(n.rows),
        births   = round(class.births),
        p_lbw    = class.mean,
        lower    = class.lwr,
        upper    = class.upr,
        eligible = class.births >= MIN.CLASS.BIRTHS,
        profile  = apply(X.matrix, 1, function(r)
                     paste(colnames(X.matrix), r, sep = "=", collapse = ", ")),
        stringsAsFactors = FALSE
      )
      risk.ranking <- risk.ranking[order(-risk.ranking$p_lbw), ]

      cat(sprintf("\n--- Empirical risk extremes ---\n"))
      cat(sprintf("Ranking metric: %s (%s)\n", risk.label, risk.desc))
      cat(sprintf("Classes: %d total, %d with >= %d births (eligible)\n",
                  n.rows, length(eligible), MIN.CLASS.BIRTHS))
      cat(sprintf("HIGHEST risk: class %d  %s=%.4f [%.4f, %.4f]  N=%s  (top in %.0f%% of replicates)\n",
                  is.high.risk.indices, risk.label,
                  class.mean[is.high.risk.indices],
                  class.lwr[is.high.risk.indices], class.upr[is.high.risk.indices],
                  format(round(class.births[is.high.risk.indices]), big.mark = ","),
                  100 * high.stability))
      cat("   ", describe.profile(empirical.high.profile), "\n")
      cat(sprintf("LOWEST  risk: class %d  %s=%.4f [%.4f, %.4f]  N=%s  (bottom in %.0f%% of replicates)\n",
                  is.low.risk.indices, risk.label,
                  class.mean[is.low.risk.indices],
                  class.lwr[is.low.risk.indices], class.upr[is.low.risk.indices],
                  format(round(class.births[is.low.risk.indices]), big.mark = ","),
                  100 * low.stability))
      cat("   ", describe.profile(empirical.low.profile), "\n")
      if (length(hi.pick$tied) > 1)
        cat(sprintf("   NOTE: %d classes tie at the maximum (same terminal node); representative = the one with most births.\n",
                    length(hi.pick$tied)))
      if (length(lo.pick$tied) > 1)
        cat(sprintf("   NOTE: %d classes tie at the minimum (same terminal node); representative = the one with most births.\n",
                    length(lo.pick$tied)))
      # rank the assumed profiles among ALL classes, not just the eligible ones -- an a priori profile can easily fall below the support threshold (the assumed high-risk one does), and reporting rank NA would read as a bug rather than as the finding it is.
      rank.all.desc <- rank(-class.mean, ties.method = "min")
      rank.all.asc  <- rank( class.mean, ties.method = "min")
      cat(sprintf("\nFor comparison, the ASSUMED profiles rank (of %d classes):\n", n.rows))
      cat(sprintf("   a priori high (class %d): %s=%.4f -> rank %d%s\n",
                  apriori.high.index, risk.label,
                  class.mean[apriori.high.index],
                  rank.all.desc[apriori.high.index],
                  if (class.births[apriori.high.index] < MIN.CLASS.BIRTHS)
                    sprintf("  [BELOW SUPPORT THRESHOLD: N=%s]",
                            format(round(class.births[apriori.high.index]), big.mark = ",")) else ""))
      cat(sprintf("   a priori low  (class %d): %s=%.4f -> rank %d from the bottom%s\n",
                  apriori.low.index, risk.label,
                  class.mean[apriori.low.index],
                  rank.all.asc[apriori.low.index],
                  if (class.births[apriori.low.index] < MIN.CLASS.BIRTHS)
                    sprintf("  [BELOW SUPPORT THRESHOLD: N=%s]",
                            format(round(class.births[apriori.low.index]), big.mark = ",")) else ""))
      cat("\nTop 5 eligible classes by P(LBW):\n")
      print(head(risk.ranking[risk.ranking$eligible, c("class","births","p_lbw","lower","upper","profile")], 5),
            row.names = FALSE)
      cat("\nBottom 5 eligible classes by P(LBW):\n")
      print(tail(risk.ranking[risk.ranking$eligible, c("class","births","p_lbw","lower","upper","profile")], 5),
            row.names = FALSE)
      # with single class per subgroup, the bootstrap CI is just the CI for that class. But we keep the code in this form so that it will also work if a future profile deliberately spans several classes.
      direct.bootstrap.calculation <- function(indices, yhat.list = Yhat) {
        B_len <- length(yhat.list)
        K_len <- ncol(yhat.list[[1]])
        prob.matrix <- matrix(0, nrow = B_len, ncol = K_len)

        for (b_i in 1:B_len) {
          if (length(indices) > 0) {
            prob.matrix[b_i, ] <- colMeans(yhat.list[[b_i]][indices, , drop = FALSE])
          }
        }

        list(
          means = colMeans(prob.matrix),
          sd    = apply(prob.matrix, 2, sd),
          lwr   = apply(prob.matrix, 2, function(x) quantile(x, 0.025)),
          upr   = apply(prob.matrix, 2, function(x) quantile(x, 0.975))
        )
      }

      high.risk.result <- direct.bootstrap.calculation(is.high.risk.indices)
      low.risk.result  <- direct.bootstrap.calculation(is.low.risk.indices)
      ap.high.result   <- direct.bootstrap.calculation(apriori.high.index)
      ap.low.result    <- direct.bootstrap.calculation(apriori.low.index)

      K_cats <- length(high.risk.result$means)
      cat.labels <- if (PROCESS.WITH.2.5KG) c(paste0("Q", 1:10), "Normal") else paste0("Q", 1:10)

      # ---- persist the risk summary table for later use in the manuscript ----
      mk.risk.block <- function(res, subgroup, definition, idx, prof) {
        data.frame(Subgroup = subgroup, Definition = definition, Category = cat.labels,
                   Mean = res$means, SE = res$sd, Lower = res$lwr, Upper = res$upr,
                   ClassIndex = idx, Births = round(class.births[idx]),
                   Profile = describe.profile(prof), stringsAsFactors = FALSE)
      }
      risk_summary <- rbind(
        mk.risk.block(high.risk.result, "High-Risk Subgroup", "empirical",
                      is.high.risk.indices, empirical.high.profile),
        mk.risk.block(low.risk.result,  "Low-Risk Subgroup",  "empirical",
                      is.low.risk.indices,  empirical.low.profile),
        mk.risk.block(ap.high.result,   "High-Risk Subgroup", "a priori",
                      apriori.high.index,   high.risk.profile),
        mk.risk.block(ap.low.result,    "Low-Risk Subgroup",  "a priori",
                      apriori.low.index,    low.risk.profile)
      )
      rownames(risk_summary) <- NULL

      # Aggregate P(LBW) per subgroup and definition
      risk_totals <- do.call(rbind, lapply(split(risk_summary, ~ Subgroup + Definition), function(d) {
        data.frame(Subgroup = d$Subgroup[1], Definition = d$Definition[1],
                   ClassIndex = d$ClassIndex[1], Births = d$Births[1],
                   P_RISK = sum(d$Mean[risk.cols]),
                   P_RISK_lower = sum(d$Lower[risk.cols]), P_RISK_upper = sum(d$Upper[risk.cols]),
                   RiskMetric = risk.label,
                   P_LBW = if (PROCESS.WITH.2.5KG) sum(d$Mean[lbw.cols]) else NA_real_,
                   P_LBW_lower = if (PROCESS.WITH.2.5KG) sum(d$Lower[lbw.cols]) else NA_real_,
                   P_LBW_upper = if (PROCESS.WITH.2.5KG) sum(d$Upper[lbw.cols]) else NA_real_,
                   P_NBW = if (PROCESS.WITH.2.5KG) d$Mean[K_cats] else NA_real_,
                   P_NBW_lower = if (PROCESS.WITH.2.5KG) d$Lower[K_cats] else NA_real_,
                   P_NBW_upper = if (PROCESS.WITH.2.5KG) d$Upper[K_cats] else NA_real_,
                   Profile = d$Profile[1], stringsAsFactors = FALSE)
      }))
      rownames(risk_totals) <- NULL
      cat(sprintf("\n%s by subgroup and definition:\n", risk.label))
      print(risk_totals[, c("Subgroup","Definition","ClassIndex","Births","P_RISK","P_RISK_lower","P_RISK_upper")],
            row.names = FALSE)

      # =======================================================================
      # subgroup bar chart -- LBW categories only
      # =======================================================================
      plot.cols <- lbw.cols
      arm.tag <- if (PROCESS.WITH.2.5KG) "full" else "lbw"
      plot_pdf <- file.path(paths$plots,
                            sprintf("high_low_risk_pred_%s_%d.pdf", arm.tag, TARGET.YEAR))
      pdf(file = plot_pdf, width = 9.5, height = 8.5, pointsize = 12)
      par(mfrow = c(2, 1), mar = c(3.2, 4.8, 4.0, 1.5), oma = c(5.2, 0, 4.2, 0))

      y.top <- max(high.risk.result$upr[plot.cols], low.risk.result$upr[plot.cols]) * 1.16
      # draw.panel is a helper function to draw a single barplot panel for either the high-risk or low-risk subgroup. It takes the result object (means, lower, upper), title, profile description, bar color, and border color as arguments. It plots the means with error bars representing the 95% bootstrap intervals.
      draw.panel <- function(res, ttl, prof, col, border) {
        bars <- barplot(res$means[plot.cols], ylab = "Probability", ylim = c(0, y.top),
                        names.arg = cat.labels[plot.cols], col = col, border = border,
                        cex.names = 1.05, cex.axis = 1.05, cex.lab = 1.15)
        title(ttl, line = 1.4, cex.main = 1.3)
        mtext(prof, side = 3, line = 0.25, cex = 0.92, col = "#374151")
        arrows(bars, res$lwr[plot.cols], bars, res$upr[plot.cols],
               angle = 90, code = 3, length = 0.045, col = "black", lwd = 1.4)
        invisible(bars)
      }
      draw.panel(high.risk.result,
                 sprintf("Highest-risk class (i = %d, N = %s)", is.high.risk.indices,
                         format(round(class.births[is.high.risk.indices]), big.mark = ",")),
                 describe.profile(empirical.high.profile), "#EF4444", "#B91C1C")
      draw.panel(low.risk.result,
                 sprintf("Lowest-risk class (i = %d, N = %s)", is.low.risk.indices,
                         format(round(class.births[is.low.risk.indices]), big.mark = ",")),
                 describe.profile(empirical.low.profile), "#3B82F6", "#1D4ED8")

      mtext(sprintf("Within-LBW Birth Weight Category Probabilities (%d)", TARGET.YEAR),
            side = 3, line = 2.1, outer = TRUE, cex = 1.35, font = 2)
      mtext(sprintf("%s | posterior means with 95%% percentile bootstrap intervals (B = %d)",
                    model.label, B),
            side = 3, line = 0.5, outer = TRUE, cex = 1.0)
      
      # exclude NBW mass
      cap1 <- sprintf("Total %s:  high-risk %.4f [%.4f, %.4f]     low-risk %.4f [%.4f, %.4f]",
                      risk.label,
                      sum(high.risk.result$means[risk.cols]), sum(high.risk.result$lwr[risk.cols]),
                      sum(high.risk.result$upr[risk.cols]),
                      sum(low.risk.result$means[risk.cols]), sum(low.risk.result$lwr[risk.cols]),
                      sum(low.risk.result$upr[risk.cols]))
      cap2 <- if (PROCESS.WITH.2.5KG) {
        sprintf("Normal (>2500 g) category omitted from the axes:  high-risk %.4f [%.4f, %.4f]     low-risk %.4f [%.4f, %.4f]",
                high.risk.result$means[K_cats], high.risk.result$lwr[K_cats], high.risk.result$upr[K_cats],
                low.risk.result$means[K_cats],  low.risk.result$lwr[K_cats],  low.risk.result$upr[K_cats])
      } else {
        "LBW-only model: no Normal birth-weight category is estimated; the ten bars are conditional on LBW."
      }
      mtext(cap1, side = 1, line = 1.6, outer = TRUE, cex = 0.95)
      mtext(cap2, side = 1, line = 3.0, outer = TRUE, cex = 0.95, col = "#4B5563")

      dev.off()
      cat("Saved subgroup plot to:", plot_pdf, "\n")

      # =======================================================================
      # the extreme top 10 classes
      # =======================================================================
      # Reporting a single argmax invites the objection that it is a lucky cell. The most stark pair rather than a real pattern. Ranking the top and bottom TEN eligible classes and reading down the covariate columns answers that directly: a covariate whose value is constant across all ten is a genuine shared feature of the extreme, whereas one that varies freely is not doing the work regardless of how it looks in a single profile.
      # `consensus` in the data frame below is exactly that: for each covariate, the modal level within the set and the share of the ten that carry it. 10/10 is a unanimous marker; 5/10 on a binary covariate is noise.

      TOP.N <- 10
      elig.rank <- risk.ranking[risk.ranking$eligible, ]           # already sorted desc
      top.n     <- head(elig.rank, TOP.N)
      bottom.n  <- tail(elig.rank, TOP.N)
      bottom.n  <- bottom.n[order(bottom.n$p_lbw), ]

      consensus.of <- function(idx) {
        do.call(rbind, lapply(colnames(X.matrix), function(v) {
          vals <- as.character(X.matrix[idx, v])
          tb   <- sort(table(vals), decreasing = TRUE)
          data.frame(variable = v, modal = names(tb)[1],
                     n = as.integer(tb[1]), of = length(idx),
                     share = as.integer(tb[1]) / length(idx),
                     unanimous = as.integer(tb[1]) == length(idx),
                     stringsAsFactors = FALSE)
        }))
      }
      top.consensus    <- consensus.of(top.n$class)
      bottom.consensus <- consensus.of(bottom.n$class)

      risk.extremes <- rbind(
        transform(top.n,    group = "highest"),
        transform(bottom.n, group = "lowest")
      )
      rownames(risk.extremes) <- NULL

      cat(sprintf("\n--- Highest %d eligible classes ---\n", TOP.N))
      print(top.n[, c("class", "births", "p_lbw", "lower", "upper", "profile")], row.names = FALSE)
      cat(sprintf("\n--- Lowest %d eligible classes ---\n", TOP.N))
      print(bottom.n[, c("class", "births", "p_lbw", "lower", "upper", "profile")], row.names = FALSE)
      cat(sprintf("\nShared features of the %d HIGHEST-risk classes:\n", TOP.N))
      print(top.consensus, row.names = FALSE)
      cat(sprintf("\nShared features of the %d LOWEST-risk classes:\n", TOP.N))
      print(bottom.consensus, row.names = FALSE)
      cat("\nUnanimous markers -- highest:",
          paste(with(top.consensus[top.consensus$unanimous, ],
                     paste(variable, modal, sep = "=")), collapse = ", "), "\n")
      cat("Unanimous markers -- lowest :",
          paste(with(bottom.consensus[bottom.consensus$unanimous, ],
                     paste(variable, modal, sep = "=")), collapse = ", "), "\n")
      
      # save as csv to copy/paste into manuscript tables
      utils::write.csv(risk.extremes,
                       file.path(paths$results, sprintf("risk_extremes_top%d_%d.csv", TOP.N, TARGET.YEAR)),
                       row.names = FALSE)

      # ---- figure: profile grid for both extreme sets ----------------------
      # Rows are profiles (ordered by risk), columns are the seven covariates,
      # each cell printed with that profile's level. Cells in a column that
      # match the set's modal level are shaded, so a fully shaded column is a
      # unanimous marker and reads at a glance. The risk with its interval is
      # drawn as a bar to the right of each grid.
      ext_pdf <- file.path(paths$article,
                           sprintf("risk_extremes_profiles_%s_%d.pdf", arm.tag, TARGET.YEAR))
      pdf(ext_pdf, width = 16, height = 9.5, pointsize = 12)
      par(mfrow = c(1, 2), mar = c(4.4, 1.0, 4.6, 1.0), oma = c(3.4, 0, 3.6, 0))
      vars <- colnames(X.matrix); nv <- length(vars)

      panel <- function(tab, cons, ttl, accent) {
        n <- nrow(tab)
        # xlim starts negative to leave room for the i=<class> labels, which
        # were being clipped at the panel edge.
        plot.new(); plot.window(xlim = c(-1.8, nv + 5.2), ylim = c(n + 0.9, 0.1))
        text((seq_len(nv)) - 0.5, 0.45, vars, srt = 32, adj = c(0, 0.5), cex = 0.95, font = 2)
        text(nv + 2.0, 0.45, sprintf("%s (95%% PI)", risk.label), cex = 0.95, font = 2)
        rng <- range(elig.rank$p_lbw)
        for (r in seq_len(n)) {
          for (cc in seq_len(nv)) {
            val <- as.character(X.matrix[tab$class[r], vars[cc]])
            hit <- val == cons$modal[cc]
            rect(cc - 1, r - 0.42, cc, r + 0.42,
                 col = if (hit) adjustcolor(accent, alpha.f = 0.20) else "white",
                 border = "#E5E7EB")
            text(cc - 0.5, r, val, cex = 0.92, font = if (hit) 2 else 1,
                 col = if (hit) accent else "#6B7280")
          }
          # risk bar, on a common scale across both panels
          x0 <- nv + 0.4; w <- 2.9
          sc <- function(p) x0 + w * (p - rng[1]) / diff(rng)
          segments(sc(tab$lower[r]), r, sc(tab$upper[r]), r, col = "#9CA3AF", lwd = 1.4)
          points(sc(tab$p_lbw[r]), r, pch = 19, col = accent, cex = 1.15)
          text(nv + 5.2, r, sprintf("%.3f", tab$p_lbw[r]), adj = c(1, 0.5), cex = 0.92)
          text(-0.15, r, sprintf("i=%d", tab$class[r]), adj = c(1, 0.5), cex = 0.88, col = "#6B7280")
        }
        title(ttl, cex.main = 1.35, font.main = 2, line = 2.4)
        mtext(paste("unanimous:",
                    paste(with(cons[cons$unanimous, ], paste(variable, modal, sep = "=")),
                          collapse = ", ")),
              side = 3, line = 0.7, cex = 0.95, col = accent)
      }
      panel(top.n,    top.consensus,    sprintf("Highest %d risk profiles", TOP.N), "#B91C1C")
      panel(bottom.n, bottom.consensus, sprintf("Lowest %d risk profiles",  TOP.N), "#1D4ED8")
      mtext(sprintf("Extreme risk profiles and their shared features (%d)", TARGET.YEAR),
            side = 3, line = 1.0, outer = TRUE, cex = 1.5, font = 2)
      mtext(sprintf("%s | %s | shaded = matches that column's modal level within the set | eligible classes only (N >= %s births)",
                    model.label, risk.label, format(MIN.CLASS.BIRTHS, big.mark = ",")),
            side = 1, line = 0.9, outer = TRUE, cex = 0.95, col = "#4B5563")
      dev.off()
      cat("Saved extreme-profile grid to:", ext_pdf, "\n")

      # =======================================================================
      # —— two-pane tree showing where each extreme class lands
      # =======================================================================
      # follow the path of each extreme class through the tree, highlighting the nodes along that path. The leaf node itself is highlighted with a darker color. The tree is drawn twice, once for the highest-risk class and once for the lowest-risk class, side by side.
      path.to <- function(fit, class.index) {
        ids  <- as.numeric(rownames(fit$frame))
        leaf <- ids[fit$where[class.index]]
        p <- leaf; while (tail(p, 1) > 1) p <- c(p, tail(p, 1) %/% 2)
        rev(p)
      }
      hi.path <- path.to(proto.tree, is.high.risk.indices)
      lo.path <- path.to(proto.tree, is.low.risk.indices)

      pane_pdf <- file.path(paths$plots, sprintf("risk_paths_tree_%d.pdf", TARGET.YEAR))
      pdf(pane_pdf, width = 20, height = 8.5, pointsize = 15)
      par(mfrow = c(1, 2), oma = c(4.2, 0, 4.0, 0))
      draw.dm.tree(proto.tree,
                   title = sprintf("Highest risk: class %d", is.high.risk.indices),
                   subtitle = sprintf("%s = %.4f [%.4f, %.4f]  |  N = %s births", risk.label,
                                      class.mean[is.high.risk.indices],
                                      class.lwr[is.high.risk.indices],
                                      class.upr[is.high.risk.indices],
                                      format(round(class.births[is.high.risk.indices]), big.mark = ",")),
                   lbw.cols = lbw.cols, alphavec = alphavec,
                   highlight = hi.path, highlight.col = "#B91C1C",
                   leaf.highlight = tail(hi.path, 1))
      draw.dm.tree(proto.tree,
                   title = sprintf("Lowest risk: class %d", is.low.risk.indices),
                   subtitle = sprintf("%s = %.4f [%.4f, %.4f]  |  N = %s births", risk.label,
                                      class.mean[is.low.risk.indices],
                                      class.lwr[is.low.risk.indices],
                                      class.upr[is.low.risk.indices],
                                      format(round(class.births[is.low.risk.indices]), big.mark = ",")),
                   lbw.cols = lbw.cols, alphavec = alphavec,
                   highlight = lo.path, highlight.col = "#1D4ED8",
                   leaf.highlight = tail(lo.path, 1))
      mtext(sprintf("Extreme-risk profiles traced through the depth-%d DM-CART tree (%d)",
                    DISPLAY.DEPTH, TARGET.YEAR),
            side = 3, line = 1.4, outer = TRUE, cex = 1.5, font = 2)
      mtext(paste(sprintf("highest: %s", describe.profile(empirical.high.profile)),
                  sprintf("     lowest: %s", describe.profile(empirical.low.profile))),
            side = 1, line = 0.6, outer = TRUE, cex = 0.95, col = "#4B5563")

      mtext(sprintf(paste0("Path traced through the depth-%d display tree; the leaf's own LBW= value is that tree's ",
                           "estimate, while the quoted %s is from the maxdepth-%d bootstrap ensemble."),
                    DISPLAY.DEPTH, risk.label, dm.control$maxdepth),
            side = 1, line = 1.9, outer = TRUE, cex = 0.85, col = "#6B7280")
      dev.off()
      cat("Saved risk-path tree to:", pane_pdf, "\n")

      #-----------------------------------------------------------
      # Save Summary Results
      #-----------------------------------------------------------
      var.names <- colnames(X.matrix)
      var.usage <- matrix(0, nrow = B, ncol = length(var.names))
      colnames(var.usage) <- var.names

      for (b in 1:B) {
        if (length(tree.structure[[b]]$variables) > 0) {
          for (var in tree.structure[[b]]$variables) {
            if (var %in% var.names) var.usage[b, var] <- 1
          }
        }
      }

      var.freq.df <- data.frame(variable = var.names, frequency = colSums(var.usage) / B)
      var.freq.df <- var.freq.df[order(-var.freq.df$frequency), ]

      top.var.df <- as.data.frame(table(top.vars.used), stringsAsFactors = FALSE)
      colnames(top.var.df) <- c("variable", "count")
      top.var.df$frequency <- top.var.df$count / B
      top.var.df <- top.var.df[order(-top.var.df$frequency), ]

      # ---- assemble the real depth distribution ----
      depth_df <- do.call(rbind, depth.rows[!vapply(depth.rows, is.null, logical(1))])
      if (is.null(depth_df) || nrow(depth_df) == 0) {
        stop("No splits recorded across any bootstrap replicate -- depth_df would be empty. ",
             "Refusing to save; the depth figures must not fall back to invented data.")
      }
      cat(sprintf("\nRecorded %d splits across %d bootstrap trees (mean depth %.2f)\n",
                  nrow(depth_df), B, mean(depth_df$Depth)))
      print(aggregate(Depth ~ Variable, depth_df, function(z)
        c(n = length(z), mean = round(mean(z), 2))))

      bootstrap.tree.results <- list(
        var.usage      = var.usage,
        var.freq.df    = var.freq.df,
        top.var.df     = top.var.df,      
        n.vars.summary = summary(n.vars.used),
        risk_summary   = risk_summary,
        risk_totals    = risk_totals,     
        risk.ranking   = risk.ranking,    # all classes, ranked
        risk.extremes  = risk.extremes,   # top/bottom 10 with profiles
        top.consensus  = top.consensus,   # shared features, highest
        bottom.consensus = bottom.consensus,
        depth_df       = depth_df,
        B              = B,
        min.class.births = MIN.CLASS.BIRTHS,
        profiles       = list(
          empirical = list(high = empirical.high.profile, low = empirical.low.profile,
                           high.index = is.high.risk.indices, low.index = is.low.risk.indices,
                           high.tied = hi.pick$tied, low.tied = lo.pick$tied,
                           high.stability = high.stability, low.stability = low.stability),
          apriori   = list(high = high.risk.profile, low = low.risk.profile,
                           high.index = apriori.high.index, low.index = apriori.low.index))
      )

      boot_out_file <- paths$bootstrap.results(TARGET.YEAR)

      save(bootstrap.tree.results, risk_summary, risk_totals, risk.ranking,
           risk.extremes, top.consensus, bottom.consensus,
           depth_df, top.var.df, file = boot_out_file)
      cat("Saved bootstrap tree results to:", boot_out_file, "\n")

      depth_out_file <- paths$tree.depth.results(TARGET.YEAR)
      save(trees, tree.sum, var.across.depths, pred.stats.compare,
           file = depth_out_file)
      cat("Saved depth comparison to:", depth_out_file, "\n")
}

cat("\nCompleted DM-CART analysis and bootstrapping for all model variants.\n\n")

