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

# Resolve target year and cohort
init.paths  <- get_paths(year = PRIOR.YEAR, with.2.5kg = TRUE)
TARGET.YEAR <- init.paths$target.year   # 2024 -- the cohort being modeled
PRIOR.YEAR  <- init.paths$prior.year    # 2023 -- cutpoints + alphavec ONLY

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

    # Category index sets. Defined HERE (rather than further down, just before
    # the bootstrap, where they used to live) because mysplit()'s many-level
    # fallback branch references lbw.cols and the primary tree is fit before
    # that point. REVISION R2.
    lbw.cols <- 1:10                                              # Q1..Q10, the LBW region
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
        # =====================================================================
        # REVISION R8: node-label callbacks were producing useless labels
        # =====================================================================
        # rpart hands these callbacks `yval` = frame$yval2, which for a user
        # method is a MATRIX (one row per node, one column per birth-weight
        # category). The old `text` callback did:
        #     total.counts <- sum(yval)          # <- sums the WHOLE matrix
        #     which.max(total.counts)            # <- which.max of a scalar = 1
        # so every leaf in every plotted tree was labelled "1", regardless of
        # its actual composition. `print` had the matching problem and emitted
        # "Deviance: NULL".
        #
        # These now handle the matrix case and report the MODAL birth-weight
        # category per node, which is what the label was always meant to be.
        # =====================================================================
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
        text = function(yval, dev, wt, ylevel, digits, n, use.n){
          m     <- if (is.matrix(yval)) yval else matrix(yval, nrow = 1)
          modal <- apply(m, 1, which.max)          # per-node modal category
          # Posterior P(LBW) for the node, same (counts + alpha)/(N + alpha0)
          # form as dm.leaf.probs(). Folded into the LABEL rather than drawn as
          # a separate text() call, so text.rpart handles multi-line centring --
          # annotating separately made it collide with the "n=" line.
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


    myeval <- function(y, wt = NULL, parms = 1) {
      counts <- colSums(y)
      dev <- -log.dm.likelihood(counts, alpha = alphavec)
      return(list(label = counts, deviance = dev))
    }

    #-----------------------------------------------------------
    # mysplit -- exhaustive search over categorical partitions
    #-----------------------------------------------------------
    # WHAT WAS WRONG
    # --------------
    # The old version returned `direction <- ux`, the sorted raw level codes.
    # rpart's user-split protocol treats `direction` as an ORDERING and then
    # only ever considers the k-1 splits that are contiguous in that ordering:
    # first level vs rest, first two vs rest, and so on. For a binary
    # predictor that is complete (k=2 -> exactly 1 possible split), which is
    # why this was invisible while every predictor was binary.
    #
    # For a k-level factor there are 2^(k-1) - 1 distinct ways to split the
    # levels into two groups. Sorted-code order reaches only k-1 of them, and
    # WHICH k-1 depends on nothing but the order the levels happen to be
    # declared in. With 5-level race that is 4 of 15; a grouping like
    # {Black, Asian} | {White, Hispanic, Other} is simply unreachable.
    # Measured on a constructed case where that grouping is optimal:
    #     sorted-code order : best goodness =  35.86
    #     exhaustive        : best goodness = 174.95   (4.9x better)
    #
    # HOW IT IS FIXED
    # ---------------
    # rpart will not accept an arbitrary subset as a split, but it does accept
    # any PERMUTATION as `direction`. So we search all 2^(k-1)-1 partitions
    # ourselves, then hand back an ordering with the winning left-group placed
    # first -- which makes the optimum reachable as "first m vs rest". rpart
    # then writes exactly that partition into fit$csplit. Verified end to end.
    #
    # Binary predictors are mathematically unaffected (k=2 -> one partition),
    # so this changes nothing about any binary-only result.
    #
    # `continuous` is now honored. Every predictor here is a factor so the
    # branch is currently unused, but the old code ignored the flag entirely
    # and would have returned a categorical-format answer for a numeric
    # predictor, which rpart would have misread.
    # =========================================================================
    MAX.EXHAUSTIVE.LEVELS <- 12   # 2^11-1 = 2047 partitions per node: still trivial

    mysplit <- function(y, wt, x, parms, continuous = FALSE) {
      ux <- sort(unique(x))
      k  <- length(ux)

      # Guard: a node can contain a single distinct x value. The old
      # `for (i in 1:(length(ux) - 1))` became `1:0` and looped over i = 1, 0,
      # indexing ux[0] and returning a goodness vector of the wrong length.
      if (k < 2) return(list(goodness = numeric(0), direction = ux))

      # ---- REVISION R14: collapse to per-level totals ONCE, up front ----------
      # Every candidate split is a partition of the k LEVELS, so the only thing
      # any evaluation needs is each level's category totals. Computing that
      # once turns each candidate evaluation from "subset an n x K matrix and
      # colSums it" into "sum a handful of rows of a k x K matrix".
      #
      # At this node size (n up to 480, k up to 5) that is roughly a 100x
      # reduction in work per candidate, which is what pays for the exhaustive
      # search being 4x larger than the old sorted-order search. Measured net
      # effect vs. the ORIGINAL pre-revision splitter: faster, not slower.
      # ------------------------------------------------------------------------
      lev.sums   <- rowsum(y, group = match(x, ux), reorder = TRUE)   # k x K
      total      <- colSums(lev.sums)
      parent.ll  <- log.dm.likelihood(total, alpha = alphavec)

      # Goodness of sending the levels flagged in `in.left` to the left child.
      # This is the same Delta-L rule as before: L_left + L_right - L_parent,
      # which equals parent.dev - child.dev.
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
      bit <- 2^(0:(k - 1))

      if (k <= MAX.EXHAUSTIVE.LEVELS) {
        # Enumerate every bipartition. Integer `m` is a BITMASK over levels:
        # bit j set means level j goes left. We stop at 2^(k-1)-1 because a
        # partition and its complement describe the same split.
        #
        # REVISION R14: results are stored BY MASK so the k-1 contiguous splits
        # scored below are looked up instead of recomputed. Previously the
        # winning partition was evaluated twice -- which for a binary predictor
        # (1 partition) meant doing double the work of the original splitter,
        # on 5 of the 7 predictors. Now a binary predictor costs exactly ONE
        # evaluation, the same as the original code.
        n.mask <- 2^(k - 1) - 1
        gvals  <- numeric(n.mask)
        for (m in seq_len(n.mask)) gvals[m] <- goodness.of(as.logical(bitwAnd(m, bit)))

        best.m   <- which.max(gvals)
        best.idx <- which(as.logical(bitwAnd(best.m, bit)))
        ord.idx  <- c(best.idx, setdiff(seq_len(k), best.idx))   # winner first

        # A prefix of `ord.idx` is itself a mask; canonicalise it against its
        # complement (same split) so it indexes into the table we just built.
        full <- 2^k - 1
        goodness <- vapply(seq_len(k - 1), function(i) {
          mask <- sum(bit[ord.idx[seq_len(i)]])
          gvals[min(mask, full - mask)]
        }, numeric(1))

      } else {
        # Too many levels to enumerate. Fall back to ordering the levels by
        # their LBW share, so the contiguous splits rpart examines are at
        # least monotone in the quantity of interest rather than arbitrary.
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




    # =========================================================================
    # REVISION R7: dm.leaf.probs -- replaces the hand-rolled mypred() traversal
    # =========================================================================
    # WHAT WAS WRONG
    # --------------
    # The old mypred() walked the tree itself for each row, deciding direction
    # with:
    #     split_val <- fit$splits[split_rows[1], "index"]
    #     go_left   <- as.numeric(cur_val) <= split_val
    # Two independent errors:
    #
    #  1. For a CATEGORICAL split, fit$splits[, "index"] is not a cutpoint --
    #     it is the ROW NUMBER of fit$csplit, the table that actually records
    #     which levels go left (1 = left, 3 = right, 2 = not present at this
    #     node). Comparing a factor's integer code against a csplit row index
    #     is meaningless. EVERY predictor in this model is a factor, so this
    #     was wrong at every split -- including the binary ones. It was never
    #     a race-only problem.
    #
    #  2. split_rows[1] takes the FIRST fit$splits row bearing that variable
    #     name, so a variable reused at several nodes was always routed by its
    #     first split's parameters regardless of which node was being visited.
    #
    # Measured against correct leaf assignment on a controlled 4-leaf tree:
    # 16 of 20 rows landed in the wrong leaf; max probability error 0.612.
    # Yhat -- and so every high/low-risk figure derived from it -- was built
    # from this.
    #
    # HOW IT IS FIXED
    # ---------------
    # No traversal is needed. X.matrix is the COMPLETE enumeration of classes
    # and every tree here is fit on exactly those rows, so rpart has already
    # done the routing: fit$where holds the frame ROW INDEX (not the node id)
    # of the leaf each training row fell into. fit$frame$yval2 holds the
    # per-node count vector returned by myeval(). Indexing one by the other
    # gives each class its leaf's counts, in the row order of X.matrix.
    #
    # Two traps deliberately avoided:
    #   * predict.rpart() on a custom method returns an n x 1 yval, not the
    #     count matrix -- useless here.
    #   * rpart:::pred.rpart() SEGFAULTS on custom-method fits (verified: R
    #     aborts). Do not reach for it.
    #
    # This function is only valid for the rows the tree was FIT on. That is
    # all this design ever needs. To score genuinely new rows you would have
    # to decode csplit per node (not per variable) -- see REVISIONS.md R7.
    # =========================================================================
    dm.leaf.probs <- function(fit, alphavec, cat.names = colnames(Y.df)) {
      stopifnot(!is.null(fit$where), is.matrix(fit$frame$yval2))

      counts <- fit$frame$yval2[fit$where, , drop = FALSE]   # leaf counts, per class row
      probs  <- t(apply(counts, 1, function(z) (z + alphavec) / sum(z + alphavec))) # posterior P(LBW) per row

      colnames(probs) <- cat.names
      rownames(probs) <- NULL
      probs
    }


    dm.method <- list(init=myinit, eval=myeval, split=mysplit, method="dm")

    # =========================================================================
    # REVISION R9: plot.dm.tree() -- prototype tree rendering that actually works
    # =========================================================================
    # Three separate things break when you try to plot a custom-method rpart
    # fit, which is why this has been painful:
    #
    #  1. rpart.plot::rpart.plot(fit) FAILS with
    #     "missing value where TRUE/FALSE needed". The package inspects
    #     fit$method to decide how to format nodes and does not understand a
    #     user-supplied method object. There is no argument that fixes this --
    #     rpart.plot simply does not support custom methods. Use base
    #     plot.rpart + text.rpart, which DO work (verified).
    #
    #  2. Split labels default to rpart's compact letter code, so a race split
    #     renders as `race=bd` -- b and d being the 2nd and 4th LEVELS. Passing
    #     `pretty = 0` to text() expands these to `race=Black,Asian`. This
    #     matters much more now that R2 lets race split into genuine subsets.
    #
    #  3. Leaf labels were all "1" -- fixed in R8 above.
    #
    # =========================================================================
    # REVISION R15: rebuilt for legibility and for BRANCH-LEVEL interpretation
    # =========================================================================
    # This no longer uses text.rpart at all. text.rpart can only print the
    # condition for going LEFT, at the node -- so a 5-level race split rendered
    # as a bare "race=Black" and there was no way to see what went RIGHT, or to
    # confirm that the multi-level partitioning from R2 was being exercised.
    #
    # Instead we draw the skeleton with plot.rpart (which handles layout) and
    # then place every label ourselves at coordinates from rpart:::rpartco():
    #
    #     at each internal node : the SPLIT VARIABLE name
    #     on the left branch    : the levels sent left
    #     on the right branch   : the levels sent right
    #     at each leaf          : modal category, P(LBW), and N births
    #
    # The left/right level sets come from labels(fit, pretty = 0,
    # collapse = FALSE), which is rpart's own decoding of csplit -- so the
    # printed groupings are exactly what the model used, not a reconstruction.
    #
    # IMPORTANT -- what "n" means. rpart's frame$n counts the DATA ROWS at a
    # node, and here a data row is a CLASS, not a birth. `use.n = TRUE` was
    # therefore printing "n=32" meaning 32 predictor classes, which reads as a
    # sample size and is off by ~5 orders of magnitude. The label now reports
    # N = sum of the node's category counts (actual births), with the class
    # count given separately as "cls".
    # =========================================================================
    # Tree-drawing helpers (dm.tree.coords / draw.dm.tree / plot.dm.tree)
    # now live in util.R -- REVISION R27 -- so visual.R can render the same
    # figures from the saved tree objects without duplicating the code.

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

    # Save Main Tree
      dm.tree.file <- paths$dm.tree.rebin(TARGET.YEAR)
      save(dm.tree, file = dm.tree.file)
      cat("Saved decision tree to:", dm.tree.file, "\n")

      #-----------------------------------------------------------
      # E. Prototype tree + depth-controlled comparison (REVISION R9/R10)
      #    Deliberately placed BEFORE the bootstrap: these are the figures you
      #    read to understand the model, and they cost seconds, whereas the
      #    bootstrap costs hours. Nothing below depends on the bootstrap.
      #-----------------------------------------------------------
      cat("\n--- Rendering prototype tree ---\n")

      # The full tree (maxdepth = 8, ~480 classes) is unreadable when drawn, so
      # the figure for the thesis is a depth-3 refit of the SAME criterion on
      # the SAME observed counts. This is a display choice, not a model choice:
      # dm.tree above is untouched and is what gets saved.
      DISPLAY.DEPTH <- 3
      proto.tree <- rpart(
        formula = Y.df ~ .,
        data    = data.frame(X.matrix),
        method  = dm.method,
        control = rpart.control(minsplit = 2, cp = 0, maxdepth = DISPLAY.DEPTH,
                                xval = 0, usesurrogate = 0)
      )

      # REVISION R27: figure count cut for the article. dm-cart.R now writes
      # ONE annotated tree per arm (the fitted maxdepth-8 model). The depth-3
      # "prototype" PDF and the four depth-2..5 PDFs are no longer written --
      # visual.R draws a single depth-comparison panel from the saved `trees`
      # object instead. `proto.tree` is still FITTED, because the risk-path
      # figure traces through it; it just no longer gets its own file.
      plot.dm.tree(
        dm.tree,
        file     = file.path(paths$plots, sprintf("full_tree_%d.pdf", TARGET.YEAR)),
        title    = sprintf("DM-CART (Full Tree, maxdepth = %d)", dm.control$maxdepth),
        subtitle = sprintf("%s | cohort %s | %d classes",
                           model.label, init.paths$cohort, nrow(X.matrix)),
        lbw.cols = lbw.cols, alphavec = alphavec
      )

      # ---- REVISION R15: split-partition report --------------------------------
      # Prints the actual grouping chosen at every internal node, so you can
      # confirm from the console (not just the figure) that the multi-level
      # factor search from R2 is being exercised -- e.g. that race is not
      # merely peeling off one level at a time, and that the 3-level dmar is
      # genuinely grouping "Unknown" with one side or the other.
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
        # Was this partition reachable by the OLD sorted-code search? That
        # search could only take a contiguous PREFIX of the levels present at
        # the node, ordered by their factor codes. Anything else is a split
        # that only the R2 exhaustive search can find. This is the concrete
        # test of whether the multi-level factor handling is doing any work.
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
        cat("  old sorted-code search, i.e. they exist only because of REVISION R2:\n")
        if (nrow(new)) {
          for (r in seq_len(nrow(new)))
            cat(sprintf("    node %-4d %-10s  %s  |  %s\n",
                        new$node[r], new$variable[r], new$left[r], new$right[r]))
        } else cat("    (none at this fit)\n")
        invisible(out)
      }
      split.report <- report.splits(dm.tree, sprintf("full tree, %s", model.label))

      # ---- Depth-controlled comparison (REVISION R10) ----
      # Reinstates the Section 3.2.4 / Figures 3.4-3.5 analysis. The previous
      # version of this lived entirely in the commented-out block after the end
      # of the model loop -- there were ZERO active lines after it, so nothing
      # produced these figures and nothing ever wrote tree.depth.results.
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

        # Posterior P(LBW) per class under this depth, via the SAME correct
        # leaf lookup as the bootstrap (REVISION R7). The old dead code called
        # predict(..., type="matrix") here, which does not work on a custom
        # method -- another reason this block never ran.
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

        # (no per-depth PDF: visual.R draws all four as one panel -- R27)
      }

      # Which variables survive at each depth (Table 3.1-style companion)
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
      # REVISION R13: B is now settable without editing the file, so the final
      # 10,000-tree run is a flag rather than a code change that can be
      # forgotten. Interactive sessions get a fast default; batch runs get the
      # publication value.
      #   DMCART_B=200 Rscript dm-cart.R 2023-2024      # smoke test
      #   DMCART_B=10000 Rscript dm-cart.R 2023-2024    # publication run
      B <- as.integer(Sys.getenv("DMCART_B",
                                 unset = if (interactive()) "200" else "10000"))
      stopifnot(is.finite(B), B >= 2)
      n.rows <- nrow(counts.df)

      # (lbw.cols / nbw.cols moved above, before the primary fit -- REVISION R2)

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
      depth.rows <- vector("list", B)   # REVISION R11: real split depths

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

        # ---- REVISION R11: record the depth of every split, per replicate ----
        # This is the input for the Fig. 3.11/3.12 depth-distribution plots.
        # Previously nothing wrote it, so visual.R fell back to
        # `sample(0:6, 5000, replace = TRUE, prob = runif(7))` -- i.e. those
        # figures were random noise, with no warning that it had happened.
        #
        # rpart names frame rows by BINARY HEAP node id (1, 2, 3, 6, 7, ...),
        # so depth = floor(log2(id)) and the root is depth 0. It must be read
        # from rownames(frame), NOT from the row position: the old dead code
        # used `floor(log2(which(frame$var == "cig_0")))`, and row position
        # equals node id only for a perfectly full tree, so that was wrong
        # wherever it appeared.
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

        # REVISION R7: was mypred(dm.tree.b, newdata = data.frame(X.matrix), ...),
        # which misrouted rows at every factor split. The bootstrap tree is fit
        # on model.data, whose predictor rows ARE X.matrix in the same order,
        # so fit$where already gives the correct leaf for every class.
        Yhat[[b]] <- dm.leaf.probs(dm.tree.b, alphavec)
        setTxtProgressBar(pb, b)
      }
      close(pb)

      #-----------------------------------------------------------
      # High/Low Risk Subgroup Prediction Plot
      #-----------------------------------------------------------
      # =======================================================================
      # REVISION R12: high/low-risk subgroups are now single CLASSES
      # =======================================================================
      # The old conditions pinned only 3 of the 7 predictors:
      #     which(X.matrix$race == "Black" & X.matrix$dmar == 0 & X.matrix$cig_0 == 1)
      # which selects 32 of the 480 rows, and the block below then AVERAGED
      # over every combination of the remaining 4 predictors. That is not the
      # single 7-dimensional profile the thesis defines (classes i=69 / i=28);
      # it produces attenuated, blended probabilities instead of the sharp
      # contrast reported in Tables 3.2/3.3.
      #
      # Each profile below therefore fixes ALL 7 predictors, and find.class()
      # asserts that exactly one class matches -- so if the design changes
      # (another race level, a DMAR.POLICY switch) this fails loudly instead of
      # silently averaging again.
      #
      # NOTE ON THE CHOSEN PROFILES: the thesis's i=69 / i=28 are indices into
      # the old 128-class binary design and do not survive the re-encoding, so
      # they cannot be carried over mechanically. These are the profiles
      # proposed in REVISION-PLAN.md §5 Q2 and are ASSUMED pending your
      # confirmation -- edit the two lists here to change them; nothing else
      # needs to move.
      # =======================================================================
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

      # =======================================================================
      # REVISION R17: find the extreme-risk profiles FROM THE DATA
      # =======================================================================
      # The two profiles above are a clinical prior -- someone's reading of
      # which covariate values ought to be worst and best. That is an
      # assumption, and the model is perfectly capable of answering the
      # question instead. Here every one of the I classes is ranked by its
      # posterior P(LBW), pooled over the bootstrap ensemble:
      #
      #   lbw.by.class[i, b] = sum_{k in LBW} Yhat[[b]][i, k]
      #
      # Yhat is the leaf posterior each class inherits (REVISION R7), so this
      # is the model's fitted LBW risk for that exact covariate combination,
      # with the bootstrap giving its sampling distribution. The extremes are
      # then argmax / argmin of the posterior MEAN, with percentile bounds.
      #
      # SUPPORT FILTER -- this matters. A class with few or no observed births
      # contributes almost nothing to its leaf, so its "risk" is really its
      # neighbours' risk plus the Dirichlet prior; letting such a class win the
      # argmax would report an artefact as a finding. Classes are therefore
      # required to have at least MIN.CLASS.BIRTHS observed births. Rank
      # stability is separately reported: how often each class is the extreme
      # across replicates, which is the honest measure of whether the winner is
      # distinguishable from its rivals at all.
      # =======================================================================
      MIN.CLASS.BIRTHS <- 1000

      # ---- the ranking metric has to differ by model arm ----------------------
      # In the FULL model, K = 11 and P(LBW) = sum over Q1..Q10 is the natural
      # risk score. In the LBW-ONLY model there is no Normal category, so
      # lbw.cols spans EVERY category and sum(P[lbw.cols]) is identically 1 for
      # every class -- ranking on it is vacuous (it produced P(LBW)=1.0000 for
      # all 480 classes, with the "highest" and "lowest" class being the same
      # one). That arm's probabilities are conditional on being LBW, so the
      # meaningful severity score is the mass in the smallest deciles.
      if (PROCESS.WITH.2.5KG) {
        risk.cols  <- lbw.cols          # Q1..Q10 out of 11
        risk.label <- "P(LBW)"
        risk.desc  <- "probability of low birth weight"
      } else {
        risk.cols  <- 1:3               # three smallest deciles, conditional on LBW
        risk.label <- "P(Q1-Q3 | LBW)"
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

      # ---- TIES ARE THE NORMAL CASE, NOT AN EDGE CASE ----------------------
      # A tree assigns every class in the same terminal node the SAME posterior,
      # so classes sharing a leaf have byte-identical P(LBW) and there is no
      # "the" highest-risk class -- only a highest-risk LEAF, whose members the
      # model cannot tell apart. Taking which.max() alone would silently pick
      # whichever tied class happened to come first in the enumeration and
      # present it as a finding.
      #
      # Instead: identify the tied set, then report the member with the most
      # observed births as the representative profile (the one describing the
      # largest real population), and carry the full tied set alongside.
      TIE.TOL <- 1e-10
      pick.extreme <- function(want.max) {
        v      <- class.mean[eligible]
        target <- if (want.max) max(v) else min(v)
        tied   <- eligible[abs(v - target) <= TIE.TOL]
        list(index = tied[which.max(class.births[tied])], tied = tied)
      }
      hi.pick <- pick.extreme(TRUE);  lo.pick <- pick.extreme(FALSE)
      is.high.risk.indices <- hi.pick$index
      is.low.risk.indices  <- lo.pick$index

      # Stability, also tie-aware: the question is how often the chosen class is
      # AMONG the extreme classes in a replicate, not whether it happens to win
      # an arbitrary tie-break. Comparing against each replicate's own max/min
      # counts a tie as a success, which is the honest reading.
      hi.col.max <- apply(lbw.by.class[eligible, , drop = FALSE], 2, max)
      lo.col.min <- apply(lbw.by.class[eligible, , drop = FALSE], 2, min)
      high.stability <- mean(abs(lbw.by.class[is.high.risk.indices, ] - hi.col.max) <= TIE.TOL)
      low.stability  <- mean(abs(lbw.by.class[is.low.risk.indices,  ] - lo.col.min) <= TIE.TOL)

      profile.of <- function(i) as.list(setNames(as.character(unlist(X.matrix[i, ])),
                                                 colnames(X.matrix)))
      empirical.high.profile <- profile.of(is.high.risk.indices)
      empirical.low.profile  <- profile.of(is.low.risk.indices)

      # Ranking table, saved so the thesis can quote the top/bottom of it.
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

      cat(sprintf("\n--- Empirical risk extremes (REVISION R17) ---\n"))
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

      # Rank the assumed profiles among ALL classes, not just the eligible
      # ones -- an a priori profile can easily fall below the support
      # threshold (the assumed high-risk one does), and reporting rank NA
      # there would read as a bug rather than as the finding it is.
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

      # With a single class per subgroup, colMeans() over one row is that row.
      # Kept in this form so the function also works if a future profile
      # deliberately spans several classes.
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

      # ---- REVISION R12/R17: persist risk_summary ----
      # visual.R looks for an object of exactly this name and, not finding it,
      # silently substituted hardcoded placeholders -- c(rep(0.016, 10), 0.840)
      # and c(rep(0.009, 10), 0.910) -- i.e. the headline subgroup figures were
      # invented constants, not model output. Nothing wrote this object before.
      #
      # Lower/Upper are PERCENTILE bootstrap bounds, carried explicitly.
      # visual.R used to reconstruct intervals as Mean +/- 1.96*SE, so the same
      # quantity had two different intervals in the same document. SE is kept
      # alongside for reference but the percentile bounds are the ones to plot.
      #
      # R17 adds a `Definition` column. "empirical" rows are the classes the
      # model actually ranks highest/lowest; "a priori" rows are the assumed
      # clinical profiles, retained so the two can be compared directly rather
      # than one silently replacing the other.
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

      # Aggregate P(LBW) per subgroup, carried alongside so the figures and the
      # text quote the same numbers.
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
      # REVISION R18: subgroup bar chart -- LBW categories only
      # =======================================================================
      # The Normal (NBW) bar is ~0.77-0.94 and dwarfs the ten LBW bars, so the
      # within-LBW shape -- which is the entire point of expanding the LBW
      # region into deciles -- was visually flattened to nothing. NBW is now
      # dropped from the axes and reported as text in the subcaption instead,
      # where the number is still available without costing the reader the
      # comparison they actually need. For the LBW-only model there is no NBW
      # category and the subcaption says so.
      plot.cols <- lbw.cols
      # REVISION R29: per-arm figures that ARE article figures go to the
      # cohort-level article/ folder, tagged with the arm so the two do not
      # collide on filename. Diagnostics stay in plots/.
      # REVISION R31: the ARTICLE version of this figure is redrawn in visual.R
      # (ggplot, matching the other article figures). What dm-cart.R writes here
      # is the base-graphics diagnostic, so it stays in plots/.
      arm.tag <- if (PROCESS.WITH.2.5KG) "full" else "lbw"
      plot_pdf <- file.path(paths$plots,
                            sprintf("high_low_risk_pred_%s_%d.pdf", arm.tag, TARGET.YEAR))
      pdf(file = plot_pdf, width = 9.5, height = 8.5, pointsize = 12)
      par(mfrow = c(2, 1), mar = c(3.2, 4.8, 4.0, 1.5), oma = c(5.2, 0, 4.2, 0))

      y.top <- max(high.risk.result$upr[plot.cols], low.risk.result$upr[plot.cols]) * 1.16
      # Title and profile are separate lines: the 7-term profile string does not
      # fit on a title line at a readable size and was running off the panel.
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

      # Subcaption: the excluded NBW mass, plus the aggregate LBW risk, so
      # nothing is lost by dropping the eleventh bar.
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
      # REVISION R9: the `par(mfrow = c(1, 1))` that used to sit here has been
      # removed. Calling par() when no device is open OPENS one -- so this line,
      # running immediately after dev.off(), silently created a stray
      # `Rplots.pdf` in the working directory on every batch run. The par()
      # settings it was trying to reset died with the device on the line above.
      cat("Saved subgroup plot to:", plot_pdf, "\n")

      # REVISION R28: the risk-profile landscape moved to visual.R, where both
      # model arms are drawn as one two-pane figure instead of two single panes.


      # =======================================================================
      # REVISION R23: the extreme TEN, not just the extreme one
      # =======================================================================
      # Reporting a single argmax invites the objection that it is a lucky cell
      # -- the most stark pair rather than a real pattern. Ranking the top and
      # bottom TEN eligible classes and reading down the covariate columns
      # answers that directly: a covariate whose value is constant across all
      # ten is a genuine shared feature of the extreme, whereas one that varies
      # freely is not doing the work regardless of how it looks in a single
      # profile.
      #
      # `consensus` below is exactly that: for each covariate, the modal level
      # within the set and the share of the ten that carry it. 10/10 is a
      # unanimous marker; 5/10 on a binary covariate is noise.
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

      cat(sprintf("\n--- Highest %d eligible classes (REVISION R23) ---\n", TOP.N))
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

      # CSV alongside the RData so the table can go straight into the thesis.
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
      # REVISION R19: two-pane tree showing where each extreme class lands
      # =======================================================================
      # Answers "what do the high/low-risk groups look like ON the tree" by
      # tracing the root-to-leaf path each class actually follows and dimming
      # everything else. The path is recovered from the leaf's node id -- rpart
      # numbers frame rows by binary-heap id, so repeatedly halving the leaf's
      # id walks straight back to the root. The leaf itself comes from
      # dm.tree$where (REVISION R7), so this is the same routing the model used
      # rather than a re-derivation that could disagree with it.
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
      # The path is traced through the depth-3 DISPLAY tree (legible), while the
      # quoted risk comes from the maxdepth-8 bootstrap ensemble -- so the leaf's
      # own LBW= value will not equal the headline number. Said outright rather
      # than left for a reader to trip over.
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

      # ---- REVISION R13: restore top.var.df ----
      # `top.vars.used[b]` (the first/root split variable of each bootstrap
      # tree) was computed in the loop and then thrown away -- it never entered
      # bootstrap.tree.results, so Table 3.1's "Initial Split Variable" row had
      # no source. The 2020-2021 results on disk DO contain a `top.var.df`,
      # so this was lost in a rewrite rather than never existing.
      top.var.df <- as.data.frame(table(top.vars.used), stringsAsFactors = FALSE)
      colnames(top.var.df) <- c("variable", "count")
      top.var.df$frequency <- top.var.df$count / B
      top.var.df <- top.var.df[order(-top.var.df$frequency), ]

      # ---- REVISION R11: assemble the real depth distribution ----
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
        top.var.df     = top.var.df,      # REVISION R13
        n.vars.summary = summary(n.vars.used),
        risk_summary   = risk_summary,    # REVISION R12/R17
        risk_totals    = risk_totals,     # REVISION R17
        risk.ranking   = risk.ranking,    # REVISION R17: all classes, ranked
        risk.extremes  = risk.extremes,   # REVISION R23: top/bottom 10 with profiles
        top.consensus  = top.consensus,   # REVISION R23: shared features, highest
        bottom.consensus = bottom.consensus,
        depth_df       = depth_df,        # REVISION R11
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
      # `risk_summary` and `depth_df` are saved as TOP-LEVEL objects as well as
      # inside the list, because visual.R loads this file and looks them up by
      # bare name. Saving them only inside bootstrap.tree.results would leave
      # visual.R's `exists("risk_summary")` check failing exactly as before.
      save(bootstrap.tree.results, risk_summary, risk_totals, risk.ranking,
           risk.extremes, top.consensus, bottom.consensus,
           depth_df, top.var.df, file = boot_out_file)
      cat("Saved bootstrap tree results to:", boot_out_file, "\n")

      # ---- REVISION R10: persist the depth-controlled comparison ----
      # Nothing wrote this file before: its only writer lived in the commented
      # -out block after the end of the model loop, so paths$tree.depth.results
      # was a path to a file that was never created.
      depth_out_file <- paths$tree.depth.results(TARGET.YEAR)
      save(trees, tree.sum, var.across.depths, pred.stats.compare,
           file = depth_out_file)
      cat("Saved depth comparison to:", depth_out_file, "\n")
}

cat("\nCompleted DM-CART analysis and bootstrapping for all model variants.\n\n")


# --- YOUR EXPLORATORY PLOTTING, PRESERVED BUT DISABLED (REVISION R9) ---------
# These four lines ran at top level, after the model loop had already ended.
# Three problems with leaving them active:
#   1. `plot()` with no open device makes Rscript dump a stray `Rplots.pdf`
#      into the working directory on every batch run.
#   2. `path.rpart()` printed to the console AFTER the "Completed" banner,
#      which reads like the script is still doing something.
#   3. `dm.tree` here is whichever model arm the loop happened to finish on
#      (the LBW-only model), not the full model the comments assume -- and
#      the comments reference `mrace15`, which no longer exists (REVISION R4/R6
#      renamed the design; race is now a 5-level factor named `race`).
#
# The working equivalent now runs INSIDE the loop, per model arm, before the
# bootstrap: see plot.dm.tree() and the "Rendering prototype tree" section.
# It writes proper PDFs, expands categorical splits to level names, and
# labels leaves with P(LBW).
#
# path.rpart() itself works fine on a custom-method fit -- it is genuinely
# useful for reading off the rule for one node. Kept here as a recipe:
#
# load(paths$dm.tree.rebin(TARGET.YEAR))       # restores `dm.tree`
# plot(dm.tree, main = "DM Tree", uniform = TRUE)
# text(dm.tree, use.n = TRUE, cex = 1, font = 3, pretty = 0)   # pretty=0 -> level names
# path.rpart(dm.tree, node = 2)   # rule for the left  child of the root
# path.rpart(dm.tree, node = 3)   # rule for the right child of the root
# ----------------------------------------------------------------------------





# # } # end of for loop over PROCESS.WITH.2.5KG


# # NEW ABOVE
# # ...
# # OLD BELOW

# # check what is loaded
# print(ls())
# print(alphavec); length(alphavec)
# Y.df <- counts.df; response.cols <- colnames(Y.df)

# #-----------------------------------------------------------
# # B. Simplified Dirichlet-Multinomial Functions
# #-----------------------------------------------------------
# log.dm.likelihood <- function(counts, alpha = 1) {
#  # cat("** counts: ", counts, "\n")
#   # 'counts': an integer vector (n_1, ..., n_K)
#   # 'alpha':  scalar Dirichlet hyperparameter (assuming alpha_k = alpha for all k)
#   #
#   # Returns the log-likelihood of the Dirichlet-Multinomial model:
#   #   log p(x | alpha) = log Gamma(alpha_0) - log Gamma(N+alpha_0)
#   #                     + sum_k [ log Gamma(n_k + alpha) - log Gamma(alpha) ]
#   # where alpha_0 = K * alpha,  N = sum(counts), and K = length(counts).
  
#   N <- sum(counts)
#   K <- length(counts)
#   #alpha <- alpha*rep(1,K)
  
#   # cat("N: ", N, "\n")
#   # cat("K: ", K, "\n")

#   # sum of alpha over all categories
#   alpha_0 <- sum(alpha)

#   # Term 1: log Gamma(alpha_0) + log Gamma(N+1) - log Gamma(N + alpha_0)
#   term1 <- lgamma(alpha_0) + 0*lgamma(N + 1) - lgamma(N + alpha_0)
  
#   # Term 2: sum over k of [log Gamma(n_k + alpha) - log Gamma(alpha) - log Gamma(n_k + 1)]
#   term2 <- sum(lgamma(counts + alpha) - lgamma(alpha) - 0*lgamma(counts + 1))
  
#   # Total log-likelihood
#   ll <- term1 + term2

#   # if numeric issues occur
#   if (is.na(ll) || is.infinite(ll)) {
#     ll <- -Inf
#   }

#   return(ll)
# }

# #-----------------------------------------------------------
# # C. Define Custom rpart Method (DM with log-likelihood deviance)
# #-----------------------------------------------------------
# #----------------- 1) init() ------------------#
# myinit <- function(y, offset, parms=NULL, wt=NULL) {
#   # y: matrix of response counts 
#   # offset: not used
#   # parms: list containing alpha parameter
#   # wt: not used
#   # cat("** myinit() called. Y dimension:", dim(y), "\n")

#   list(
#     method="dm",
#     y = y,
#     parms = list(alpha=1),
#     numresp = ncol(y),
#     numy = ncol(y),
#     summary = function(yval, dev, wt, ylevel, digits) {
#       cat("DM Model\n")
#       cat("Deviance: ", format(dev,digits=digits), "\n")
#       cat("Counts: ", paste(format(yval, digits=digits), collapse=", "), "\n")
#     },
#     print = function(yval, dev, wt, ylevel, digits) {
#       cat("Deviance:", format(dev, digits = digits), "\n")
#     },
#     text = function(yval, dev, wt, ylevel, digits, n, use.n){
#       # yval is the colsum from myeval(); a vector of category counts.
#       # keep labels short.
#       total.counts <- sum(yval)
#       lbl <- paste0(format(which.max(total.counts), digits=4))
#       if (use.n) lbl <- paste(lbl, "\nn", n)
#       return(lbl)
#     }
#   )
# }

# #----------------- 2) eval() ------------------#
# myeval <- function(y, wt=NULL, parms=1) {
#   # y: response matrix subset for current node
#   # wt: weights (not used here)
#   if (nrow(y)==1) {counts = y}else{
#     counts <- colSums(y)}
#   # dev <- dm.deviance(counts,alpha=1) # return -log.dm.likelihood
#   dev <- -log.dm.likelihood(counts, alpha=alphavec) # counts is a vector!
#   # 'label' can be something to print for the node. We'll just store colSums  
#   return(list(label=counts, deviance=dev))
# }

# #----------------- 3) spit() ------------------#
# mysplit <- function(y, wt, x, parms, continuous = FALSE) {
  
#   if (nrow(y)==1) {counts = y}else{
#     counts <- colSums(y)}

#   parent.dev <- -log.dm.likelihood(counts, alpha=alphavec)
  
#   # Suppose 'x' is a factor or a binary 0/1 variable
#   ux <- sort(unique(x)) # unique values of x
#   goodness <- numeric(length(ux) - 1)# store improvement in deviance
#   direction <- ux
  
#   for(i in 1:(length(ux) - 1)) {
#     split.val <- ux[i]
#     left.idx <- x == split.val
#     right.idx <- !left.idx
    
#     if (sum(left.idx)>1){
#          left.counts  <- colSums(y[left.idx,  , drop=FALSE])}else{left.counts <- y[left.idx, , drop = FALSE]}
         
         
#          if (sum(right.idx)>1){
#            right.counts  <- colSums(y[right.idx,  , drop=FALSE])}else{right.counts <- y[right.idx, , drop = FALSE]}
         
    
    
#     child.dev <- -log.dm.likelihood(left.counts,alpha=alphavec) - log.dm.likelihood(right.counts,alpha=alphavec)

#     goodness[i] <- parent.dev - child.dev
#   }
#   return(list(goodness = pmax(goodness, 0), direction = direction))

# }
 
# #----------------- 4) pred() ------------------#
# mypred <- function(fit, newdata = NULL, type = c("vector", "prob", "class", "matrix")) {
#   type <- match.arg(type)
  
#   # Create a modified copy of the tree to avoid changing the original
#   fit_copy <- fit
  
#   # Check and repair the var field for the root node if needed
#   root_row <- which(rownames(fit_copy$frame) == "1")
#   if (length(root_row) > 0) {
#     if (is.na(fit_copy$frame$var[root_row]) || fit_copy$frame$var[root_row] == "<leaf>") {
#       cat("WARNING: Root node var is NA or leaf - attempting to fix\n")
      
#       # Try to infer the root variable from splits
#       split_vars <- rownames(fit_copy$splits)
#       if (length(split_vars) > 0) {
#         # Find the first variable used in splits
#         root_var <- split_vars[1]
#         fit_copy$frame$var[root_row] <- root_var
#         cat("Setting root node var to:", root_var, "\n")
#       } else {
#         cat("ERROR: Could not determine root variable from splits\n")
#       }
#     }
#   }
  
  
#   # Determine terminal nodes for each observation
#   if (is.null(newdata)) {
#     where <- fit_copy$where
#   } else {
#     newdata <- as.data.frame(newdata)
#     where <- integer(nrow(newdata))
    
#     for (i in 1:nrow(newdata)) {
#       # Start at root
#       node <- 1
      
#       while (TRUE) {
#         node_str <- as.character(node)
#         if (!node_str %in% rownames(fit_copy$frame)) break
        
#         node_row <- which(rownames(fit_copy$frame) == node_str)
#         varname <- fit_copy$frame$var[node_row]
        
#         if (is.na(varname) || varname == "<leaf>") break
        
#         # Find split information for this variable
#         split_rows <- which(rownames(fit_copy$splits) == varname)
#         if (length(split_rows) == 0) break
        
#         # Get the split value
#         split_val <- fit_copy$splits[split_rows[1], "index"]
        
#         # Get current observation's value for the split variable
#         if (!varname %in% colnames(newdata)) {
#           cat("WARNING: Variable", varname, "not found in newdata for observation", i, "\n")
#           break
#         }
        
#         cur_val <- newdata[i, varname]
        
#         # Determine direction based on split
#         go_left <- !is.na(cur_val) && as.numeric(cur_val) <= split_val
        
#         # Move to child node or stop if not available
#         child_node <- if (go_left) 2*node else 2*node + 1
#         child_str <- as.character(child_node)
        
#         if (child_str %in% rownames(fit_copy$frame)) {
#           node <- child_node
#         } else {
#           break
#         }
#       }
      
#       where[i] <- node
#     }
#   }
  
#   # Create result matrix
#   result_matrix <- matrix(0, nrow=length(where), ncol=ncol(Y.df))
#   colnames(result_matrix) <- colnames(Y.df)
  
#   # Extract counts for each observation
#   unique_nodes <- unique(where)
#   # cat("Number of unique terminal nodes:", length(unique_nodes), "\n")
  
#   for (i in 1:length(where)) {
#     node_str <- as.character(where[i])
#     row_idx <- which(rownames(fit_copy$frame) == node_str)
    
#     if (length(row_idx) > 0) {
#       if (!is.null(fit_copy$frame$yval2)) {
#         # Try to get counts from yval2
#         if (is.matrix(fit_copy$frame$yval2) && ncol(fit_copy$frame$yval2) == ncol(Y.df)) {
#           result_matrix[i,] <- fit_copy$frame$yval2[row_idx, ]
#         }
#       } else if ("label" %in% names(fit_copy$frame)) {
#         # Try to get label from frame
#         node_label <- fit_copy$frame$label[row_idx]
#         if (is.list(node_label)) node_label <- unlist(node_label)
#         if (length(node_label) == ncol(Y.df)) {
#           result_matrix[i, ] <- node_label
#         }
#       }
#     }
#   }
  
#   # Apply Dirichlet smoothing and return appropriate format
#   if (type == "matrix") {
#     return(result_matrix)
#   }
  
#   probs <- t(apply(result_matrix, 1, function(counts) {
#     (counts + alphavec) / sum(counts + alphavec)
#   }))
#   colnames(probs) <- colnames(Y.df)
  
#   unique_probs <- nrow(unique(round(probs, 8)))
#   # cat("Number of unique probability patterns:", unique_probs, "\n")
  
#   if (type == "vector" || type == "prob") {
#     return(probs)
#   } else if (type == "class") {
#     return(apply(probs, 1, which.max))
#   }
# }

# #-----------------------------------------------
# # original
# #-----------------------------------------------
# #----------------- 4) pred() ------------------#
# # mypred <- function(fit, newdata = NULL, type = c("vector", "prob", "class", "matrix")) {
# #   # 1. argument setup
# #   type <- match.arg(type) # ensure type is valid
# #   
# #   # 2. determine terminal nodes i.e. "where" for each observation 
# #   if (is.null(newdata)){
# #     # no new data provided, use terminal node assignment from fitting
# #     where <- fit$where
# #   } else {
# #     newdata <- as.data.frame(newdata) 
# #     where <- integer(nrow(newdata)) # initialize node assignments
# #     
# #     for (i in 1:nrow(newdata)){
# #       # start at root
# #       node <- 1
# #       
# #       while (TRUE){ # loop until leaf
# #         if (fit$frame$var[node] == "<leaf>") break # reached leaf
# #         
# #         split.var <- as.character(fit$frame$var[node]) # split on this variable
# #         split.rows <- which(rownames(fit$splits) == split.var) # find split info
# #         
# #         if(length(split.rows) == 0) break # no split info: treat as leaf
# #         
# #         cur.val <- newdata[i, split.var] # get current value 
# #         split.row <- split.rows[1] # first split info
# #         split.val <- fit$splits[split.row, "index"] # get split value
# #         
# #         
# #         # decide direction based on split
# #         go.left <- !is.na(cur.val) && as.numeric(cur.val) <= split.val
# #         
# #         # traverse to child node of stop if not available 
# #         if (go.left){
# #           left.child <- 2*node
# #           if (left.child %in% as.numeric(rownames(fit$frame))){
# #             node <- left.child 
# #           } else break # left child doesn't exist
# #         } else {
# #           right.child <- 2*node + 1
# #           if (right.child %in% as.numeric(rownames(fit$frame))){
# #             node <- right.child
# #           } else break # right child doesn't exist
# #         }
# #       }
# #       
# #       where[i] <- node # assign final node
# #     }
# #   }
# #   
# #   # 3. construct result matrix w/ label counts from terminal nodes
# #   frame <- fit$frame
# #   result.matrix <- matrix(0, nrow=length(where), ncol=ncol(Y.df))
# #   colnames(result.matrix) <- colnames(Y.df)
# #   
# #   for (i in 1:length(where)) { # for each observation 
# #     node.idx <- where[i]
# #     
# #     if ("label" %in% names(frame)){ # try to get label from frame
# #       
# #       node.label <- frame$label[node.idx]
# #       
# #       if (is.list(node.label)) node.label <- unlist(node.label) # unlist
# #       
# #       if (length(node.label) == ncol(Y.df)) {
# #         result.matrix[i, ] <- node.label # assign label counts
# #       }
# #       
# #     } else if (!is.null(frame$yval2) && ncol(frame$yval2) == ncol(Y.df)) { # try to get yval2 for counts
# #       # if no label, try yval2
# #       result.matrix[i, ] <- frame$yval2[node.idx, ]
# #     }
# #   }
# #   
# #   # 4. convert counts to probabilities or predicted class
# #   if (type == "matrix") {
# #     return(result.matrix) # return raw counts
# #   }
# #   
# #   # use dirichlet smoothing to convert counts to probabilities
# #   probs <- t(apply(result.matrix, 1, function(counts){
# #     (counts + alphavec) / sum(counts + alphavec)
# #   }))
# #   colnames(probs) <- colnames(Y.df)
# #   
# #   # 5. return based on type 
# #   if (type == "vector" || type == "prob") {
# #     return(probs) # return smoothed probabilities
# #   } else if (type == "class") {
# #     result <- apply(probs, 1, which.max) # return index of highest prob birth weight category
# #     return(result)
# #   }
# #   
# #   
# # }

# dm.method <- list(init=myinit, eval=myeval, split=mysplit, pred=mypred, method="dm")
# # pred <- mypred(dm.tree, newdata = data.frame(X.matrix), type = "prob")


# #-----------------------------------------------------------
# # D. Fit rpart Tree Using Our Custom DM Method
# #-----------------------------------------------------------
# cat("\n--- Fitting the DM-based rpart tree ---\n")

# dm.control <- rpart.control(minsplit=2, cp=0, maxdepth=8, xval=0, usesurrogate = 0)
# dm.tree <- rpart(
#   formula = Y.df~.,
#   data = data.frame(X.matrix),
#   method = dm.method,
#   control = dm.control
# )


# # preds <- mypred(dm.tree, newdata = data.frame(X.matrix), type = "prob")
# # head(preds, 10)

# # write.csv(preds, file = sprintf("/Users/adamkurth/Downloads/preds_%d.csv", year), row.names = FALSE)
# # dim(mypreds) # 128 x 11
# #   - hits "where <- fit$where" and outputs pred. prob. vector for 128 rows
# #   - retrieves terminal node indices for each row
# #   - "result.matrix" builds one row per terminal node (i.e. 128 total)
# #   - "frame$label[node.idx]" contains label counts at leaf
# #   - each i-th row in result.matrix represents raw counts of label occurrences at term. node that i-th obs falls in.
# #    -- i.e. node 12 has 3 examples w/ labels [1,0,0,...][1,1,0,...][0,1,0,...] then counts might look like [2,2,0,...]
# #   - "probs" smooths the counts into prob. using alphavec, sum to 1
# #   - sum of each row = 1, each column = category


# #-----------------------------------------------------------
# # E. Save and Inspect Results
# #-----------------------------------------------------------
# save(dm.tree, file=sprintf("%s/dm_tree_rebin_%d.RData", data.pwd, year))
# print(dm.tree)
# printcp(dm.tree)

# plot(dm.tree, main="DM Tree",uniform=TRUE)
# text(dm.tree,use.n=TRUE,cex=1,font=3)

# path.rpart(dm.tree, node = 2)     # left child  => mrace15 = 0
# path.rpart(dm.tree, node = 3)     # right child => mrace15 = 1

# #-----------------------------------------------------------
# # F. Multinomial Parametric Bootstrap Sampling
# #-----------------------------------------------------------

# B <- 10000  # Number of bootstrap samples
# # B <- 100
# n.rows <- nrow(counts.df)

# # category definitions
# get.category.cols <- function(type = "lbw"){

#     if (type == "lbw") {
#       return(1:10)  # Indices of LBW categories 
#     } else if (type == "nbw") {
#       return(11:ncol(Y.df))  # Adjust based on your data structure
#     } else if (type == "all") {
#       return(1:ncol(Y.df))  # All categories
#     } else {
#       stop("Invalid category type. Use 'lbw', 'nbw', or 'all'")
#     }

# }

# # define category cols 
# lbw.cols <- get.category.cols(type = "lbw")
# nbw.cols <- get.category.cols(type = "nbw")

# # create one structure to store all probability vectors 
# # each element in list is bootstrap sample, of matrix n.rows x ncol(Y.df)
# Yhat <- vector("list", length = B)


# # Dirichlet-smoothed cell probabilities
# multinomial.probs <- matrix(0, nrow = n.rows, ncol = ncol(Y.df))

# for(i in 1:n.rows) {
  
#   # Get counts for this cell/row
#   cell.counts <- as.numeric(Y.df[i,])
  
#   # Calculate probability with Dirichlet smoothing
#   cell.probs <- (cell.counts + alphavec) / (sum(cell.counts) + sum(alphavec))
  
#   multinomial.probs[i, ] <- cell.probs
# }


# # row-probability vector for total-count bootstrap
# row.probs <- rowMeans(counts.df) # positive weights
# row.probs <- row.probs / sum(row.probs)  # Normalize to sum to 1

# barplot(row.probs, main="Row Probabilities for Bootstrap Sampling", 
#         xlab="Row Index", ylab="Probability")

# # initialize list for tree structure information
# tree.structure <- list()
# top.vars.used <- character(B)
# n.vars.used <- integer(B)

# cat("Multinomial parametric bootstrap     (inter-/intra- resampling)\n")
# pb <- txtProgressBar(0, B, style = 3)


# #-----------------------------------------------------------
# # G. Main Loop over B Bootstrap Samples
# #-----------------------------------------------------------
# for (b in seq_len(B)) {
#     # Note: Y.df == counts.df
#     #---------------------------------------
#     # (A) 2-stage parametric bootstrap of the 128×11 count table
#     #---------------------------------------
#     # 1st stage: generate new counts from multinom

#     n.star <- as.vector( rmultinom(1, size = sum(Y.df[]), prob = row.probs) )
    
#     # 2nd stage: generate counts for each row using multinom
#     boot.counts <- t(
#       sapply(1:n.rows, function(j) rmultinom(1, size = n.star[j], prob = multinomial.probs[j, ]))
#     )
    
#     colnames(boot.counts) <- colnames(Y.df)
#     boot.df <- as.data.frame(boot.counts)

    
#     #---------------------------------------------------------
#     # (B) fit the full DM tree on the bootstrap counts
#     #---------------------------------------------------------
#     model.data <- cbind(boot.df, X.matrix)
#     r.formula <- as.formula( paste0("cbind(", paste(colnames(boot.df), collapse = ","), ") ~ .") )
    
#     # fit the tree on the bootstrap sample
#     dm.tree.b <- rpart(
#       formula = r.formula,
#       data = model.data,
#       method  = dm.method,
#       control = dm.control
#     )
    
#     #---------------------------------------------------------
#     # (C) record tree statistics
#     #--------------------------------------------------------
#     vars.used <- unique(dm.tree.b$frame$var[dm.tree.b$frame$var != "<leaf>"])
#     n.vars.used[b] <- length(vars.used)
#     top.vars.used[b] <- if (length(vars.used)) vars.used[1] else "none" # store first split variable
    
#     tree.structure[[b]] <-list(
#       variables = vars.used,
#       frame = dm.tree.b$frame,
#       splits = dm.tree.b$splits
#     )


#     #---------------------------------------------------------
#     # (D) Get predicted probability vectors
#     #--------------------------------------------------------
#     # custom predict function to get probability vectors
#     all.preds <- mypred(dm.tree.b, newdata = data.frame(X.matrix), type = "prob")
    
#     # Store the complete probability matrix for this bootstrap sample
#     Yhat[[b]] <- all.preds
  
    
#     setTxtProgressBar(pb, b)
# }
# close(pb)
# cat("\n--- Bootstrap sampling completed ---\n")


# #-----------------------------------------------------------
# # Visualize the bootstrap results
# #-----------------------------------------------------------
# is.high.risk.indices <- which(
#   X.matrix$mrace15 == 1 &
#     X.matrix$dmar == 0 & 
#     X.matrix$cig_0 == 1 &
#     X.matrix$sex == 0 &  # female 
#     X.matrix$mager == 0 &  # young
#     X.matrix$precare5 == 0 &  # not adequate prenatal
#     X.matrix$meduc == 0  # not high school grad
# )

# is.low.risk.indices <- which(
#   X.matrix$mrace15 == 0 &
#     X.matrix$dmar == 1 & 
#     X.matrix$cig_0 == 0 &
#     X.matrix$sex == 1 &  # male 
#     X.matrix$mager == 1 &  # older
#     X.matrix$precare5 == 0 &  # adequate prenatal
#     X.matrix$meduc == 1  # high school grad
# )


# direct.bootstrap.calculation <- function(indicies, yhat.list = Yhat) {
#   B <- length(yhat.list)
#   K <- ncol(yhat.list[[1]])
  
#   # for each b and category k, calculate the mean probability across B
#   prob.matrix <- matrix(0, nrow = B, ncol = K)
  
#   for (b in 1:B) {
#     if (length(indicies) > 0) {
#       # Get probabilities for the selected indices in this bootstrap sample
#       prob.matrix[b, ] <- colMeans(yhat.list[[b]][indicies, , drop = FALSE])
#     }
#   }
  
#   # mean probabilities by category
#   mean.probs <- colMeans(prob.matrix)
  
#   # 95% percentile bootstrap conf int
#   ci.lwr <- apply(prob.matrix, 2, function(x) quantile(x, 0.025))
#   ci.upr <- apply(prob.matrix, 2, function(x) quantile(x, 0.975))

#   return(list(
#     means = mean.probs,
#     lwr = ci.lwr,
#     upr = ci.upr
#   ))
# }


# high.risk.result <- direct.bootstrap.calculation(indicies = is.high.risk.indices, yhat.list = Yhat)
# high.risk.probs <- high.risk.result$means
# high.risk.lwr <- high.risk.result$lwr
# high.risk.upr <- high.risk.result$upr

# low.risk.result <- direct.bootstrap.calculation(indicies = is.low.risk.indices, yhat.list = Yhat)
# low.risk.probs <- low.risk.result$means
# low.risk.lwr <- low.risk.result$lwr
# low.risk.upr <- low.risk.result$upr

# K <- length(high.risk.probs)



# # 
# if (type == 1) {
#   result.tab <- read.csv(sprintf("%s/result_table_1.csv", results.pwd))
#   result.tab[-1]
# } else {
#   result.tab <- read.csv(sprintf("%s/result_table_2.csv", results.pwd))
#   result.tab[-1]
# }

# # 
# result.tab.1 <- read.csv("~/Downloads/result_table_1.csv")[-1]
# result.tab.2 <- read.csv("~/Downloads/result_table_2.csv")[-1]
# diff.pi.hat <- result.tab.1$high.risk.prob[1:10] - result.tab.1$low.risk.prob[1:10]



# cat.labels <- paste0("C", 1:K)
# pdf(file = sprintf("%s/high_low_risk_pred_%d.pdf", boot.pwd, year), width = 8.5, height = 8)
# par(mfrow = c(2, 1), 
#     mar = c(3, 4.5, 3, 2),    # More space around individual plots
#     oma = c(2, 0, 4, 0))      # Larger outer margin at top for title and subtitle

# # top
# high.bars <- barplot(high.risk.probs, 
#                      main = "",
#                      xlab = "", 
#                      ylab = "Probability",
#                      ylim = c(0, max(high.risk.upr) * 1.1),
#                      names.arg = cat.labels,  
#                      col = "red",
#                      cex.names = 1.0,
#                      cex.axis = 1.0,
#                      border = "darkred")  # Add border for better definition
# title("High Risk Subgroup", line = 1, cex.main = 1.2)

# arrows(high.bars, high.risk.lwr, 
#        high.bars, high.risk.upr, 
#        angle = 90, code = 3, length = 0.05, 
#        col = "black", lwd = 1.5)


# # bottom
# low.bars <- barplot(low.risk.probs,
#                     main = "",  # Remove individual title
#                     xlab = "", 
#                     ylab = "Probability",
#                     ylim = c(0, max(high.risk.upr) * 1.1),
#                     names.arg = cat.labels,
#                     col = "blue",
#                     cex.names = 1.0,
#                     cex.axis = 1.0,
#                     border = "darkblue")  # Add border

# title("Low Risk Subgroup", line = 1, cex.main = 1.2)

# arrows(low.bars, low.risk.lwr, 
#        low.bars, low.risk.upr, 
#        angle = 90, code = 3, length = 0.05, 
#        col = "black", lwd = 1.5)

# mtext(sprintf("Birth Weight Category Probabilities (%d)", year), 
#       side = 3, line = 2, outer = TRUE, cex = 1.4, font = 2)

# # Add subtitle
# mtext("Predicted probabilities with 95% confidence intervals", 
#       side = 3, line = 0.5, outer = TRUE, cex = 1.1)

# # Add note about categories if needed
# mtext("Categories represent different birth weight ranges", 
#       side = 1, line = 0, outer = TRUE, cex = 0.9, col = "darkgray")

# dev.off()
# par(mfrow = c(1, 1))




# #-----------------------------------------------------------
# # Quick check for pred. prob. vectors 
# #-----------------------------------------------------------
# analyze.bootstrap.samples <- function(Y.hat) {
#   B <- length(Y.hat)
  
#   lbw.cols <- get.category.cols(type = "lbw")
#   nbw.cols <- get.category.cols(type = "nbw")
  
#   # --- Category-level variation ---
#   category.stats <- lapply(1:ncol(Y.hat[[1]]), function(i) {
#     probs <- sapply(Y.hat, function(b) b[, i])
#     list(
#       mean = rowMeans(probs),
#       sd = apply(probs, 1, sd),
#       cv = apply(probs, 1, sd) / rowMeans(probs)
#     )
#   })
  
#   # --- LBW probability analysis ---
#   lbw.probs <- sapply(Y.hat, function(b) rowSums(b[, lbw.cols]))
#   lbw.means <- rowMeans(lbw.probs)
#   lbw.sds <- apply(lbw.probs, 1, sd)
#   lbw.cvs <- lbw.sds / lbw.means
  
#   # # --- NBW probability analysis ---
#   # nbw.probs <- rows)
#   # nbw.means <- rowMeans(nbw.probs)
#   # nbw.sds <- apply(nbw.probs, 1, sd)
#   # nbw.cvs <- nbw.sds / nbw.means
  
#   cat("LBW mean prob:", mean(lbw.means),
#       "| SD:", sd(lbw.means),
#       "| Mean CV:", mean(lbw.cvs), "\n")
  
#   # cat("NBW mean prob:", mean(nbw.means),
#   #     "| SD:", sd(nbw.means),
#   #     "| Mean CV:", mean(nbw.cvs), "\n")
  
#   # --- Category proportion consistency ---
#   category.props <- sapply(Y.hat, function(b) colMeans(b))
#   category.means <- rowMeans(category.props)
#   category.sds <- apply(category.props, 1, sd)
#   category.cvs <- category.sds / category.means
  
#   proportion.summary <- data.frame(
#     # category = colnames(Y.hat[[1]]),
#     mean.prob = category.means,
#     sd = category.sds,
#     cv = category.cvs
#   )
#   print(proportion.summary)
  
#   # --- Terminal node pattern consistency ---
#   pattern.counts <- sapply(Y.hat, function(b) {
#     length(unique(apply(round(b, 6), 1, paste, collapse = ",")))
#   })
#   cat("Unique terminal node patterns per sample:\n")
#   print(summary(pattern.counts))
  
#   list(
#     category.stats = category.stats,
#     lbw.stats = list(means = lbw.means, sds = lbw.sds, cvs = lbw.cvs),
#     category.proportions = proportion.summary,
#     pattern.counts = pattern.counts
#   )
# }

# compare.category.predictions <- function(Y.hat, category.idx) {
#   B <- length(Y.hat)
#   category.preds <- sapply(Y.hat, function(b) b[, category.idx])
  
#   means <- rowMeans(category.preds)
#   sds <- apply(category.preds, 1, sd)
#   cvs <- sds / means
  
#   par(mfrow = c(3, 1))
  
#   hist(means,
#        main = paste("Distribution of Mean Probabilities -", colnames(Y.hat[[1]])[category.idx]),
#        xlab = "Mean Probability", col = "lightblue", border = "white")
  
#   hist(sds,
#        main = "Distribution of Standard Deviations",
#        xlab = "Standard Deviation", col = "lightblue", border = "white")
  
#   hist(cvs,
#        main = "Distribution of Coefficients of Variation",
#        xlab = "CV", col = "lightblue", border = "white")
#   # 
#   # plot(means, sds,
#   #      main = "Standard Deviation vs. Mean",
#   #      xlab = "Mean Probability", ylab = "Standard Deviation",
#   #      pch = 19, col = "blue")
#   # 
#   par(mfrow = c(1, 1))
  
#   list(
#     means = means,
#     sds = sds,
#     cvs = cvs,
#     overall.mean = mean(means),
#     overall.sd = sd(means),
#     mean.cv = mean(cvs)
#   )
# }

# validate.full.probability.vectors <- function(Yhat) {
#   validations <- lapply(seq_along(Yhat), function(i) {
#     probs <- Yhat[[i]]
    
#     row.sums <- rowSums(probs)
#     all.rows.sum.to.1 <- all(abs(row.sums - 1) < .Machine$double.eps^0.5)
#     all.values.in.01 <- all(probs >= 0 & probs <= 1)
    
#     if (!all.rows.sum.to.1) {
#       cat(sprintf("Warning (bootstrap %d): Not all rows sum to 1. Min: %.5f | Max: %.5f\n",
#                   i, min(row.sums), max(row.sums)))
#     }
#     if (!all.values.in.01) {
#       cat(sprintf("Warning (bootstrap %d): Some values outside [0,1] range.\n", i))
#     }
    
#     faulty.rows <- which(abs(row.sums - 1) > .Machine$double.eps^0.5 |
#                            probs < 0 | probs > 1, arr.ind = TRUE)
    
#     list(
#       bootstrap.index = i,
#       all.rows.sum.to.1 = all.rows.sum.to.1,
#       all.values.in.01 = all.values.in.01,
#       faulty.rows = faulty.rows
#     )
#   })
  
#   return(validations)
# }

# bootstrap.analysis <- analyze.bootstrap.samples(Yhat)
# print(bootstrap.analysis$category.proportions)

# cat.1.analysis <- compare.category.predictions(Yhat, 1)
# cat.6.analysis <- compare.category.predictions(Yhat, 6)
# cat.10.analysis <- compare.category.predictions(Yhat, 10)
# # cat.11.analysis <- compare.category.predictions(Yhat, 11)

# # category.comparison <- data.frame(
# #   category = c(colnames(Yhat[[1]])[1], colnames(Yhat[[1]])[6], colnames(Yhat[[1]])[11]),
# #   mean.prob = c(cat.1.analysis$overall.mean, cat.6.analysis$overall.mean, cat.11.analysis$overall.mean),
# #   sd.prob = c(cat.1.analysis$overall.sd, cat.6.analysis$overall.sd, cat.11.analysis$overall.sd),
# #   mean.cv = c(cat.1.analysis$mean.cv, cat.6.analysis$mean.cv, cat.11.analysis$mean.cv)
# # )
# category.comparison <- data.frame(
#   category = c(colnames(Yhat[[100]])[1], colnames(Yhat[[100]])[6], colnames(Yhat[[100]])[10]),
#   mean.prob = c(cat.1.analysis$overall.mean, cat.6.analysis$overall.mean, cat.10.analysis$overall.mean),
#   sd.prob = c(cat.1.analysis$overall.sd, cat.6.analysis$overall.sd, cat.10.analysis$overall.sd),
#   mean.cv = c(cat.1.analysis$mean.cv, cat.6.analysis$mean.cv, cat.10.analysis$mean.cv)
# )

# print(category.comparison)

# hist(bootstrap.analysis$category.proportions$mean.prob,
#      main = "Distribution of Mean Predicted Probabilities",
#      xlab = "Mean Probability", col = "lightblue")

# vector.validations <- validate.full.probability.vectors(Yhat)
# print(vector.validations[1:5])





# compare.probability.matrices <- function(Yhat, alphavec, return.full.diff = FALSE) {
#   comparisons <- lapply(seq_along(Yhat), function(i) {
#     probs <- Yhat[[i]]
    
#     if (!all(dim(probs) == dim(alphavec))) {
#       stop(sprintf("Dimension mismatch in bootstrap %d: Yhat has dim %s but alphavec has dim %s",
#                    i,
#                    paste(dim(probs), collapse = "x"),
#                    paste(dim(alphavec), collapse = "x")))
#     }
    
#     diff.matrix <- probs - alphavec
#     max.abs.diff <- max(abs(diff.matrix))
#     mean.abs.diff <- mean(abs(diff.matrix))
    
#     list(
#       bootstrap.index = i,
#       max.abs.diff = max.abs.diff,
#       mean.abs.diff = mean.abs.diff,
#       diff.matrix = if (return.full.diff) diff.matrix else NULL
#     )
#   })
  
#   return(comparisons)
# }


# comparison.results <- compare.probability.matrices(Yhat, alphavec)
# str(comparison.results)




# #-----------------------------------------------------------
# # Depth distributions visualization
# #-----------------------------------------------------------
# library(reshape2)
# library(dplyr)
# library(ggplot2)

# # get all possible predictor variables
# all.vars <- colnames(X.matrix)
# n.vars <- length(all.vars)
# co.occur.mat <- matrix(0, nrow = n.vars, ncol = n.vars, dimnames = list(all.vars, all.vars))
# depth.sums <- setNames(numeric(n.vars), all.vars)
# depth.counts <- setNames(numeric(n.vars), all.vars)
# depth.values <- vector("list", length = n.vars)
# names(depth.values) <- all.vars

# for (b in seq_len(B)){ 
#   tree <- tree.structure[[b]]
#   frame <- tree$frame
#   vars.used <- tree$variables
  
#   # update co-coccurence matrix
#   if (length(vars.used) >= 1) {
#       idx <- which(all.vars %in% vars.used)
#       co.occur.mat[idx, idx] <- co.occur.mat[idx, idx] + 1
#   }
  
#   # compute depth per node
#   node.ids <- as.numeric(rownames(frame))
#   node.depths <- floor(log2(node.ids))
#   vars <- frame$var
#   internal.nodes <- which(vars != "<leaf>")
  
#     for (i in internal.nodes) {
#       var <- vars[i]
#       depth <- node.depths[i]
#         if (var %in% names(depth.sums)) {
#             depth.sums[var] <- depth.sums[var] + depth
#             depth.counts[var] <- depth.counts[var] + 1
#             depth.values[[var]] <- c(depth.values[[var]], depth)
#         }
#     }
# }    

# co.occur.prop <- co.occur.mat / B
# avg.depths <- depth.sums / depth.counts
# depth.df <- data.frame(
#   Variable = names(avg.depths),
#   AvgDepth = avg.depths
# )

# depth.long <- bind_rows(
#   lapply(names(depth.values), function(var) {
#     depths <- depth.values[[var]]
#     if (length(depths) > 0) {
#       data.frame(Variable = var, Depth = depths)
#     }
#   })
# )
# depth.long$Variable <- factor(depth.long$Variable, levels = depth.df$Variable[order(depth.df$AvgDepth)])

# means <- depth.long %>%
#   group_by(Variable) %>%
#   summarise(mean_depth = mean(Depth))

# ggplot(depth.long, aes(x = Depth, fill = Variable)) +
#   geom_histogram(binwidth = 1, alpha = 0.6) +
#   geom_vline(data = means, aes(xintercept = mean_depth),
#              color = "black", linetype = "dashed", size = 0.5) +
#   facet_wrap(~ Variable, scales = "free_y") +
#   labs(title = "Distribution of Depths for Each Predictor Variable",
#        x = "Depth",
#        y = "Number of Splits") +
#   theme_minimal()

# ggsave(sprintf("%s/depth_distributions_%d.png", boot.pwd, type), width = 12, height = 8, dpi = 400)

# #-----------------------------------------------------------
# # H. Decision Tree Depth Analysis for cig_0 Predictor
# #-----------------------------------------------------------
# cat("\n--- Analyzing cig_0 predictor across different tree depths ---\n")
# depths <- c(2, 3, 4, 5)
# trees <- list()
# tree.sum <- list()
# var.mats <- list()

# cat("\n--- Comparing Decision Trees at Different Depths ---\n")

# for (depth in depths) {
    
#     cat(sprintf("\nFitting tree with max depth = %d...\n", depth))
    
#     # control parameters with specific depths
#     depth.control <- rpart.control(
#       minsplit = 2,     # min obs in node for split
#       cp = 0,           # complexity param
#       maxdepth = depth, # max depth 
#       xval = 0,         # no cross-validation
#       usesurrogate = 0  # no surrogate splits
#     )
    
#     depth.tree <- rpart(
#       formula = r.formula,
#       data = model.data,
#       method = dm.method,
#       control = depth.control
#     )
  
#     # store tree 
#     trees[[as.character(depth)]] <- depth.tree
    
#     # extract info about variables used
#     vars.used <- unique(depth.tree$frame$var[depth.tree$frame$var != "<leaf>"])
    
#     tree.sum[[as.character(depth)]] <- list(
#       variables = vars.used,
#       n.variables = length(vars.used),
#       n.terminal.nodes = sum(depth.tree$frame$var == "<leaf>"),
#       first.split = ifelse(length(vars.used) > 0, 
#                           as.character(depth.tree$frame$var[1]), 
#                           "none")
#     )
    
#     # check if "cig_0" is used in this tree
#     cig.used <- "cig_0" %in% vars.used
    
#     if (cig.used){
        
#         # find where cig0 used
#         cig.nodes <- which(depth.tree$frame$var == "cig_0")
#         node.depths <- floor(log2(cig.nodes))
#         cat("cig_0 appears at tree level(s):", paste(unique(node.depths), collapse = ", "), "\n")
        
#     } else {

#         cat("cig_0 not used in this tree.\n")
#     }
    
#       # create presence/absence matrix for variables used
#       # visualize which variables used at different depths
#       var.names <- colnames(X.matrix)
#       var.mat <- matrix(0, nrow = length(var.names), ncol = depth)
#       rownames(var.mat) <- var.names  
#       colnames(var.mat) <- paste0("level_", 1:depth)
      
#       # for each node in tree 
#       for (i in 1:nrow(depth.tree$frame)) {
        
#           var.name <- depth.tree$frame$var[i]
          
#           if (var.name != "<leaf>"){
#               # calculate the level of this node (approximately)
#               node.level <- min(floor(log2(i)) + 1, depth)
#               var.mat[var.name, node.level] <- 1
#           }
#       }
      
#       var.mats[[as.character(depth)]] <- var.mat
      
#       # Print a summary of tree structure
#       cat("\nTree Summary (Depth =", depth, "):\n")
#       cat("Number of terminal nodes:", sum(depth.tree$frame$var == "<leaf>"), "\n")
#       cat("Number of variables used:", length(vars.used), "\n")
#       cat("First split on variable:", ifelse(length(vars.used) > 0, as.character(depth.tree$frame$var[1]), "none"), "\n")
  
#        # print top 3 splits
#       if (length(vars.used) >= 3) {
#           top.splits <- as.character(depth.tree$frame$var[1:min(3, length(vars.used) + 1)])
#           top.splits <- top.splits[top.splits != "<leaf>"]
#           cat("Top 3 splits: ", paste(top.splits, collapse = ", "), "\n")
#       }
      
#       # Predict LBW probabilities with this tree
#       preds <- predict(depth.tree, newdata = data.frame(X.matrix), type = "matrix")
#       lbw.prob <- rowSums(preds[, lbw.cols]) / rowSums(preds)
      
#       # Store prediction statistics directly in tree summary (not in a nested pred.stats list)
#       tree.sum[[as.character(depth)]]$mean.lbw.prob = mean(lbw.prob)
#       tree.sum[[as.character(depth)]]$median.lbw.prob = median(lbw.prob)
#       tree.sum[[as.character(depth)]]$min.lbw.prob = min(lbw.prob)
#       tree.sum[[as.character(depth)]]$max.lbw.prob = max(lbw.prob)
#       tree.sum[[as.character(depth)]]$sd.lbw.prob = sd(lbw.prob)
#       # These are already stored earlier, but we'll store them again for consistency
#       tree.sum[[as.character(depth)]]$n.terminal.nodes = sum(depth.tree$frame$var == "<leaf>")
#       tree.sum[[as.character(depth)]]$n.variables = length(vars.used)
      
#       cat("LBW Probability - Mean:", round(mean(lbw.prob), 4), 
#           "Min:", round(min(lbw.prob), 4), 
#           "Max:", round(max(lbw.prob), 4), "\n")
# }


# #-----------------------------------------------------------
# # J) Comparative Analysis of Trees
# #-----------------------------------------------------------
# cat("\n--- Comparative Analysis of Trees at Different Depths ---\n")

# # Compare variable usage across depths
# var.across.depths <- data.frame(
#   variable = colnames(X.matrix)
# )

# for(depth in depths) {
  
#   depth.str <- as.character(depth)
  
#   var.across.depths[[paste0("depth_", depth)]] <- 
#     ifelse(var.across.depths$variable %in% tree.sum[[depth.str]]$variables, 1, 0)
  
# }

# # Calculate total usage across all depths
# var.across.depths$total.usage <- rowSums(var.across.depths[, paste0("depth_", depths)])

# # Sort by total usage
# var.across.depths <- var.across.depths[order(-var.across.depths$total.usage), ]

# cat("\nVariable Usage Across Different Tree Depths:\n")
# print(var.across.depths)

# cat("\nFocus on cig_0 Predictor:\n")
# cig.usage <- var.across.depths[var.across.depths$variable == "cig_0", ]
# print(cig.usage)

# pred.stats.compare <- data.frame(
#   depth = depths,
#   mean.lbw.prob = sapply(as.character(depths), function(d) tree.sum[[d]]$mean.lbw.prob),
#   min.lbw.prob = sapply(as.character(depths), function(d) tree.sum[[d]]$min.lbw.prob),
#   max.lbw.prob = sapply(as.character(depths), function(d) tree.sum[[d]]$max.lbw.prob),
#   sd.lbw.prob = sapply(as.character(depths), function(d) tree.sum[[d]]$sd.lbw.prob),
#   n.terminal.nodes = sapply(as.character(depths), function(d) tree.sum[[d]]$n.terminal.nodes),
#   n.variables = sapply(as.character(depths), function(d) tree.sum[[d]]$n.variables)
# )


# cat("\nPrediction and Tree Structure Comparison:\n")
# print(pred.stats.compare)

# # check consistency
# first.splits <- sapply(as.character(depths), function(d) tree.sum[[d]]$first.split)
# cat("\nFirst Split Variable at Each Depth:", paste(first.splits, collapse = ", "), "\n")

# if (length(unique(first.splits)) == 1){
#     cat("First split is consistent across all depths\n")
# } else {
#     cat("First split varies across different depths\n")
# }

# save(trees, tree.sum, var.mats, var.across.depths, pred.stats.compare,
#      file = sprintf("%s/tree_depth_comparison_%d.RData", results.pwd, year))

# #-----------------------------------------------------------
# # K) Analyze Bootstrap Tree Structure Results
# #-----------------------------------------------------------
# # Summarize variable usage across bootstrap samples
# var.names <- colnames(X.matrix) 
# var.usage <- matrix(0, nrow = B, ncol=length(var.names))
# colnames(var.usage) <- var.names

# for ( b in 1:B ){
  
#   if(length(tree.structure[[b]]$variables) > 0 ){
    
#     for (var in tree.structure[[b]]$variables){
      
#       if(var %in% var.names) {
#         var.usage[b,var] <- 1
#       }
#     }
#   }
# }

# # calclate frequency table of variables used
# var.freq <- colSums(var.usage) / B
# var.freq.df <- data.frame(
#   variable = var.names,
#   frequency = var.freq
# )
# var.freq.df <- var.freq.df[order(-var.freq.df$frequency),]


# # top split variable 
# top.var.freq <- table(top.vars.used) / B 
# top.var.df <- data.frame(
#   variable = names(top.var.freq),
#   frequency = as.numeric(top.var.freq)
# )
# top.var.freq <- top.var.df[order(-top.var.df$frequency),]


# # analyze number of variables used 
# n.vars.summary <- data.frame(
#   mean = mean(n.vars.used),
#   median = median(n.vars.used),
#   min = min(n.vars.used),
#   max = max(n.vars.used), 
#   sd = sd(n.vars.used)
# )


# bootstrap.tree.results <- list(
#   var.usage = var.usage,
#   var.freq.df = var.freq.df,
#   top.var.df = top.var.df,
#   n.vars.summary = n.vars.summary
# )

# #-----------------------------------------------------------
# cat("\n--- Bootstrap Tree Structure Analysis ---\n")
# cat("\nFrequency of Variables Used in Trees:\n")
# print(var.freq.df)

# cat("\nFrequency of Top Split Variables:\n")
# print(top.var.freq)

# cat("\nSummary of Number of Variables Used:\n")
# print(n.vars.summary)

# save(bootstrap.tree.results,file = sprintf("%s/bootstrap_tree_results_%d.RData", results.pwd, year))