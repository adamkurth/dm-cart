library(rprojroot)

# ---- Path resolver ----
# year        : the natality data year being modeled (2021, 2023, 2024, ...)
# prior.year  : year used to build the informed Dirichlet prior / quantile
#               cutpoints, avoiding double-dipping. Defaults to year - 1,
#               matching the thesis methodology (prior derived from the
#               *previous* year's data). Unused by quantile.R itself --
#               see note in the quantile.R preamble example below.
# with.2.5kg  : TRUE  -> "full model" (11 categories, includes NBW)
#               FALSE -> "LBW-only model" (10 categories)
get_paths <- function(year, prior.year = NULL, with.2.5kg = TRUE) {
    
    # Fallback to global variable if with.2.5kg is omitted, else default to TRUE
    if (is.null(with.2.5kg)) {
        with.2.5kg <- if (exists("PROCESS.WITH.2.5KG", envir = .GlobalEnv)) {
        get("PROCESS.WITH.2.5KG", envir = .GlobalEnv)
        } else {
        TRUE
        }
    }
    
    
    requested.year <- year
    # Automatically resolve cohort mapping
    if (is.null(prior.year)) {
        if (year %in% c(2020, 2021)) {
        prior.year  <- 2020
        target.year <- 2021
        cohort.dir  <- "2020-2021"
        } else if (year %in% c(2023, 2024)) {
        prior.year  <- 2023
        target.year <- 2024
        cohort.dir  <- "2023-2024"
        } else {
        stop(sprintf("Unsupported year %d. Supported years are 2020, 2021 or 2023, 2024.", year))
        }
    } else {
        cohort.dir <- sprintf("%d-%d", prior.year, year)
        valid.cohorts <- c("2020-2021", "2023-2024")
        if (!cohort.dir %in% valid.cohorts) {
        stop(sprintf("Invalid cohort '%s'. Supported cohorts are strictly: %s", 
                    cohort.dir, paste(valid.cohorts, collapse = ", ")))
        }
        target.year <- year
    }

    if (prior.year >= target.year) {
        stop(sprintf(paste0("Cohort resolution failed: prior.year (%d) must be strictly ",
                            "earlier than target.year (%d). The prior/cutpoints must come ",
                            "from a different year than the counts."),
                     prior.year, target.year))
    }

    root         <- rprojroot::find_root(rprojroot::is_git_root)
    data.root    <- file.path(root, "data")
    results.root <- file.path(root, "results")

    # Rebin path: data/rebin/<cohort>/ or data/rebin_without_2.5kg/<cohort>/
    rebin.dir  <- if (with.2.5kg) "rebin" else "rebin_without_2.5kg"
    data.rebin <- file.path(data.root, rebin.dir, cohort.dir)

    # Results path: results/<cohort>/fullmodel/ or results/<cohort>/lbwmodel/
    model.dir   <- if (with.2.5kg) "fullmodel" else "lbwmodel"
    results.dir <- file.path(results.root, cohort.dir, model.dir)
    plots.dir   <- file.path(results.dir, "plots")
    
    # `article` holds only what is meant to be published, so it is always the same path regardless of model type.
    article.dir <- file.path(results.root, cohort.dir, "article")

    list(
    root       = root,
    data.root  = data.root,
    data.rebin = data.rebin,
    results    = results.dir,
    plots      = plots.dir,
    article    = article.dir,

    # these four fields are the fix. `year` used to be the
    # requested year; it is now the resolved TARGET year (the year whose births
    # are modeled). `requested.year` preserves whatever the caller passed in.
    year           = target.year,
    target.year    = target.year,
    prior.year     = prior.year,
    requested.year = requested.year,
    cohort         = cohort.dir,

    with.2.5kg = with.2.5kg,

    # is the TARGET year
    raw.csv             = function(yr = target.year) file.path(data.root, sprintf("natality%dus-original.csv", yr)),
    quantile.cutpoints  = function(yr)         file.path(data.rebin, sprintf("quantile_cutpoints_%d.RData", yr)),
    informed.prior      = function(yr)         file.path(data.rebin, sprintf("informed_prior_%d.RData", yr)),
    rebin.data          = function(yr = target.year) file.path(data.rebin, sprintf("birthweight_data_rebin_%d.RData", yr)),
    final.data.csv      = function(yr = target.year) file.path(data.rebin, sprintf("final_data_%d.csv", yr)),
    x.matrix.csv        = function(yr = target.year) file.path(data.rebin, sprintf("X_matrix_%d.csv", yr)),
    y.matrix.csv        = function(yr = target.year) file.path(data.rebin, sprintf("Y_matrix_%d.csv", yr)),
    counts.df.csv       = function(yr = target.year) file.path(data.rebin, sprintf("counts_df_%d.csv", yr)),
    dm.tree             = function(yr = target.year) file.path(data.rebin, sprintf("dm_tree_%d.RData", yr)),
    dm.tree.rebin       = function(yr = target.year) file.path(data.rebin, sprintf("dm_tree_rebin_%d.RData", yr)),
    bootstrap.results   = function(yr = target.year) file.path(results.dir, sprintf("bootstrap_tree_results_%d.RData", yr)),
    tree.depth.results  = function(yr = target.year) file.path(results.dir, sprintf("tree_depth_comparison_%d.RData", yr))
  )
}

# ---- Directory bootstrapping ----
# Creates any missing output directories (results/, results/plots/, etc.)
ensure_dirs <- function(paths) {
  dirs <- c(paths$data.rebin, paths$results, paths$plots, paths$article)
  for (d in dirs) {
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE)
      cat("Created directory:", d, "\n")
    }
  }
}

# ==========================================================================
# Tree rendering
# ==========================================================================
# These live in util.R rather than dm-cart.R so visual.R can redraw the same
# figures from the SAVED tree objects without dm-cart.R having to re-run.  
# dm-cart.R fits and saves; visual.R draws

dm.tree.coords <- function(fit, branch = 0.30) {
  tmp <- tempfile(fileext = ".pdf")
  pdf(tmp, width = 7, height = 7); on.exit({dev.off(); unlink(tmp)}, add = TRUE)
  xy <- plot(fit, uniform = TRUE, branch = branch, margin = 0.02, compress = FALSE)
  list(x = xy$x, y = xy$y)
}

# Draws into the CURRENT device 
# `highlight` is a vector of node ids; when supplied, branches on those
# paths are drawn bold and colored and everything else is dimmed.
draw.dm.tree <- function(fit, title = NULL, subtitle = NULL,
                         lbw.cols = 1:10, alphavec, cex = NULL,
                         branch = 0.30, highlight = NULL,
                         highlight.col = "#B91C1C", leaf.highlight = NULL,
                         leaf.detail = c("full", "compact", "none"),
                         branch.labels = c("both", "left"),
                         var.labels = NULL, cex.max = 1.15,
                         title.cex = 1.65, subtitle.cex = 1.05) {

  # -------`leaf.detail` and `branch.labels` exist for small multiples.
  #   leaf.detail   "full"    k / LBW / N / cls        -- the standalone figure
  #                 "compact" LBW only                 -- panelled figures
  #                 "none"    no leaf text            -- structure-only panels
  #   branch.labels "both"    levels on both branches  -- the standalone figure
  #                 "left"    left branch only; the right child is its
  #                           complement, so printing both doubles the text for
  #                           no information. At panel size that is the
  #                           difference between legible and not.
  leaf.detail   <- match.arg(leaf.detail)
  branch.labels <- match.arg(branch.labels)
  # `var.labels` is an optional named vector of SHORT display names for the
  # split variables (e.g. precare1st -> "care")
  var.name <- function(v) {
    v <- as.character(v)
    if (is.null(var.labels)) v else ifelse(v %in% names(var.labels), var.labels[v], v)
  }

  frame    <- fit$frame
  is.leaf  <- frame$var == "<leaf>"
  n.leaves <- sum(is.leaf)
  node.id  <- as.numeric(rownames(frame))
  depth.max <- max(floor(log2(node.id)))
  xy  <- dm.tree.coords(fit, branch = branch)
  lab <- labels(fit, pretty = 0, collapse = FALSE)   # [,1] left  [,2] right

  # ------- wrap long level-set labels.
  wrap.levels <- function(z, per.line = 1) {
    vapply(z, function(one) {
      parts <- strsplit(one, ",", fixed = TRUE)[[1]]
      if (length(parts) < 3) return(one)
      grp <- split(parts, ceiling(seq_along(parts) / per.line))
      paste(vapply(grp, paste, character(1), collapse = ","), collapse = ",\n")
    }, character(1), USE.NAMES = FALSE)
  }
  lab[, 1] <- wrap.levels(lab[, 1])
  lab[, 2] <- wrap.levels(lab[, 2])
  births <- rowSums(frame$yval2)
  par(mar = c(0.3, 0.3, if (is.null(title)) 1.6 else 4.0, 0.3), xpd = NA)

  # ---- how large can the type actually be? -------------------------------
  padx <- 0.055 * diff(range(xy$x))
  xlim <- range(xy$x) + c(-padx, padx)
  plot.new()
  plot.window(xlim = xlim, ylim = range(xy$y))   # provisional, for measuring

  leaf.x  <- sort(xy$x[is.leaf])
  min.gap <- if (n.leaves > 1) min(diff(leaf.x)) else diff(xlim)

  rotate  <- TRUE
  stagger <- FALSE

  pin <- par("pin")
  in.per.ux <- pin[1] / diff(xlim)              # inches per user x-unit
  gap.in    <- min.gap * in.per.ux

  if (is.null(cex)) {
    line.in <- 1.30 * par("ps") / 72
    char.in <- 0.60 * par("ps") / 72
    if (leaf.detail == "none") {
      # No leaf text: the tightest thing on the page is the node label
      max.chr <- max(nchar(var.name(frame$var[!is.leaf])), 3L)
      cex <- min(cex.max, max(0.45, gap.in / (max.chr * char.in)))
    } else if (leaf.detail == "compact") {
      # the compact label is a single short number ("0.179") horizontally
      cex <- min(cex.max, max(0.42, (2 * gap.in) / (8.5 * char.in)))
    } else {
      cex <- min(cex.max, max(0.42, gap.in / (4.3 * line.in)))
    }
  }

  # ---- vertical room for the leaf labels ---------------------------------
  y.gap    <- diff(range(xy$y)) / max(depth.max, 1)
  char.in.1 <- 0.60 * par("ps") / 72     # width of one character at cex = 1
  in.per.x  <- pin[1] / diff(xlim)
  depth.of  <- floor(log2(node.id))

  local.gap <- function(idx, step = 1) {
    # distance to the nearest same-depth neighbor `step` positions away
    vapply(idx, function(i) {
      sib <- idx[depth.of[idx] == depth.of[i]]
      xs  <- sort(xy$x[sib])
      k   <- match(xy$x[i], xs)
      cand <- c(if (k > step) xs[k] - xs[k - step] else Inf,
                if (k + step <= length(xs)) xs[k + step] - xs[k] else Inf)
      m <- min(cand)
      if (is.finite(m)) m * in.per.x else diff(xlim) * in.per.x
    }, numeric(1))
  }

  int.idx  <- which(!is.leaf)
  leaf.idx.all <- which(is.leaf)
  node.cex <- rep(cex * 1.05, nrow(frame))
  if (length(int.idx)) {
    g <- local.gap(int.idx, 1)
    w <- nchar(var.name(frame$var[int.idx])) * char.in.1
    node.cex[int.idx] <- pmin(cex * 1.05, pmax(0.50, 0.92 * g / w))
  }
  stagger.rows <- 1L
  if (leaf.detail == "compact" && length(leaf.idx.all) > 1) {
    lx0 <- sort(xy$x[leaf.idx.all])
    min.gap.in <- min(diff(lx0)) * pin[1] / diff(xlim)
    need <- (5.6 * 0.60 * par("ps") * cex / 72) / max(min.gap.in, 1e-6)
    stagger.rows <- max(1L, min(3L, as.integer(ceiling(need))))
  }

  leaf.cex <- rep(cex, nrow(frame))
  if (leaf.detail == "compact" && length(leaf.idx.all)) {
    o  <- order(xy$x[leaf.idx.all])
    xs <- xy$x[leaf.idx.all][o]
    n  <- length(xs)
    st <- stagger.rows
    g  <- vapply(seq_len(n), function(k) {
      cand <- c(if (k > st) xs[k] - xs[k - st] else Inf,
                if (k + st <= n) xs[k + st] - xs[k] else Inf)
      m <- min(cand)
      if (is.finite(m)) m * in.per.x else diff(xlim) * in.per.x
    }, numeric(1))
    leaf.cex[leaf.idx.all[o]] <- pmin(cex, pmax(0.45, 0.92 * g / (8.5 * char.in.1)))
  }


  block.in <- switch(leaf.detail,
    full    = 9.6 * 0.60 * par("ps") * cex / 72,
    compact = (2.6 + 1.35 * (stagger.rows + 2)) * 1.28 * par("ps") * cex / 72,
    none    = 0)
  base.rng <- diff(range(xy$y)); top <- 0.10 * base.rng
  pady <- block.in * (base.rng + top) / max(pin[2] - block.in, 0.15 * pin[2])
  pady <- pady + 0.22 * y.gap                      # room for the leader line

  plot.window(xlim = xlim, ylim = c(min(xy$y) - pady, max(xy$y) + top))

  dim.on <- !is.null(highlight)
  col.of <- function(i) if (!dim.on) "black" else if (node.id[i] %in% highlight) highlight.col else "#C7CBD1"

  # ---- branches: short vertical stub off the parent, then a diagonal ----
  for (i in which(!is.leaf)) {
    for (side in 1:2) {
      ci <- match(2 * node.id[i] + (side - 1), node.id)
      if (is.na(ci)) next
      on.path <- !dim.on || (node.id[i] %in% highlight && node.id[ci] %in% highlight)
      cc  <- if (!dim.on) "black" else if (on.path) highlight.col else "#D7DBE0"
      lwd <- if (!dim.on) 1 else if (on.path) 3.2 else 1
      ystub <- xy$y[i] - branch * (xy$y[i] - xy$y[ci])
      segments(xy$x[i], xy$y[i], xy$x[i], ystub, col = cc, lwd = lwd)
      segments(xy$x[i], ystub, xy$x[ci], xy$y[ci], col = cc, lwd = lwd)
    }
  }

  # ---- keep a register of what has been drawn ---------------
  placed <- list()
  box.of <- function(x, y, txt, cx, adj) {
    w <- max(strwidth(strsplit(txt, "\n")[[1]], cex = cx))
    h <- strheight(txt, cex = cx)
    c(x - adj[1] * w, x + (1 - adj[1]) * w,
      y - adj[2] * h, y + (1 - adj[2]) * h)
  }
  hits <- function(b, pad.x = 0, pad.y = 0) {
    b <- b + c(-pad.x, pad.x, -pad.y, pad.y)
    for (q in placed)
      if (b[1] < q[2] && b[2] > q[1] && b[3] < q[4] && b[4] > q[3]) return(TRUE)
    FALSE
  }

  # ---- internal nodes: variable name + the level sets on each branch ----
  for (i in which(!is.leaf)) {
    nl <- var.name(frame$var[i]); ny <- xy$y[i] + 0.155 * y.gap
    text(xy$x[i], ny, nl, cex = node.cex[i], font = 2, adj = c(0.5, 0),
         col = col.of(i))
    placed[[length(placed) + 1]] <- box.of(xy$x[i], ny, nl, node.cex[i], c(0.5, 0))
    # for a BINARY predictor the branch label is just "0"
    v.lev <- attr(fit, "xlevels")[[as.character(frame$var[i])]]
    if (!is.null(v.lev) && length(v.lev) <= 2) next

  }

  # ---- leaves: modal category, posterior P(LBW), births, classes -------
  # "N" is BIRTHS (sum of the node's category counts); "cls" is the number
  # of predictor classes routed here. rpart's frame$n is the latter
  # 
  leaf.i  <- which(is.leaf)
  ord     <- order(xy$x[leaf.i])
  lead.top <- 0.06 * y.gap
  lead.bot <- 0.22 * y.gap
  # the staggered rows must be separated by a TEXT height
  line.y <- 1.30 * par("ps") * cex / 72 * diff(par("usr")[3:4]) / par("pin")[2]
  for (j in seq_along(leaf.i)) {
    i    <- leaf.i[ord[j]]
    cnt  <- frame$yval2[i, ]
    post <- (cnt + alphavec) / sum(cnt + alphavec)
    hot  <- !is.null(leaf.highlight) && node.id[i] %in% leaf.highlight
    lc   <- if (hot) highlight.col else col.of(i)
    txt  <- switch(leaf.detail,
      full    = sprintf("k%d\nLBW=%.3f\nN=%s\ncls=%d", which.max(cnt), sum(post[lbw.cols]),
                        format(round(births[i]), big.mark = ","), frame$n[i]),
      compact = sprintf("%.3f\n(n=%s)", sum(post[lbw.cols]),
                        format(round(births[i]), big.mark = ",")),
      none    = "")

    if (leaf.detail != "none")
    lead.col <- if (hot) highlight.col else if (dim.on) "#DDE1E6" else "#9CA3AF"
    if (leaf.detail == "compact") {
      base.drop <- 0.10 * y.gap + ((j - 1) %% stagger.rows) * 1.35 * line.y
      drop <- base.drop
      for (extra in 0:4) {
        drop <- base.drop + extra * 1.35 * line.y
        if (!hits(box.of(xy$x[i], xy$y[i] - drop, txt, leaf.cex[i], c(0.5, 1)),
                  0.20 * strwidth("M", cex = leaf.cex[i]), 0.10 * line.y)) break
      }
      segments(xy$x[i], xy$y[i] - lead.top, xy$x[i], xy$y[i] - drop,
               col = lead.col, lty = 3, lwd = if (hot) 1.2 else 0.7)
      text(xy$x[i], xy$y[i] - drop, txt,
           cex = leaf.cex[i], adj = c(0.5, 1), col = lc, font = if (hot) 2 else 1)
      placed[[length(placed) + 1]] <- box.of(xy$x[i], xy$y[i] - drop, txt,
                                             leaf.cex[i], c(0.5, 1))
    } else if (leaf.detail != "none") {
      segments(xy$x[i], xy$y[i] - lead.top, xy$x[i], xy$y[i] - lead.bot,
               col = lead.col, lty = 3, lwd = if (hot) 1.2 else 0.7)
      text(xy$x[i], xy$y[i] - lead.bot, txt, srt = 90,
           cex = cex * 0.95, adj = c(1, 0.5), col = lc, font = if (hot) 2 else 1)
    }

    if (hot)
      points(xy$x[i], xy$y[i], pch = 21, cex = cex * 1.5,
             bg = highlight.col, col = highlight.col)
  }

  # ---- branch labels, placed last and de-collided ------------------------
  for (i in which(!is.leaf)) {
    v.lev <- attr(fit, "xlevels")[[as.character(frame$var[i])]]
    if (!is.null(v.lev) && length(v.lev) <= 2) next   # binary: convention covers it

    for (side in 1:2) {
      if (branch.labels == "left" && side == 2) next
      ci <- match(2 * node.id[i] + (side - 1), node.id)
      if (is.na(ci)) next
      on.path <- !dim.on || (node.id[i] %in% highlight && node.id[ci] %in% highlight)
      lcol <- if (!dim.on) "#1D4ED8" else if (on.path) highlight.col else "#C7CBD1"
      lab.i <- lab[i, side]
      adj.i <- c(if (side == 1) -0.06 else 1.06, 0.5)

      bx <- xy$x[i] + 0.36 * (xy$x[ci] - xy$x[i])
      by <- xy$y[i] + 0.36 * (xy$y[ci] - xy$y[i])

      # Step outward along the branch normal until the label is clear. The
      # normal points away from the subtree, i.e. up-and-out, which is where
      # the free space is.
      nx  <- if (side == 1) -1 else 1
      padx <- 0.30 * strwidth("M", cex = cex)
      pady <- 0.18 * line.y
      # Start at step 1, never 0: a label anchored exactly on the branch has the
      # line drawn through it
      fx <- bx; fy <- by; moved <- 1
      for (step in 1:6) {
        fx <- bx + nx * step * 0.55 * strwidth("M", cex = cex)
        fy <- by + step * 0.42 * line.y
        if (!hits(box.of(fx, fy, lab.i, cex * 0.95, adj.i), padx, pady)) { moved <- step; break }
        moved <- step
      }

      if (moved > 1)   # leader back to the branch, so the move is unambiguous
        segments(bx, by, fx, fy, col = lcol, lty = 3, lwd = 0.6)
      text(fx, fy, lab.i, cex = cex * 0.95, col = lcol,
           font = if (on.path && dim.on) 2 else 1, adj = adj.i)
      placed[[length(placed) + 1]] <- box.of(fx, fy, lab.i, cex * 0.95, adj.i)
    }
  }
  if (!is.null(title))    title(main = title, cex.main = title.cex, font.main = 2, line = 2.1)
  if (!is.null(subtitle)) mtext(subtitle, side = 3, line = 0.5, cex = subtitle.cex, col = "#4B5563")
  invisible(xy)
}

# PDF wrapper. Canvas grows with the tree, and `pointsize` grows with it too, so the type reads larger on the page than it would if the same tree were drawn on a fixed-size canvas with a fixed pointsize. 
plot.dm.tree <- function(fit, file, title, subtitle = NULL,
                         lbw.cols = 1:10, alphavec,
                         cex = NULL, width = NULL, height = NULL, ...) {

  frame    <- fit$frame
  n.leaves <- sum(frame$var == "<leaf>")
  if (nrow(frame) <= 1) {
    cat("  [skip]", basename(file), "-- tree has no splits\n")
    return(invisible(NULL))
  }
  depth.max <- max(floor(log2(as.numeric(rownames(frame)))))
  if (is.null(width))  width  <- max(11, min(34, 1.15 * n.leaves + 4))
  if (is.null(height)) height <- max(7.5, 1.75 * (depth.max + 1) + 2.6)
  base.pt <- 13.5 * width / 11

  pdf(file, width = width, height = height, pointsize = base.pt)
  on.exit(dev.off(), add = TRUE)
  draw.dm.tree(fit, title = title, subtitle = subtitle, lbw.cols = lbw.cols,
               alphavec = alphavec, cex = cex, ...)
  cat("  Saved tree plot:", file, "\n")
  invisible(fit)
}
