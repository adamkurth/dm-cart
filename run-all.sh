#!/usr/bin/env bash
# =============================================================================
#  run-all.sh -- reproduce the DM-CART analysis end to end.
#
#    ./run-all.sh                     full pipeline, B = 10000, skips finished stages
#    ./run-all.sh --b 500             smoke run at a smaller ensemble size
#    ./run-all.sh --figures-only      redraw article outputs from saved results
#    ./run-all.sh --force             redo every stage, ignoring saved outputs
#    ./run-all.sh --cohort 2020-2021  a different cohort
#
#  Stages, in order (each consumes the previous one's output):
#    1 quantile.R  prior-year file  -> decile cutpoints + Dirichlet prior
#    2 data.R      target-year file -> class design matrix + counts
#    3 dm-cart.R   counts           -> fitted trees + bootstrap ensemble
#    4 visual.R    saved results    -> figures and tables in results/<cohort>/article/
#    5 verify-numbers.R             -> checks every number quoted in main.tex
#
#  Only stage 3 is expensive. Stages 1-2 read the raw natality CSVs (~2 GB each);
#  stage 4 reads nothing but saved .RData and is safe to re-run at any time.
# =============================================================================
set -euo pipefail

COHORT="2023-2024"; B=10000; FORCE=0; FIGURES_ONLY=0
while [ $# -gt 0 ]; do
  case "$1" in
    --b)            B="$2"; shift 2 ;;
    --cohort)       COHORT="$2"; shift 2 ;;
    --force)        FORCE=1; shift ;;
    --figures-only) FIGURES_ONLY=1; shift ;;
    -h|--help)      sed -n '2,20p' "$0"; exit 0 ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
PRIOR="${COHORT%%-*}"; TARGET="${COHORT##*-}"
REBIN="data/rebin/$COHORT"; ART="results/$COHORT/article"

say()  { printf '\n\033[1m== %s\033[0m\n' "$*"; }
skip() { printf '   (skip) %s\n' "$*"; }

# --- stage 1: cutpoints and prior from the PRIOR year ------------------------
say "1/5  quantile.R  --  $PRIOR cutpoints and Dirichlet prior"
if [ $FIGURES_ONLY -eq 1 ] || { [ $FORCE -eq 0 ] && [ -f "$REBIN/informed_prior_$PRIOR.RData" ]; }; then
  skip "informed_prior_$PRIOR.RData exists"
else
  Rscript quantile.R "$COHORT"
fi

# --- stage 2: class design and counts from the TARGET year -------------------
say "2/5  data.R  --  $TARGET counts binned on $PRIOR cutpoints"
if [ $FIGURES_ONLY -eq 1 ] || { [ $FORCE -eq 0 ] && [ -f "$REBIN/birthweight_data_rebin_$TARGET.RData" ]; }; then
  skip "birthweight_data_rebin_$TARGET.RData exists"
else
  Rscript data.R "$COHORT"
fi

# --- stage 3: trees and bootstrap -------------------------------------------
say "3/5  dm-cart.R  --  fit + bootstrap (B = $B)"
BOOT="results/$COHORT/fullmodel/bootstrap_tree_results_$TARGET.RData"
EXISTING=$(Rscript -e '
  f <- commandArgs(TRUE)[1]
  if (file.exists(f)) { e <- new.env(); load(f, envir = e)
    cat(e$bootstrap.tree.results$B) } else cat(0)' "$BOOT" 2>/dev/null || echo 0)

if [ $FIGURES_ONLY -eq 1 ]; then
  skip "figures-only: keeping the saved ensemble (B = $EXISTING)"
elif [ $FORCE -eq 0 ] && [ "$EXISTING" -gt "$B" ]; then
  # Guard: a smoke run must never silently replace a publication run.
  echo "   REFUSING: saved results use B = $EXISTING, larger than the requested B = $B."
  echo "   Re-run with --force to overwrite, or --figures-only to keep them."
  exit 1
elif [ $FORCE -eq 0 ] && [ "$EXISTING" -eq "$B" ]; then
  skip "ensemble already at B = $B"
else
  DMCART_B="$B" Rscript dm-cart.R "$COHORT"
fi

# --- stage 4: article figures and tables ------------------------------------
say "4/5  visual.R  --  figures and tables -> $ART"
Rscript visual.R "$COHORT"
Rscript visual.R "$COHORT" lbw

FINAL=$(Rscript -e '
  f <- commandArgs(TRUE)[1]; e <- new.env(); load(f, envir = e)
  cat(e$bootstrap.tree.results$B)' "$BOOT")
say "done  --  outputs built from B = $FINAL"
ls -1 "$ART"

# --- stage 5: check the manuscript's numbers against what was just produced --
# A re-fit moves interval endpoints in the third decimal and selection
# percentages by a few tenths of a point. Without this the drift is silent.
say "5/5  verify-numbers.R  --  manuscript claims vs saved output"
Rscript verify-numbers.R || {
  echo "   Numbers in tex/article/main.tex no longer match the saved results."
  echo "   Update the manuscript, then update the expected values in verify-numbers.R."
  exit 1
}
