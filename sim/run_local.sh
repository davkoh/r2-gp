#!/usr/bin/env bash
# Local runner. Defaults to the smoke config, which takes a few minutes and
# exercises every rule: compile -> data -> fit -> aggregate -> plot -> prior.
#
#   bash sim/run_local.sh                       # smoke test
#   bash sim/run_local.sh config.yaml 8         # the real grid, 8 cores
#
# Local versus cluster is only the absence of --profile slurm_profile/.

set -euo pipefail
cd "$(dirname "$0")"

CONFIG="${1:-config_smoke.yaml}"
CORES="${2:-4}"

snakemake --configfile "${CONFIG}" --cores "${CORES}" --unlock 2>/dev/null || true

# Targets come first: --configfile and --rerun-triggers both take lists, so
# either would otherwise swallow the target names as values.
# --rerun-triggers mtime input: "mtime" alone stops Snakemake noticing that the
# SET of input files changed, so adding a model or a condition silently produces
# no work even though its result files are missing. Adding "input" restores that
# check while still ignoring cosmetic edits to the R scripts, which would
# otherwise retrigger every fit.
snakemake -p all prior_all \
  --configfile "${CONFIG}" \
  --cores "${CORES}" \
  --rerun-triggers mtime input
