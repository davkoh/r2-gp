#!/bin/bash
#SBATCH --job-name=hsgp_calib
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=5G
#SBATCH --output=calibrate-%j.out
#SBATCH --error=calibrate-%j.err

# Measure what a fit actually costs, before the full grid is launched.
#
#   sbatch run_calibrate.sh
#
# 24 fits at the expensive corner of the grid (Matern 3/2, rho = 0.3, D = 5,
# 30 and 100) with the SAME MCMC settings as the real run.
#
# Two things differ from run_snakemake.sh on purpose:
#
#   1. Grouping is OFF (--group-components fit_model=1). Each fit becomes its
#      own Slurm job, so sacct reports Elapsed and MaxRSS PER FIT. A grouped
#      job would report the total for 20 fits and tell you nothing useful.
#   2. The limits are deliberately generous (12 h, 32 GB). A fit that is killed
#      at its limit gives no measurement, only a failure.
#
# The run ends by printing the numbers to put in slurm_profile/config.yaml.

set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  cd "${SLURM_SUBMIT_DIR}"
  if [[ ! -f Snakefile && -f sim/Snakefile ]]; then cd sim; fi
else
  cd "$(dirname "$0")"
fi

SIM_ENV_PREFIX="${SIM_ENV_PREFIX:-./env}"
if [[ ! -d "${SIM_ENV_PREFIX}/conda-meta" && -d "${HOME}/.conda/envs/hsgp_r2_sim/conda-meta" ]]; then
  SIM_ENV_PREFIX="${HOME}/.conda/envs/hsgp_r2_sim"
fi

module load triton/2024.1-gcc
module load mamba
source activate "${SIM_ENV_PREFIX}"
unset R_LIBS R_LIBS_USER R_LIBS_SITE
module load gcc gmake binutils
module use /scratch/cs/bayes_ave/modulefiles
module load cmdstan

if [[ -z "${CMDSTAN:-}" || "${CMDSTAN}" == "${CONDA_PREFIX}"* ]]; then
  echo "ERROR: CMDSTAN must point to the module CmdStan, but is '${CMDSTAN:-unset}'." >&2
  exit 1
fi
echo "Using module CmdStan: CMDSTAN=${CMDSTAN}"

# Record the start time so sacct can be filtered to this run only.
START_TIME=$(date +%Y-%m-%dT%H:%M:%S)
echo "calibration started at ${START_TIME}"

snakemake -p fit_all aggregate_results \
  --configfile config_calibrate.yaml \
  --rerun-triggers mtime input \
  --cores 4 \
  --profile slurm_profile/ \
  --group-components fit_model=1 \
  --set-resources fit_model:runtime=720 fit_model:mem_mb=32000 \
  "$@"

echo
echo "=============================================================="
echo " Slurm accounting for this run (per fit)"
echo "=============================================================="
sacct --starttime="${START_TIME}" --user="${USER}" \
      --format=JobID%20,JobName%18,Elapsed,TotalCPU,MaxRSS,State \
      --units=M 2>/dev/null || echo "sacct unavailable"

echo
Rscript calibration_report.R summaries_calib
