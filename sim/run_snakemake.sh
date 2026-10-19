#!/bin/bash
#SBATCH --job-name=hsgp_r2_sim
#SBATCH --time=48:00:00            # Must exceed the whole workflow's wall time.
#SBATCH --cpus-per-task=4          # For the Snakemake driver itself.
#SBATCH --mem=5G                   # For the Snakemake driver itself.
#SBATCH --output=snakemake_slurm-%j.out
#SBATCH --error=snakemake_slurm-%j.err

# Cluster driver for the additive-HSGP simulation (Aalto Triton).
#
# Submit with:   sbatch sim/run_snakemake.sh
#
# The toolchain is split deliberately: CmdStan and the compilers come from
# cluster modules, while Snakemake AND R come from the conda env built from
# create_py_env/env.yml. We do NOT `module load r` -- that would shadow the
# conda R and its bundled packages.

set -euo pipefail

# Snakemake must run from sim/, whether submitted from the repo root or sim/.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  cd "${SLURM_SUBMIT_DIR}"
  if [[ ! -f Snakefile && -f sim/Snakefile ]]; then
    cd sim
  fi
else
  cd "$(dirname "$0")"
fi

# Override SIM_ENV_PREFIX if /scratch has hit its inode quota and the env had
# to be created under $HOME instead.
SIM_ENV_PREFIX="${SIM_ENV_PREFIX:-./env}"
if [[ ! -d "${SIM_ENV_PREFIX}/conda-meta" && -d "${HOME}/.conda/envs/hsgp_r2_sim/conda-meta" ]]; then
  SIM_ENV_PREFIX="${HOME}/.conda/envs/hsgp_r2_sim"
fi

module load triton/2024.1-gcc
module load mamba
source activate "${SIM_ENV_PREFIX}"

# Keep the site R library out: a stale system-wide package here silently
# shadows the conda one and the failure surfaces much later.
unset R_LIBS R_LIBS_USER R_LIBS_SITE

module load gcc gmake binutils
module use /scratch/cs/bayes_ave/modulefiles
module load cmdstan

# cmdstanr will happily use a CmdStan bundled inside the conda env, which is
# built with the wrong toolchain here. Fail loudly rather than compile against it.
if [[ -z "${CMDSTAN:-}" || "${CMDSTAN}" == "${CONDA_PREFIX}"* ]]; then
  echo "ERROR: CMDSTAN must point to the module CmdStan, but is '${CMDSTAN:-unset}'." >&2
  exit 1
fi
echo "Using module CmdStan: CMDSTAN=${CMDSTAN}"

# Check the boundary and basis rules before spending days on fits.
snakemake -p basis_adequacy --cores 4

# --rerun-triggers mtime input: "mtime" alone stops Snakemake noticing that the
# SET of input files changed, so adding a model or a condition silently produces
# no work even though its result files are missing. Adding "input" restores that
# check while still ignoring cosmetic edits to the R scripts, which would
# otherwise retrigger every fit.
# --cores 4 is not just local parallelism: it is the `_cores` constraint that
# Snakemake uses to pack a group job into layers. With fit_model at 4 threads,
# --cores 4 runs the group's fits ONE AT A TIME, so the allocation stays at
# 4 cpus and the walltime becomes components x per-fit runtime. Raising it
# makes fits run in parallel inside one allocation instead, multiplying cpus
# and memory. See slurm_profile/config.yaml.
snakemake -p all prior_all \
  --rerun-triggers mtime input \
  --cores 4 \
  --profile slurm_profile/ \
  "$@"
