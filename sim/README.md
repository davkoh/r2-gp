# Additive-HSGP simulation pipeline

A Snakemake pipeline comparing variance priors for additive Hilbert-space
approximate Gaussian processes, across dimension, length-scale and kernel.

## What a run computes

For every cell of the grid and every prior model it fits the model to simulated
data with a held-out test set, and produces two figures and one table, all in
the manuscript's own style:

| Output | What | Styled after |
|---|---|---|
| `figures/elpd_difference.pdf` | paired held-out ELPD against the reference model | `paper/figures/fig04_elpd_diff.R` |
| `figures/posterior_r2.pdf` | posterior R2, in sample and out of sample | `paper/figures/fig03_r2_recovery.R` |
| `summaries/parameter_recovery.tex` | RMSE, bias and coverage for `alpha`, `rho`, `sigma`, plus Rhat, ESS, divergences and runtime | booktabs longtable |

Parameter recovery is a table rather than a figure: the numbers are read, not
compared visually, and the sampling diagnostics belong beside them.

The out-of-sample R2 panel uses `R2_oos`, the residual-based Bayesian R2 of
Gelman et al. (2019). The manuscript's figure 3 used the model-based quantity
at the test design, which is what the `\todo` in
`paper/sections/04-simulation.tex` refers to.

Two further figures depend only on the config, not on any fit, and are built by
the `paper_figures` rule:

| Output | What |
|---|---|
| `figures/fig01_variance_gap.pdf` | pointwise versus sample variance of a GP realisation, and the normaliser `Q(rho)` |
| `figures/fig02_prior_predictive.pdf` | implied prior CDF of the empirical R2 under each prior on the component scales, as D grows |

Both are ports of the manuscript scripts, unchanged in appearance. What changes
is that the priors, spectral densities and basis settings now come from
`config.yaml` rather than being restated in the figure script, so a retuned
prior cannot leave a figure stale.

## Quick start

```bash
bash run_local.sh                    # smoke test: ~1 min, exercises every rule
bash run_local.sh config.yaml 8      # the real grid on 8 local cores
sbatch run_snakemake.sh              # the real grid on Slurm
```

Run the pre-flight check before committing to a long run. It costs seconds and
catches an undersized boundary that would otherwise surface after every fit:

```bash
snakemake basis_adequacy --cores 4
```

## Layout

| Path | What |
|---|---|
| `Snakefile` | 7 rules; expands the grid, composes all paths |
| `config.yaml` | the grid, the models, priors, MCMC settings |
| `config_smoke.yaml` | cut-down version for the smoke test |
| `R/` | one file per rule, plus four shared modules |
| `slurm_profile/` | resources and job grouping; the only cluster-specific file |
| `create_py_env/env.yml` | conda environment (Snakemake + R) |
| `data/ results/ summaries/ figures/ models/` | generated, git-ignored |

Shared R modules, sourced rather than run:

| File | What |
|---|---|
| `R/scaling.R` | columnwise input scaling to `[-1, 1]` |
| `R/conditions.R` | grid expansion, and the boundary / basis-count rules |
| `R/dgp.R` | the data-generating process and the spectral densities |
| `R/stan_data.R` | the Stan data interface |
| `R/metrics.R` | RMSE, pooled ELPD, diagnostics, failure stubs |
| `R/theme_hsgp.R` | the shared palette, `theme_hsgp()` and `theme_paper()` |

`theme_paper()` reproduces the manuscript's figure style exactly; every figure
destined for the paper uses it, so the pipeline and `paper/figures/` cannot
drift apart. The repo-root plotting scripts and `paper/figures/fig03`/`fig04`
all source `R/theme_hsgp.R` for their palette.

## Local versus cluster

The only difference is `--profile slurm_profile/`. The Snakefile contains no
Slurm-specific logic; resources, threads and job grouping live in
`slurm_profile/config.yaml`.

On Triton, `run_snakemake.sh` takes CmdStan and the compilers from cluster
modules while Snakemake and R come from the conda env. It deliberately does not
`module load r` (that would shadow the conda R), unsets `R_LIBS*` so the site
library cannot leak in, and fails loudly if `CMDSTAN` still points inside the
conda prefix.

### First run on the cluster

Do not launch the full grid first. Work up in four steps. Each one fails fast
and cheaply if the step before it was wrong.

```bash
# 1. Does the environment work at all?  (~10 min)
mamba env create -f create_py_env/env.yml -p ./env
sbatch run_snakemake.sh --configfile config_smoke.yaml

# 2. Is the grid itself sane?  (seconds, runs on the login node)
snakemake basis_adequacy --cores 4

# 3. What does a fit actually cost?  (~1 h)
sbatch run_calibrate.sh

# 4. Put the measured numbers in slurm_profile/config.yaml, then:
sbatch run_snakemake.sh
```

Step 1 proves the modules, the conda environment, the `CMDSTAN` guard and the
Slurm plugin. Those fail for reasons that have nothing to do with the science,
so find out in ten minutes rather than after a day in the queue.

Step 3 runs 24 fits at the most expensive corner of the grid (Matern 3/2,
rho = 0.3, D = 5, 30 and 100) with the **same MCMC settings as the real run**.
Grouping is off, so `sacct` reports `Elapsed` and `MaxRSS` per fit. It ends by
printing the numbers to paste into the profile.

### Job grouping and Slurm

`group-components: fit_model: 20` submits **one** Slurm job that runs 20 fits.
Slurm sees 90 submissions for 1800 fits, not 1800. Short jobs annoy schedulers,
so this is deliberate.

The allocation is computed, not declared. Snakemake packs the group into layers
that fit the driver's `--cores`, sums cpus and memory inside a layer, and **sums
runtime across layers**. With `--cores 4` and `fit_model` at 4 threads, one fit
fills a layer, so:

```
cpus    = 4                    (one fit at a time)
mem_mb  = 8000                 (one fit at a time)
runtime = 20 x 30 min = 10 h   (SUMMED)
```

**The runtime multiplies.** `group-components x per-job runtime` must stay below
the partition time limit. A per-fit runtime of 8h with 20 components asks for
160 hours and is rejected.

Calibrate `runtime` before the full run. `summaries/parameter_recovery.csv` has
a `seconds` column. Take the slowest cell and double it.

Create the environment once:

```bash
mamba env create -f create_py_env/env.yml -p ./env
```

## Input scaling and the boundary

Every input column is mapped to `[-1, 1]` by `R/scaling.R::scale_inputs()`,
using the **training** design only. The training half-range is then exactly 1,
so Soham's rule `L_d = c * half-range(x_d)` collapses to `L_d = c` for every
dimension: one boundary factor governs the whole design and is guaranteed to
contain it. Test points may land slightly outside `[-1, 1]`, which is the
honest situation for held-out data and still far inside `[-L, L]`.

Both HSGP rules come from Riutort-Mayol et al. (2023), with the kernel
constants from `../varhsgpfns.R`:

```
boundary   c >= 4.5 * rho / S,  S = 1 after scaling   ->  L = c
basis      M >= const * c / rho,  const = 1.75 (SE), 2.65 (M5/2), 3.42 (M3/2)
```

`basis.size_from` decides which length-scale feeds those rules:

- `truth` (default) sizes them at the generating length-scale. Cheap, and the
  approximation is accurate where the data actually live (kernel error 0.0004
  to 0.027 across the grid), but it degrades in the tails of the length-scale
  prior.
- `prior` sizes them for the 5th-95th percentile range of the length-scale
  priors. Conservative and valid anywhere the posterior can wander, but roughly
  6x the cost.

Either way `summaries/basis_adequacy.csv` reports the actual HSGP-vs-exact
kernel error across the prior, so the cost of the cheaper choice is measured
rather than assumed.

## The grid

48 conditions: `D x length-scale x kernel`, minus the cells dropped by
`grid.max_dimension_by_kernel`. The Matérns need about twice the basis
functions of the squared exponential, so `D * M` is largest exactly at
`D = 200` for them; that corner is dropped rather than paid for.

```
D            5, 10, 30, 50, 100, 200   (Matérns capped at 100)
length-scale 0.3, 0.5, 0.7
kernel       se, m32, m52
```

Four models, one `(alpha prior, rho prior)` pair each. `normal_centre: true` puts
the normal prior's mean at the generating length-scale, making it a
well-specified informative prior rather than a misspecified one.

The inverse-gamma is a containment prior with 1% tails at 0.15 and 1.5. It
replaces IG(5,4), whose median was 0.86 and whose 95th percentile was 2.03 — a
length-scale twice the domain width — and which put only 0.3% of its mass below
the shortest length-scale in the grid.

## The data-generating process

Truth is drawn from the same basis the model uses, but with
`basis_gen_multiplier` times as many terms (default 4). That keeps `alpha_d` in
the model's exact normalised convention — the expected in-sample sample
variance of component `d`, as derived in `../R2D2_HSGP_NORMALIZATION.md` — while
making the truth close enough to an exact GP that a fit with too few basis
functions is detectable rather than correct by construction.

## Adding to the pipeline

- **A kernel**: add its spectral density to `spd()` in
  `../Stanmodels/aHSGP_oos.stan` and to `spectral_density()` in `R/dgp.R`, then
  register it in `KERNELS` and `KERNEL_BASIS_CONSTANT` in `R/conditions.R`.
  The two implementations are checked against each other to machine precision.
- **A model**: add an entry to `models:` in `config.yaml`. RHS and normal-gamma
  priors on `alpha` are already implemented in Stan and need only a config line.
- **A grid axis**: extend the lists under `grid:`; the Snakefile expands the
  product and `R/conditions.R::build_conditions()` mirrors it.

## Notes

- Seeding: data uses `seed + rep`, sampling uses `seed + 1000 * rep`, so a
  replicate is reproducible independently of how jobs are grouped.
- A failed fit degrades to a stub with the same shape and `NA` metrics, so
  aggregation and plotting survive and the cell renders as missing rather than
  vanishing from the grid.
- Snakemake tracks code provenance by default and will rerun steps after a
  purely cosmetic edit. Both runner scripts pass `--rerun-triggers mtime`.
- `parameter_recovery.tex` is a `longtable` and needs **two** LaTeX passes to
  settle its column widths. After a single pass the header rules appear
  narrower than the body; that is the rerun warning, not a broken table.
- The Stan source is the tracked artifact; its compiled binary and everything
  under `data/ results/ summaries/ figures/ models/` are git-ignored.
