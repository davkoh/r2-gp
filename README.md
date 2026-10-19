# R2-GP: Variance-Calibrated Priors for Additive Gaussian Processes

Priors on component variances for additive Hilbert-space approximate Gaussian processes (HSGPs). The repo compares the **R2-GP** prior with the regularised horseshoe (RHS), normal-gamma (NG) and half-normal (Normal) priors on the component scales α_d. It includes a simulation study and a real-data case study.

## Structure

| Path | Contents |
|---|---|
| `Stanmodels/aHSGP_oos.stan` | Additive HSGP model with out-of-sample prediction, used by both studies |
| `sim/` | Simulation study: a Snakemake pipeline over dimension, length-scale and kernel (see `sim/README.md`) |
| `Case study/` | BloodBrain case study (208 compounds, 101 descriptors, response logBBB), in plain R |
| `R2D2_HSGP_NORMALIZATION.md` | Derivation of the design-based variance normaliser the Stan model uses |

## Case study

| File | Role |
|---|---|
| `vardecomp_functions.R` | Shared helpers: Stan data, ELPD, plotting, Slurm submission |
| `case_study_bloodbrain.R` | 10 random 80/20 splits × 4 priors × 2 kernels (SE, Matérn 3/2) |
| `case_study_bloodbrain_effects.R` | Full-data fits; per-covariate effect curves and variance shares |
| `final_outputs.R` | Builds the paper's table and figures from saved results |

```bash
Rscript "Case study/case_study_bloodbrain.R"           # fit locally
Rscript "Case study/case_study_bloodbrain.R" submit    # or on Slurm
EFFECTS_KERNEL=se  Rscript "Case study/case_study_bloodbrain_effects.R"
EFFECTS_KERNEL=m32 Rscript "Case study/case_study_bloodbrain_effects.R"
Rscript "Case study/final_outputs.R"                   # paper outputs → vardecomp_casestudy/final/
```

## Reproducing the paper

Every figure and table in the paper is produced by one of the commands below.
Nothing in this repository is a stored result: all outputs are regenerated.

| Paper output | Built by | Lands in |
|---|---|---|
| Figure 1, variance gap | `cd sim && snakemake paper_figures --cores 4` | `sim/figures/fig01_variance_gap.pdf` |
| Figure 2, prior predictive | same rule | `sim/figures/fig02_prior_predictive.pdf` |
| Figure 3, posterior R2 recovery | `cd sim && snakemake all` | `sim/figures/posterior_r2.pdf` |
| Figure 4, ELPD difference | same rule | `sim/figures/elpd_difference.pdf` |
| Parameter recovery table | same rule | `sim/summaries/parameter_recovery.tex` |
| Case-study table and figures | `Rscript "Case study/final_outputs.R"` | `Case study/vardecomp_casestudy/final/` |

Figures 1 and 2 need only the config, so they build in minutes. Figures 3 and 4 and the table need the full simulation grid. Read `sim/README.md` first: it sets out a four-step order for the first cluster run, and the resource settings in  sim/slurm_profile/config.yaml` are estimates that should be measured with `sim/run_calibrate.sh` before the full grid is launched.

## Requirements

R with `cmdstanr` (and CmdStan), `posterior`, `dplyr`, `ggplot2`. The simulation also needs Snakemake.

## License

The code is released under the MIT License (see `LICENSE`). `Case study/data/bloodbrain.rds` is a copy of the `BloodBrain` data from the R package `caret` and is distributed under that package's license (GPL).
