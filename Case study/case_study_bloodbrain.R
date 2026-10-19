## Case study: additive HSGP on the BloodBrain data (caret::BloodBrain).
##
## 208 compounds, 134 molecular descriptors (bbbDescr), response logBBB, the
## log ratio of a compound's concentration in the brain to that in the blood.
## D is close to N, which is where the prior on the component variances
## should matter most.
##
## This script sets the case study's settings and prepares its data; the
## workflow itself (splits, fits, figures, table, Slurm submission) is
## run_case_study() in vardecomp_functions.R. Each of 10 random 80/20
## train/test splits plays the role of one simulation replication, and the
## held-out ELPD is the simulation's (pointwise_elpd(), not PSIS-LOO).
##
## Needed: this script, vardecomp_functions.R, data/bloodbrain.rds (a copy of
## caret's BloodBrain data) and ../Stanmodels/aHSGP_oos.stan, the model the
## simulation study fits. Only R is needed: no conda, no Snakemake, no sim/
## folder, no caret.
##
## Outputs, in vardecomp_casestudy/bloodbrain/:
##   summaries/covariate_screening.csv  which descriptors were kept, and why
##   fits/                              one summary file per fit
##   figures/posterior_r2.pdf           Fig. 3 style
##   figures/elpd_difference.pdf        Fig. 4 style
##   summaries/fit_summary.tex          table: posterior and sampling summary
##   summaries/variance_decomposition.csv  each descriptor's variance share
##   logs/                              Slurm job scripts and logs
##
## Modes (first command-line argument):
##   all      (default, and what source() in RStudio runs) prepare, run every
##            fit in parallel on this machine, then report
##   submit   on a Slurm cluster: submit prepare, one array
##            task per fit, and report as three dependent jobs
##   prepare  screen the descriptors and compile the model
##   fit i    run fit number i of the grid (or $SLURM_ARRAY_TASK_ID)
##   report   collect the finished fits into the figures and the table
##
## Local, from the repo root or from Case study/:
##   Rscript "Case study/case_study_bloodbrain.R"
##
## Slurm cluster, one-time setup on the login node (load a C++ compiler
## and R first if they are environment modules):
##   Rscript -e 'install.packages(c("cmdstanr", "posterior", "dplyr",
##     "ggplot2"), repos = c("https://stan-dev.r-universe.dev",
##     "https://cloud.r-project.org"))'
##   Rscript -e 'cmdstanr::install_cmdstan()'     # ~10 min; skipped if a
##                                                # CmdStan module is used
## then, to run:
##   Rscript "Case study/case_study_bloodbrain.R" submit
## The job scripts load the modules in settings$slurm$modules.

## vardecomp_functions.R sits next to this script. It is found whether the
## script was started with Rscript, source()d from RStudio, or run from the
## repo root or from Case study/.
local({
  directory <- NULL
  file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE),
                        value = TRUE)
  if (length(file_argument) > 0) {
    ## Rscript encodes spaces in the path ("Case study") as ~+~.
    directory <- dirname(gsub("~+~", " ", sub("^--file=", "", file_argument[1]),
                              fixed = TRUE))
  } else {
    for (frame in rev(seq_len(sys.nframe()))) {
      if (!is.null(sys.frame(frame)$ofile)) {
        directory <- dirname(sys.frame(frame)$ofile)
        break
      }
    }
  }
  candidates <- c(directory, getwd(), file.path(getwd(), "Case study"))
  functions_file <- file.path(candidates, "vardecomp_functions.R")
  functions_file <- functions_file[file.exists(functions_file)]
  if (length(functions_file) == 0) {
    stop("vardecomp_functions.R not found; run from the repo root or from ",
         "inside Case study/.", call. = FALSE)
  }
  source(functions_file[1], local = globalenv())
})
case_study_dir <- case_study_directory()

## ---------------------------------------------------------------------------
## Settings
## ---------------------------------------------------------------------------

settings <- list(
  ## Outputs go to vardecomp_casestudy/bloodbrain/ next to this script, or to
  ## $VARDECOMP_OUT/bloodbrain/ if that is set.
  output_root = case_study_output_root(case_study_dir),
  ## The simulation study's model, so the two cannot drift apart.
  stan_file = file.path(dirname(case_study_dir), "Stanmodels", "aHSGP_oos.stan"),
  seed = 20261002,
  refit = FALSE,
  ## Also keep every fitted model with all its draws in fits/models/ (tens of
  ## MB per fit). Summaries are always kept.
  save_fits = FALSE,

  ## -- Descriptor screen ----------------------------------------------------
  ## 16 descriptors are integer counts or flags with at most 10 values (e.g.
  ## nbasic, alert, hbdnr); min_unique = 20 removes them, the remaining counts
  ## (rotatablebonds, a_aro) and the coarsest binned surface areas: 29 in all.
  ## Such a descriptor gives a GP component whose length-scale is not
  ## identified.
  min_unique = 20,
  ## Catches descriptors dominated by one value. At 0.5 it drops only
  ## n_sp2 (51% at one value); o_sp2 and o_sp3 sit just under, at 48-49%.
  max_top_share = 0.5,
  ## Only near-duplicates are removed: mw (= weight), fnsa1 (= fpsa1) and
  ## volume (r = 0.992 with surface_area). Strong but imperfect correlation is
  ## kept on purpose, e.g. tpsa and tpsa.1 at r = 0.986: how each prior splits
  ## variance among correlated descriptors is part of what the case study
  ## tests. Lower it (e.g. 0.9, which leaves 71) for a smaller, easier problem.
  ## At the defaults 101 of the 134 descriptors are kept.
  max_correlation = 0.99,
  ## Descriptors to keep or drop regardless of the screen.
  force_keep = character(0),
  force_drop = character(0),
  ## Descriptors to log-transform (log1p, as some are zero) before scaling.
  log_covariates = character(0),

  ## -- Model grid -----------------------------------------------------------
  ## Every kernel is crossed with every prior on alpha and every split.
  ## Kernels: se, m32, m52, exp. Priors on alpha: R2D2, RHS, Independent, NG.
  kernels = c("se", "m32"),
  alpha_priors = c("R2D2", "RHS", "Independent", "NG"),
  ## Prior on the length-scales, shared by every component:
  ##   ig      inverse-gamma containment prior, as in the simulation
  ##   vague   lognormal(0, rho_vague_sd)
  ##   normal  normal(rho_normal_mean, rho_normal_sd), truncated at zero
  ## Inputs are scaled to [-1, 1], so rho is relative to a unit half-range.
  rho_prior = "normal",
  rho_normal_mean = 0.5,
  rho_normal_sd = 0.1,
  ## Fig. 4 compares every prior against this one, within the same kernel.
  reference_prior = "R2D2",
  ## Regularized horseshoe prior guess at the number of active descriptors.
  ## NULL would use D / 2 = 50 as in the simulation, which is implausible for
  ## this many correlated descriptors.
  hs_p0 = 10,

  ## -- Splits ---------------------------------------------------------------
  ## 10 splits, matching the simulation's 10 replications. 20% is ~42 test
  ## compounds per split.
  n_splits = 10,
  test_fraction = 0.2,

  ## -- HSGP approximation ---------------------------------------------------
  ## Inputs are scaled to [-1, 1] over all 208 compounds (no response
  ## information is used), so test compounds never fall outside the training
  ## range and L = boundary for every component. M follows
  ## basis_count() in vardecomp_functions.R: M = 21 for se at the defaults.
  boundary = 2,
  rho_min = 0.2,

  ## Priors, matching the simulation (sim/config.yaml).
  priors = list(
    r2_mean = 0.3333,
    r2_precision = 3,
    dirichlet_concentration = 1,
    independent_alpha_scale = 1,
    sigma_scale = 1,
    intercept_scale = 2,
    lengthscale_containment = c(0.15, 1.5),
    lengthscale_containment_tail = 0.01,
    rho_vague_sd = 3,
    hs_df = 3,
    hs_df_global = 1,
    hs_df_slab = 4,
    hs_scale_slab = 2
  ),

  ## -- MCMC and parallelism -------------------------------------------------
  mcmc = list(
    chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    ## The same for every prior, so that the sampling statistics compare the
    ## priors under identical sampler settings. 0.99 rather than Stan's 0.8 or
    ## the simulation's 0.95, because RHS and NG have a global-local funnel
    ## that needs a small step size.
    adapt_delta = 0.99,
    ## Per-prior overrides, e.g. c(RHS = 0.99). Leave empty as numeric(0)
    ## (not NULL) for none.
    adapt_delta_by_prior = numeric(0),
    max_treedepth = 10,
    refresh = 0
  ),
  ## Fits run at the same time in "all" mode, each with its chains in
  ## parallel. NULL uses (physical cores) / chains. Under Slurm each array task
  ## is one fit and uses $SLURM_CPUS_PER_TASK cores for its chains.
  parallel_fits = NULL,

  ## -- Slurm (submit mode) --------------------------------------------------
  slurm = list(
    ## Environment modules loaded at the start of every job, in this order.
    ## CmdStan needs a C++ compiler to build the model, so gcc is loaded
    ## first; pin the version (e.g. "gcc/12.3.0") if the default one does not
    ## match the compiler CmdStan was installed with.
    modules = c("gcc", "r"),
    ## Extra module paths, e.g. "/path/to/modulefiles" to use a shared
    ## CmdStan module (then add "cmdstan" to modules). Without it,
    ## cmdstanr uses the CmdStan from cmdstanr::install_cmdstan().
    module_use = character(0),
    partition = "batch",
    ## Per fit. Measured on our cluster (SE): R2D2 23-51 min at adapt_delta 0.95,
    ## RHS up to 108 min at 0.99; 0.99 makes every prior slower. Check
    ## `seconds` in summaries/fits.csv after a run and adjust.
    fit_time = "06:00:00",
    fit_mem = "8G",
    ## Array tasks running at the same time.
    max_parallel = 30
  )
)

## ---------------------------------------------------------------------------
## Data and descriptor screen
## ---------------------------------------------------------------------------

skewness <- function(v) {
  centred <- v - mean(v)
  mean(centred^3) / mean(centred^2)^1.5
}

## Returns the screening table; the kept descriptors are those with kept = TRUE.
screen_descriptors <- function(descriptors) {
  screening <- data.frame(
    covariate = names(descriptors),
    n_unique = vapply(descriptors, function(v) length(unique(v)), integer(1)),
    top_share = vapply(descriptors, function(v) max(table(v)) / length(v),
                       numeric(1)),
    integer_valued = vapply(descriptors, function(v) all(v == round(v)),
                            logical(1)),
    skewness = vapply(descriptors, skewness, numeric(1)),
    row.names = NULL
  )

  screening$reason <- with(screening, ifelse(
    covariate %in% settings$force_drop, "forced drop",
    ifelse(n_unique < settings$min_unique,
           sprintf("fewer than %d distinct values", settings$min_unique),
           ifelse(top_share > settings$max_top_share,
                  sprintf("one value shared by more than %.0f%%",
                          100 * settings$max_top_share),
                  NA_character_))
  ))
  screening$reason[screening$covariate %in% settings$force_keep] <- NA_character_

  ## Near-duplicates are resolved among the descriptors that survived the
  ## checks above, keeping the first of each pair in the data's column order so
  ## that the result is deterministic.
  candidates <- screening$covariate[is.na(screening$reason)]
  correlation <- abs(cor(descriptors[candidates]))
  kept <- character(0)
  screening$max_correlation_with_kept <- NA_real_
  for (covariate in candidates) {
    strongest <- if (length(kept) > 0) max(correlation[covariate, kept]) else 0
    partner <- if (length(kept) > 0) {
      kept[which.max(correlation[covariate, kept])]
    }
    row <- screening$covariate == covariate
    screening$max_correlation_with_kept[row] <- strongest
    if (strongest > settings$max_correlation &&
        !covariate %in% settings$force_keep) {
      screening$reason[row] <- sprintf("|r| = %.3f with %s", strongest, partner)
    } else {
      kept <- c(kept, covariate)
    }
  }
  screening$kept <- is.na(screening$reason)
  screening$reason[screening$kept] <- "kept"
  screening
}

## The data are read from data/bloodbrain.rds, an unchanged copy of caret's
## data(BloodBrain) (bbbDescr and logBBB), so that caret, whose dependencies
## need compiling, is not needed on the cluster. caret is only a fallback.
data_file <- file.path(case_study_dir, "data", "bloodbrain.rds")
if (file.exists(data_file)) {
  bloodbrain <- readRDS(data_file)
  bbbDescr <- bloodbrain$bbbDescr
  logBBB <- bloodbrain$logBBB
} else if (requireNamespace("caret", quietly = TRUE)) {
  data("BloodBrain", package = "caret", envir = environment())
} else {
  stop("BloodBrain data not found: copy Case study/data/bloodbrain.rds, or ",
       "install the caret package.", call. = FALSE)
}

screening <- screen_descriptors(bbbDescr)
covariates <- screening$covariate[screening$kept]
x <- as.matrix(bbbDescr[covariates])
for (column in intersect(settings$log_covariates, covariates)) {
  x[, column] <- log1p(x[, column])
}
## Scaled to [-1, 1] over all 208 compounds (no response information is used),
## so test compounds never fall outside the training range.
x <- scale_inputs(x)$train

## Written by the prepare step: the screen, and which descriptors it dropped.
write_screening <- function(output_path) {
  write.csv(screening, output_path("summaries", "covariate_screening.csv"),
            row.names = FALSE)
  dropped <- screening[!screening$kept, ]
  reason_group <- ifelse(
    startsWith(dropped$reason, "|r|"),
    sprintf("near-duplicate, |r| > %g", settings$max_correlation),
    dropped$reason
  )
  for (reason in unique(reason_group)) {
    matches <- dropped$covariate[reason_group == reason]
    message(sprintf("  dropped, %s (%d): %s", reason, length(matches),
                    paste(matches, collapse = ", ")))
  }
}

## ---------------------------------------------------------------------------
## Run
## ---------------------------------------------------------------------------

run_case_study(list(
  name = "bloodbrain",
  script = file.path(case_study_dir, "case_study_bloodbrain.R"),
  case_study_dir = case_study_dir,
  settings = settings,
  x = x,
  y = logBBB,
  covariates = covariates,
  text = list(data = "the BloodBrain data", unit = "compounds",
              covariate_word = "descriptors", response_tex = "logBBB",
              label = "tab:bloodbrain"),
  prepare_extra = write_screening
))
