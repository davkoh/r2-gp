## Shared functions for the case studies.
##
## Everything a case study needs besides its own settings, so that the
## Case study folder runs on its own (plain R, no conda, no Snakemake, no
## sim/ folder). Contents:
##
##   1. HSGP rules: kernels, basis count, the inverse-gamma length-scale prior,
##      input scaling
##   2. Stan data assembly for aHSGP_oos.stan
##   3. Held-out ELPD
##   4. Figure style
##   5. Paths and Slurm submission
##   6. The case-study workflow: splits, fitting, report, and the modes that
##      run it locally or on a Slurm cluster (run_case_study())
##   7. Paper figures and tables from saved fits, per study or side by side
##
## A case-study script only sets its settings, prepares its data and calls
## run_case_study(); see case_study_bloodbrain.R.
##
## Sections 1-4 are copied from the simulation pipeline (sim/R/conditions.R,
## scaling.R, stan_data.R, metrics.R and theme_hsgp.R) with the functions the
## simulation alone needs left out, so the case studies use exactly the
## simulation's priors, basis rules, ELPD and figure style. If those change in
## sim/R/, change them here too.

`%||%` <- function(x, y) if (is.null(x)) y else x

## ===========================================================================
## 1. HSGP rules (sim/R/conditions.R, sim/R/scaling.R)
## ===========================================================================

## Kernel identifiers. The integer is the covfn passed to Stan; see
## aHSGP_oos.stan.
KERNELS <- c(se = 0L, m32 = 1L, m52 = 2L, exp = 3L)

## Minimum-M constants from the HSGP basis-size rules:
## m >= const * c / (rho / S). Smoother kernels have lighter spectral tails
## and so need fewer terms; the Materns need roughly twice what the squared
## exponential does.
KERNEL_BASIS_CONSTANT <- c(se = 1.75, m32 = 3.42, m52 = 2.65, exp = 3.42)

check_kernel <- function(kernel) {
  if (!kernel %in% names(KERNELS)) {
    stop("Unknown kernel: ", kernel, ". Expected one of ",
         paste(names(KERNELS), collapse = ", "), call. = FALSE)
  }
  kernel
}

## The inverse-gamma length-scale prior, from the interval it must contain:
## `tail` of the prior mass falls below the lower limit and `tail` above the
## upper limit. An explicit `lengthscale_inv_gamma: [shape, scale]` wins.
.ig_cache <- new.env(parent = emptyenv())

lengthscale_inv_gamma <- function(priors) {
  if (!is.null(priors$lengthscale_inv_gamma)) {
    return(as.numeric(priors$lengthscale_inv_gamma))
  }
  limits <- priors$lengthscale_containment
  if (is.null(limits)) {
    stop("Set priors$lengthscale_containment (or lengthscale_inv_gamma).",
         call. = FALSE)
  }
  lower <- as.numeric(limits[[1]])
  upper <- as.numeric(limits[[2]])
  tail <- as.numeric(priors$lengthscale_containment_tail %||% 0.01)
  stopifnot(lower > 0, upper > lower, tail > 0, tail < 0.5)

  key <- sprintf("%.10g_%.10g_%.10g", lower, upper, tail)
  if (!is.null(.ig_cache[[key]])) {
    return(.ig_cache[[key]])
  }

  ## If X ~ Gamma(shape, rate = 1) then scale / X ~ InvGamma(shape, scale), so
  ## the inverse-gamma tail probabilities are gamma tail probabilities.
  objective <- function(par) {
    shape <- exp(par[1])
    scale <- exp(par[2])
    below_lower <- pgamma(scale / lower, shape, lower.tail = FALSE)
    below_upper <- pgamma(scale / upper, shape, lower.tail = FALSE)
    (below_lower - tail)^2 + (below_upper - (1 - tail))^2
  }
  fit <- optim(c(log(5), log(4)), objective, method = "Nelder-Mead",
               control = list(reltol = 1e-14, maxit = 5000))
  parameters <- c(shape = exp(fit$par[1]), scale = exp(fit$par[2]))

  ## The solve is numerical, so check it rather than trust it.
  achieved <- c(
    pgamma(parameters[["scale"]] / lower, parameters[["shape"]],
           lower.tail = FALSE),
    1 - pgamma(parameters[["scale"]] / upper, parameters[["shape"]],
               lower.tail = FALSE)
  )
  if (max(abs(achieved - tail)) > 1e-4) {
    stop(sprintf(
      paste("Could not solve the inverse-gamma for [%g, %g] at tail %g;",
            "achieved tails %.4g and %.4g."),
      lower, upper, tail, achieved[1], achieved[2]
    ), call. = FALSE)
  }
  .ig_cache[[key]] <- unname(parameters)
  unname(parameters)
}

## Basis count per component. M = ceil(safety * const * c / rho), where rho is
## the SHORTEST length-scale that must be resolved.
basis_count <- function(rho, kernel, boundary, safety = 1.2,
                        m_min = 8L, m_max = 128L) {
  check_kernel(kernel)
  const <- KERNEL_BASIS_CONSTANT[[kernel]]
  m <- ceiling(safety * const * boundary / rho)
  if (m > m_max) {
    warning(sprintf(
      paste("basis_count: %s at rho = %.3f, L = %.2f needs M = %d but m_max",
            "is %d; the HSGP basis-size rule is not satisfied."),
      kernel, rho, boundary, m, m_max
    ), call. = FALSE)
  }
  as.integer(min(max(m, m_min), m_max))
}

## Map each column of x_train to [-1, 1] and apply the same affine map to
## x_test, so that L_d = c for every component.
scale_inputs <- function(x_train, x_test = NULL) {
  x_train <- as.matrix(x_train)

  column_min <- apply(x_train, 2, min)
  column_max <- apply(x_train, 2, max)
  center <- 0.5 * (column_min + column_max)
  halfwidth <- 0.5 * (column_max - column_min)
  ## A constant column has no range to scale; leave it centred at zero rather
  ## than dividing by zero.
  halfwidth[halfwidth == 0] <- 1

  apply_map <- function(x) {
    if (is.null(x)) {
      return(NULL)
    }
    x <- as.matrix(x)
    sweep(sweep(x, 2, center, "-"), 2, halfwidth, "/")
  }

  list(
    train = apply_map(x_train),
    test = apply_map(x_test),
    center = center,
    halfwidth = halfwidth
  )
}

## ===========================================================================
## 2. Stan data (sim/R/stan_data.R)
## ===========================================================================

## Prior models. The integer is the prior_type passed to Stan.
PRIOR_TYPES <- c(R2D2 = 1L, Independent = 2L, RHS = 3L, NG = 4L)

## Length-scale priors. The integer is rho_prior_type.
RHO_PRIOR_TYPES <- c(ig = 1L, vague = 2L, normal = 3L)

## `condition` carries dimension, n_basis, boundary, covfn and rho; `model`
## carries alpha_prior, rho_prior and optionally normal_mean / normal_sd.
make_stan_data <- function(dataset, condition, model, config, prior_only = FALSE) {
  priors <- config$priors
  dimension <- condition$dimension

  if (!model$alpha_prior %in% names(PRIOR_TYPES)) {
    stop("Unknown alpha_prior: ", model$alpha_prior, call. = FALSE)
  }
  if (!model$rho_prior %in% names(RHO_PRIOR_TYPES)) {
    stop("Unknown rho_prior: ", model$rho_prior, call. = FALSE)
  }

  ## Stan expects length-D vectors even when D is 1.
  dimension_vector <- function(value) array(rep(value, dimension), dim = dimension)

  rho_normal_mean <- if (isTRUE(model$normal_centre)) {
    condition$rho
  } else {
    model$normal_mean %||% condition$rho
  }

  ## The inputs span [-1, 1] after scale_inputs(), so the midpoint Stan
  ## subtracts is zero and L = boundary for every dimension.
  list(
    N = nrow(dataset$inputs),
    D = dimension,
    M = condition$n_basis,
    inputs = dataset$inputs,
    y = dataset$y,
    N_test = nrow(dataset$inputs_test),
    inputs_test = dataset$inputs_test,
    y_test = dataset$y_test,

    L = dimension_vector(condition$boundary),
    x_min = dimension_vector(-1),
    x_max = dimension_vector(1),
    is_centered = 1L,
    covfn = condition$covfn,

    prior_only = as.integer(isTRUE(prior_only)),
    prior_type = PRIOR_TYPES[[model$alpha_prior]],
    normalize_variance = as.integer(!identical(model$normalize_variance, FALSE)),

    rho_prior_type = RHO_PRIOR_TYPES[[model$rho_prior]],
    rho_vague_sd = priors$rho_vague_sd,
    rho_normal = c(rho_normal_mean, model$normal_sd %||% priors$rho_normal_sd),
    estimate_rho = as.integer(!identical(config$fit$estimate_rho, FALSE)),
    lengthscale_inv_gamma = lengthscale_inv_gamma(priors),
    rho_fixed = dimension_vector(condition$rho),

    r2_mean = priors$r2_mean,
    r2_precision = priors$r2_precision,
    dirichlet_concentration = priors$dirichlet_concentration,
    independent_alpha_scale = priors$independent_alpha_scale,
    sigma_scale = priors$sigma_scale,
    intercept_scale = priors$intercept_scale,

    ## Regularized horseshoe. p0 = D/2 by default.
    hs_df = priors$hs_df,
    hs_df_global = priors$hs_df_global,
    hs_df_slab = priors$hs_df_slab,
    hs_scale_slab = priors$hs_scale_slab,
    hs_p0 = dimension / 2,

    ## Normal-gamma, matched to the R2D2 prior mean total signal variance.
    ng_nu_scale = (priors$r2_mean / (1 - priors$r2_mean)) / dimension
  )
}

## ===========================================================================
## 3. Held-out ELPD (sim/R/metrics.R)
## ===========================================================================

log_sum_exp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

## Pointwise held-out log predictive density: log of the posterior-mean
## likelihood at each test point. log_lik_test is draws x test points with the
## chains pooled. Not PSIS-LOO: the test points were never used in fitting.
pointwise_elpd <- function(log_lik_test) {
  log_lik_test <- as.matrix(log_lik_test)
  n_draws <- nrow(log_lik_test)
  apply(log_lik_test, 2, log_sum_exp) - log(n_draws)
}

## ===========================================================================
## 4. Figure style (sim/R/theme_hsgp.R)
## ===========================================================================

library(ggplot2)

## Okabe-Ito, colour carrying the prior on alpha.
alpha_prior_colours <- c(
  R2D2 = "#0072B2", Independent = "#D55E00", RHS = "#009E73",
  NG = "#CC79A7", Oracle = "grey35"
)
alpha_prior_levels <- c("R2D2", "RHS", "NG", "Independent", "Oracle")

kernel_labels <- c(
  se = "Squared exponential", m32 = "Matern 3/2",
  m52 = "Matern 5/2", exp = "Exponential"
)

order_alpha_priors <- function(x) {
  factor(x, levels = alpha_prior_levels[alpha_prior_levels %in% unique(x)])
}

## The manuscript's figure style, matched to paper/figures/fig03 and fig04.
theme_paper <- function(base_size = 15.4, grid_y = FALSE) {
  theme_minimal(base_family = "sans", base_size = base_size) +
    theme(
      axis.title = element_text(size = 14.3),
      axis.text = element_text(size = 12.7),
      axis.line = element_line(colour = "grey30", linewidth = 0.3),
      axis.ticks = element_line(colour = "grey30", linewidth = 0.3),
      strip.text = element_text(size = 14.3),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = if (grid_y) {
        element_line(colour = "grey92", linewidth = 0.3)
      } else {
        element_blank()
      },
      legend.position = "bottom",
      legend.text = element_text(size = 13),
      plot.caption = element_text(hjust = 0, size = 9, colour = "grey35")
    )
}

## Manuscript figures are vector PDFs at the paper's own dimensions.
## device = cairo_pdf renders non-ASCII labels such as "Mat\u00e9rn" correctly
## whatever the locale; the default "pdf" device needs a UTF-8 locale for them.
save_paper_plot <- function(path, plot, width, height, device = "pdf") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggsave(path, plot = plot, width = width, height = height, device = device)
  invisible(path)
}

## ===========================================================================
## 5. Paths and Slurm submission
## ===========================================================================

## The folder holding the running script, whether it was started with Rscript,
## source()d (e.g. from RStudio), or run from the repo root or this folder.
case_study_directory <- function(marker = "vardecomp_functions.R") {
  file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE),
                        value = TRUE)
  if (length(file_argument) > 0) {
    ## Rscript encodes spaces in the path ("Case study") as ~+~.
    path <- gsub("~+~", " ", sub("^--file=", "", file_argument[1]),
                 fixed = TRUE)
    return(dirname(normalizePath(path)))
  }
  for (frame in rev(seq_len(sys.nframe()))) {
    path <- sys.frame(frame)$ofile
    if (!is.null(path)) {
      return(dirname(normalizePath(path)))
    }
  }
  candidates <- c(getwd(), file.path(getwd(), "Case study"))
  hit <- candidates[file.exists(file.path(candidates, marker))]
  if (length(hit) > 0) {
    return(normalizePath(hit[1]))
  }
  stop("Cannot locate the Case study folder; run the script from the repo ",
       "root or from inside Case study/.", call. = FALSE)
}

## Write one Slurm batch script and submit it; returns the job id.
##   resources  named list of #SBATCH options, e.g. list(time = "01:00:00")
##   setup      shell lines run before the command (module loads)
submit_slurm_job <- function(name, command, resources, setup, log_dir,
                             working_dir, dependency = NULL, array = NULL) {
  options <- c(
    sprintf("--job-name=%s", name),
    sprintf("--%s=%s", names(resources), unlist(resources)),
    if (!is.null(array)) sprintf("--array=%s", array),
    if (!is.null(dependency)) sprintf("--dependency=%s", dependency)
  )
  ## Given on the sbatch command line rather than as #SBATCH, where a space in
  ## the path ("Case study") would split the option.
  log_file <- file.path(
    log_dir, if (is.null(array)) paste0(name, "-%j.out")
    else paste0(name, "-%A_%a.out")
  )
  script <- c(
    "#!/bin/bash",
    paste("#SBATCH", options),
    "set -euo pipefail",
    setup,
    sprintf("cd %s", shQuote(working_dir)),
    "echo \"$(date) on $(hostname)\"",
    command
  )
  script_file <- file.path(log_dir, paste0(name, ".sbatch"))
  writeLines(script, script_file)
  output <- system2(
    "sbatch",
    c("--parsable", shQuote(paste0("--output=", log_file)), shQuote(script_file)),
    stdout = TRUE, stderr = TRUE
  )
  status <- attr(output, "status")
  if (!is.null(status) && status != 0) {
    stop("sbatch failed for ", name, ": ", paste(output, collapse = "\n"),
         call. = FALSE)
  }
  ## --parsable prints "jobid" or "jobid;cluster".
  sub(";.*", "", trimws(output[length(output)]))
}

## ===========================================================================
## 6. The case-study workflow
## ===========================================================================
##
## The design mirrors the simulation study. Each of `n_splits` random
## train/test splits plays the role of one simulation replication; every prior
## is fitted to the same splits, so the held-out ELPD differences are paired
## per test point, and Fig. 4 averages them over splits with their Monte Carlo
## standard error, as the simulation does over replications.
##
## Modes (first command-line argument of the case-study script):
##   all      (default, and what source() in RStudio runs) prepare, run every
##            fit in parallel on this machine, then report
##   submit   on a Slurm cluster: submit prepare, one array
##            task per fit, and report as three dependent jobs
##   prepare  write the data summaries and compile the model
##   fit i    run fit number i of the grid (or $SLURM_ARRAY_TASK_ID)
##   report   collect the finished fits into the figures and the table

## Where a case study writes: one folder per study inside vardecomp_casestudy/
## (next to the scripts, or $VARDECOMP_OUT if set).
case_study_output_root <- function(case_study_dir) {
  Sys.getenv("VARDECOMP_OUT", file.path(case_study_dir, "vardecomp_casestudy"))
}

## The test-set indices of each random split, from the seed. Shared by the
## fits and by any exploration, so that both see the same splits.
case_study_splits <- function(n_obs, settings) {
  set.seed(settings$seed)
  lapply(seq_len(settings$n_splits), function(split) {
    sort(sample(n_obs, round(settings$test_fraction * n_obs)))
  })
}

## `study` holds:
##   name            short name, used for the output folder and the job names
##   script          path of the case-study script, rerun by the Slurm jobs
##   case_study_dir  the Case study folder
##   settings        the case study's settings list
##   x               covariate matrix, already scaled to [-1, 1]
##   y               the modelled response (e.g. logBBB), unstandardised
##   covariates      column names of x
##   text            list(data, unit, covariate_word, response_tex, label) for
##                   the table caption
##   prepare_extra   optional function(output_path) writing study-specific
##                   summaries in the prepare step (e.g. a covariate screen)
run_case_study <- function(study) {
  ## Lets another script source a case-study script for its settings and
  ## prepared data without running anything (e.g. final_outputs.R).
  if (isTRUE(getOption("vardecomp.load_only"))) {
    return(invisible(study))
  }
  suppressPackageStartupMessages({
    library(cmdstanr)
    library(posterior)
    library(dplyr)
  })

  settings <- study$settings
  x <- study$x
  y <- study$y
  covariates <- study$covariates
  n_obs <- nrow(x)
  dimension <- ncol(x)
  output_dir <- file.path(settings$output_root, study$name)
  output_path <- function(...) file.path(output_dir, ...)

  ## Fixed once, from the seed, so every prior sees the same splits and the
  ## ELPD differences are paired.
  splits <- case_study_splits(n_obs, settings)

  ## The response is standardised with the training mean and sd, so that the
  ## unit-scale priors on sigma, the intercept and alpha mean the same thing
  ## here as in the simulation.
  make_dataset <- function(split) {
    test <- splits[[split]]
    y_train <- y[-test]
    centre <- mean(y_train)
    spread <- sd(y_train)
    list(
      inputs = x[-test, , drop = FALSE],
      y = (y_train - centre) / spread,
      inputs_test = x[test, , drop = FALSE],
      y_test = (y[test] - centre) / spread,
      y_centre = centre,
      y_scale = spread
    )
  }

  ## make_stan_data() takes a simulation condition and a model spec; a
  ## case-study fit is described in the same terms so the priors cannot drift
  ## from the simulation's.
  config <- list(priors = settings$priors, fit = list(estimate_rho = TRUE))

  make_condition <- function(kernel) {
    check_kernel(kernel)
    list(
      id = kernel, kernel = kernel, covfn = KERNELS[[kernel]],
      dimension = dimension,
      boundary = settings$boundary,
      n_basis = basis_count(settings$rho_min, kernel, settings$boundary),
      ## Used only if estimate_rho were FALSE.
      rho = settings$rho_min
    )
  }

  grid <- expand.grid(
    split = seq_len(settings$n_splits),
    alpha_prior = settings$alpha_priors,
    kernel = settings$kernels,
    stringsAsFactors = FALSE
  )[, c("kernel", "alpha_prior", "split")]

  ## The length-scale prior is part of the file name so that switching it
  ## cannot silently reuse fits made under another one.
  fit_stem <- function(i) {
    cell <- grid[i, ]
    rho_tag <- if (settings$rho_prior == "normal") {
      sprintf("normal%g-%g", settings$rho_normal_mean, settings$rho_normal_sd)
    } else {
      settings$rho_prior
    }
    sprintf("%s_%s_%s_split%02d", cell$kernel, cell$alpha_prior, rho_tag,
            cell$split)
  }
  fit_file <- function(i) output_path("fits", paste0(fit_stem(i), ".rds"))
  model_file <- function(i) {
    if (isTRUE(settings$save_fits)) {
      output_path("fits", "models", paste0(fit_stem(i), "_model.rds"))
    }
  }
  ## A fit is done when its summary is on disk, and also its full model when
  ## save_fits asks for one, so turning save_fits on later fills them in.
  fit_done <- function(i) {
    file.exists(fit_file(i)) &&
      (is.null(model_file(i)) || file.exists(model_file(i)))
  }

  ## The executable is built into the output folder rather than next to the
  ## .stan source.
  load_model <- function() {
    dir.create(output_path("stan"), recursive = TRUE, showWarnings = FALSE)
    cmdstan_model(settings$stan_file, dir = output_path("stan"))
  }

  summarise_draws_90 <- function(draws) {
    draws <- as.numeric(draws)
    c(mean = mean(draws), q05 = unname(quantile(draws, 0.05)),
      q95 = unname(quantile(draws, 0.95)))
  }

  fit_one <- function(compiled, i, cores) {
    cell <- grid[i, ]
    condition <- make_condition(cell$kernel)
    model <- list(id = cell$alpha_prior, alpha_prior = cell$alpha_prior,
                  rho_prior = settings$rho_prior,
                  normal_mean = settings$rho_normal_mean,
                  normal_sd = settings$rho_normal_sd)
    dataset <- make_dataset(cell$split)
    stan_data <- make_stan_data(dataset, condition, model, config)
    if (!is.null(settings$hs_p0)) {
      stan_data$hs_p0 <- settings$hs_p0
    }

    mcmc <- settings$mcmc
    adapt_delta <- mcmc$adapt_delta_by_prior[cell$alpha_prior]
    if (is.na(adapt_delta)) {
      adapt_delta <- mcmc$adapt_delta
    }

    started <- proc.time()[["elapsed"]]
    fit <- compiled$sample(
      data = stan_data,
      seed = settings$seed + cell$split,
      chains = mcmc$chains,
      parallel_chains = min(mcmc$chains, cores),
      iter_warmup = mcmc$iter_warmup,
      iter_sampling = mcmc$iter_sampling,
      adapt_delta = unname(adapt_delta),
      max_treedepth = mcmc$max_treedepth,
      refresh = mcmc$refresh,
      show_messages = FALSE
    )
    seconds <- proc.time()[["elapsed"]] - started

    if (!is.null(model_file(i))) {
      ## The same reads that cmdstanr's $save_object() does, so the object no
      ## longer needs the temporary CSV files; stored with what is needed to
      ## reproduce or reuse the fit.
      fit$draws()
      try(fit$sampler_diagnostics(), silent = TRUE)
      try(fit$init(), silent = TRUE)
      try(fit$profiles(), silent = TRUE)
      dir.create(dirname(model_file(i)), recursive = TRUE, showWarnings = FALSE)
      saveRDS(
        list(
          fit = fit, stan_data = stan_data,
          kernel = cell$kernel, alpha_prior = cell$alpha_prior,
          split = cell$split,
          rho_prior = model[c("rho_prior", "normal_mean", "normal_sd")],
          covariates = covariates, test_index = splits[[cell$split]],
          ## To put predictions back on the modelled response scale:
          ## y = y_centre + y_scale * f.
          y_centre = dataset$y_centre, y_scale = dataset$y_scale,
          settings = settings
        ),
        model_file(i)
      )
    }

    summary <- fit$summary(
      variables = c("alpha", "rho", "sigma", "R2_empirical", "lp__"),
      "rhat", "ess_bulk"
    )
    sampler <- fit$diagnostic_summary(quiet = TRUE)

    ## Held-out RMSE of the posterior mean prediction, on the scale of the
    ## modelled response, not the standardised one.
    f_test <- colMeans(as_draws_matrix(fit$draws("f_test")))
    prediction <- dataset$y_centre + dataset$y_scale * f_test
    observed <- dataset$y_centre + dataset$y_scale * dataset$y_test

    ## Share of the total signal variance per covariate, alpha_d^2 / sum
    ## alpha^2: the quantity the R2D2 prior places its Dirichlet on, and
    ## comparable across priors.
    alpha <- as_draws_matrix(fit$draws("alpha"))
    share <- alpha^2 / rowSums(alpha^2)
    rho <- as_draws_matrix(fit$draws("rho"))
    components <- data.frame(
      covariate = covariates,
      share_mean = colMeans(share),
      share_q05 = apply(share, 2, quantile, 0.05),
      share_q95 = apply(share, 2, quantile, 0.95),
      alpha_mean = colMeans(alpha),
      rho_mean = colMeans(rho),
      row.names = NULL
    )

    sigma <- summarise_draws_90(fit$draws("sigma"))
    r2 <- summarise_draws_90(fit$draws("R2_empirical"))
    r2_oos <- summarise_draws_90(fit$draws("R2_oos"))
    pointwise <- pointwise_elpd(as_draws_matrix(fit$draws("log_lik_test")))

    list(
      meta = list(kernel = cell$kernel, alpha_prior = cell$alpha_prior,
                  split = cell$split, n_basis = condition$n_basis,
                  adapt_delta = unname(adapt_delta)),
      metrics = list(
        sigma_mean = sigma[["mean"]], sigma_q05 = sigma[["q05"]],
        sigma_q95 = sigma[["q95"]],
        r2_mean = r2[["mean"]], r2_q05 = r2[["q05"]], r2_q95 = r2[["q95"]],
        r2_oos_mean = r2_oos[["mean"]], r2_oos_q05 = r2_oos[["q05"]],
        r2_oos_q95 = r2_oos[["q95"]],
        rmse_test = sqrt(mean((prediction - observed)^2)),
        elpd_test = sum(pointwise),
        se_elpd_test = sqrt(length(pointwise) * var(pointwise))
      ),
      diagnostics = list(
        max_rhat = max(summary$rhat, na.rm = TRUE),
        min_ess_bulk = min(summary$ess_bulk, na.rm = TRUE),
        divergences = sum(sampler$num_divergent),
        max_treedepth_hits = sum(sampler$num_max_treedepth),
        sampling_seconds = seconds
      ),
      components = components,
      pointwise_elpd_test = pointwise
    )
  }

  ## One grid cell. Skips a fit that is already done, and writes through a
  ## temporary file so that report never reads a half-written result.
  run_task <- function(i, compiled, cores) {
    if (fit_done(i) && !isTRUE(settings$refit)) {
      message(sprintf("Fit %d already done: %s", i, fit_stem(i)))
      return(invisible(fit_file(i)))
    }
    cell <- grid[i, ]
    message(sprintf("Fit %d/%d: %s / %s / split %d", i, nrow(grid),
                    cell$kernel, cell$alpha_prior, cell$split))
    result <- fit_one(compiled, i, cores)
    file <- fit_file(i)
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    partial <- paste0(file, ".partial")
    saveRDS(result, partial)
    file.rename(partial, file)
    message(sprintf("Fit %d done: max Rhat %.3f, %d divergences, %.0fs", i,
                    result$diagnostics$max_rhat, result$diagnostics$divergences,
                    result$diagnostics$sampling_seconds))
    invisible(file)
  }

  prepare <- function() {
    for (folder in c("fits", "figures", "summaries", "stan")) {
      dir.create(output_path(folder), recursive = TRUE, showWarnings = FALSE)
    }
    message(sprintf(
      "prepare %s: N = %d %s, D = %d %s, %d fits in the grid",
      study$name, n_obs, study$text$unit, dimension, study$text$covariate_word,
      nrow(grid)
    ))
    if (!is.null(study$prepare_extra)) {
      study$prepare_extra(output_path)
    }
    ## Compiled once here, so that parallel fits never race to compile.
    invisible(load_model())
  }

  fit_task <- function(i) {
    if (is.na(i) || i < 1 || i > nrow(grid)) {
      stop(sprintf("Fit number must be between 1 and %d.", nrow(grid)),
           call. = FALSE)
    }
    cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK",
                                   as.character(settings$mcmc$chains)))
    run_task(i, load_model(), cores)
  }

  fit_all_local <- function() {
    compiled <- load_model()
    chains <- settings$mcmc$chains
    workers <- settings$parallel_fits %||%
      max(1L, parallel::detectCores(logical = FALSE) %/% chains)
    pending <- which(isTRUE(settings$refit) |
                       !vapply(seq_len(nrow(grid)), fit_done, logical(1)))
    message(sprintf("Running %d of %d fits, %d at a time", length(pending),
                    nrow(grid), workers))
    ## Forked workers: each runs whole fits, and each fit runs its chains in
    ## parallel, so the machine uses about workers x chains cores.
    outcome <- parallel::mclapply(pending, function(i) {
      tryCatch(run_task(i, compiled, chains),
               error = function(e) {
                 structure(conditionMessage(e), class = "failed")
               })
    }, mc.cores = workers, mc.preschedule = FALSE)
    failed <- vapply(outcome, inherits, logical(1), "failed")
    for (k in which(failed)) {
      warning(sprintf("Fit %d failed: %s", pending[k], outcome[[k]]),
              call. = FALSE)
    }
  }

  report <- function() {
    files <- vapply(seq_len(nrow(grid)), fit_file, character(1))
    done <- file.exists(files)
    if (!any(done)) {
      stop("No finished fits in ", output_path("fits"), call. = FALSE)
    }
    if (!all(done)) {
      warning(sprintf("%d of %d fits are missing; reporting the rest.",
                      sum(!done), length(done)), call. = FALSE)
    }
    results <- lapply(files[done], readRDS)
    finished <- grid[done, ]

    bind_component <- function(component) {
      do.call(rbind, lapply(results, function(result) {
        data.frame(result$meta, result[[component]])
      }))
    }
    metrics <- bind_component("metrics")
    diagnostics <- bind_component("diagnostics")
    components <- bind_component("components")

    ## Paired per test point within a split: both priors saw identical data,
    ## so only the predictive difference remains. Positive favours the
    ## competitor, matching the axis label "ELPD - ELPD(reference)".
    pointwise_of <- function(kernel, alpha_prior, split) {
      k <- which(finished$kernel == kernel &
                   finished$alpha_prior == alpha_prior &
                   finished$split == split)
      if (length(k) == 0L) NULL else results[[k]]$pointwise_elpd_test
    }
    elpd <- do.call(rbind, lapply(seq_len(nrow(finished)), function(i) {
      cell <- finished[i, ]
      if (cell$alpha_prior == settings$reference_prior) {
        return(NULL)
      }
      reference <- pointwise_of(cell$kernel, settings$reference_prior,
                                cell$split)
      competitor <- pointwise_of(cell$kernel, cell$alpha_prior, cell$split)
      if (is.null(reference) || is.null(competitor)) {
        return(NULL)
      }
      difference <- competitor - reference
      data.frame(cell, elpd_difference = sum(difference),
                 se_difference = sqrt(length(difference) * var(difference)))
    }))

    write.csv(merge(metrics, diagnostics), output_path("summaries", "fits.csv"),
              row.names = FALSE)
    write.csv(elpd %||% data.frame(),
              output_path("summaries", "elpd_difference.csv"),
              row.names = FALSE)

    ## Averaged over splits, for reading which covariates carry the signal.
    component_summary <- components %>%
      group_by(kernel, alpha_prior, covariate) %>%
      summarise(across(c(share_mean, share_q05, share_q95, alpha_mean,
                         rho_mean), mean), .groups = "drop") %>%
      arrange(kernel, alpha_prior, desc(share_mean))
    write.csv(component_summary,
              output_path("summaries", "variance_decomposition.csv"),
              row.names = FALSE)

    kernel_factor <- function(k) {
      factor(kernel_labels[k], levels = kernel_labels[settings$kernels])
    }

    ## -- Fig. 3 style: posterior R2, in sample and out of sample ------------
    r2_panels <- bind_rows(
      metrics %>% transmute(kernel, alpha_prior, split, estimate = r2_mean,
                            lo = r2_q05, hi = r2_q95, panel = "In sample"),
      ## Residual-based Bayesian R2, as in the simulation study.
      metrics %>% transmute(kernel, alpha_prior, split,
                            estimate = r2_oos_mean, lo = r2_oos_q05,
                            hi = r2_oos_q95, panel = "Out of sample")
    ) %>%
      mutate(
        panel = factor(panel, levels = c("In sample", "Out of sample")),
        kernel_label = kernel_factor(kernel),
        alpha_prior = order_alpha_priors(alpha_prior)
      )

    ## One point per split, dodged within prior, so split-to-split variation
    ## is visible in the same way replications are in the simulation figure.
    dodge <- position_dodge(width = 0.75)
    r2_plot <- ggplot(
      r2_panels,
      aes(kernel_label, estimate, colour = alpha_prior,
          group = interaction(alpha_prior, split))
    ) +
      geom_linerange(aes(ymin = lo, ymax = hi), position = dodge,
                     linewidth = 0.35, alpha = 0.7) +
      geom_point(position = dodge, size = 1.5, alpha = 0.85) +
      scale_colour_manual(values = alpha_prior_colours, name = NULL) +
      coord_cartesian(ylim = c(0, 1)) +
      labs(x = NULL, y = expression(R[emp]^2)) +
      facet_wrap(~panel) +
      theme_paper(grid_y = TRUE)
    save_paper_plot(output_path("figures", "posterior_r2.pdf"), r2_plot,
                    width = 8.0, height = 4.2)

    ## -- Fig. 4 style: held-out ELPD difference against the reference -------
    if (!is.null(elpd)) {
      ## Mean over splits with its Monte Carlo standard error, as fig04
      ## computes them over replications.
      elpd_summary <- elpd %>%
        group_by(kernel, alpha_prior) %>%
        summarise(elpd_diff = mean(elpd_difference),
                  se = sd(elpd_difference) / sqrt(n()),
                  n_splits = n(), .groups = "drop") %>%
        mutate(kernel_label = kernel_factor(kernel),
               alpha_prior = order_alpha_priors(alpha_prior))
      write.csv(elpd_summary, output_path("figures", "elpd_difference.csv"),
                row.names = FALSE)

      ## The window always includes the conventional +-4 band.
      bounds <- c(elpd_summary$elpd_diff - elpd_summary$se,
                  elpd_summary$elpd_diff + elpd_summary$se,
                  elpd_summary$elpd_diff)
      y_lo <- min(-4, bounds, na.rm = TRUE)
      y_hi <- max(4, bounds, na.rm = TRUE)
      dodge_elpd <- position_dodge(width = 0.4)

      elpd_plot <- ggplot(
        elpd_summary, aes(kernel_label, elpd_diff, colour = alpha_prior)
      ) +
        geom_hline(yintercept = 0, linetype = "dotted", colour = "grey45",
                   linewidth = 0.3) +
        geom_hline(yintercept = c(-4, 4), linetype = "dashed",
                   colour = "grey45", linewidth = 0.3) +
        geom_errorbar(aes(ymin = elpd_diff - se, ymax = elpd_diff + se),
                      position = dodge_elpd, width = 0.15, linewidth = 0.55,
                      na.rm = TRUE) +
        geom_point(position = dodge_elpd, size = 2) +
        scale_y_continuous(labels = function(v) sprintf("%.0f", v)) +
        scale_colour_manual(values = alpha_prior_colours, name = NULL) +
        coord_cartesian(ylim = c(y_lo, y_hi)) +
        labs(x = NULL,
             y = sprintf("ELPD - ELPD(%s)", settings$reference_prior)) +
        theme_paper()
      save_paper_plot(output_path("figures", "elpd_difference.pdf"),
                      elpd_plot, width = 6.6, height = 4.0)
    }

    ## -- Table: posterior summaries and sampling statistics -----------------
    table_data <- merge(metrics, diagnostics) %>%
      group_by(kernel, alpha_prior) %>%
      summarise(
        n_splits = n(),
        sigma = mean(sigma_mean), r2 = mean(r2_mean),
        r2_oos = mean(r2_oos_mean), rmse = mean(rmse_test),
        elpd = mean(elpd_test),
        ## Worst case over splits, so one failed chain is not averaged away.
        max_rhat = max(max_rhat), min_ess = min(min_ess_bulk),
        divergences = mean(divergences), seconds = median(sampling_seconds),
        .groups = "drop"
      ) %>%
      mutate(kernel = factor(kernel, levels = settings$kernels),
             alpha_prior = factor(alpha_prior,
                                  levels = settings$alpha_priors)) %>%
      arrange(kernel, alpha_prior)
    write.csv(table_data, output_path("summaries", "fit_summary.csv"),
              row.names = FALSE)

    fmt <- function(v, digits) formatC(v, format = "f", digits = digits)
    kernel_tex <- c(se = "Squared exp.", m32 = "Mat\\'ern $3/2$",
                    m52 = "Mat\\'ern $5/2$", exp = "Exponential")
    alpha_prior_tex <- c(R2D2 = "R2D2", Independent = "Normal", RHS = "RHS",
                         NG = "NG")

    kernel_column <- kernel_tex[as.character(table_data$kernel)]
    kernel_column[duplicated(table_data$kernel)] <- ""
    body <- paste0(paste(
      kernel_column,
      alpha_prior_tex[as.character(table_data$alpha_prior)],
      fmt(table_data$sigma, 3), fmt(table_data$r2, 3),
      fmt(table_data$r2_oos, 3), fmt(table_data$rmse, 3),
      fmt(table_data$elpd, 1), fmt(table_data$max_rhat, 3),
      fmt(table_data$min_ess, 0), fmt(table_data$divergences, 1),
      fmt(table_data$seconds, 0),
      sep = " & "
    ), " \\\\")
    new_kernel <- c(FALSE, diff(as.integer(table_data$kernel)) != 0)
    body[new_kernel] <- paste0("\\midrule\n", body[new_kernel])

    text <- study$text
    caption <- sprintf(
      paste(
        "Additive HSGP fits to %s ($N = %d$ %s, $D = %d$ %s), averaged over",
        "%d random %d/%d train/test splits. $\\sigma$ and $R^2$ are",
        "posterior means on the standardised response; out-of-sample $R^2$",
        "is the residual-based Bayesian $R^2$ on the test set. RMSE is of",
        "the posterior mean prediction of %s on the test set, and ELPD is",
        "the summed held-out log predictive density. $\\widehat{R}$ and ESS",
        "are the worst case over all $\\alpha_d$, $\\rho_d$, $\\sigma$ and",
        "$\\mathrm{lp}$ and over splits; divergences are the mean per fit and",
        "time is the median sampling time in seconds."
      ),
      text$data, n_obs, text$unit, dimension, text$covariate_word,
      settings$n_splits, round(100 * (1 - settings$test_fraction)),
      round(100 * settings$test_fraction), text$response_tex
    )

    header <- c(
      "% Requires \\usepackage{booktabs}.",
      sprintf("%% Generated by Case study/%s -- do not edit by hand.",
              basename(study$script)),
      "\\begin{table}",
      "\\centering\\small",
      "\\setlength{\\tabcolsep}{4pt}",
      sprintf("\\caption{%s}\\label{%s}", caption, text$label),
      "\\begin{tabular}{ll rrr rr rrrr}",
      "\\toprule",
      paste("Kernel & Prior &",
            "\\multicolumn{3}{c}{Posterior} &",
            "\\multicolumn{2}{c}{Held out} &",
            "\\multicolumn{4}{c}{Sampling} \\\\"),
      "\\cmidrule(lr){3-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-11}",
      paste("& & $\\sigma$ & $R^2$ & $R^2_{\\mathrm{oos}}$ & RMSE & ELPD &",
            "$\\widehat{R}$ & ESS & Div. & s \\\\"),
      "\\midrule"
    )
    writeLines(
      c(header, body, "\\bottomrule", "\\end{tabular}", "\\end{table}"),
      output_path("summaries", "fit_summary.tex")
    )

    message(sprintf("report %s: %d of %d fits; outputs in %s", study$name,
                    sum(done), length(done), output_dir))
  }

  ## Three dependent Slurm jobs: prepare (compile once, so the array tasks
  ## never race to compile), one array task per fit with a core per chain,
  ## and report, which runs even if some fits failed (afterany).
  submit <- function() {
    if (!nzchar(Sys.which("sbatch"))) {
      stop("sbatch not found: submit mode is for a Slurm cluster. ",
           "Locally, run without an argument.", call. = FALSE)
    }
    slurm <- settings$slurm
    log_dir <- output_path("logs")
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

    setup <- c(
      if (length(slurm$module_use) > 0) paste("module use", slurm$module_use),
      if (length(slurm$modules) > 0) {
        paste("module load", paste(slurm$modules, collapse = " "))
      },
      ## The jobs write where this submission would have written.
      sprintf("export VARDECOMP_OUT=%s", shQuote(settings$output_root)),
      ## Without TZ, R asks timedatectl for the time zone, which can fail on
      ## a cluster (an OpenSSL version clash) and puts a warning in every
      ## log.
      sprintf("export TZ=\"${TZ:-%s}\"", settings$slurm$time_zone %||% "UTC")
    )
    script <- shQuote(study$script)
    job <- function(step, ...) {
      submit_slurm_job(paste(study$name, step, sep = "_"), ..., setup = setup,
                       log_dir = log_dir, working_dir = study$case_study_dir)
    }

    prepare_id <- job(
      "prepare", paste("Rscript", script, "prepare"),
      resources = list(partition = slurm$partition, time = "01:00:00",
                       `cpus-per-task` = 1, mem = "4G")
    )
    fit_id <- job(
      "fit", paste("Rscript", script, "fit"),
      resources = list(partition = slurm$partition, time = slurm$fit_time,
                       `cpus-per-task` = settings$mcmc$chains,
                       mem = slurm$fit_mem),
      array = sprintf("1-%d%%%d", nrow(grid), slurm$max_parallel),
      dependency = paste0("afterok:", prepare_id)
    )
    report_id <- job(
      "report", paste("Rscript", script, "report"),
      resources = list(partition = slurm$partition, time = "00:30:00",
                       `cpus-per-task` = 1, mem = "4G"),
      dependency = paste0("afterany:", fit_id)
    )
    message(sprintf(
      "Submitted %s: prepare %s, %d fits as array %s, report %s\nLogs: %s",
      study$name, prepare_id, nrow(grid), fit_id, report_id, log_dir
    ))
  }

  arguments <- commandArgs(trailingOnly = TRUE)
  mode <- if (length(arguments) > 0) arguments[[1]] else "all"
  switch(
    mode,
    submit = submit(),
    prepare = prepare(),
    fit = fit_task(as.integer(
      if (length(arguments) > 1) arguments[[2]] else
        Sys.getenv("SLURM_ARRAY_TASK_ID", NA)
    )),
    report = report(),
    all = {
      prepare()
      fit_all_local()
      report()
    },
    stop("Unknown mode '", mode, "'. Use all, submit, prepare, fit or report.",
         call. = FALSE)
  )
  invisible(NULL)
}

## ===========================================================================
## 7. Paper figures and tables, from saved fits
## ===========================================================================
##
## Nothing here fits anything: the per-fit result files are read directly, so
## the output always reflects exactly the fits on disk. Building blocks:
##   read_case_study_fits()   one study's fits and paired ELPD differences
##   plot_elpd_differences()  ELPD of each prior minus the reference prior's
##   plot_posterior_r2()      posterior R2 per prior, in versus out of sample
##   write_latex_table()      a booktabs table, with an optional PDF preview
## Both plots facet by a `dataset` column when the data hold more than one
## study. final_outputs.R builds the paper's figures and table from these.

## One study's fits, filtered by the parts of their file names
## (<kernel>_<prior>_<rho prior>_splitNN.rds), and the held-out ELPD
## differences against the reference prior, paired per test point within each
## split. Positive favours the competitor, matching the axis label.
read_case_study_fits <- function(study_dir, kernel, rho_prior_tag,
                                 reference_prior, label = basename(study_dir)) {
  pattern <- sprintf("^%s_(.+)_%s_split[0-9]+\\.rds$", kernel, rho_prior_tag)
  files <- list.files(file.path(study_dir, "fits"), pattern = pattern,
                      full.names = TRUE)
  if (length(files) == 0) {
    stop("No fits matching ", pattern, " in ", file.path(study_dir, "fits"),
         call. = FALSE)
  }
  results <- lapply(files, readRDS)
  fits <- do.call(rbind, lapply(results, function(result) {
    data.frame(result$meta, result$metrics, result$diagnostics)
  }))

  pointwise <- setNames(lapply(results, `[[`, "pointwise_elpd_test"),
                        paste(fits$alpha_prior, fits$split))
  elpd <- do.call(rbind, lapply(seq_len(nrow(fits)), function(i) {
    if (fits$alpha_prior[i] == reference_prior) {
      return(NULL)
    }
    reference <- pointwise[[paste(reference_prior, fits$split[i])]]
    if (is.null(reference)) {
      return(NULL)
    }
    data.frame(alpha_prior = fits$alpha_prior[i], split = fits$split[i],
               elpd_difference = sum(pointwise[[i]] - reference))
  }))

  ## Say what went in, so a missing or extra fit, or a prior that ran with
  ## different sampler settings, is noticed.
  counts <- table(fits$alpha_prior)
  adapt_delta <- tapply(fits$adapt_delta, fits$alpha_prior, function(v) {
    paste(unique(v), collapse = "/")
  })
  message(label, ": fits per prior: ",
          paste(names(counts), counts, collapse = ", "),
          "; adapt_delta: ", paste(names(adapt_delta), adapt_delta,
                                   collapse = ", "))
  list(fits = fits, elpd = elpd)
}

## theme_bw() at a LaTeX-like font size: with the figure saved at the width it
## is shown at in the document (e.g. \includegraphics[width=\textwidth] and a
## figure width equal to \textwidth), text prints at font_size points, tick
## labels at 0.8 of it, matching an 11pt document at the default of 11.
theme_case_study <- function(figure_settings) {
  theme_bw(base_size = figure_settings$font_size %||% 11) +
    theme(legend.position = "bottom",
          legend.margin = margin(0, 0, 0, 0),
          panel.grid.minor = element_blank())
}

## Factor of display labels in the configured order, keeping only the priors
## present.
prior_label_factor <- function(alpha_prior, figure_settings) {
  present <- intersect(figure_settings$prior_order, unique(alpha_prior))
  factor(figure_settings$prior_labels[alpha_prior],
         levels = figure_settings$prior_labels[present])
}

## Facet by dataset when there is more than one, with each panel's own y
## range for the ELPD (the scales differ with the number of test points).
## An optional `title` labels the whole plot with a facet strip, e.g.
## "Predictive performance"; with several data sets it sits above the data-set
## name in each strip. The layers' data need a `panel_title` column for it
## (add_panel_title()).
facet_datasets <- function(data, scales, title = NULL) {
  several <- "dataset" %in% names(data) && length(unique(data$dataset)) > 1
  if (!is.null(title) && several) {
    facet_wrap(~panel_title + dataset, nrow = 1, scales = scales)
  } else if (!is.null(title)) {
    facet_wrap(~panel_title, scales = scales)
  } else if (several) {
    facet_wrap(~dataset, nrow = 1, scales = scales)
  }
}

add_panel_title <- function(data, title) {
  if (is.null(title)) data else mutate(data, panel_title = title)
}

## Mean over splits of the per-split differences, with the Monte Carlo
## standard error of that mean, as the simulation's Fig. 4.
summarise_elpd_differences <- function(elpd, figure_settings) {
  elpd %>%
    filter(alpha_prior %in% figure_settings$prior_order) %>%
    group_by(across(any_of("dataset")), alpha_prior) %>%
    summarise(elpd_diff = mean(elpd_difference),
              se = sd(elpd_difference) / sqrt(n()),
              n_splits = n(), .groups = "drop")
}

plot_elpd_differences <- function(elpd, figure_settings, title = NULL) {
  elpd <- elpd %>%
    filter(alpha_prior %in% figure_settings$prior_order) %>%
    mutate(prior = prior_label_factor(alpha_prior, figure_settings))
  elpd_summary <- summarise_elpd_differences(elpd, figure_settings) %>%
    mutate(prior = prior_label_factor(alpha_prior, figure_settings)) %>%
    add_panel_title(title)
  elpd <- add_panel_title(elpd, title)

  ## Every panel always shows the conventional +-4 band; invisible points
  ## stretch each panel's y range to include it.
  band <- if ("dataset" %in% names(elpd)) {
    expand.grid(dataset = unique(elpd$dataset), y = c(-4, 4))
  } else {
    data.frame(y = c(-4, 4))
  }
  band$prior <- factor(levels(elpd_summary$prior)[1],
                       levels = levels(elpd_summary$prior))
  band <- add_panel_title(band, title)

  ## Each prior in its own colour, or all in one (figure_settings$elpd_colour).
  prior_colours <- if (is.null(figure_settings$elpd_colour)) {
    alpha_prior_colours
  } else {
    setNames(rep(figure_settings$elpd_colour, length(alpha_prior_colours)),
             names(alpha_prior_colours))
  }

  elpd_plot <- ggplot(elpd_summary, aes(prior, elpd_diff,
                                        colour = alpha_prior)) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey45",
               linewidth = 0.5) +
    geom_hline(yintercept = c(-4, 4), linetype = "dashed", colour = "grey45",
               linewidth = 0.5) +
    geom_blank(data = band, aes(prior, y), inherit.aes = FALSE)
  if (figure_settings$plot_style == "violin") {
    elpd_plot <- elpd_plot +
      geom_violin(data = elpd, aes(prior, elpd_difference, fill = alpha_prior),
                  alpha = 0.25, linewidth = 0.4, scale = "width", width = 0.6,
                  trim = TRUE) +
      geom_point(data = elpd, aes(prior, elpd_difference),
                 position = position_jitter(width = 0.06, height = 0, seed = 1),
                 size = 0.9, alpha = 0.7) +
      geom_point(colour = "black", size = 1.2)
  } else if (figure_settings$plot_style == "box") {
    ## Box: median and interquartile range over splits, whiskers to the most
    ## extreme split within 1.5 IQR. Every split is drawn as a point, so the
    ## box's own outlier marks are switched off; the mean is the black point.
    elpd_plot <- elpd_plot +
      geom_boxplot(data = elpd, aes(prior, elpd_difference, fill = alpha_prior),
                   alpha = 0.25, linewidth = 0.4, width = 0.45,
                   outlier.shape = NA) +
      geom_point(data = elpd, aes(prior, elpd_difference),
                 position = position_jitter(width = 0.08, height = 0, seed = 1),
                 size = 0.9, alpha = 0.7) +
      geom_point(colour = "black", size = 1.2)
  } else {
    elpd_plot <- elpd_plot +
      geom_errorbar(aes(ymin = elpd_diff - se, ymax = elpd_diff + se),
                    width = 0.15, linewidth = 0.55) +
      geom_point(size = 2)
  }
  elpd_plot +
    scale_colour_manual(values = prior_colours, guide = "none") +
    scale_fill_manual(values = prior_colours, guide = "none") +
    ## Whole-number breaks: rounding the labels of half-step breaks (2.5 shown
    ## as "2") would misplace the +-4 band visually.
    scale_y_continuous(breaks = function(limits) {
      breaks <- pretty(limits)
      if (any(breaks != round(breaks))) pretty(limits, n = 3) else breaks
    }) +
    facet_datasets(elpd, scales = "free_y", title = title) +
    ## figure_settings$elpd_axis_title overrides the default "ELPD - ELPD(<ref>)".
    labs(x = NULL,
         y = figure_settings$elpd_axis_title %||%
           sprintf("ELPD - ELPD(%s)",
                   figure_settings$prior_labels[[figure_settings$reference_prior]])) +
    theme_case_study(figure_settings)
}

## In-sample (model-based) and out-of-sample (residual-based) R2, one row
## per fit and sample.
r2_long_format <- function(fits, figure_settings) {
  bind_rows(
    fits %>% transmute(across(any_of("dataset")), alpha_prior, split,
                       estimate = r2_mean, lo = r2_q05, hi = r2_q95,
                       sample = "In sample"),
    fits %>% transmute(across(any_of("dataset")), alpha_prior, split,
                       estimate = r2_oos_mean, lo = r2_oos_q05,
                       hi = r2_oos_q95, sample = "Out of sample")
  ) %>%
    filter(alpha_prior %in% figure_settings$prior_order) %>%
    mutate(prior = prior_label_factor(alpha_prior, figure_settings),
           sample = factor(sample, levels = c("In sample", "Out of sample")))
}

summarise_r2 <- function(r2_long, figure_settings) {
  posterior <- figure_settings$r2_interval == "posterior"
  r2_long %>%
    group_by(across(any_of("dataset")), alpha_prior, prior, sample) %>%
    summarise(
      mean = mean(estimate),
      lo = if (posterior) mean(lo) else mean(estimate) - sd(estimate) / sqrt(n()),
      hi = if (posterior) mean(hi) else mean(estimate) + sd(estimate) / sqrt(n()),
      n_splits = n(), .groups = "drop"
    )
}

plot_posterior_r2 <- function(fits, figure_settings, title = NULL) {
  r2_long <- r2_long_format(fits, figure_settings)
  r2_summary <- summarise_r2(r2_long, figure_settings) %>%
    add_panel_title(title)
  r2_long <- add_panel_title(r2_long, title)
  if (figure_settings$plot_style == "violin") {
    dodge <- position_dodge(width = 0.8)
    r2_plot <- ggplot(r2_long, aes(prior, estimate, colour = sample,
                                   fill = sample)) +
      geom_violin(position = dodge, alpha = 0.25, linewidth = 0.4,
                  scale = "width", width = 0.7, trim = TRUE) +
      geom_point(position = position_jitterdodge(jitter.width = 0.08,
                                                 dodge.width = 0.8, seed = 1),
                 size = 0.9, alpha = 0.7) +
      geom_point(data = r2_summary, aes(prior, mean, group = sample),
                 position = dodge, colour = "black", size = 1.2,
                 show.legend = FALSE)
  } else if (figure_settings$plot_style == "box") {
    ## As for the ELPD: box over the per-split posterior means, every split as
    ## a point, the mean as the black point.
    dodge <- position_dodge(width = 0.8)
    r2_plot <- ggplot(r2_long, aes(prior, estimate, colour = sample,
                                   fill = sample)) +
      geom_boxplot(position = dodge, alpha = 0.25, linewidth = 0.4,
                   width = 0.6, outlier.shape = NA) +
      geom_point(position = position_jitterdodge(jitter.width = 0.1,
                                                 dodge.width = 0.8, seed = 1),
                 size = 0.9, alpha = 0.7) +
      geom_point(data = r2_summary, aes(prior, mean, group = sample),
                 position = dodge, colour = "black", size = 1.2,
                 show.legend = FALSE)
  } else {
    dodge <- position_dodge(width = 0.45)
    r2_plot <- ggplot(r2_summary, aes(prior, mean, colour = sample)) +
      geom_errorbar(aes(ymin = lo, ymax = hi), position = dodge, width = 0.18,
                    linewidth = 0.55) +
      geom_point(position = dodge, size = 2)
  }
  r2_plot +
    scale_colour_manual(values = figure_settings$r2_colours, name = NULL) +
    scale_fill_manual(values = figure_settings$r2_colours, name = NULL) +
    facet_datasets(r2_long, scales = "fixed", title = title) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = NULL, y = expression(R^2)) +
    theme_case_study(figure_settings)
}

## Write a booktabs table fragment to \input{} in the manuscript, and, if
## asked and pdflatex is available, compile a one-page PDF preview next to it.
write_latex_table <- function(lines, path, preview = TRUE) {
  writeLines(lines, path)
  if (!isTRUE(preview)) {
    return(invisible(path))
  }
  if (!nzchar(Sys.which("pdflatex"))) {
    message("pdflatex not found; table preview skipped.")
    return(invisible(path))
  }
  stem <- sub("\\.tex$", "", basename(path))
  preview_stem <- paste0(stem, "_preview")
  directory <- dirname(path)
  writeLines(c(
    "\\documentclass[11pt]{article}",
    "\\usepackage[margin=2.5cm]{geometry}",
    "\\usepackage{booktabs}",
    "\\pagestyle{empty}",
    "\\begin{document}",
    sprintf("\\input{%s}", basename(path)),
    "\\end{document}"
  ), file.path(directory, paste0(preview_stem, ".tex")))
  previous <- setwd(directory)
  status <- system2("pdflatex", c("-interaction=nonstopmode", "-halt-on-error",
                                  paste0(preview_stem, ".tex")),
                    stdout = FALSE, stderr = FALSE)
  setwd(previous)
  unlink(file.path(directory, paste0(preview_stem, c(".aux", ".log", ".tex"))))
  if (status == 0) {
    message("Table preview: ", file.path(directory, paste0(preview_stem, ".pdf")))
  } else {
    warning("pdflatex could not compile the preview of ", basename(path),
            call. = FALSE)
  }
  invisible(path)
}

## Per prior: worst R-hat and smallest ESS over all alpha_d, rho_d, sigma and
## lp and over splits (so one failed chain is not averaged away), mean
## divergences per fit and median sampling time.
summarise_sampling <- function(fits, figure_settings) {
  fits %>%
    filter(alpha_prior %in% figure_settings$prior_order) %>%
    group_by(across(any_of("dataset")), alpha_prior) %>%
    summarise(max_rhat = max(max_rhat), min_ess = min(min_ess_bulk),
              divergences = mean(divergences),
              seconds = median(sampling_seconds), n_splits = n(),
              .groups = "drop")
}

latex_number <- function(v, digits) formatC(v, format = "f", digits = digits)

## ===========================================================================
## 8. Performance and sampling table
## ===========================================================================
##
## One row per prior; a column group per value of `fits$dataset` (data sets,
## or kernels), each with a Predictive part (mean held-out ELPD over splits
## +- its standard error over splits; posterior R2 in and out of sample, each
## the posterior mean with its 90% interval, averaged over splits) and a
## Sampling part (worst R-hat and smallest ESS over all alpha_d, rho_d, sigma
## and lp and over splits, mean divergences per fit, median time). `lead` is
## the caption's first sentence; the definitions are appended. Writes the
## booktabs fragment to `path`, a CSV next to it, and optionally a preview.
## layout = "columns" puts the groups side by side (one column group each);
## "rows" stacks them, each under its own heading, which keeps a table of
## several groups within the text width.
write_performance_table <- function(fits, figure_settings, lead, path, label,
                                    preview = TRUE,
                                    layout = c("columns", "rows")) {
  layout <- match.arg(layout)
  groups <- levels(factor(fits$dataset))
  predictive <- fits %>%
    filter(alpha_prior %in% figure_settings$prior_order) %>%
    group_by(dataset, alpha_prior) %>%
    summarise(elpd = mean(elpd_test), elpd_se = sd(elpd_test) / sqrt(n()),
              r2 = mean(r2_mean), r2_lo = mean(r2_q05), r2_hi = mean(r2_q95),
              r2_oos = mean(r2_oos_mean), r2_oos_lo = mean(r2_oos_q05),
              r2_oos_hi = mean(r2_oos_q95), .groups = "drop")
  table_data <- summarise_sampling(fits, figure_settings) %>%
    left_join(predictive, by = c("dataset", "alpha_prior"))
  write.csv(table_data, sub("\\.tex$", ".csv", path), row.names = FALSE)

  n_columns <- 7
  interval <- function(estimate, lo, hi) {
    sprintf("%s [%s, %s]", latex_number(estimate, 2), latex_number(lo, 2),
            latex_number(hi, 2))
  }
  cells <- function(prior, group) {
    row <- table_data[table_data$alpha_prior == prior &
                        table_data$dataset == group, ]
    if (nrow(row) == 0) {
      return(rep("---", n_columns))
    }
    c(sprintf("$%s \\pm %s$", latex_number(row$elpd, 1),
              latex_number(row$elpd_se, 1)),
      interval(row$r2, row$r2_lo, row$r2_hi),
      interval(row$r2_oos, row$r2_oos_lo, row$r2_oos_hi),
      latex_number(row$max_rhat, 3), latex_number(row$min_ess, 0),
      latex_number(row$divergences, 1), latex_number(row$seconds, 0))
  }
  priors <- intersect(figure_settings$prior_order,
                      unique(table_data$alpha_prior))
  row_of <- function(prior, in_groups) {
    paste0(paste(c(figure_settings$prior_labels[[prior]],
                   unlist(lapply(in_groups, cells, prior = prior))),
                 collapse = " & "), " \\\\")
  }
  if (layout == "columns") {
    body <- vapply(priors, row_of, character(1), in_groups = groups)
    n_groups <- length(groups)
  } else {
    ## A heading row per group, then its priors; a rule between groups.
    body <- unlist(lapply(seq_along(groups), function(i) {
      c(if (i > 1) "\\midrule",
        sprintf("\\multicolumn{%d}{l}{\\textit{%s}} \\\\", n_columns + 1,
                groups[i]),
        vapply(priors, row_of, character(1), in_groups = groups[i]))
    }))
    n_groups <- 1
  }
  first_columns <- 2 + n_columns * (seq_len(n_groups) - 1)
  caption <- paste(
    lead,
    "ELPD is the held-out log predictive density summed over the test set,",
    "as the mean over splits $\\pm$ its standard error over splits. $R^2$",
    "(in sample) and $R^2_{\\mathrm{oos}}$ (residual-based, on the test set)",
    "are posterior means with 90\\% posterior intervals, averaged over",
    "splits. $\\widehat{R}$ and ESS are the worst case over all $\\alpha_d$,",
    "$\\rho_d$, $\\sigma$ and $\\mathrm{lp}$ and over splits; divergences are",
    "the mean per fit and time is the median sampling time in seconds."
  )
  group_rules <- unlist(lapply(first_columns, function(first) {
    sprintf("\\cmidrule(lr){%d-%d} \\cmidrule(lr){%d-%d}",
            first, first + 2, first + 3, first + 6)
  }))
  write_latex_table(c(
    "% Requires \\usepackage{booktabs}.",
    "% Generated by the Case study scripts -- do not edit by hand.",
    "\\begin{table}",
    "\\centering\\small",
    "\\setlength{\\tabcolsep}{3pt}",
    sprintf("\\caption{%s}\\label{%s}", caption, label),
    sprintf("\\begin{tabular}{l %s}",
            paste(rep("ccc rrrr", n_groups), collapse = " ")),
    "\\toprule",
    if (n_groups > 1) {
      c(paste0("& ", paste(sprintf("\\multicolumn{%d}{c}{%s}", n_columns,
                                    groups), collapse = " & "), " \\\\"),
        paste(sprintf("\\cmidrule(lr){%d-%d}", first_columns,
                      first_columns + n_columns - 1), collapse = " "))
    },
    paste0("& ", paste(rep(paste("\\multicolumn{3}{c}{Predictive} &",
                                 "\\multicolumn{4}{c}{Sampling}"), n_groups),
                       collapse = " & "), " \\\\"),
    paste(group_rules, collapse = " "),
    paste0("Prior & ", paste(rep(paste(
      "ELPD & $R^2$ & $R^2_{\\mathrm{oos}}$ &",
      "$\\widehat{R}$ & ESS & Div. & Time (s)"), n_groups), collapse = " & "),
      " \\\\"),
    "\\midrule",
    body,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  ), path, preview = preview)
  invisible(table_data)
}
