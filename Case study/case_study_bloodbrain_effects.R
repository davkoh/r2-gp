## BloodBrain: which descriptors matter? Posterior component functions.
##
## The additive HSGP writes logBBB as an intercept plus one smooth function per
## descriptor, f_d(x_d). This script fits the model to ALL 208 compounds (no
## train/test split: the question is inference, not prediction), keeps the
## posterior draws, rebuilds every f_d from them and summarises how much each
## descriptor contributes:
##
##   share        the descriptor's share of the variance of the fitted function
##                over the compounds, var_x(f_d) / sum_d var_x(f_d), posterior
##                mean and 90% interval
##   P(share > 1/D)  posterior probability that it carries more than an even
##                share of the signal
##   band_excl_0  fraction of the observed descriptor range where the 90%
##                credible band of the centred f_d excludes zero
##   range        posterior mean of max - min of f_d over the observed range,
##                on the logBBB scale
##
## The f_d are rebuilt in R from the draws of z, alpha and rho with the same
## basis functions, spectral densities and variance normaliser as
## aHSGP_oos.stan, so the Stan model is unchanged. The rebuild is checked
## against Stan's own fitted values before anything is summarised.
##
## Everything else is the main case study's: case_study_bloodbrain.R is
## sourced without running, so the 101 screened descriptors, their scaling,
## the priors, the length-scale prior, the HSGP basis and the sampler settings
## are identical.
##
## With heavily correlated descriptors, individual f_d are only weakly
## identified: two correlated descriptors can share an effect in many ways.
## Read the shares of a block of correlated descriptors together.
##
## Outputs in vardecomp_casestudy/bloodbrain_effects/:
##   fits/<prior>_draws.rds          thinned posterior draws of each fit
##   summaries/covariate_relevance.csv  all 101 descriptors, every prior
##   summaries/effect_curves.csv     posterior mean and 90% band on a grid
##   figures/effects.pdf             f_d for the strongest descriptors
##   figures/relevance.pdf           variance shares with 90% intervals
##
## Modes (first command-line argument), as for the main case study:
##   all      (default, and what source() in RStudio runs) compile, run every
##            fit in parallel on this machine, then report
##   submit   on a Slurm cluster: compile, one array task per
##            prior, then report, as three dependent jobs
##   prepare  compile the model
##   fit i    fit prior number i of effects_settings$alpha_priors
##   report   rebuild the f_d from the saved draws and write the outputs
##
##   Rscript "Case study/case_study_bloodbrain_effects.R"           # locally
##   Rscript "Case study/case_study_bloodbrain_effects.R" submit    # Slurm

## ---------------------------------------------------------------------------
## The main case study's settings and prepared data, without running it
## ---------------------------------------------------------------------------

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
  script <- file.path(candidates, "case_study_bloodbrain.R")
  script <- script[file.exists(script)]
  if (length(script) == 0) {
    stop("case_study_bloodbrain.R not found; run from the repo root or from ",
         "inside Case study/.", call. = FALSE)
  }
  options(vardecomp.load_only = TRUE)
  on.exit(options(vardecomp.load_only = NULL))
  source(script[1], local = globalenv())
})
## Now available: settings (the main case study's), x (208 x 101, scaled to
## [-1, 1]), covariates, logBBB, bbbDescr, case_study_dir, and the functions
## in vardecomp_functions.R.

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
})

## ---------------------------------------------------------------------------
## Settings of this analysis
## ---------------------------------------------------------------------------

effects_settings <- list(
  output_dir = file.path(settings$output_root, "bloodbrain_effects"),
  ## Covariance function of these fits: "se" or "m32" (or "m52", "exp"). Set
  ## here, or per run with the EFFECTS_KERNEL environment variable, e.g.
  ##   EFFECTS_KERNEL=m32 Rscript "Case study/case_study_bloodbrain_effects.R"
  ## Fits are named by kernel and each kernel's figures and summaries go to
  ## its own folder, bloodbrain_effects/<kernel>/, so both can coexist.
  kernel = Sys.getenv("EFFECTS_KERNEL", "se"),
  ## Priors to fit, one full-data fit each. Independent does not converge
  ## (R-hat 2-3); it is kept to show how its effects differ, not to be read
  ## on its own.
  alpha_priors = c("R2D2", "RHS", "NG", "Independent"),
  ## Display labels, as in the paper figures.
  prior_labels = c(R2D2 = "R2-GP", RHS = "RHS", NG = "NG",
                   Independent = "Normal"),
  ## Descriptors are ranked by this prior's posterior mean share; the same
  ## top descriptors are then shown for every prior.
  rank_by = "R2D2",
  ## Number of descriptors in the effects figure, and in the relevance figure.
  n_effects = 12,
  ## Each effects panel on its own y range (TRUE), which shows the shape of
  ## small effects; or one shared range (FALSE), which keeps their sizes
  ## comparable across descriptors.
  effects_free_y = TRUE,
  n_relevance = 25,
  ## Posterior draws kept per fit (thinned evenly from 4000), and grid points
  ## per descriptor for the effect curves.
  n_draws = 1000,
  n_grid = 100,
  ## Credible level of the bands and intervals.
  level = 0.90,
  ## Fits run at the same time locally; NULL uses (physical cores) / chains.
  parallel_fits = NULL,
  ## Slurm, as in the main case study. A full-data fit has 208 rather than 166
  ## compounds; the main case study's fits took up to ~2 h.
  slurm = modifyList(settings$slurm, list(fit_time = "08:00:00")),
  ## Figures.
  font_size = 11,
  effects_size = c(width = 6.5, height = 5.5),
  ## The separate figure of the strongest descriptors (R2-GP only):
  ## n_top_effects panels in top_effects_ncol columns, effects_r2gp_top6.pdf.
  n_top_effects = 6,
  top_effects_ncol = 2,
  top_effects_size = c(width = 6.5, height = 7.0),
  relevance_size = c(width = 6.5, height = 5.0)
)

output_path <- function(...) file.path(effects_settings$output_dir, ...)
kernel <- check_kernel(effects_settings$kernel)
## Figures and summaries of this kernel; fits, logs and the model are shared.
kernel_path <- function(...) output_path(kernel, ...)
for (folder in c("fits", "stan", "logs")) {
  dir.create(output_path(folder), recursive = TRUE, showWarnings = FALSE)
}
for (folder in c("figures", "summaries")) {
  dir.create(kernel_path(folder), recursive = TRUE, showWarnings = FALSE)
}
n_obs <- nrow(x)
dimension <- ncol(x)

## logBBB standardised on all compounds; the effects are reported back on the
## logBBB scale.
y_centre <- mean(logBBB)
y_scale <- sd(logBBB)
y_std <- (logBBB - y_centre) / y_scale

## The descriptors on their original scale, for the figure axes: the same
## transform and [-1, 1] map as the case study's data preparation.
x_original <- as.matrix(bbbDescr[covariates])
for (column in intersect(settings$log_covariates, covariates)) {
  x_original[, column] <- log1p(x_original[, column])
}
scaling <- scale_inputs(x_original)

## ---------------------------------------------------------------------------
## Fitting
## ---------------------------------------------------------------------------

load_model <- function() {
  cmdstan_model(settings$stan_file, dir = output_path("stan"))
}

## The length-scale prior is part of the name, as in the main case study, so
## that switching it cannot silently reuse draws made under another one.
rho_tag <- if (settings$rho_prior == "normal") {
  sprintf("normal%g-%g", settings$rho_normal_mean, settings$rho_normal_sd)
} else {
  settings$rho_prior
}
draws_file <- function(prior) {
  output_path("fits", sprintf("%s_%s_%s_draws.rds", kernel, prior, rho_tag))
}

fit_prior <- function(prior, compiled, cores) {
  condition <- list(
    id = kernel, kernel = kernel, covfn = KERNELS[[kernel]],
    dimension = dimension, boundary = settings$boundary,
    n_basis = basis_count(settings$rho_min, kernel, settings$boundary),
    rho = settings$rho_min
  )
  model <- list(id = prior, alpha_prior = prior,
                rho_prior = settings$rho_prior,
                normal_mean = settings$rho_normal_mean,
                normal_sd = settings$rho_normal_sd)
  ## The model needs a test set; two training compounds stand in, and their
  ## held-out quantities are simply not used.
  dataset <- list(inputs = x, y = y_std,
                  inputs_test = x[1:2, , drop = FALSE], y_test = y_std[1:2])
  config <- list(priors = settings$priors, fit = list(estimate_rho = TRUE))
  stan_data <- make_stan_data(dataset, condition, model, config)
  if (!is.null(settings$hs_p0)) {
    stan_data$hs_p0 <- settings$hs_p0
  }

  mcmc <- settings$mcmc
  adapt_delta <- mcmc$adapt_delta_by_prior[prior]
  if (is.na(adapt_delta)) {
    adapt_delta <- mcmc$adapt_delta
  }
  started <- proc.time()[["elapsed"]]
  fit <- compiled$sample(
    data = stan_data, seed = settings$seed,
    chains = mcmc$chains, parallel_chains = min(mcmc$chains, cores),
    iter_warmup = mcmc$iter_warmup, iter_sampling = mcmc$iter_sampling,
    adapt_delta = unname(adapt_delta), max_treedepth = mcmc$max_treedepth,
    refresh = mcmc$refresh, show_messages = FALSE
  )
  seconds <- proc.time()[["elapsed"]] - started

  ## Thinned evenly so that the saved draws stay small; f is kept only to
  ## check the rebuild of the components against Stan's own fitted values.
  draws <- as_draws_matrix(fit$draws(
    c("z", "alpha", "rho", "sigma", "intercept", "R2_empirical", "f")
  ))
  keep <- unique(round(seq(1, nrow(draws),
                           length.out = effects_settings$n_draws)))
  summary <- fit$summary(c("alpha", "rho", "sigma", "R2_empirical", "lp__"),
                         "rhat", "ess_bulk")
  sampler <- fit$diagnostic_summary(quiet = TRUE)
  ## Stored as a plain numeric matrix: a column of a draws_matrix keeps that
  ## class, which breaks the matrix algebra of the rebuild.
  draws <- unclass(draws[keep, ])
  attributes(draws) <- attributes(draws)[c("dim", "dimnames")]
  list(
    prior = prior, draws = draws, stan_data = stan_data,
    y_centre = y_centre, y_scale = y_scale, covariates = covariates,
    adapt_delta = unname(adapt_delta),
    diagnostics = list(
      max_rhat = max(summary$rhat, na.rm = TRUE),
      min_ess_bulk = min(summary$ess_bulk, na.rm = TRUE),
      divergences = sum(sampler$num_divergent),
      max_treedepth_hits = sum(sampler$num_max_treedepth),
      sampling_seconds = seconds
    )
  )
}

run_fit <- function(i, compiled, cores) {
  prior <- effects_settings$alpha_priors[[i]]
  file <- draws_file(prior)
  if (file.exists(file) && !isTRUE(settings$refit)) {
    message("Fit ", prior, " already done")
    return(invisible(file))
  }
  message("Fitting ", prior, " to all ", n_obs, " compounds")
  result <- fit_prior(prior, compiled, cores)
  partial <- paste0(file, ".partial")
  saveRDS(result, partial)
  file.rename(partial, file)
  message(sprintf("Fit %s done: max Rhat %.3f, %d divergences, %.0fs", prior,
                  result$diagnostics$max_rhat, result$diagnostics$divergences,
                  result$diagnostics$sampling_seconds))
  invisible(file)
}

## ---------------------------------------------------------------------------
## Rebuilding the component functions
## ---------------------------------------------------------------------------

## As in aHSGP_oos.stan: basis functions on [-L, L] and the spectral densities
## without the alpha^2 factor.
basis_function <- function(L, m, x) (1 / sqrt(L)) * sin(m * pi / (2 * L) * (x + L))
spectral_density <- function(rho, w) {
  switch(
    kernel,
    se = sqrt(2 * pi) * rho * exp(-0.5 * (rho * w)^2),
    m32 = (4 * 3^1.5 / rho^3) * (3 / rho^2 + w^2)^(-2),
    m52 = (4 * 5^2.5 / (0.75 * rho^5)) * (5 / rho^2 + w^2)^(-3),
    exp = (2 * rho) / (1 + (rho * w)^2)
  )
}

## Draws x points matrix of f_d at x_new (scaled inputs), for descriptor d.
## The normaliser uses the TRAINING design, as in Stan, so that alpha_d^2 is
## the expected sample variance of the component over the compounds.
component_draws <- function(draws, stan_data, d, x_new) {
  L <- stan_data$L[d]
  M <- stan_data$M
  m <- seq_len(M)
  omega <- m * pi / (2 * L)
  basis_train <- sapply(m, basis_function, L = L, x = stan_data$inputs[, d])
  basis_new <- sapply(m, basis_function, L = L, x = x_new)
  basis_variance <- apply(basis_train, 2, var)
  alpha <- draws[, sprintf("alpha[%d]", d)]
  rho <- draws[, sprintf("rho[%d]", d)]
  z <- draws[, sprintf("z[%d,%d]", d, m), drop = FALSE]
  density <- outer(rho, omega, spectral_density)        # draws x M
  normaliser <- as.vector(density %*% basis_variance)
  weights <- alpha * sqrt(density / normaliser) * z      # draws x M
  weights %*% t(basis_new)                               # draws x points
}

## Checks that intercept + sum_d f_d reproduces Stan's fitted f at the
## compounds, so the rebuilt components can be trusted.
check_rebuild <- function(result) {
  draws <- result$draws
  total <- matrix(draws[, "intercept"], nrow(draws), n_obs)
  for (d in seq_len(dimension)) {
    total <- total + component_draws(draws, result$stan_data, d, x[, d])
  }
  stan_f <- draws[, sprintf("f[%d]", seq_len(n_obs)), drop = FALSE]
  max(abs(total - stan_f))
}

summarise_prior <- function(result) {
  draws <- result$draws
  lower <- (1 - effects_settings$level) / 2
  upper <- 1 - lower

  ## Variance of each component over the compounds, per draw.
  component_variance <- sapply(seq_len(dimension), function(d) {
    apply(component_draws(draws, result$stan_data, d, x[, d]), 1, var)
  })
  share <- component_variance / rowSums(component_variance)

  curves <- list()
  relevance <- list()
  for (d in seq_len(dimension)) {
    observed <- range(x[, d])
    grid <- seq(observed[1], observed[2], length.out = effects_settings$n_grid)
    ## Centred so that each f_d has mean zero over the compounds, which takes
    ## out the constant the intercept absorbs; on the logBBB scale.
    at_compounds <- component_draws(draws, result$stan_data, d, x[, d])
    on_grid <- component_draws(draws, result$stan_data, d, grid)
    on_grid <- (on_grid - rowMeans(at_compounds)) * result$y_scale
    band <- apply(on_grid, 2, quantile, c(lower, upper))
    excludes_zero <- band[1, ] > 0 | band[2, ] < 0
    curves[[d]] <- data.frame(
      covariate = covariates[d],
      x_scaled = grid,
      x = grid * scaling$halfwidth[d] + scaling$center[d],
      mean = colMeans(on_grid), lo = band[1, ], hi = band[2, ]
    )
    relevance[[d]] <- data.frame(
      covariate = covariates[d],
      share_mean = mean(share[, d]),
      share_lo = unname(quantile(share[, d], lower)),
      share_hi = unname(quantile(share[, d], upper)),
      p_above_even = mean(share[, d] > 1 / dimension),
      band_excl_0 = mean(excludes_zero),
      range_mean = mean(apply(on_grid, 1, function(v) diff(range(v)))),
      alpha_mean = mean(draws[, sprintf("alpha[%d]", d)]),
      rho_mean = mean(draws[, sprintf("rho[%d]", d)])
    )
  }
  label <- effects_settings$prior_labels[[result$prior]]
  list(
    curves = cbind(prior = result$prior, prior_label = label,
                   do.call(rbind, curves)),
    relevance = cbind(prior = result$prior, prior_label = label,
                      do.call(rbind, relevance))
  )
}

## ---------------------------------------------------------------------------
## Report
## ---------------------------------------------------------------------------

report <- function() {
  priors <- effects_settings$alpha_priors
  done <- file.exists(vapply(priors, draws_file, character(1)))
  if (!any(done)) {
    stop("No finished fits in ", output_path("fits"), call. = FALSE)
  }
  if (!all(done)) {
    warning("Missing fits, reported without them: ",
            paste(priors[!done], collapse = ", "), call. = FALSE)
  }
  results <- lapply(priors[done], function(p) {
    result <- readRDS(draws_file(p))
    ## Plain numeric matrix, also for draws saved before this was enforced.
    draws <- unclass(result$draws)
    attributes(draws) <- attributes(draws)[c("dim", "dimnames")]
    result$draws <- draws
    result
  })

  for (result in results) {
    error <- check_rebuild(result)
    message(sprintf(
      "%s: max Rhat %.3f, min ESS %.0f, %d divergences | rebuild vs Stan f: max |difference| %.1e",
      result$prior, result$diagnostics$max_rhat, result$diagnostics$min_ess_bulk,
      result$diagnostics$divergences, error
    ))
    ## CmdStan writes draws to 6 significant digits, so the rebuild matches to
    ## about that precision; anything larger would be a real mismatch.
    if (error > 1e-4) {
      stop("The rebuilt components do not reproduce Stan's f for ",
           result$prior, "; not reporting.", call. = FALSE)
    }
  }

  summaries <- lapply(results, summarise_prior)
  curves <- do.call(rbind, lapply(summaries, `[[`, "curves"))
  relevance <- do.call(rbind, lapply(summaries, `[[`, "relevance"))

  ## Rank by the chosen prior (or the first one fitted) and keep that order
  ## for every prior.
  rank_prior <- if (effects_settings$rank_by %in% relevance$prior) {
    effects_settings$rank_by
  } else {
    relevance$prior[1]
  }
  ranking <- relevance %>%
    filter(prior == rank_prior) %>%
    arrange(desc(share_mean)) %>%
    pull(covariate)
  relevance <- relevance %>%
    mutate(rank = match(covariate, ranking)) %>%
    arrange(rank, prior)
  write.csv(relevance, kernel_path("summaries", "covariate_relevance.csv"),
            row.names = FALSE)
  write.csv(curves, kernel_path("summaries", "effect_curves.csv"),
            row.names = FALSE)

  labels <- effects_settings$prior_labels[priors[done]]
  colours <- setNames(alpha_prior_colours[priors[done]], labels)
  plot_theme <- theme_case_study(list(font_size = effects_settings$font_size))
  prior_factor <- function(v) factor(v, levels = labels)

  ## -- Effects of the strongest descriptors ---------------------------------
  top <- head(ranking, effects_settings$n_effects)
  effect_data <- curves %>%
    filter(covariate %in% top) %>%
    mutate(covariate = factor(covariate, levels = top),
           prior_label = prior_factor(prior_label))
  rug <- data.frame(
    covariate = factor(rep(top, each = n_obs), levels = top),
    x = as.vector(x_original[, top])
  )
  ## One panel per descriptor in `data`, in `ncol` columns; the rug is limited
  ## to the same descriptors so that no empty panels appear.
  plot_effects <- function(data, ncol = 4) {
    shown <- levels(droplevels(data$covariate))
    data <- mutate(data, covariate = factor(covariate, levels = shown))
    rug_shown <- rug %>%
      filter(covariate %in% shown) %>%
      mutate(covariate = factor(covariate, levels = shown))
    ggplot(data, aes(x, mean)) +
      geom_hline(yintercept = 0, linetype = "dotted", colour = "grey45") +
      geom_ribbon(aes(ymin = lo, ymax = hi, fill = prior_label), alpha = 0.2) +
      geom_line(aes(colour = prior_label), linewidth = 0.6) +
      geom_rug(data = rug_shown, aes(x = x), inherit.aes = FALSE, sides = "b",
               alpha = 0.3, length = unit(0.03, "npc")) +
      facet_wrap(~covariate, ncol = ncol,
                 scales = if (isTRUE(effects_settings$effects_free_y)) "free" else "free_x") +
      scale_colour_manual(values = colours, name = NULL) +
      scale_fill_manual(values = colours, name = NULL) +
      labs(x = "Covariate value",
           y = expression("Effect on logBBB, " * f[d](x) - bar(f)[d])) +
      plot_theme
  }
  ## All fitted priors overlaid.
  save_paper_plot(kernel_path("figures", "effects.pdf"),
                  plot_effects(effect_data),
                  width = effects_settings$effects_size[["width"]],
                  height = effects_settings$effects_size[["height"]])
  ## R2-GP alone, so that each panel is scaled to its own band rather than to
  ## the much wider bands of the other priors. Same descriptors and order.
  single_label <- effects_settings$prior_labels[["R2D2"]]
  if (single_label %in% effect_data$prior_label) {
    single <- filter(effect_data, prior_label == single_label)
    save_paper_plot(kernel_path("figures", "effects_r2gp.pdf"),
                    plot_effects(single) + theme(legend.position = "none"),
                    width = effects_settings$effects_size[["width"]],
                    height = effects_settings$effects_size[["height"]])
    ## The strongest few on their own, as a 3 x 2 grid (3 rows, 2 columns).
    top_few <- head(top, effects_settings$n_top_effects)
    save_paper_plot(kernel_path("figures", "effects_r2gp_top6.pdf"),
                    plot_effects(filter(single, covariate %in% top_few),
                                 ncol = effects_settings$top_effects_ncol) +
                      theme(legend.position = "none"),
                    width = effects_settings$top_effects_size[["width"]],
                    height = effects_settings$top_effects_size[["height"]])
  }

  ## -- Variance shares with intervals ---------------------------------------
  shown <- head(ranking, effects_settings$n_relevance)
  relevance_data <- relevance %>%
    filter(covariate %in% shown) %>%
    mutate(covariate = factor(covariate, levels = rev(shown)),
           prior_label = prior_factor(prior_label))
  dodge <- position_dodge(width = 0.6)
  relevance_plot <- ggplot(relevance_data,
                           aes(share_mean, covariate, colour = prior_label)) +
    geom_vline(xintercept = 1 / dimension, linetype = "dashed",
               colour = "grey45") +
    geom_errorbarh(aes(xmin = share_lo, xmax = share_hi), position = dodge,
                   height = 0, linewidth = 0.5) +
    geom_point(position = dodge, size = 1.5) +
    scale_colour_manual(values = colours, name = NULL) +
    labs(x = "Share of the variance of the fitted function",
         y = NULL,
         caption = sprintf(
           "Dashed line: an even share, 1/D = 1/%d. Intervals: %d%% posterior.",
           dimension, round(100 * effects_settings$level)
         )) +
    plot_theme
  save_paper_plot(kernel_path("figures", "relevance.pdf"), relevance_plot,
                  width = effects_settings$relevance_size[["width"]],
                  height = effects_settings$relevance_size[["height"]])

  ## -- Console summary ------------------------------------------------------
  message("\nStrongest descriptors by ", rank_prior,
          " (share of fitted variance, mean and ",
          round(100 * effects_settings$level), "% interval):")
  print(relevance %>%
          filter(prior == rank_prior, rank <= effects_settings$n_effects) %>%
          transmute(rank, covariate, share = round(share_mean, 3),
                    lo = round(share_lo, 3), hi = round(share_hi, 3),
                    p_above_even = round(p_above_even, 2),
                    band_excl_0 = round(band_excl_0, 2),
                    range = round(range_mean, 2)),
        row.names = FALSE)
  message("Outputs in ", effects_settings$output_dir)
}

## ---------------------------------------------------------------------------
## Modes
## ---------------------------------------------------------------------------

fit_all_local <- function() {
  compiled <- load_model()
  chains <- settings$mcmc$chains
  workers <- effects_settings$parallel_fits %||%
    max(1L, parallel::detectCores(logical = FALSE) %/% chains)
  outcome <- parallel::mclapply(
    seq_along(effects_settings$alpha_priors),
    function(i) tryCatch(run_fit(i, compiled, chains),
                         error = function(e) {
                           structure(conditionMessage(e), class = "failed")
                         }),
    mc.cores = min(workers, length(effects_settings$alpha_priors)),
    mc.preschedule = FALSE
  )
  for (k in which(vapply(outcome, inherits, logical(1), "failed"))) {
    warning("Fit ", effects_settings$alpha_priors[k], " failed: ",
            outcome[[k]], call. = FALSE)
  }
}

submit <- function() {
  if (!nzchar(Sys.which("sbatch"))) {
    stop("sbatch not found: submit mode is for a Slurm cluster. ",
         "Locally, run without an argument.", call. = FALSE)
  }
  slurm <- effects_settings$slurm
  log_dir <- output_path("logs")
  setup <- c(
    if (length(slurm$module_use) > 0) paste("module use", slurm$module_use),
    if (length(slurm$modules) > 0) {
      paste("module load", paste(slurm$modules, collapse = " "))
    },
    sprintf("export VARDECOMP_OUT=%s", shQuote(settings$output_root)),
    sprintf("export EFFECTS_KERNEL=%s", kernel),
    sprintf("export TZ=\"${TZ:-%s}\"", slurm$time_zone %||% "UTC")
  )
  script <- shQuote(file.path(case_study_dir,
                              "case_study_bloodbrain_effects.R"))
  job <- function(step, ...) {
    submit_slurm_job(paste0("bb_effects_", step), ..., setup = setup,
                     log_dir = log_dir, working_dir = case_study_dir)
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
    array = sprintf("1-%d", length(effects_settings$alpha_priors)),
    dependency = paste0("afterok:", prepare_id)
  )
  report_id <- job(
    "report", paste("Rscript", script, "report"),
    resources = list(partition = slurm$partition, time = "01:00:00",
                     `cpus-per-task` = 1, mem = "8G"),
    dependency = paste0("afterany:", fit_id)
  )
  message(sprintf("Submitted: prepare %s, %d fits as array %s, report %s",
                  prepare_id, length(effects_settings$alpha_priors), fit_id,
                  report_id))
}

arguments <- commandArgs(trailingOnly = TRUE)
mode <- if (length(arguments) > 0) arguments[[1]] else "all"
switch(
  mode,
  prepare = invisible(load_model()),
  fit = {
    i <- as.integer(if (length(arguments) > 1) arguments[[2]] else
      Sys.getenv("SLURM_ARRAY_TASK_ID", NA))
    if (is.na(i) || i < 1 || i > length(effects_settings$alpha_priors)) {
      stop("Fit number must be between 1 and ",
           length(effects_settings$alpha_priors), call. = FALSE)
    }
    cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK",
                                   as.character(settings$mcmc$chains)))
    run_fit(i, load_model(), cores)
  },
  report = report(),
  submit = submit(),
  all = {
    fit_all_local()
    report()
  },
  stop("Unknown mode '", mode, "'. Use all, submit, prepare, fit or report.",
       call. = FALSE)
)
