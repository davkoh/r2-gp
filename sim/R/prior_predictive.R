## Stan prior_only cross-check of the variance normalisation.
##
## Figure 2 (R/paper_figures.R) is the presentational prior-predictive output.
## This script exists for a different purpose: to confirm that Stan's own
## implied_R2 and R2_empirical agree, and therefore that alpha[d]^2 really is
## delivering the expected sample variance that
## ../../R2D2_HSGP_NORMALIZATION.md derives. If the R-side and Stan-side
## normalisations ever drift apart, this is where it shows up.
##
## Run on a subset of conditions: it is an implementation check, not a sweep.

suppressPackageStartupMessages(library(posterior))
library(cmdstanr)

source("R/dgp.R")
source("R/stan_data.R")

config <- yaml::read_yaml(snakemake@params[["config_file"]])
conditions <- build_conditions(config)
settings <- config$prior_predictive %||% list()
stan_conditions <- settings$stan_conditions %||% character(0)

rows <- list()
if (length(stan_conditions) > 0) {
  compiled <- cmdstan_model(snakemake@input[["stan_file"]])
  for (condition_id in stan_conditions) {
    condition <- find_condition(conditions, condition_id)
    dataset <- simulate_dataset(condition, config$design,
                                seed = as.integer(config$seed))
    for (model in config$models) {
      fit <- compiled$sample(
        data = make_stan_data(dataset, condition, model, config,
                              prior_only = TRUE),
        seed = as.integer(config$seed), chains = 2, parallel_chains = 2,
        iter_warmup = as.integer(settings$stan_warmup %||% 500L),
        iter_sampling = as.integer(settings$stan_sampling %||% 500L),
        refresh = 0
      )
      ## implied_R2 uses dot_self(alpha) - the variance the prior INTENDS.
      ## R2_empirical uses variance(f) - the variance actually realised on the
      ## design. The normalisation is what makes the two agree.
      implied <- as.numeric(as_draws_matrix(fit$draws("implied_R2")))
      empirical <- as.numeric(as_draws_matrix(fit$draws("R2_empirical")))
      rows[[length(rows) + 1L]] <- data.frame(
        condition_id = condition_id, dimension = condition$dimension,
        kernel = condition$kernel, rho = condition$rho, model_id = model$id,
        alpha_prior = model$alpha_prior, rho_prior = model$rho_prior,
        implied_r2_median = median(implied),
        empirical_r2_median = median(empirical),
        median_gap = median(empirical) - median(implied),
        stringsAsFactors = FALSE
      )
    }
  }
}

output <- if (length(rows) > 0) do.call(rbind, rows) else data.frame()
dir.create(dirname(snakemake@output[["stan_csv"]]), recursive = TRUE,
           showWarnings = FALSE)
write.csv(output, snakemake@output[["stan_csv"]], row.names = FALSE)

if (nrow(output) > 0) {
  worst <- max(abs(output$median_gap))
  message(sprintf(
    "prior_predictive: %d Stan cells, largest |R2_empirical - implied_R2| = %.3f",
    nrow(output), worst
  ))
  print(output[c("condition_id", "model_id", "implied_r2_median",
                 "empirical_r2_median", "median_gap")], row.names = FALSE)
} else {
  message("prior_predictive: no stan_conditions configured; nothing to check.")
}
