## Fit one model to one replicate of one condition.
##
## Only summaries leave this script. The alpha and rho draws are ~13 MB per fit
## at D = 200, so the recovery metrics are computed here and the draws
## discarded; the one array that must survive is the pointwise held-out log
## predictive density, which the paired ELPD differences need at aggregation.

library(cmdstanr)
library(posterior)

source("R/stan_data.R")
source("R/metrics.R")

config <- yaml::read_yaml(snakemake@params[["config_file"]])
condition <- find_condition(
  build_conditions(config), snakemake@params[["condition_id"]]
)
model <- find_model(config, snakemake@params[["model_id"]])
replication <- as.integer(snakemake@params[["replication"]])
dataset <- readRDS(snakemake@input[["dataset"]])
output_file <- snakemake@output[["result"]]

mcmc <- config$mcmc
## Priors with a global-local funnel (RHS, NG) need a smaller step size to
## avoid divergences. Under-tuning a model would handicap it in the comparison,
## so each gets what it needs and the setting is recorded in the output.
adapt_delta <- model$adapt_delta %||% mcmc$adapt_delta

fit_one <- function() {
  compiled <- cmdstan_model(snakemake@input[["stan_file"]])
  stan_data <- make_stan_data(dataset, condition, model, config)

  started <- proc.time()[["elapsed"]]
  fit <- compiled$sample(
    data = stan_data,
    seed = as.integer(config$seed) + 1000L * replication,
    chains = as.integer(mcmc$chains),
    parallel_chains = min(as.integer(mcmc$parallel_chains), snakemake@threads),
    iter_warmup = as.integer(mcmc$iter_warmup),
    iter_sampling = as.integer(mcmc$iter_sampling),
    adapt_delta = as.numeric(adapt_delta),
    max_treedepth = as.integer(mcmc$max_treedepth),
    refresh = as.integer(mcmc$refresh %||% 0L)
  )
  sampling_seconds <- proc.time()[["elapsed"]] - started

  ## One summary call covers both the recovery metrics and the diagnostics.
  summary <- fit$summary(
    variables = c("alpha", "rho", "sigma", "R2_empirical", "lp__"),
    "mean", "sd", "rhat", "ess_bulk", "ess_tail"
  )

  alpha_draws <- as_draws_matrix(fit$draws("alpha"))
  rho_draws <- as_draws_matrix(fit$draws("rho"))
  sigma_draws <- as_draws_matrix(fit$draws("sigma"))
  alpha <- recovery_metrics(alpha_draws, dataset$alpha_true)
  rho <- recovery_metrics(rho_draws, dataset$rho_true)
  sigma <- recovery_metrics(sigma_draws, dataset$sigma_true)

  r2 <- scalar_summary(fit$draws("R2_empirical"))
  r2_oos <- scalar_summary(fit$draws("R2_oos"))
  pointwise <- pointwise_elpd(as_draws_matrix(fit$draws("log_lik_test")))

  alpha_block <- block_diagnostics(summary, "alpha")
  rho_block <- block_diagnostics(summary, "rho")
  sampler <- fit$diagnostic_summary(quiet = TRUE)

  list(
    meta = list(
      condition_id = condition$id, dimension = condition$dimension,
      rho = condition$rho, kernel = condition$kernel,
      boundary = condition$boundary, n_basis = condition$n_basis,
      model_id = model$id, model_label = model$label,
      alpha_prior = model$alpha_prior, rho_prior = model$rho_prior,
      replication = replication, adapt_delta = adapt_delta,
      ok = TRUE, error = NA_character_
    ),
    metrics = list(
      rmse_alpha_point = alpha$rmse_point, rmse_alpha_draws = alpha$rmse_draws,
      bias_alpha = alpha$bias, coverage_alpha = alpha$coverage,
      rmse_rho_point = rho$rmse_point, rmse_rho_draws = rho$rmse_draws,
      bias_rho = rho$bias, coverage_rho = rho$coverage,
      rmse_sigma_point = sigma$rmse_point, rmse_sigma_draws = sigma$rmse_draws,
      bias_sigma = sigma$bias, coverage_sigma = sigma$coverage,
      ## In-sample R2 is model-based; out-of-sample R2 is residual-based
      ## (Gelman et al. 2019), so a model that predicts badly is penalised.
      r2_mean = r2$mean, r2_median = r2$median,
      r2_q05 = r2$q05, r2_q95 = r2$q95,
      r2_oos_mean = r2_oos$mean, r2_oos_median = r2_oos$median,
      r2_oos_q05 = r2_oos$q05, r2_oos_q95 = r2_oos$q95,
      true_r2 = dataset$true_r2,
      elpd_test = sum(pointwise),
      mlpd_test = mean(pointwise),
      se_elpd_test = sqrt(length(pointwise) * var(pointwise)),
      oracle_mlpd_test = dataset$oracle_elpd_test
    ),
    diagnostics = list(
      max_rhat = max(summary$rhat, na.rm = TRUE),
      min_ess_bulk = min(summary$ess_bulk, na.rm = TRUE),
      max_rhat_alpha = alpha_block$max_rhat,
      min_ess_bulk_alpha = alpha_block$min_ess_bulk,
      max_rhat_rho = rho_block$max_rhat,
      min_ess_bulk_rho = rho_block$min_ess_bulk,
      divergences = sum(sampler$num_divergent),
      max_treedepth_hits = sum(sampler$num_max_treedepth),
      sampling_seconds = sampling_seconds
    ),
    pointwise_elpd_test = pointwise
  )
}

result <- tryCatch(fit_one(), error = function(e) {
  warning("fit failed for ", condition$id, " / ", model$id, " / rep ",
          replication, ": ", conditionMessage(e))
  failed_result(conditionMessage(e), condition, model, replication,
                n_test = nrow(dataset$inputs_test))
})

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(result, output_file)

message(sprintf(
  "fit_model: %s / %s / rep %d - %s",
  condition$id, model$id, replication,
  if (isTRUE(result$meta$ok)) {
    sprintf("max Rhat %.3f, %.0fs", result$diagnostics$max_rhat,
            result$diagnostics$sampling_seconds)
  } else {
    "FAILED"
  }
))
