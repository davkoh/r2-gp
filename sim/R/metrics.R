## Evaluation metrics for one fit.
##
## Everything here is computed inside the fitting job and only summaries are
## written to disk: at D = 200 the alpha and rho draws alone are ~13 MB per
## fit, which is not worth keeping 1920 times over.

## Numerically stable log-sum-exp, as in ../../aHSGP_oos_simulation.R.
log_sum_exp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

## Draw-level recovery statistics, ported off rstan from Soham's
## ../../varhsgpfns.R:73-76. `draws` is draws x parameters, `truth` is the
## matching vector of true values.
##
## Two RMSEs are reported because they answer different questions:
##   *_point  RMSE of the posterior mean across components - how well the
##            parameter is recovered.
##   *_draws  RMSE across draws and components - includes posterior spread, so
##            a diffuse but unbiased posterior is penalised.
recovery_metrics <- function(draws, truth, probs = c(0.05, 0.95)) {
  draws <- as.matrix(draws)
  truth <- as.numeric(truth)
  stopifnot(ncol(draws) == length(truth))

  posterior_mean <- colMeans(draws)
  errors <- posterior_mean - truth
  interval <- apply(draws, 2, quantile, probs = probs, names = FALSE)

  list(
    rmse_point = sqrt(mean(errors^2)),
    rmse_draws = sqrt(mean(sweep(draws, 2, truth, "-")^2)),
    bias = mean(errors),
    mae = mean(abs(errors)),
    coverage = mean(truth >= interval[1, ] & truth <= interval[2, ])
  )
}

## Pooled pointwise held-out log predictive density.
##
## log_lik_test arrives as draws x test points with the chains already pooled
## by posterior::as_draws_matrix, so the pointwise elpd is a plain log-mean-exp
## down each column. Kept as a named function because the paired ELPD
## differences at aggregation need the pointwise vector, not just its sum.
pointwise_elpd <- function(log_lik_test) {
  log_lik_test <- as.matrix(log_lik_test)
  n_draws <- nrow(log_lik_test)
  apply(log_lik_test, 2, log_sum_exp) - log(n_draws)
}

## Collapse a posterior::summarise_draws() frame over a block of parameters
## (e.g. all D entries of alpha) into the worst-case diagnostics.
block_diagnostics <- function(summary, prefix) {
  rows <- summary[startsWith(summary$variable, prefix), , drop = FALSE]
  if (nrow(rows) == 0L) {
    return(list(max_rhat = NA_real_, min_ess_bulk = NA_real_,
                min_ess_tail = NA_real_))
  }
  list(
    max_rhat = max(rows$rhat, na.rm = TRUE),
    min_ess_bulk = min(rows$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(rows$ess_tail, na.rm = TRUE)
  )
}

## Posterior summary of a scalar generated quantity.
scalar_summary <- function(draws, probs = c(0.05, 0.95)) {
  draws <- as.numeric(draws)
  quantiles <- quantile(draws, probs = probs, names = FALSE)
  list(
    mean = mean(draws), median = median(draws),
    q05 = quantiles[1], q95 = quantiles[2]
  )
}

## A crashed or infeasible fit degrades to a stub with the same component
## shapes, so that aggregation and plotting still run and the cell renders as
## NA rather than silently disappearing from the grid.
failed_result <- function(message, condition, model, replication, n_test) {
  list(
    meta = list(
      condition_id = condition$id, dimension = condition$dimension,
      rho = condition$rho, kernel = condition$kernel,
      boundary = condition$boundary, n_basis = condition$n_basis,
      model_id = model$id, model_label = model$label,
      alpha_prior = model$alpha_prior, rho_prior = model$rho_prior,
      replication = replication, adapt_delta = NA_real_,
      ok = FALSE, error = message
    ),
    metrics = list(
      rmse_alpha_point = NA_real_, rmse_alpha_draws = NA_real_,
      bias_alpha = NA_real_, coverage_alpha = NA_real_,
      rmse_rho_point = NA_real_, rmse_rho_draws = NA_real_,
      bias_rho = NA_real_, coverage_rho = NA_real_,
      rmse_sigma_point = NA_real_, rmse_sigma_draws = NA_real_,
      bias_sigma = NA_real_, coverage_sigma = NA_real_,
      r2_mean = NA_real_, r2_median = NA_real_, r2_q05 = NA_real_,
      r2_q95 = NA_real_, r2_oos_mean = NA_real_, r2_oos_median = NA_real_,
      r2_oos_q05 = NA_real_, r2_oos_q95 = NA_real_,
      true_r2 = NA_real_, elpd_test = NA_real_, mlpd_test = NA_real_,
      se_elpd_test = NA_real_, oracle_mlpd_test = NA_real_
    ),
    diagnostics = list(
      max_rhat = NA_real_, min_ess_bulk = NA_real_,
      max_rhat_alpha = NA_real_, min_ess_bulk_alpha = NA_real_,
      max_rhat_rho = NA_real_, min_ess_bulk_rho = NA_real_,
      divergences = NA_integer_, max_treedepth_hits = NA_integer_,
      sampling_seconds = NA_real_
    ),
    ## Length-matched so rbind-style aggregation of the pointwise matrix works.
    pointwise_elpd_test = rep(NA_real_, n_test)
  )
}
