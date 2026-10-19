## Data-generating process for the additive-HSGP simulation.
##
## Ported from ../../aHSGP_oos_common.R::simulate_dataset() with three changes:
##
##   1. the spectral density follows the condition's kernel, mirroring spd()
##      in ../../Stanmodels/aHSGP_oos.stan formula for formula;
##   2. generation uses a finer basis than fitting (n_basis_gen = 4 * n_basis),
##      so the truth is effectively an exact GP under the model's own
##      normalisation convention and a fit with too few basis functions is
##      detectable rather than exactly correct by construction; and
##   3. the inputs pass through scale_inputs(), so the boundary is always
##      derived from a unit half-range.
##
## Everything else is kept from the original: equal variance allocation across
## components, spectral coefficients shared between the training and test
## designs, the realized (not nominal) true R2, and the oracle log score.

source("R/conditions.R")
source("R/scaling.R")

## HSGP eigenfunctions on [-L, L]. Matches phi() in aHSGP_oos.stan.
phi_matrix <- function(x, n_basis, boundary) {
  vapply(
    seq_len(n_basis),
    function(m) sin(m * pi * (x + boundary) / (2 * boundary)) / sqrt(boundary),
    numeric(length(x))
  )
}

## Spectral densities with alpha^2 factored out. These must stay numerically
## identical to spd_se / spd_m32 / spd_m52 / spd_exponential in
## ../../Stanmodels/aHSGP_oos.stan; tests/test_spectral_density.R checks that.
spectral_density <- function(kernel, rho, omega) {
  check_kernel(kernel)
  switch(
    kernel,
    se = sqrt(2 * pi) * rho * exp(-0.5 * (rho * omega)^2),
    m32 = (4 * 3^1.5 / rho^3) * (3 / rho^2 + omega^2)^(-2),
    m52 = (4 * 5^2.5 / (0.75 * rho^5)) * (5 / rho^2 + omega^2)^(-3),
    exp = (2 * rho) / (1 + (rho * omega)^2)
  )
}

## The exact stationary kernels corresponding to those spectral densities, with
## unit marginal variance. Kept here so the exact/approximate pair cannot drift:
## basis_adequacy.R measures the gap between them, and paper_figures.R draws
## from the exact process.
exact_kernel <- function(kernel, r, rho) {
  check_kernel(kernel)
  r <- abs(r)
  switch(
    kernel,
    se = exp(-0.5 * (r / rho)^2),
    m32 = (1 + sqrt(3) * r / rho) * exp(-sqrt(3) * r / rho),
    m52 = (1 + sqrt(5) * r / rho + 5 * r^2 / (3 * rho^2)) *
      exp(-sqrt(5) * r / rho),
    exp = exp(-r / rho)
  )
}

## Angular frequencies of the first n_basis eigenfunctions on [-L, L].
basis_frequencies <- function(n_basis, boundary) {
  seq_len(n_basis) * pi / (2 * boundary)
}

## Draw one replicate for a condition.
##
## The returned list carries the true alpha, rho and sigma so that fit_model.R
## can compute recovery metrics without re-deriving the DGP.
simulate_dataset <- function(condition, design, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }

  dimension <- condition$dimension
  n_obs <- as.integer(design$n_obs)
  n_test <- as.integer(design$n_test)

  draw_inputs <- function(n) {
    matrix(
      runif(n * dimension, min = design$x_min, max = design$x_max),
      nrow = n
    )
  }
  ## Drawn on the raw axis and then scaled, so the transform is exercised on
  ## every run rather than bypassed because the raw range happens to be [-1, 1].
  scaling <- scale_inputs(draw_inputs(n_obs), draw_inputs(n_test))
  inputs <- scaling$train
  inputs_test <- scaling$test

  n_basis_gen <- condition$n_basis_gen
  boundary <- condition$boundary
  omega <- basis_frequencies(n_basis_gen, boundary)
  density <- spectral_density(condition$kernel, condition$rho, omega)

  ## Equal allocation: every component carries the same share of the target
  ## signal variance, so alpha_d^2 = sigma^2 R2/(1-R2) / D.
  global_variance <- design$sigma^2 * design$target_r2 / (1 - design$target_r2)
  alpha <- sqrt(global_variance / dimension)

  components <- matrix(0, nrow = n_obs, ncol = dimension)
  components_test <- matrix(0, nrow = n_test, ncol = dimension)

  for (d in seq_len(dimension)) {
    basis <- phi_matrix(inputs[, d], n_basis_gen, boundary)
    basis_test <- phi_matrix(inputs_test[, d], n_basis_gen, boundary)
    ## The normalizer is a training-design quantity, exactly as in Stan, so
    ## alpha_d^2 is the expected sample variance of component d in-sample.
    normalizer <- sum(density * apply(basis, 2, var))
    coefficients <- alpha * sqrt(density / normalizer) * rnorm(n_basis_gen)
    components[, d] <- basis %*% coefficients
    components_test[, d] <- basis_test %*% coefficients
  }

  latent_mean <- rowSums(components)
  latent_mean_test <- rowSums(components_test)
  y <- rnorm(n_obs, mean = latent_mean, sd = design$sigma)
  y_test <- rnorm(n_test, mean = latent_mean_test, sd = design$sigma)
  realized_signal_variance <- var(latent_mean)

  list(
    condition_id = condition$id,
    inputs = inputs,
    y = y,
    inputs_test = inputs_test,
    y_test = y_test,
    latent_mean = latent_mean,
    latent_mean_test = latent_mean_test,
    scaling = scaling,
    ## Truth, in the parameterisation the model uses: alpha_d is the expected
    ## in-sample sample standard deviation of component d.
    alpha_true = rep(alpha, dimension),
    rho_true = rep(condition$rho, dimension),
    sigma_true = design$sigma,
    ## Realized rather than nominal: each replicate uses one draw of the
    ## spectral coefficients, so the achieved R2 scatters around the target.
    true_r2 = realized_signal_variance /
      (realized_signal_variance + design$sigma^2),
    ## The best log score any model can achieve on this test set.
    oracle_elpd_test = mean(dnorm(
      y_test, mean = latent_mean_test, sd = design$sigma, log = TRUE
    ))
  )
}
