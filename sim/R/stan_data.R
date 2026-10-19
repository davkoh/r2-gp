## Stan data assembly for ../../Stanmodels/aHSGP_oos.stan.
##
## One function so that the fitting path and the prior-predictive path cannot
## drift apart, which is the same reason ../../aHSGP_oos_common.R exists.

source("R/conditions.R")
source("R/scaling.R")

## Prior models. The integer is the prior_type passed to Stan.
PRIOR_TYPES <- c(R2D2 = 1L, Independent = 2L, RHS = 3L, NG = 4L)

## Length-scale priors. The integer is rho_prior_type.
##   ig     inverse-gamma containment prior, the existing default
##   vague  lognormal(0, rho_vague_sd), the control used when studying the
##          alpha-rho ridge
##   normal normal(mu, sd) truncated at zero
RHO_PRIOR_TYPES <- c(ig = 1L, vague = 2L, normal = 3L)

make_stan_data <- function(dataset, condition, model, config, prior_only = FALSE) {
  design <- config$design
  priors <- config$priors
  dimension <- condition$dimension

  if (!model$alpha_prior %in% names(PRIOR_TYPES)) {
    stop("Unknown alpha_prior: ", model$alpha_prior, call. = FALSE)
  }
  if (!model$rho_prior %in% names(RHO_PRIOR_TYPES)) {
    stop("Unknown rho_prior: ", model$rho_prior, call. = FALSE)
  }

  ## Stan expects length-D vectors even when D is 1, so every scalar
  ## hyperparameter that varies by dimension is wrapped in a 1-d array.
  dimension_vector <- function(value) array(rep(value, dimension), dim = dimension)

  ## normal_centre puts the prior mean at the generating length-scale, making
  ## it a well-specified informative prior rather than a misspecified one.
  rho_normal_mean <- if (isTRUE(model$normal_centre)) {
    condition$rho
  } else {
    model$normal_mean %||% condition$rho
  }

  ## The training design spans exactly [-1, 1] after scale_inputs(), so the
  ## midpoint Stan subtracts is zero and L = boundary for every dimension.
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

    ## Regularized horseshoe, Piironen & Vehtari (2017). p0 = D/2 is an
    ## agnostic guess at how many components are active; the DGP here is
    ## dense, so this is deliberately not tuned in the horseshoe's favour.
    hs_df = priors$hs_df,
    hs_df_global = priors$hs_df_global,
    hs_df_slab = priors$hs_df_slab,
    hs_scale_slab = priors$hs_scale_slab,
    hs_p0 = dimension / 2,

    ## Normal-gamma. E[alpha_d^2] = sigma^2 nu and E[nu] = ng_nu_scale, so this
    ## matches the R2D2 model's implied prior mean total signal variance at every
    ## D; without it the comparison would turn on a scale mismatch.
    ng_nu_scale = (priors$r2_mean / (1 - priors$r2_mean)) / dimension
  )
}
