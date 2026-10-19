## Simulation conditions, and the HSGP boundary / basis-count rules.
##
## A condition is one cell of the grid: (D, length-scale, kernel). The boundary
## factor c and the basis count M are derived from the length-scale and the
## kernel rather than hardcoded.
## Both rules come from Riutort-Mayol, Solin, Vehtari et al. (2023), and the
## kernel constants are the ones in ../../varhsgpfns.R:106.

source("R/scaling.R")

`%||%` <- function(x, y) if (is.null(x)) y else x

## Kernel identifiers. The integer is the covfn passed to Stan; see
## ../../Stanmodels/aHSGP_oos.stan. "exp" is kept so the pipeline can also
## reproduce the original exponential-kernel results.
KERNELS <- c(se = 0L, m32 = 1L, m52 = 2L, exp = 3L)

## Minimum-M constants from Riutort-Mayol et al. (2023), eqs. 18-23:
## m >= const * c / (rho / S). Smoother kernels have lighter spectral tails
## and so need fewer terms; the Materns need roughly twice what the squared
## exponential does, which is the "more basis functions for the Materns" note.
KERNEL_BASIS_CONSTANT <- c(se = 1.75, m32 = 3.42, m52 = 2.65, exp = 3.42)

KERNEL_LABELS <- c(
  se = "Squared exponential", m32 = "Matern 3/2",
  m52 = "Matern 5/2", exp = "Exponential"
)

check_kernel <- function(kernel) {
  if (!kernel %in% names(KERNELS)) {
    stop("Unknown kernel: ", kernel, ". Expected one of ",
         paste(names(KERNELS), collapse = ", "), call. = FALSE)
  }
  kernel
}

## Boundary factor. Riutort-Mayol et al. require c >= 4.5 * rho / S for the
## approximation to be accurate away from the boundary. With the inputs scaled
## to [-1, 1] the half-range S is 1, so the rule reads c >= 4.5 * rho and L = c.
## The floor keeps a sensible boundary at very short length-scales, where the
## binding constraint is M rather than c.
boundary_factor <- function(rho, c_mult = 4.5, c_min = 1.5) {
  pmax(c_min, c_mult * rho)
}

## The inverse-gamma length-scale prior, from the interval it must contain.
##
## A shape and a scale are not quantities anyone knows. An interval is: "this
## design cannot resolve a length-scale below 0.15, and a length-scale above
## 1.5 is longer than the domain". So the config states the interval, and the
## shape and scale are solved for here, such that `tail` of the prior mass
## falls below the lower limit and `tail` above the upper limit.
##
## An explicit `lengthscale_inv_gamma: [shape, scale]` in the config still
## wins, for reproducing an older run.
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

## Quantiles of the length-scale prior a model uses. Needed because the
## approximation has to hold wherever the POSTERIOR can go, which the prior
## bounds; sizing the basis from the true length-scale alone silently assumes
## the model already knows the answer.
rho_prior_quantiles <- function(model, rho_true, priors, probs = c(0.05, 0.95)) {
  switch(
    model$rho_prior,
    ## If X ~ Gamma(a, b) then b / X ~ InvGamma(a, b), so the quantiles invert.
    ig = {
      parameters <- lengthscale_inv_gamma(priors)
      parameters[2] / qgamma(1 - probs, shape = parameters[1])
    },
    vague = qlnorm(probs, meanlog = 0, sdlog = priors$rho_vague_sd),
    normal = {
      mean <- if (isTRUE(model$normal_centre)) rho_true else model$normal_mean
      sd <- model$normal_sd %||% priors$rho_normal_sd
      ## Truncated at zero, matching Stan's <lower=0> on rho_parameter.
      lower <- pnorm(0, mean, sd)
      qnorm(lower + probs * (1 - lower), mean, sd)
    },
    stop("Unknown rho_prior: ", model$rho_prior, call. = FALSE)
  )
}

## The length-scale range the basis must support for a condition, taken as the
## union over every model. All models share one boundary and one basis count so
## that the comparison between them is not confounded with differences in
## approximation quality: they must be fitted on an identical basis.
prior_lengthscale_range <- function(config, rho_true, probs = c(0.05, 0.95)) {
  quantiles <- vapply(
    config$models,
    function(model) rho_prior_quantiles(model, rho_true, config$priors, probs),
    numeric(2)
  )
  ## Include the truth: a prior concentrated away from it must still be fitted
  ## on a basis that can represent the data-generating function.
  c(lower = min(quantiles[1, ], rho_true), upper = max(quantiles[2, ], rho_true))
}

## Basis count per component. M = ceil(safety * const * c / rho).
##
## `rho` here is the SHORTEST length-scale that must be resolved, which is
## generally not the generating one: see prior_lengthscale_range(). The kernel
## constant is what makes the Materns need roughly twice the terms the squared
## exponential does at the same boundary and length-scale.
basis_count <- function(rho, kernel, boundary, safety = 1.2,
                        m_min = 8L, m_max = 128L) {
  check_kernel(kernel)
  const <- KERNEL_BASIS_CONSTANT[[kernel]]
  m <- ceiling(safety * const * boundary / rho)
  if (m > m_max) {
    ## Capping is a silent downgrade of the approximation, so say so. Either
    ## raise m_max and pay for it, or tighten the length-scale prior so the
    ## basis does not have to cover such a wide range.
    warning(sprintf(
      paste("basis_count: %s at rho = %.3f, L = %.2f needs M = %d but m_max",
            "is %d; the Riutort-Mayol rule is not satisfied."),
      kernel, rho, boundary, m, m_max
    ), call. = FALSE)
  }
  as.integer(min(max(m, m_min), m_max))
}

## Expand the config grid into a list of condition specs.
##
## grid$max_dimension_by_kernel drops the most expensive corner of the grid
## (the Materns at large D need ~2x the basis functions of the SE, so D * M
## grows fastest exactly there).
build_conditions <- function(config) {
  grid <- config$grid
  basis <- config$basis
  max_by_kernel <- grid$max_dimension_by_kernel %||% list()

  conditions <- list()
  for (kernel in grid$kernels) {
    check_kernel(kernel)
    max_dimension <- max_by_kernel[[kernel]] %||% Inf
    for (rho in grid$lengthscales) {
      for (dimension in grid$dimensions) {
        if (dimension > max_dimension) {
          next
        }
        ## size_from = "prior" (default) sizes the basis for the whole range
        ## the length-scale priors allow, so the approximation stays valid
        ## wherever the posterior goes. "truth" sizes it at the generating
        ## length-scale only, which is cheaper but assumes the model already
        ## knows the answer; basis_adequacy will flag the difference.
        range <- if (identical(basis$size_from %||% "prior", "truth")) {
          c(lower = rho, upper = rho)
        } else {
          prior_lengthscale_range(config, rho)
        }
        boundary <- basis$boundary %||% boundary_factor(
          range[["upper"]],
          c_mult = basis$boundary_multiplier %||% 4.5,
          c_min = basis$boundary_min %||% 1.5
        )
        ## The basis count is set by the SHORTEST length-scale that has to be
        ## resolved: short length-scales need high frequencies.
        n_basis <- basis$n_basis %||% basis_count(
          range[["lower"]], kernel, boundary,
          safety = basis$basis_safety %||% 1.2,
          m_min = as.integer(basis$m_min %||% 8L),
          m_max = as.integer(basis$m_max %||% 128L)
        )
        conditions[[length(conditions) + 1L]] <- list(
          id = condition_id(dimension, rho, kernel),
          dimension = as.integer(dimension),
          rho = as.numeric(rho),
          kernel = kernel,
          covfn = KERNELS[[kernel]],
          boundary = as.numeric(boundary),
          n_basis = as.integer(n_basis),
          rho_lower = as.numeric(range[["lower"]]),
          rho_upper = as.numeric(range[["upper"]]),
          ## The data-generating basis is deliberately finer than the fitting
          ## basis, so that the truth is close to an exact GP under the same
          ## normalisation convention and an under-resolved fit is detectable.
          n_basis_gen = as.integer(
            n_basis * (basis$basis_gen_multiplier %||% 4L)
          )
        )
      }
    }
  }
  conditions
}

## Format a number for a directory name: 0.5 -> 0p5, -2 -> m2.
num_token <- function(value) {
  text <- if (is.numeric(value)) formatC(value, format = "g") else as.character(value)
  gsub("-", "m", gsub(".", "p", text, fixed = TRUE), fixed = TRUE)
}

## Path-safe identifier for one cell of the grid. Mirrored by condition_id()
## in the Snakefile, which is the single place that composes paths.
condition_id <- function(dimension, rho, kernel) {
  sprintf("D%s_rho%s_%s", num_token(dimension), num_token(rho), kernel)
}

find_condition <- function(conditions, id) {
  match <- Filter(function(x) identical(x$id, id), conditions)
  if (length(match) == 0L) {
    stop("Unknown condition id: ", id, call. = FALSE)
  }
  match[[1]]
}

find_model <- function(config, id) {
  match <- Filter(function(x) identical(x$id, id), config$models)
  if (length(match) == 0L) {
    stop("Unknown model id: ", id, call. = FALSE)
  }
  match[[1]]
}
