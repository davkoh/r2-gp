## Is the boundary large enough and the basis fine enough, for every condition?
##
## This is the cheap check that answers the L question in
## ../../Meeting_notes.md, and it should be run before committing to a long
## cluster run: it costs seconds, whereas discovering an undersized boundary
## after 1920 fits costs days.
##
## Two rules are checked, both from Riutort-Mayol et al. (2023), and crucially
## they are checked at the PRIOR quantiles of rho as well as at the truth. The
## model does not know the true length-scale; if the posterior can wander to
## rho values the boundary cannot support, the approximation degrades there.
##
##   boundary rule   c >= 4.5 * rho / S,   S = 1 after scaling
##   basis rule      M >= const * c / rho
##
## On top of the rules it measures the thing the rules are proxies for: the
## largest absolute discrepancy between the HSGP kernel and the exact kernel
## over the design.

source("R/dgp.R")
source("R/theme_hsgp.R")

library(dplyr)

config <- yaml::read_yaml(snakemake@params[["config_file"]])
conditions <- build_conditions(config)
priors <- config$priors

## The HSGP approximation k(x, x') = sum_m S(omega_m) phi_m(x) phi_m(x').
approx_kernel_matrix <- function(kernel, x, rho, boundary, n_basis) {
  basis <- phi_matrix(x, n_basis, boundary)
  omega <- basis_frequencies(n_basis, boundary)
  density <- spectral_density(kernel, rho, omega)
  basis %*% diag(density, nrow = n_basis) %*% t(basis)
}

## Largest absolute error over a dense grid on the scaled design [-1, 1].
kernel_error <- function(kernel, rho, boundary, n_basis, n_grid = 60L) {
  x <- seq(-1, 1, length.out = n_grid)
  approx <- approx_kernel_matrix(kernel, x, rho, boundary, n_basis)
  exact <- outer(x, x, function(a, b) exact_kernel(kernel, a - b, rho))
  max(abs(approx - exact))
}

## Prior quantiles of rho, from conditions.R so that this check and the sizing
## rule cannot disagree about what the prior says.
rho_quantiles <- function(model, condition, probs = c(0.05, 0.5, 0.95)) {
  rho_prior_quantiles(model, condition$rho, priors, probs)
}

rows <- list()
for (condition in conditions) {
  const <- KERNEL_BASIS_CONSTANT[[condition$kernel]]
  for (model in config$models) {
    quantiles <- rho_quantiles(model, condition)
    names(quantiles) <- c("q05", "q50", "q95")
    for (label in names(quantiles)) {
      rho <- quantiles[[label]]
      rows[[length(rows) + 1L]] <- data.frame(
        condition_id = condition$id, dimension = condition$dimension,
        rho_true = condition$rho, kernel = condition$kernel,
        model_id = model$id, rho_prior = model$rho_prior,
        quantile = label, rho = rho,
        boundary = condition$boundary, n_basis = condition$n_basis,
        required_boundary = 4.5 * rho,
        required_basis = ceiling(const * condition$boundary / rho),
        boundary_ok = condition$boundary >= 4.5 * rho,
        basis_ok = condition$n_basis >= ceiling(const * condition$boundary / rho),
        max_kernel_error = kernel_error(
          condition$kernel, rho, condition$boundary, condition$n_basis
        ),
        stringsAsFactors = FALSE
      )
    }
  }
}
adequacy <- do.call(rbind, rows)

dir.create(dirname(snakemake@output[["csv"]]), recursive = TRUE,
           showWarnings = FALSE)
write.csv(adequacy, snakemake@output[["csv"]], row.names = FALSE)

## A condition passes only if it holds across the prior, not just at the truth.
flagged <- adequacy %>% filter(!boundary_ok | !basis_ok)
if (nrow(flagged) > 0) {
  message(sprintf(
    "basis_adequacy: %d of %d (condition, model, quantile) cells violate a rule.",
    nrow(flagged), nrow(adequacy)
  ))
  print(utils::head(
    flagged[c("condition_id", "model_id", "quantile", "rho", "boundary",
              "required_boundary", "n_basis", "required_basis",
              "max_kernel_error")],
    20
  ))
} else {
  message("basis_adequacy: every condition satisfies both rules across the ",
          "5th-95th prior percentile range of rho.")
}

## The figure shows the quantity the rules stand in for, so an amber cell can
## be judged on its actual approximation error rather than on the rule alone.
adequacy_plot <- adequacy %>%
  mutate(kernel_label = order_kernels(kernel),
         rule = if_else(boundary_ok & basis_ok, "satisfied", "violated")) %>%
  ggplot(aes(x = rho, y = max_kernel_error, colour = kernel_label,
             shape = rule)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey40") +
  geom_point(alpha = 0.8, size = 2) +
  facet_wrap(~quantile, nrow = 1) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  scale_colour_manual(values = setNames(kernel_colours[names(kernel_labels)],
                                        kernel_labels),
                      name = "Kernel") +
  scale_shape_manual(values = c(satisfied = 19, violated = 4),
                     name = "Riutort-Mayol rules") +
  theme_hsgp() +
  labs(
    x = "Length-scale at that prior quantile (log scale)",
    y = "Max |HSGP kernel - exact kernel| (log scale)",
    title = "HSGP approximation quality across the prior range of rho",
    subtitle = paste("Each point is one (condition, model) at the 5th, 50th and",
                     "95th prior percentile of rho.\nDashed line at an",
                     "absolute kernel error of 0.05."),
    caption = paste("Crosses violate c >= 4.5 rho or M >= const c / rho.",
                    "A cross with a small error is tolerable; a cross with a",
                    "large one is not.")
  )
save_hsgp_plot(snakemake@output[["figure"]], adequacy_plot, width = 11,
               height = 5)
