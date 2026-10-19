## The two manuscript figures, in the paper's own style.
##
##   elpd_difference.pdf  styled exactly as paper/figures/fig04_elpd_diff.R:
##                        dotted line at zero, dashed lines at the conventional
##                        +-4 threshold, the mean difference as a solid line
##                        with points, its Monte Carlo standard error as dotted
##                        bounds, no panel grid.
##   posterior_r2.pdf     styled exactly as paper/figures/fig03_r2_recovery.R:
##                        dodged points with 90% intervals, truth as crosses,
##                        in-sample and out-of-sample panels.
##
## The one deliberate departure from fig03 is that its out-of-sample panel
## showed the MODEL-based R2 at the test design (posterior_mean_test), which is
## the mistake flagged by the \todo in paper/sections/04-simulation.tex. This
## uses R2_oos, the residual-based Bayesian R2 of Gelman et al. (2019), which is
## the genuine out-of-sample quantity and now carries an interval too.
##
## Both figures gain facets the paper versions did not need: the manuscript
## swept D alone, whereas this grid also varies the kernel and the length-scale.
## Everything else - geometry, theme, sizes, colours - is unchanged.
##
## Parameter recovery and sampling statistics are NOT plotted; they go to the
## LaTeX table produced by R/make_table.R.

library(dplyr)
library(tidyr)

source("R/theme_hsgp.R")

results <- readRDS(snakemake@input[["rds"]])
metrics <- results$metrics
diagnostics <- results$diagnostics
elpd <- results$elpd
config <- results$config

figures <- snakemake@params[["figure_dir"]]
dir.create(figures, recursive = TRUE, showWarnings = FALSE)
figure <- function(name) file.path(figures, name)

dimensions <- sort(unique(metrics$dimension))
single_kernel <- length(unique(metrics$kernel)) == 1L
single_rho <- length(unique(metrics$rho)) == 1L

## Facet only over the factors that actually vary, so a single-cell run renders
## as the paper's own single panel rather than a one-panel facet grid. The
## length-scale strips are plotmath, hence the labeller.
rho_labeller <- labeller(rho_label = label_parsed)

add_facets <- function(plot) {
  if (single_kernel && single_rho) {
    plot
  } else if (single_rho) {
    plot + facet_wrap(~kernel_label, nrow = 1)
  } else if (single_kernel) {
    plot + facet_wrap(~rho_label, nrow = 1, labeller = rho_labeller)
  } else {
    plot + facet_grid(kernel_label ~ rho_label, labeller = rho_labeller)
  }
}

label_factors <- function(data) {
  data %>% mutate(
    kernel_label = order_kernels(kernel),
    rho_label = factor(sprintf("rho == %.1f", rho),
                       levels = sprintf("rho == %.1f", sort(unique(rho))))
  )
}

## The paper labels lines by the prior on alpha. Where a second length-scale
## prior is also in play it becomes the linetype, which keeps the single-factor
## case identical to fig04.
single_rho_prior <- length(unique(metrics$rho_prior)) == 1L

## ---------------------------------------------------------------------------
## Figure 4 style: ELPD difference against the reference model
## ---------------------------------------------------------------------------

reference_label <- Filter(
  function(a) identical(a$id, config$reference_model), config$models
)[[1]]$label

## Mean over replications with the Monte Carlo standard error of that mean,
## exactly as fig04_elpd_diff.R computes them.
elpd_summary <- elpd %>%
  filter(!is.na(elpd_difference)) %>%
  group_by(dimension, kernel, rho, model_id, alpha_prior, rho_prior) %>%
  summarise(
    elpd_diff = mean(elpd_difference),
    se = sd(elpd_difference) / sqrt(n()),
    n_rep = n(),
    .groups = "drop"
  ) %>%
  mutate(alpha_prior = order_alpha_priors(alpha_prior)) %>%
  label_factors()

## fig04 fixes the window to include the +-4 threshold lines so that the
## conventional band is always visible, however small the differences are.
y_lo <- min(-4, min(elpd_summary$elpd_diff - elpd_summary$se, na.rm = TRUE))
y_hi <- max(4, max(elpd_summary$elpd_diff + elpd_summary$se, na.rm = TRUE))

elpd_plot <- ggplot(
  elpd_summary,
  aes(dimension, elpd_diff, colour = alpha_prior, group = model_id)
) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey45",
             linewidth = 0.3) +
  geom_hline(yintercept = c(-4, 4), linetype = "dashed", colour = "grey45",
             linewidth = 0.3) +
  geom_line(aes(y = elpd_diff + se), linetype = "dotted", linewidth = 0.55) +
  geom_line(aes(y = elpd_diff - se), linetype = "dotted", linewidth = 0.55) +
  {
    if (single_rho_prior) {
      geom_line(linewidth = 0.75)
    } else {
      geom_line(aes(linetype = rho_prior), linewidth = 0.75)
    }
  } +
  geom_point(size = 2) +
  scale_x_continuous(breaks = dimensions) +
  scale_y_continuous(labels = function(x) sprintf("%.0f", x)) +
  scale_colour_manual(values = alpha_prior_colours, name = NULL) +
  coord_cartesian(ylim = c(y_lo, y_hi)) +
  labs(x = "Number of additive components  D",
       y = sprintf("ELPD - ELPD(%s)", reference_label)) +
  theme_paper()
if (!single_rho_prior) {
  elpd_plot <- elpd_plot + scale_linetype_rho_prior(name = NULL)
}
elpd_plot <- add_facets(elpd_plot)

## Panel size follows the paper's 6.6 x 4.0 for a single panel and grows with
## the facet grid rather than squeezing it.
n_col <- if (single_rho) 1L else length(unique(elpd_summary$rho))
n_row <- if (single_kernel) 1L else length(unique(elpd_summary$kernel))
save_paper_plot(figure("elpd_difference.pdf"), elpd_plot,
                width = 6.6 * max(1, n_col * 0.72),
                height = 4.0 * max(1, n_row * 0.68))
write.csv(elpd_summary, file.path(figures, "elpd_difference.csv"),
          row.names = FALSE)

## ---------------------------------------------------------------------------
## Figure 3 style: posterior R2 recovery, in sample and out of sample
## ---------------------------------------------------------------------------

r2_panels <- bind_rows(
  metrics %>% transmute(
    dimension, kernel, rho, alpha_prior, rho_prior, model_id, replication,
    truth = true_r2, estimate = r2_mean, lo = r2_q05, hi = r2_q95,
    panel = "In sample"
  ),
  ## Residual-based, so a model that predicts badly is penalised. fig03 used
  ## the model-based quantity here, which is the error its \todo refers to.
  metrics %>% transmute(
    dimension, kernel, rho, alpha_prior, rho_prior, model_id, replication,
    truth = true_r2, estimate = r2_oos_mean, lo = r2_oos_q05, hi = r2_oos_q95,
    panel = "Out of sample"
  )
) %>%
  filter(!is.na(estimate)) %>%
  mutate(
    panel = factor(panel, levels = c("In sample", "Out of sample")),
    alpha_prior = order_alpha_priors(alpha_prior),
    D = factor(dimension, levels = dimensions)
  ) %>%
  label_factors()

dodge <- position_dodge(width = 0.6)

r2_plot <- ggplot(r2_panels, aes(D, estimate, colour = alpha_prior,
                                 group = model_id)) +
  geom_linerange(aes(ymin = lo, ymax = hi), position = dodge,
                 linewidth = 0.35, na.rm = TRUE, alpha = 0.7) +
  geom_point(position = dodge, size = 1.5, alpha = 0.85) +
  geom_point(
    data = r2_panels %>% distinct(panel, kernel_label, rho_label, D,
                                  replication, truth),
    aes(x = D, y = truth), inherit.aes = FALSE,
    colour = "grey20", shape = 4, size = 1.9, stroke = 0.6
  ) +
  scale_colour_manual(values = alpha_prior_colours, name = NULL) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Number of additive components  D",
       y = expression(R[emp]^2)) +
  theme_paper(grid_y = TRUE)

## The paper facets on the panel alone; here the panel is crossed with whichever
## grid factors vary.
if (single_kernel && single_rho) {
  r2_plot <- r2_plot + facet_wrap(~panel)
} else if (single_rho) {
  r2_plot <- r2_plot + facet_grid(kernel_label ~ panel)
} else {
  r2_plot <- r2_plot + facet_grid(
    kernel_label + rho_label ~ panel,
    labeller = labeller(rho_label = label_parsed)
  )
}

n_row_r2 <- max(1L, length(unique(r2_panels$kernel)) *
                  length(unique(r2_panels$rho)))
save_paper_plot(figure("posterior_r2.pdf"), r2_plot,
                width = 8.0, height = max(4.2, 2.0 * n_row_r2))

message(sprintf("plot_results: wrote elpd_difference.pdf and posterior_r2.pdf to %s",
                figures))
