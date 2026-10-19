## Shared plotting style for every figure in the project.
##
## Before this file there were two incompatible palettes in use: Okabe-Ito in
## the repo-root scripts (R2D2 = #0072B2) and the ggplot2 default hue wheel in
## paper/figures/fig03 and fig04 (R2D2 = #00BFC4). Okabe-Ito wins because it is
## colour-blind safe and already dominant. Source this file from any script
## that draws, including the ones under paper/figures/.

library(ggplot2)

## Okabe-Ito. A model is an (alpha prior x rho prior) pair, so rather than
## giving four models four arbitrary colours, the two factors get two aesthetics:
## colour carries the alpha prior, which is what the paper is about, and
## linetype carries the length-scale prior.
alpha_prior_colours <- c(
  R2D2 = "#0072B2", Independent = "#D55E00", RHS = "#009E73",
  NG = "#CC79A7", Oracle = "grey35"
)
alpha_prior_levels <- c("R2D2", "RHS", "NG", "Independent", "Oracle")

rho_prior_linetypes <- c(ig = "solid", normal = "22", vague = "42")
## Kept generic: the inverse-gamma's parameters are a config setting, so
## hardcoding them here would go stale the moment the prior is retuned.
rho_prior_labels <- c(
  ig = "Inverse-gamma", normal = "Normal", vague = "Lognormal (vague)"
)

## Ablation models in the identifiability study are variants of the two main
## priors, so they shade within each family rather than taking unrelated hues:
## R2D2 stays in the blue family, Independent in the orange one.
ablation_colours <- c(
  "R2D2, containment rho" = "#0072B2",
  "R2D2, vague rho" = "#56B4E9",
  "Independent, vague rho, normalised" = "#D55E00",
  "Independent, vague rho, unnormalised" = "#E69F00"
)

kernel_colours <- c(se = "#0072B2", m32 = "#E69F00", m52 = "#009E73",
                    exp = "#CC79A7")
kernel_labels <- c(
  se = "Squared exponential", m32 = "Matern 3/2",
  m52 = "Matern 5/2", exp = "Exponential"
)

## Sequential, because the length-scales are ordered.
lengthscale_colours <- c(`0.3` = "#56B4E9", `0.5` = "#0072B2", `0.7` = "#002F5B")

## Filled marks a comparison in which every model involved converged; open marks
## one that did not. This is the existing convention in aHSGP_oos_plots.R.
converged_shapes <- c(`TRUE` = 19, `FALSE` = 21)

## Fits whose worst Rhat exceeds this are not treated as reliable. The figures
## say so in the caption rather than dropping the points.
RHAT_THRESHOLD <- 1.05

order_alpha_priors <- function(x) {
  factor(x, levels = alpha_prior_levels[alpha_prior_levels %in% unique(x)])
}

order_kernels <- function(x) {
  factor(x, levels = names(kernel_labels)[names(kernel_labels) %in% unique(x)],
         labels = kernel_labels[names(kernel_labels) %in% unique(x)])
}

## The manuscript's figure style, matched exactly to paper/figures/fig03 and
## fig04: minimal, no panel grid (fig04) or horizontal rules only (fig03), with
## drawn axis lines and ticks in grey30. Every pipeline figure destined for the
## paper uses this rather than theme_hsgp(), so the two cannot drift apart.
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
save_paper_plot <- function(path, plot, width, height) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggsave(path, plot = plot, width = width, height = height, device = "pdf")
  invisible(path)
}

theme_hsgp <- function(base_size = 12) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = base_size + 2),
      plot.title = element_text(size = base_size + 4, face = "bold"),
      plot.subtitle = element_text(colour = "grey35"),
      plot.caption = element_text(hjust = 0, size = base_size - 3,
                                  colour = "grey25"),
      strip.text = element_text(face = "bold"),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

## The two scales used together wherever models are compared.
scale_colour_alpha_prior <- function(name = "Prior on alpha", ...) {
  scale_colour_manual(values = alpha_prior_colours, name = name, drop = TRUE, ...)
}

scale_fill_alpha_prior <- function(...) {
  scale_fill_manual(values = alpha_prior_colours, guide = "none", drop = TRUE, ...)
}

scale_linetype_rho_prior <- function(name = "Prior on rho", ...) {
  scale_linetype_manual(values = rho_prior_linetypes, labels = rho_prior_labels,
                        name = name, drop = TRUE, ...)
}

scale_shape_converged <- function(name = NULL) {
  scale_shape_manual(
    values = converged_shapes, name = name,
    labels = c(`TRUE` = sprintf("Rhat <= %.2f", RHAT_THRESHOLD),
               `FALSE` = sprintf("Rhat > %.2f", RHAT_THRESHOLD))
  )
}

save_hsgp_plot <- function(path, plot, width = 9, height = 6) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggsave(path, plot = plot, width = width, height = height, dpi = 200,
         bg = "white")
  invisible(path)
}
