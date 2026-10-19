## Final figures and table of the BloodBrain case study.
##
##   final/effects_top2.pdf       the 2 strongest descriptors (by SE R2-GP),
##                                all priors overlaid; SE top row, Matern 3/2
##                                bottom row
##   final/performance_table.tex  mean held-out ELPD and posterior R2 in and
##                                out of sample per prior; SE and Matern 3/2
##                                as column groups
##   final/se/effects_top12.pdf   the 12 strongest descriptors, all priors
##   final/m32/effects_top12.pdf  the same 12 descriptors under Matern 3/2
##   final/{se,m32}/variance_explained.pdf  each descriptor's share of the
##                                fitted variance, all priors
##
## Nothing is fitted. The table comes from the main case study's split fits
## (vardecomp_casestudy/bloodbrain/fits/); the effects figures from the
## summaries that case_study_bloodbrain_effects.R writes per kernel
## (vardecomp_casestudy/bloodbrain_effects/<kernel>/summaries/), which hold the
## priors listed in that script's alpha_priors. A kernel whose summaries are
## missing is skipped with a message saying how to make them.
##
##   Rscript "Case study/final_outputs.R"


## ---------------------------------------------------------------------------
## The case study's settings and data, without running it
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
## Now available: settings, x, covariates, logBBB, bbbDescr and the functions
## in vardecomp_functions.R.

suppressPackageStartupMessages(library(dplyr))


## ---------------------------------------------------------------------------
## Settings
## ---------------------------------------------------------------------------

final_settings <- list(
  output_dir = file.path(settings$output_root, "final"),
  main_dir = file.path(settings$output_root, "bloodbrain"),
  effects_dir = file.path(settings$output_root, "bloodbrain_effects"),
  ## Which split fits to use, by the length-scale prior part of their names.
  rho_prior_tag = "normal0.5-0.1",
  rho_prior_text = "a normal(0.5, 0.1) prior on the length-scales",
  ## Kernels, in this order (rows of the top-2 figure, column groups of the
  ## table), with their labels.
  kernels = c(se = "SE", m32 = "Matérn 3/2"),
  kernel_tex = c(se = "SE", m32 = "Mat\\'ern 3/2"),
  kernel_text = c(se = "squared exponential (SE)", m32 = "Mat\\'ern 3/2"),
  data_text = "the BloodBrain data ($N = 208$, $D = 101$)",
  ## Descriptors are ranked by this prior's mean share under SE; the same
  ## descriptors are shown for both kernels.
  rank_prior = "R2D2",
  rank_kernel = "se",
  n_top2 = 2,
  n_top12 = 12,
  n_variance = 25,
  prior_order = c("R2D2", "RHS", "NG", "Independent"),
  prior_labels = c(R2D2 = "R2-GP", RHS = "RHS", NG = "NG",
                   Independent = "Normal"),
  reference_prior = "R2D2",
  font_size = 11,
  ## One column of a two-column page: set the width to the template's
  ## \columnwidth in inches (\the\columnwidth prints it in pt; / 72.27).
  top2_size = c(width = 3.2, height = 3.6),
  top12_size = c(width = 6.5, height = 5.5),
  variance_size = c(width = 6.5, height = 5.0),
  compile_table_preview = TRUE,
  ## Also write effects_top2.tex, the top-2 figure as TikZ code (tikzDevice).
  tikz = TRUE,
  ## The standard PDF device (Cairo needs XQuartz on macOS); a UTF-8 locale is
  ## set below so that it can encode the accent in "Matern".
  device = "pdf"
)
## Rscript often runs in the "C" locale, where the PDF device cannot encode
## non-ASCII labels; a UTF-8 locale lets it map them to its Latin-1 fonts.
if (!isTRUE(l10n_info()$`UTF-8`)) {
  invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", "en_US.UTF-8")))
}
final_path <- function(...) file.path(final_settings$output_dir, ...)
plot_theme <- theme_case_study(list(font_size = final_settings$font_size))
kernel_ids <- names(final_settings$kernels)

## Descriptors on their original scale, for the rug of observed values.
x_original <- as.matrix(bbbDescr[covariates])
for (column in intersect(settings$log_covariates, covariates)) {
  x_original[, column] <- log1p(x_original[, column])
}

## ---------------------------------------------------------------------------
## Pieces
## ---------------------------------------------------------------------------

effects_summaries <- function(kernel) {
  directory <- file.path(final_settings$effects_dir, kernel, "summaries")
  files <- file.path(directory, c("effect_curves.csv",
                                  "covariate_relevance.csv"))
  if (!all(file.exists(files))) {
    message(sprintf(
      paste("No effects summaries for %s in %s; skipping its effects.",
            "Make them with:\n  EFFECTS_KERNEL=%s Rscript",
            "\"Case study/case_study_bloodbrain_effects.R\""),
      kernel, directory, kernel
    ))
    return(NULL)
  }
  summaries <- list(curves = read.csv(files[1]),
                    relevance = read.csv(files[2]))
  message(kernel, " effects: priors ",
          paste(intersect(final_settings$prior_order,
                          unique(summaries$curves$prior)), collapse = ", "))
  summaries
}

prior_factor <- function(prior) {
  present <- intersect(final_settings$prior_order, unique(prior))
  factor(final_settings$prior_labels[prior],
         levels = final_settings$prior_labels[present])
}
prior_colours <- setNames(alpha_prior_colours[final_settings$prior_order],
                          final_settings$prior_labels[final_settings$prior_order])

## Effect curves (posterior mean, pointwise 90% band) of every prior, one
## panel per descriptor (and per kernel, if `by_kernel`), each panel on its
## own y range.
## for_tikz = TRUE writes every label as LaTeX (for tikzDevice): the y-axis as
## LaTeX math and underscores in descriptor names escaped.
## short_y = TRUE drops "Effect on logBBB, " from the y-axis title, for small
## figures (say so in the caption instead).
plot_curves <- function(curves, shown, by_kernel = FALSE, ncol = 4,
                        for_tikz = FALSE, short_y = FALSE) {
  covariate_label <- if (for_tikz) {
    function(v) gsub("_", "\\_", v, fixed = TRUE)   # one backslash: \_
  } else {
    identity
  }
  data <- curves %>%
    filter(covariate %in% shown) %>%
    mutate(covariate = factor(covariate, levels = shown),
           prior = prior_factor(prior))
  rug <- data.frame(covariate = factor(rep(shown, each = nrow(x_original)),
                                       levels = shown),
                    x = as.vector(x_original[, shown]))
  if (by_kernel) {
    rug <- merge(rug, data.frame(kernel = unique(data$kernel)))
  }
  facets <- if (by_kernel) {
    facet_grid(kernel ~ covariate, scales = "free",
               labeller = labeller(covariate = as_labeller(covariate_label)))
  } else {
    facet_wrap(~covariate, ncol = ncol, scales = "free",
               labeller = labeller(covariate = as_labeller(covariate_label)))
  }
  plot <- ggplot(data, aes(x, mean)) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey45") +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = prior), alpha = 0.15) +
    geom_line(aes(colour = prior), linewidth = 0.6) +
    geom_rug(data = rug, aes(x = x), inherit.aes = FALSE, sides = "b",
             alpha = 0.3, length = unit(0.03, "npc")) +
    facets +
    ## About 3 ticks per panel, so that long values do not run together.
    scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
    scale_colour_manual(values = prior_colours, name = NULL) +
    scale_fill_manual(values = prior_colours, name = NULL) +
    labs(x = "Covariate value",
         y = if (for_tikz) {
           paste0(if (!short_y) "Effect on logBBB, ", "$f_d(x) - \\bar{f}_d$")
         } else if (short_y) {
           expression(f[d](x) - bar(f)[d])
         } else {
           expression("Effect on logBBB, " * f[d](x) - bar(f)[d])
         }) +
    plot_theme
  if (for_tikz) {
    ## Tick labels in math mode (a proper minus sign), and a little more room
    ## between the y-axis title and the tick labels, which LaTeX sets wider.
    math <- function(v) ifelse(is.na(v), NA, paste0("$", v, "$"))
    plot <- plot +
      scale_x_continuous(breaks = scales::breaks_pretty(n = 3), labels = math) +
      scale_y_continuous(labels = math) +
      theme(axis.title.y = element_text(margin = margin(r = 8)))
  }
  plot
}

plot_variance <- function(relevance, shown) {
  data <- relevance %>%
    filter(covariate %in% shown, prior %in% final_settings$prior_order) %>%
    mutate(covariate = factor(covariate, levels = rev(shown)),
           prior = prior_factor(prior))
  ggplot(data, aes(share_mean, covariate, colour = prior)) +
    geom_vline(xintercept = 1 / length(covariates), linetype = "dashed",
               colour = "grey45") +
    geom_errorbarh(aes(xmin = share_lo, xmax = share_hi),
                   position = position_dodge(width = 0.7), height = 0,
                   linewidth = 0.45) +
    geom_point(position = position_dodge(width = 0.7), size = 1.3) +
    scale_colour_manual(values = prior_colours, name = NULL) +
    labs(x = "Share of the variance of the fitted function", y = NULL) +
    plot_theme
}

## ---------------------------------------------------------------------------
## Table: mean ELPD and R2, and sampling diagnostics; SE and Matern 3/2 as
## column groups
## ---------------------------------------------------------------------------

dir.create(final_path(), recursive = TRUE, showWarnings = FALSE)
fits <- bind_rows(lapply(kernel_ids, function(kernel) {
  data <- read_case_study_fits(final_settings$main_dir, kernel,
                               final_settings$rho_prior_tag,
                               final_settings$reference_prior,
                               label = final_settings$kernels[[kernel]])
  data$fits
}))
## Means over the splits, and per prior the worst R-hat and smallest bulk ESS
## over all alpha_d, rho_d, sigma and lp and over splits (so one bad split is
## not averaged away), mean divergences per fit and median time per fit.
table_data <- fits %>%
  filter(alpha_prior %in% final_settings$prior_order) %>%
  group_by(kernel, alpha_prior) %>%
  summarise(elpd = mean(elpd_test), r2 = mean(r2_mean),
            r2_oos = mean(r2_oos_mean),
            max_rhat = max(max_rhat), min_ess = min(min_ess_bulk),
            divergences = mean(divergences),
            minutes = median(sampling_seconds) / 60,
            n_splits = n(), .groups = "drop")
write.csv(table_data, final_path("performance_table.csv"), row.names = FALSE)

n_columns <- 7
priors <- intersect(final_settings$prior_order, unique(table_data$alpha_prior))
cells <- function(prior, kernel) {
  row <- table_data[table_data$alpha_prior == prior &
                      table_data$kernel == kernel, ]
  if (nrow(row) == 0) return(rep("---", n_columns))
  c(sprintf("$%s$", latex_number(row$elpd, 1)), latex_number(row$r2, 2),
    latex_number(row$r2_oos, 2), latex_number(row$max_rhat, 3),
    latex_number(row$min_ess, 0), latex_number(row$divergences, 1),
    latex_number(row$minutes, 0))
}
body <- vapply(priors, function(prior) {
  paste0(paste(c(final_settings$prior_labels[[prior]],
                 unlist(lapply(kernel_ids, cells, prior = prior))),
               collapse = " & "), " \\\\")
}, character(1))
first_columns <- 2 + n_columns * (seq_along(kernel_ids) - 1)
caption <- sprintf(
  paste(
    "Predictive performance and sampling behaviour of the additive HSGP with",
    "a %s kernel on %s, with %s, over 10 random train/test splits with the",
    "same splits and sampler settings for every prior. ELPD (the held-out",
    "log predictive density summed over the test set), $R^2$ (the in-sample",
    "posterior mean) and $R^2_{\\mathrm{oos}}$ (the residual-based posterior",
    "mean on the test set) are means over the splits. $\\widehat{R}$ and",
    "ESS are the worst case over all $\\alpha_d$, $\\rho_d$, $\\sigma$ and",
    "$\\mathrm{lp}$ and over splits; Div.\\ is the mean number of divergent",
    "transitions per fit and time the median sampling time per fit in",
    "minutes."
  ),
  paste(final_settings$kernel_text[kernel_ids], collapse = " and a "),
  final_settings$data_text, final_settings$rho_prior_text
)
write_latex_table(c(
  "% Requires \\usepackage{booktabs}.",
  "% Generated by Case study/final_outputs.R -- do not edit by hand.",
  "\\begin{table}",
  "\\centering\\footnotesize",
  "\\setlength{\\tabcolsep}{3pt}",
  sprintf("\\caption{%s}\\label{tab:bloodbrain-performance}", caption),
  sprintf("\\begin{tabular}{l %s}",
          paste(rep("rcc rrrr", length(kernel_ids)), collapse = " ")),
  "\\toprule",
  paste0("& ", paste(sprintf("\\multicolumn{%d}{c}{%s}", n_columns,
                             final_settings$kernel_tex[kernel_ids]),
                     collapse = " & "), " \\\\"),
  paste(sprintf("\\cmidrule(lr){%d-%d}", first_columns,
                first_columns + n_columns - 1), collapse = " "),
  paste0("Prior & ", paste(rep(paste(
    "ELPD & $R^2$ & $R^2_{\\mathrm{oos}}$ &",
    "$\\widehat{R}$ & ESS & Div. & Time"), length(kernel_ids)),
    collapse = " & "), " \\\\"),
  "\\midrule",
  body,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
), final_path("performance_table.tex"),
preview = final_settings$compile_table_preview)

## ---------------------------------------------------------------------------
## Effects
## ---------------------------------------------------------------------------

effects <- Filter(Negate(is.null), setNames(lapply(kernel_ids, effects_summaries),
                                            kernel_ids))
if (final_settings$rank_kernel %in% names(effects)) {
  ranking <- effects[[final_settings$rank_kernel]]$relevance %>%
    filter(prior == final_settings$rank_prior) %>%
    arrange(desc(share_mean)) %>%
    pull(covariate)

  ## Per kernel: the 12 strongest descriptors and their variance shares.
  for (kernel in names(effects)) {
    dir.create(final_path(kernel), recursive = TRUE, showWarnings = FALSE)
    top12 <- head(ranking, final_settings$n_top12)
    save_paper_plot(final_path(kernel, "effects_top12.pdf"),
                    plot_curves(effects[[kernel]]$curves, top12),
                    width = final_settings$top12_size[["width"]],
                    height = final_settings$top12_size[["height"]],
                    device = final_settings$device)
    ## The same figure as TikZ code, <kernel>/effects_top12.tex (all labels
    ## as LaTeX, underscores in covariate names escaped). Full text width:
    ## include it with figure* in a two-column document.
    if (isTRUE(final_settings$tikz) &&
        requireNamespace("tikzDevice", quietly = TRUE)) {
      ggsave(final_path(kernel, "effects_top12.tex"),
             plot = plot_curves(effects[[kernel]]$curves, top12,
                                for_tikz = TRUE),
             device = tikzDevice::tikz,
             width = final_settings$top12_size[["width"]],
             height = final_settings$top12_size[["height"]])
      message("TikZ version: ", final_path(kernel, "effects_top12.tex"))
    }
    save_paper_plot(final_path(kernel, "variance_explained.pdf"),
                    plot_variance(effects[[kernel]]$relevance,
                                  head(ranking, final_settings$n_variance)),
                    width = final_settings$variance_size[["width"]],
                    height = final_settings$variance_size[["height"]],
                    device = final_settings$device)
  }

  ## The 2 strongest: kernels as rows, descriptors as columns, all priors.
  curves <- bind_rows(lapply(names(effects), function(kernel) {
    effects[[kernel]]$curves %>%
      mutate(kernel = factor(final_settings$kernels[[kernel]],
                             levels = final_settings$kernels[kernel_ids]))
  }))
  top2 <- head(ranking, final_settings$n_top2)
  top2_plot <- plot_curves(curves, top2, by_kernel = TRUE, short_y = TRUE)
  save_paper_plot(final_path("effects_top2.pdf"), top2_plot,
                  width = final_settings$top2_size[["width"]],
                  height = final_settings$top2_size[["height"]],
                  device = final_settings$device)

  ## The same figure as TikZ code, effects_top2.tex, for \input{} in the
  ## manuscript: LaTeX typesets its text in the document's own font. All
  ## labels are LaTeX there, so the kernels use their LaTeX names.
  if (isTRUE(final_settings$tikz)) {
    if (requireNamespace("tikzDevice", quietly = TRUE)) {
      curves_tex <- curves %>%
        mutate(kernel = factor(final_settings$kernel_tex[match(
          as.character(kernel), final_settings$kernels)],
          levels = final_settings$kernel_tex[kernel_ids]))
      ggsave(final_path("effects_top2.tex"),
             plot = plot_curves(curves_tex, top2, by_kernel = TRUE,
                                for_tikz = TRUE, short_y = TRUE),
             device = tikzDevice::tikz,
             width = final_settings$top2_size[["width"]],
             height = final_settings$top2_size[["height"]])
      message("TikZ version: ", final_path("effects_top2.tex"))
    } else {
      message("tikzDevice not installed; effects_top2.tex skipped.")
    }
  }
} else {
  message("Effects figures skipped: they need the ", final_settings$rank_kernel,
          " effects summaries for the ranking.")
}

message("final_outputs: written to ", final_settings$output_dir)
