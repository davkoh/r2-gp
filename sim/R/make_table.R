## Publication-ready LaTeX table: parameter recovery and sampling statistics.
##
## Replaces the parameter-RMSE figures. One row per (kernel, length-scale, D,
## prior model), averaged over replications, with two column groups:
##
##   Parameter recovery   RMSE of the posterior mean for alpha, rho and sigma,
##                        plus 90% interval coverage for alpha.
##   Sampling             worst Rhat, smallest bulk ESS, mean divergences per
##                        fit, and median wall time.
##
## Requires \usepackage{booktabs} and, for the long form, \usepackage{longtable}.
## Written by hand rather than through knitr::kable so that the grouped column
## headers and the blanking of repeated row labels come out right, and so the
## pipeline gains no extra package dependency.

library(dplyr)

source("R/conditions.R")

results <- readRDS(snakemake@input[["rds"]])
metrics <- results$metrics
diagnostics <- results$diagnostics
config <- results$config
output_file <- snakemake@output[["tex"]]

## Model order follows the config, so the table reads in the order the paper
## introduces the priors rather than alphabetically.
model_order <- vapply(config$models, `[[`, character(1), "id")

## A model is a (prior on alpha, prior on rho) pair, so it gets two narrow
## columns rather than one wide concatenated label. Spelling out
## "Independent, inverse-gamma" in every row pushed the table ~100pt past the
## text block in an 11pt article; these abbreviations bring it back inside.
alpha_prior_tex <- c(R2D2 = "R2D2", Independent = "Indep.",
                     RHS = "RHS", NG = "NG")
rho_prior_tex <- function(model) {
  switch(
    model$rho_prior,
    ig = "IG",
    vague = "LN",
    normal = sprintf("$N(\\rho, %s)$",
                     formatC(model$normal_sd %||% config$priors$rho_normal_sd,
                             format = "g"))
  )
}
model_alpha_labels <- setNames(
  alpha_prior_tex[vapply(config$models, `[[`, character(1), "alpha_prior")],
  model_order
)
model_rho_labels <- setNames(
  vapply(config$models, rho_prior_tex, character(1)), model_order
)
kernel_order <- intersect(names(KERNEL_LABELS), unique(metrics$kernel))

summary_table <- metrics %>%
  left_join(
    diagnostics %>% select(condition_id, model_id, replication, max_rhat,
                           min_ess_bulk, divergences, sampling_seconds),
    by = c("condition_id", "model_id", "replication")
  ) %>%
  group_by(kernel, rho, dimension, model_id) %>%
  summarise(
    n_fits = n(),
    n_ok = sum(ok),
    ## Averaged over replications; NA-tolerant so a failed fit removes that
    ## replicate from the average rather than the whole cell.
    rmse_alpha = mean(rmse_alpha_point, na.rm = TRUE),
    rmse_rho = mean(rmse_rho_point, na.rm = TRUE),
    rmse_sigma = mean(rmse_sigma_point, na.rm = TRUE),
    coverage_alpha = mean(coverage_alpha, na.rm = TRUE),
    ## Worst case across replications: a table that averaged Rhat would hide
    ## the one chain that failed.
    max_rhat = suppressWarnings(max(max_rhat, na.rm = TRUE)),
    min_ess = suppressWarnings(min(min_ess_bulk, na.rm = TRUE)),
    divergences = mean(divergences, na.rm = TRUE),
    seconds = median(sampling_seconds, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    kernel = factor(kernel, levels = kernel_order),
    model_id = factor(model_id, levels = model_order)
  ) %>%
  arrange(kernel, rho, dimension, model_id)

## ---------------------------------------------------------------------------
## Formatting
## ---------------------------------------------------------------------------

## A cell with no successful fit prints as an em dash rather than NaN.
fmt <- function(x, digits) {
  out <- formatC(x, format = "f", digits = digits)
  out[is.na(x) | is.nan(x) | is.infinite(x)] <- "---"
  out
}

latex_escape <- function(x) gsub("([&%$#_{}])", "\\\\\\1", x)

## Matern 3/2 reads better as a fraction in maths mode.
kernel_tex <- c(
  se = "Squared exp.", m32 = "Mat\\'ern $3/2$",
  m52 = "Mat\\'ern $5/2$", exp = "Exponential"
)

rows <- summary_table %>%
  transmute(
    kernel_tex = kernel_tex[as.character(kernel)],
    rho_tex = fmt(rho, 1),
    dimension_tex = as.character(dimension),
    alpha_tex = model_alpha_labels[as.character(model_id)],
    rho_prior_tex = model_rho_labels[as.character(model_id)],
    rmse_alpha = fmt(rmse_alpha, 3),
    rmse_rho = fmt(rmse_rho, 3),
    rmse_sigma = fmt(rmse_sigma, 3),
    coverage_alpha = fmt(100 * coverage_alpha, 0),
    max_rhat = fmt(max_rhat, 3),
    min_ess = fmt(min_ess, 0),
    divergences = fmt(divergences, 1),
    seconds = fmt(seconds, 0),
    ## Carried for the blanking and rule logic below.
    .kernel = as.character(kernel), .rho = rho, .dimension = dimension
  )

## Blank a label when it repeats, so the eye groups the rows. Done on the
## formatted strings, after sorting, so the first row of each block keeps it.
blank_repeats <- function(x, group) {
  ifelse(duplicated(group), "", x)
}
kernel_group <- rows$.kernel
rho_group <- paste(rows$.kernel, rows$.rho)
dimension_group <- paste(rows$.kernel, rows$.rho, rows$.dimension)

rows$kernel_tex <- blank_repeats(rows$kernel_tex, kernel_group)
rows$rho_tex <- blank_repeats(rows$rho_tex, rho_group)
rows$dimension_tex <- blank_repeats(rows$dimension_tex, dimension_group)

columns <- c("kernel_tex", "rho_tex", "dimension_tex", "alpha_tex",
             "rho_prior_tex", "rmse_alpha", "rmse_rho", "rmse_sigma",
             "coverage_alpha", "max_rhat", "min_ess", "divergences", "seconds")
body <- apply(rows[columns], 1, paste, collapse = " & ")
body <- paste0(body, " \\\\")

## A rule between kernels, a little air between length-scale blocks.
new_kernel <- c(FALSE, kernel_group[-1] != kernel_group[-length(kernel_group)])
new_rho <- c(FALSE, rho_group[-1] != rho_group[-length(rho_group)])
body[new_rho & !new_kernel] <- paste0("\\addlinespace\n", body[new_rho & !new_kernel])
body[new_kernel] <- paste0("\\midrule\n", body[new_kernel])

## ---------------------------------------------------------------------------
## Assembly
## ---------------------------------------------------------------------------

n_rep <- config$sim_reps
## Escaping: these R literals hold two backslashes, which reach sprintf as
## one escaped token and are emitted as \% and \alpha. Doubling again
## would leave a stray backslash in the LaTeX.
caption <- sprintf(
  paste(
    "Parameter recovery and sampling behaviour for the additive HSGP, over",
    "%d replications per cell. RMSE is of the posterior mean, taken across the",
    "$D$ components; coverage is of the nominal $90\\%%$ posterior interval",
    "for $\\alpha$. $\\widehat{R}$ and ESS are the worst case over all",
    "$\\alpha_d$, $\\rho_d$, $\\sigma$ and $\\mathrm{lp}$, taken over",
    "replications; divergences are the mean per fit and time is the median",
    "sampling time in seconds. Inputs are scaled to $[-1,1]$, so $\\rho$ is",
    "expressed relative to a unit half-range."
  ),
  n_rep
)

header <- c(
  "% Requires \\usepackage{booktabs} and \\usepackage{longtable}.",
  "% longtable needs TWO LaTeX passes to settle the column widths; after a",
  "% single pass the header rules appear narrower than the body.",
  "% Generated by sim/R/make_table.R -- do not edit by hand.",
  ## \small keeps the numeric columns compact; the table is wide even after
  ## the prior labels were split into two narrow columns.
  "\\begingroup\\small",
  "\\setlength{\\tabcolsep}{4pt}",
  "\\begin{longtable}{lllll rrrr rrrr}",
  sprintf("\\caption{%s}\\label{tab:parameter-recovery}\\\\", caption),
  "\\toprule",
  paste("Kernel & $\\rho$ & $D$ &",
        "\\multicolumn{2}{c}{Prior} &",
        "\\multicolumn{4}{c}{Parameter recovery} &",
        "\\multicolumn{4}{c}{Sampling} \\\\"),
  "\\cmidrule(lr){4-5} \\cmidrule(lr){6-9} \\cmidrule(lr){10-13}",
  paste("& & & $\\alpha$ & $\\rho$ &",
        "$\\alpha$ & $\\rho$ & $\\sigma$ & Cov.\\ \\% &",
        "$\\widehat{R}$ & ESS & Div. & s \\\\"),
  "\\midrule",
  "\\endfirsthead",
  ## Repeated header on continuation pages, which is the whole point of
  ## longtable for a table this long.
  "\\toprule",
  paste("Kernel & $\\rho$ & $D$ &",
        "\\multicolumn{2}{c}{Prior} &",
        "\\multicolumn{4}{c}{Parameter recovery} &",
        "\\multicolumn{4}{c}{Sampling} \\\\"),
  "\\cmidrule(lr){4-5} \\cmidrule(lr){6-9} \\cmidrule(lr){10-13}",
  paste("& & & $\\alpha$ & $\\rho$ &",
        "$\\alpha$ & $\\rho$ & $\\sigma$ & Cov.\\ \\% &",
        "$\\widehat{R}$ & ESS & Div. & s \\\\"),
  "\\midrule",
  "\\endhead",
  "\\bottomrule",
  "\\endfoot"
)

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
writeLines(c(header, body, "\\end{longtable}", "\\endgroup"), output_file)

## The same numbers as CSV, for inspection and for pulling single values into
## the manuscript text.
write.csv(summary_table, snakemake@output[["csv"]], row.names = FALSE)

n_incomplete <- sum(summary_table$n_ok < summary_table$n_fits)
message(sprintf(
  "make_table: %d rows written to %s%s",
  nrow(summary_table), output_file,
  if (n_incomplete > 0) sprintf(" (%d cells have failed fits)", n_incomplete) else ""
))
