## Turn a calibration run into the numbers that go in slurm_profile/config.yaml.
##
##   Rscript calibration_report.R summaries_calib
##
## Reads the per-fit sampling times the pipeline already records, fits the way
## cost grows with the number of basis coefficients, and prints a recommended
## runtime. Memory must come from sacct MaxRSS: an R process cannot see its own
## peak resident size reliably.

args <- commandArgs(trailingOnly = TRUE)
summaries <- if (length(args) > 0) args[1] else "summaries_calib"

path <- file.path(summaries, "diagnostics.csv")
if (!file.exists(path)) {
  stop("No ", path, ". Run the calibration first.", call. = FALSE)
}
diagnostics <- read.csv(path)
metrics <- read.csv(file.path(summaries, "metrics.csv"))

fits <- merge(
  diagnostics[c("condition_id", "model_id", "replication", "sampling_seconds",
                "max_rhat", "divergences")],
  unique(metrics[c("condition_id", "dimension", "kernel", "rho", "n_basis")]),
  by = "condition_id"
)
fits$coefficients <- fits$dimension * fits$n_basis

cat("\n== measured sampling time per fit ==\n")
by_cell <- aggregate(sampling_seconds ~ dimension + kernel + rho + n_basis,
                     fits, function(x) c(median = median(x), max = max(x)))
out <- data.frame(
  by_cell[c("kernel", "rho", "dimension", "n_basis")],
  coefficients = by_cell$dimension * by_cell$n_basis,
  median_s = round(by_cell$sampling_seconds[, "median"]),
  max_s = round(by_cell$sampling_seconds[, "max"])
)
print(out[order(out$coefficients), ], row.names = FALSE)

## Cost is modelled as a power law in the number of basis coefficients. With
## three D values this is a crude fit, but it is enough to extrapolate to a
## cell the calibration did not run.
slowest <- max(fits$sampling_seconds, na.rm = TRUE)
if (length(unique(fits$coefficients)) >= 2) {
  model <- lm(log(sampling_seconds) ~ log(coefficients), data = fits)
  exponent <- unname(coef(model)[2])
  cat(sprintf("\ncost scales as (D x M)^%.2f\n", exponent))
} else {
  exponent <- NA_real_
}

## Double the slowest observed fit. The margin covers the cells the calibration
## did not run and the variation between replications.
recommended_minutes <- ceiling(2 * slowest / 60)
cat(sprintf("\nslowest fit observed: %.0f s (%.1f min)\n", slowest, slowest / 60))
cat(sprintf("recommended runtime  : %d min  (2x the slowest)\n",
            recommended_minutes))

cat("\n== put this in slurm_profile/config.yaml ==\n")
cat(sprintf("  fit_model:\n    runtime: %d\n", recommended_minutes))
cat("    mem_mb: <1.5 x the MaxRSS that sacct reported above>\n")

cat("\n== choose group-components ==\n")
cat("A group job's walltime is components x runtime. Keep it under the\n")
cat("partition limit, and prefer jobs of a few hours over jobs of minutes.\n\n")
for (components in c(5, 10, 20, 40)) {
  hours <- components * recommended_minutes / 60
  submissions <- ceiling(1800 / components)
  flag <- if (hours > 120) "  <- over a 5 day limit" else ""
  cat(sprintf("  fit_model: %-3d -> %6.1f h per job, %4d submissions%s\n",
              components, hours, submissions, flag))
}

cat("\n== sampler health at these settings ==\n")
cat(sprintf("  worst Rhat        : %.3f\n", max(fits$max_rhat, na.rm = TRUE)))
cat(sprintf("  total divergences : %d\n", sum(fits$divergences, na.rm = TRUE)))
cat("A poor Rhat here means the MCMC settings need attention before the",
    "full run,\nnot more wall time.\n")
