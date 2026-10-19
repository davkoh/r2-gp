## Collect every fit into tidy tables.
##
## Two outputs: flat CSVs, one row per (condition, model, replication), for
## inspection and for the plotting rules; and an .rds carrying the same tables
## plus the settings, so figures can be redrawn without re-reading 1920 files.

source("R/conditions.R")

config <- yaml::read_yaml(snakemake@params[["config_file"]])
result_files <- unlist(snakemake@input[["results"]], use.names = FALSE)
results <- lapply(result_files, readRDS)

## One row per fit, flattening meta + the requested component.
bind_component <- function(component) {
  do.call(rbind, lapply(results, function(result) {
    data.frame(
      c(result$meta[c("condition_id", "dimension", "rho", "kernel",
                      "boundary", "n_basis", "model_id", "model_label",
                      "alpha_prior", "rho_prior", "replication", "ok")],
        result[[component]]),
      stringsAsFactors = FALSE
    )
  }))
}

metrics <- bind_component("metrics")
diagnostics <- bind_component("diagnostics")

## Labels are presentation, not data, but they are stored in each fit result at
## fit time. Re-resolve them from the config so that renaming a model does not
## require re-running 1920 fits.
model_labels <- vapply(config$models, `[[`, character(1), "label")
names(model_labels) <- vapply(config$models, `[[`, character(1), "id")
metrics$model_label <- unname(model_labels[metrics$model_id])
diagnostics$model_label <- unname(model_labels[diagnostics$model_id])

## ---------------------------------------------------------------------------
## Paired held-out ELPD differences against the reference model
## ---------------------------------------------------------------------------
##
## Paired per test point within a (condition, replication), which is what makes
## the standard error meaningful: the two models saw identical data, so the
## replicate-level noise cancels and only the predictive difference remains.

reference_model <- config$reference_model
pointwise_by_fit <- setNames(
  lapply(results, `[[`, "pointwise_elpd_test"),
  vapply(results, function(r) {
    paste(r$meta$condition_id, r$meta$replication, r$meta$model_id, sep = "\r")
  }, character(1))
)

elpd_rows <- list()
fits <- unique(metrics[c("condition_id", "dimension", "rho", "kernel",
                         "replication")])
for (i in seq_len(nrow(fits))) {
  cell <- fits[i, ]
  key <- function(model_id) {
    paste(cell$condition_id, cell$replication, model_id, sep = "\r")
  }
  reference <- pointwise_by_fit[[key(reference_model)]]
  if (is.null(reference)) {
    next
  }
  for (model in config$models) {
    if (identical(model$id, reference_model)) {
      next
    }
    competitor <- pointwise_by_fit[[key(model$id)]]
    if (is.null(competitor)) {
      next
    }
    ## Positive favours the reference model.
    difference <- reference - competitor
    elpd_rows[[length(elpd_rows) + 1L]] <- data.frame(
      cell,
      reference = reference_model,
      model_id = model$id,
      model_label = model$label,
      alpha_prior = model$alpha_prior,
      rho_prior = model$rho_prior,
      elpd_difference = sum(difference),
      se_difference = sqrt(length(difference) * var(difference)),
      mlpd_difference = mean(difference),
      se_mlpd_difference = sd(difference) / sqrt(length(difference)),
      stringsAsFactors = FALSE
    )
  }
}
elpd <- if (length(elpd_rows) > 0L) do.call(rbind, elpd_rows) else NULL

## ---------------------------------------------------------------------------
## Long-form parameter RMSE, so plotting can facet over the parameter
## ---------------------------------------------------------------------------

parameter_rmse <- do.call(rbind, lapply(
  c("alpha", "rho", "sigma"),
  function(parameter) {
    data.frame(
      metrics[c("condition_id", "dimension", "rho", "kernel", "model_id",
                "model_label", "alpha_prior", "rho_prior", "replication", "ok")],
      parameter = parameter,
      rmse_point = metrics[[paste0("rmse_", parameter, "_point")]],
      rmse_draws = metrics[[paste0("rmse_", parameter, "_draws")]],
      bias = metrics[[paste0("bias_", parameter)]],
      coverage = metrics[[paste0("coverage_", parameter)]],
      stringsAsFactors = FALSE
    )
  }
))

summaries <- dirname(snakemake@output[["metrics"]])
dir.create(summaries, recursive = TRUE, showWarnings = FALSE)
write.csv(metrics, snakemake@output[["metrics"]], row.names = FALSE)
write.csv(diagnostics, snakemake@output[["diagnostics"]], row.names = FALSE)
write.csv(parameter_rmse, snakemake@output[["parameter_rmse"]], row.names = FALSE)
write.csv(elpd %||% data.frame(), snakemake@output[["elpd"]], row.names = FALSE)
saveRDS(
  list(metrics = metrics, diagnostics = diagnostics, elpd = elpd,
       parameter_rmse = parameter_rmse, config = config),
  snakemake@output[["rds"]]
)

n_failed <- sum(!metrics$ok)
message(sprintf(
  "aggregate_results: %d fits (%d failed), %d ELPD comparisons vs %s",
  nrow(metrics), n_failed, if (is.null(elpd)) 0L else nrow(elpd), reference_model
))
