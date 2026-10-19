## Generate every replicate for one condition.

source("R/dgp.R")

config <- yaml::read_yaml(snakemake@params[["config_file"]])
condition <- find_condition(
  build_conditions(config), snakemake@params[["condition_id"]]
)
output_files <- unlist(snakemake@output[["datasets"]], use.names = FALSE)

## The replicate id comes from the filename, never from loop order, so a
## partial rerun writes the dataset it was asked for.
rep_ids <- suppressWarnings(
  as.integer(sub("^rep_([0-9]+)\\.rds$", "\\1", basename(output_files)))
)
if (anyNA(rep_ids)) {
  stop("Could not parse a replicate id from: ",
       paste(basename(output_files)[is.na(rep_ids)], collapse = ", "),
       call. = FALSE)
}

for (directory in unique(dirname(output_files))) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

for (i in seq_along(output_files)) {
  dataset <- simulate_dataset(
    condition, config$design,
    ## Data uses seed + rep; sampling uses seed + 1000 * rep (fit_model.R).
    ## Keeping the two streams apart makes a replicate reproducible
    ## independently of how jobs are grouped.
    seed = as.integer(config$seed) + rep_ids[i]
  )
  saveRDS(dataset, output_files[i])
}

## The boundary rule only holds if the scaling did what it claims, so assert it
## once per condition rather than trusting it.
last <- readRDS(output_files[length(output_files)])
observed_half_range <- half_range(last$inputs)
stopifnot(all(abs(observed_half_range - 1) < 1e-8))

message(sprintf(
  "generate_data: %s - %d reps, L = %.2f, M = %d (generated with M = %d)",
  condition$id, length(output_files), condition$boundary,
  condition$n_basis, condition$n_basis_gen
))
