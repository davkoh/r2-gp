## Compile the Stan model once, ahead of the fitting jobs.
##
## The output is a stamp, not the binary: cmdstanr writes the executable next
## to the .stan source and finds it again by mtime. The stamp exists so that
## Snakemake can order compilation before fitting and so that N concurrent fit
## jobs do not race to compile the same model.

library(cmdstanr)

stan_file <- snakemake@input[["stan_file"]]
output_file <- snakemake@output[["stamp"]]

model <- cmdstan_model(stan_file)

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(
  list(
    stan_file = normalizePath(stan_file),
    exe_file = model$exe_file(),
    compiled_at = Sys.time()
  ),
  output_file
)

message("compile_model: ", model$exe_file())
