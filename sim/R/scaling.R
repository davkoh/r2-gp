## Columnwise input scaling for the HSGP boundary.
##
## The HSGP eigenfunctions are defined on [-L, L] and the approximation
## degrades near the boundary, so L has to be comfortably larger than the
## range of the design. Soham's choose_L_multi() in ../varhsgpfns.R derives
## L_d = c * half-range(x_d) from the empirical range of each column. If the
## columns are first mapped to [-1, 1] that half-range is exactly 1, so
##
##     L_d = c   for every d,
##
## a single boundary factor governs every dimension, and rho is always
## expressed relative to a unit half-range. That is what makes the
## Riutort-Mayol et al. (2023) rules in conditions.R apply uniformly across a
## grid that varies D, the length-scale, and the kernel.
##
## This map is computed from the TRAINING design only. Test points may
## then fall slightly outside [-1, 1], which is the honest situation for
## held-out data and is still far inside [-L, L] for any c >= 1.5.

## Map each column of x_train to [-1, 1] and apply the same affine map to
## x_test. Returns the scaled matrices together with the map itself, so that
## real datasets can be transformed back.
scale_inputs <- function(x_train, x_test = NULL) {
  x_train <- as.matrix(x_train)

  column_min <- apply(x_train, 2, min)
  column_max <- apply(x_train, 2, max)
  center <- 0.5 * (column_min + column_max)
  halfwidth <- 0.5 * (column_max - column_min)
  ## A constant column has no range to scale; leave it centred at zero rather
  ## than dividing by zero.
  halfwidth[halfwidth == 0] <- 1

  apply_map <- function(x) {
    if (is.null(x)) {
      return(NULL)
    }
    x <- as.matrix(x)
    sweep(sweep(x, 2, center, "-"), 2, halfwidth, "/")
  }

  list(
    train = apply_map(x_train),
    test = apply_map(x_test),
    center = center,
    halfwidth = halfwidth
  )
}

## Invert scale_inputs() for one matrix on the scaled axis.
unscale_inputs <- function(x_scaled, scaling) {
  x_scaled <- as.matrix(x_scaled)
  sweep(sweep(x_scaled, 2, scaling$halfwidth, "*"), 2, scaling$center, "+")
}

## The empirical half-range per column, i.e. the S in L_d = c * S_d. After
## scale_inputs() this is 1 for every training column by construction; the
## pipeline asserts that rather than assuming it.
half_range <- function(x) {
  x <- as.matrix(x)
  apply(x, 2, function(column) {
    range_d <- diff(range(column, na.rm = TRUE))
    if (range_d == 0) 1 else range_d / 2
  })
}
