## Manuscript figures 1 and 2, reproduced inside the pipeline.
##
##   fig01_variance_gap.pdf     the gap between the pointwise and the sample
##                              variance of a GP realisation, and the
##                              normaliser Q(rho) that closes it.
##   fig02_prior_predictive.pdf the implied prior on the empirical R2 of an
##                              additive normalised HSGP under each prior on
##                              the component scales, as D grows.
##
## Both are ports of paper/figures/fig01_variance_gap.R and
## fig02_prior_predictive.R, unchanged in appearance. What changes is that the
## quantities they depend on now come from the pipeline rather than from
## constants copied into the figure script: figure 1 uses exact_kernel() and
## figure 2 uses the config's priors, spectral densities and basis settings. A
## figure that restated its own settings would silently disagree with the
## simulation the moment either was retuned.
##
## Figure 1 deliberately uses the EXACT kernel, not the HSGP basis: the
## introduction motivates the problem for the Gaussian process itself, before
## the approximation is introduced.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

source("R/dgp.R")
source("R/theme_hsgp.R")

config <- yaml::read_yaml(snakemake@params[["config_file"]])
settings <- config$paper_figures %||% list()
priors <- config$priors
design <- config$design

## ---------------------------------------------------------------------------
## Figure 1: pointwise variance versus sample variance
## ---------------------------------------------------------------------------

set.seed(as.integer(settings$fig01_seed %||% 20260908L))

fig01_file <- snakemake@output[["fig01"]]
dir.create(dirname(fig01_file), recursive = TRUE, showWarnings = FALSE)

n_grid <- as.integer(settings$fig01_n_grid %||% 400L)
n_show <- as.integer(settings$fig01_n_show %||% 3L)
n_mc <- as.integer(settings$fig01_n_mc %||% 4000L)
fig01_kernel <- settings$fig01_kernel %||% "exp"
x <- seq(-1, 1, length.out = n_grid)
rho_set <- as.numeric(settings$fig01_lengthscales %||% c(0.1, 0.5, 2))
panel_lab <- c("(a)", "(b)", "(c)")

kernel_matrix <- function(x, rho) {
  outer(x, x, function(a, b) exact_kernel(fig01_kernel, a - b, rho))
}

## Q(rho) = E[sample variance of a unit-amplitude realisation]
##        = tr(P_N K) / (N - 1), with P_N the centring matrix. This is the
## quantity derived in ../../R2D2_HSGP_NORMALIZATION.md and computed in Stan as
## dot_product(spectral_density, basis_sample_variance).
q_normalizer <- function(x, rho) {
  n <- length(x)
  k <- kernel_matrix(x, rho)
  (sum(diag(k)) - sum(k) / n) / (n - 1)
}

draw_gp <- function(x, rho, n_draw) {
  k <- kernel_matrix(x, rho)
  root <- chol(k + diag(1e-9, length(x)))
  crossprod(root, matrix(rnorm(length(x) * n_draw), nrow = length(x)))
}

## Base graphics, as in the manuscript: four narrow panels with an L-shaped box.
pal <- c("#1b3a6b", "#c2571a", "#4c8055")

pdf(fig01_file, width = 7.2, height = 2.8, pointsize = 13)
op <- par(mfrow = c(1, 4), mar = c(3.3, 3.2, 3.2, 0.6), mgp = c(2.0, 0.6, 0),
          las = 1, bty = "l", cex.axis = 0.9, cex.lab = 1.05)

report <- data.frame()
for (j in seq_along(rho_set)) {
  rho <- rho_set[j]
  q <- q_normalizer(x, rho)
  f <- draw_gp(x, rho, n_show)
  sample_sd <- mean(apply(draw_gp(x, rho, n_mc), 2, sd))

  plot(NA, xlim = range(x), ylim = c(-3.1, 3.1), xlab = "x", ylab = "f(x)",
       main = "")
  ## Shaded band: plus or minus one pointwise standard deviation.
  rect(par("usr")[1], -1, par("usr")[2], 1, col = "grey92", border = NA)
  abline(h = c(-1, 1), lty = 2, col = "grey40")
  abline(h = 0, col = "grey70")
  for (i in seq_len(n_show)) lines(x, f[, i], col = pal[i], lwd = 1.6)
  box(bty = "l")
  title(main = bquote(.(panel_lab[j]) ~ ~ rho == .(format(rho, nsmall = 1))),
        line = 2.0)
  mtext(sprintf("sample sd %.2f    Q = %.2f", sample_sd, q),
        side = 3, line = 0.6, cex = 0.78, col = "grey25")
  report <- rbind(report, data.frame(rho = rho, Q = q, sample_sd = sample_sd,
                                     sqrt_Q = sqrt(q)))
}

## Panel (d): Q against the length scale, on the same input window.
rho_grid <- exp(seq(log(0.02), log(10), length.out = 120))
q_grid <- vapply(rho_grid, function(r) q_normalizer(x, r), numeric(1))

plot(rho_grid, q_grid, type = "l", log = "x", lwd = 2, col = "black",
     ylim = c(0, 1.02), xlab = expression(rho), ylab = expression(Q(rho)),
     main = "")
title(main = "(d)  normaliser", line = 2.0)
mtext(expression(Q(rho) == E * "[sample var]"), side = 3, line = 0.6,
      cex = 0.78, col = "grey25")
abline(h = 1, lty = 3, col = "grey55")
points(rho_set, vapply(rho_set, function(r) q_normalizer(x, r), numeric(1)),
       pch = 21, bg = pal, col = "black", cex = 1.3)
text(rho_set, vapply(rho_set, function(r) q_normalizer(x, r), numeric(1)),
     labels = c("a", "b", "c"), pos = c(1, 4, 4), offset = 0.7, cex = 0.85)

par(op)
invisible(dev.off())

report$rho_over_range <- report$rho / diff(range(x))
write.csv(report, snakemake@output[["fig01_csv"]], row.names = FALSE)

## ---------------------------------------------------------------------------
## Figure 2: implied prior on the empirical R2 as D grows
## ---------------------------------------------------------------------------

set.seed(as.integer(settings$fig02_seed %||% 20260915L))

n_draws <- as.integer(settings$fig02_n_draws %||% 2000L)
n_obs <- as.integer(settings$fig02_n_obs %||% 200L)
n_basis <- as.integer(settings$fig02_n_basis %||% 100L)
dims <- as.integer(settings$fig02_dimensions %||% c(1, 2, 5, 10, 30, 100))
boundary <- as.numeric(settings$fig02_boundary %||% 1.5)
fig02_kernel <- settings$fig02_kernel %||% "exp"
fig02_rho <- settings$fig02_lengthscale %||% "prior"
sigma <- design$sigma

freq <- basis_frequencies(n_basis, boundary)
inputs <- matrix(runif(n_obs * max(dims), -1, 1), nrow = n_obs)
rinvgamma <- function(n, shape, scale) 1 / rgamma(n, shape = shape, rate = scale)
half_t <- function(n, df, scale = 1) abs(rt(n, df)) * scale

## The D component scales under each prior, mirroring aHSGP_oos.stan exactly.
alpha_draws <- function(prior, D) {
  if (prior == "R2D2") {
    r2 <- rbeta(n_draws, priors$r2_mean * priors$r2_precision,
                (1 - priors$r2_mean) * priors$r2_precision)
    g <- matrix(rgamma(n_draws * D, priors$dirichlet_concentration),
                nrow = n_draws)
    sqrt(g / rowSums(g) * (sigma^2 * r2 / (1 - r2)))
  } else if (prior == "Independent") {
    matrix(abs(rnorm(n_draws * D, sd = priors$independent_alpha_scale)),
           nrow = n_draws)
  } else if (prior == "Horseshoe") {
    ## Piironen & Vehtari (2017); p0 = D/2 so tau0 = sigma / sqrt(N).
    tau0 <- (D / 2) / (D - D / 2) * sigma / sqrt(n_obs)
    tau <- half_t(n_draws, priors$hs_df_global, tau0)
    slab <- rinvgamma(n_draws, 0.5 * priors$hs_df_slab, 0.5 * priors$hs_df_slab)
    c2 <- priors$hs_scale_slab^2 * slab
    lambda <- matrix(half_t(n_draws * D, priors$hs_df), nrow = n_draws)
    lambda2 <- lambda^2
    scale <- sqrt(c2 * lambda2 / (c2 + tau^2 * lambda2)) * tau
    matrix(abs(rnorm(n_draws * D)), nrow = n_draws) * scale
  } else {
    ## Griffin & Brown (2010); nu_scale set so E[total signal variance]
    ## matches the R2D2 model at every D.
    lambda <- rexp(n_draws)
    nu <- rinvgamma(n_draws, 2, (priors$r2_mean / (1 - priors$r2_mean)) / D)
    psi2 <- matrix(rgamma(n_draws * D, shape = lambda, rate = lambda / nu),
                   nrow = n_draws)
    sigma * matrix(abs(rnorm(n_draws * D)), nrow = n_draws) * sqrt(psi2)
  }
}

prior_names <- settings$fig02_priors %||%
  c("R2D2", "Independent", "Horseshoe", "Normal-gamma")

## All priors are applied to the SAME normalised latent components, drawn once
## per (D, d) and scaled each way, so the panels differ only through the prior
## on alpha and not through Monte Carlo noise.
draw_lengthscales <- function() {
  if (identical(fig02_rho, "prior")) {
    do.call(rinvgamma, c(list(n_draws), as.list(lengthscale_inv_gamma(priors))))
  } else {
    rep(as.numeric(fig02_rho), n_draws)
  }
}

r2_rows <- list()
quantile_rows <- list()
for (dimension in dims) {
  alpha <- lapply(prior_names, alpha_draws, D = dimension)
  names(alpha) <- prior_names
  f <- lapply(prior_names, function(p) matrix(0, n_obs, n_draws))
  names(f) <- prior_names

  for (d in seq_len(dimension)) {
    basis <- phi_matrix(inputs[, d], n_basis, boundary)
    basis_variance <- apply(basis, 2, var)
    rho <- draw_lengthscales()
    ## One spectral density column per draw, then the design normaliser Q.
    density <- vapply(rho, function(r) spectral_density(fig02_kernel, r, freq),
                      numeric(n_basis))
    q <- as.vector(crossprod(basis_variance, density))
    coefficients <- sqrt(density / rep(q, each = n_basis)) *
      matrix(rnorm(n_basis * n_draws), n_basis, n_draws)
    ## Unit-signal-variance component, shared across the priors.
    component <- basis %*% coefficients
    for (p in prior_names) {
      f[[p]] <- f[[p]] + component * rep(alpha[[p]][, d], each = n_obs)
    }
  }

  for (p in prior_names) {
    variance <- apply(f[[p]], 2, var)
    r2 <- variance / (variance + sigma^2)
    r2_rows[[length(r2_rows) + 1L]] <- data.frame(prior = p, D = dimension,
                                                  r2 = r2)
    quantile_rows[[length(quantile_rows) + 1L]] <- data.frame(
      prior = p, D = dimension, q05 = quantile(r2, 0.05, names = FALSE),
      median = median(r2), q95 = quantile(r2, 0.95, names = FALSE),
      mean = mean(r2)
    )
  }
  message("paper_figures: fig02 D = ", dimension, " done")
}

r2_draws <- bind_rows(r2_rows) %>%
  mutate(prior = factor(prior, levels = prior_names))
r2_quantiles <- bind_rows(quantile_rows) %>%
  mutate(prior = factor(prior, levels = prior_names))
write.csv(r2_quantiles, snakemake@output[["fig02_csv"]], row.names = FALSE)

## The competitor priors place large mass at the boundaries, which makes a
## density plot unreadable. The empirical CDF shows the whole distribution,
## boundary mass included, and quantiles are read off as horizontal cuts.
target <- data.frame(r2 = seq(0, 1, length.out = 512))
target$cdf <- pbeta(target$r2, priors$r2_mean * priors$r2_precision,
                    (1 - priors$r2_mean) * priors$r2_precision)

fig02 <- ggplot(r2_draws, aes(r2, colour = factor(D))) +
  stat_ecdf(geom = "step", linewidth = 0.55) +
  geom_line(data = target, aes(r2, y = cdf), inherit.aes = FALSE,
            linetype = "dashed", linewidth = 0.75, colour = "black") +
  facet_wrap(~prior, nrow = 1) +
  scale_colour_viridis_d(option = "plasma", end = 0.88,
                         name = "Components  D") +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = expression(R[emp]^2), y = "Prior CDF") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey93", colour = NA),
        legend.key.width = unit(1.1, "lines"))

ggsave(snakemake@output[["fig02"]], fig02, width = 7.8, height = 3.1,
       device = "pdf")

message("paper_figures: wrote ", fig01_file, " and ",
        snakemake@output[["fig02"]])
