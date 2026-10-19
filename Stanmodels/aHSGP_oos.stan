// Additive univariate HSGPs with a held-out test set, supporting four priors
// on the D component scales alpha[d].
//
//   prior_type = 1  R2D2         alpha[d]^2 = sigma^2 R2/(1-R2) * psi[d],
//                                psi ~ Dirichlet  (joint variance decomposition)
//   prior_type = 2  Independent  alpha[d] ~ half-normal(0, s)   (no pooling)
//   prior_type = 3  RHS          regularized horseshoe on alpha, Piironen &
//                                Vehtari (2017), tau0 = p0/(D-p0) sigma/sqrt(N)
//   prior_type = 4  NG           normal-gamma on alpha, Griffin & Brown (2010)
//
// Types 3 and 4 are the non-Dirichlet competitors: they shrink the component
// scales marginally and independently rather than partitioning a fixed total.
//
// The covariance function is a data switch (covfn): squared exponential,
// Matern 3/2, Matern 5/2, or exponential. covfn = 3
// (exponential) is the default and reproduces the original behaviour.
//
// Identical prior/likelihood structure to aHSGP_R2D2.stan. The additions are
//   * a test design (inputs_test, y_test) evaluated on the same basis, and
//   * log_lik_test / f_test in generated quantities.
// The spectral normalizer is still computed from the TRAINING design, so
// alpha[d]^2 keeps its interpretation as the expected in-sample variance of
// component d; test predictions reuse the same spectral scales.
//
// f_component and the spectral scales are kept local so that only N-sized
// quantities are stored, which matters once D reaches 100.
// See ../R2D2_HSGP_NORMALIZATION.md for the derivation and data interface.

functions {
  vector phi(real L, int m, vector x) {
    return (1 / sqrt(L)) * sin((m * pi()) / (2 * L) * (x + L));
  }

  // Spectral densities, with alpha^2 factored out: the amplitude enters
  // through spectral_scale below. Under normalize_variance = 1 any constant
  // factor in S cancels in alpha * sqrt(S / dot(S, v)), so these constants
  // matter only for the unnormalised control arm. Lifted from
  // Stanmodels/sHSGP.stan, which carries the same formulas with alpha^2 in.
  real spd_se(real rho, real w) {
    return sqrt(2 * pi()) * rho * exp(-0.5 * square(rho * w));
  }

  real spd_m32(real rho, real w) {
    return (4 * pow(3, 1.5) / pow(rho, 3))
           * pow(3.0 / square(rho) + square(w), -2);
  }

  real spd_m52(real rho, real w) {
    return (4 * pow(5, 2.5) / (0.75 * pow(rho, 5)))
           * pow(5.0 / square(rho) + square(w), -3);
  }

  real spd_exponential(real rho, real omega) {
    return (2 * rho) / (1 + square(rho * omega));
  }

  // covfn: 0 = squared exponential, 1 = Matern 3/2, 2 = Matern 5/2,
  // 3 = exponential (Matern 1/2). 3 is the default so that callers predating
  // the kernel switch keep their original behaviour.
  real spd(int covfn, real rho, real w) {
    if (covfn == 0) {
      return spd_se(rho, w);
    } else if (covfn == 1) {
      return spd_m32(rho, w);
    } else if (covfn == 2) {
      return spd_m52(rho, w);
    }
    return spd_exponential(rho, w);
  }

  // Regularized-horseshoe scales: the slab c2 caps how large a component
  // scale the heavy-tailed local parameter can produce.
  vector regularized_horseshoe_scale(vector local_scale, real tau, real c2) {
    vector[rows(local_scale)] local2 = square(local_scale);
    return sqrt(c2 * local2 ./ (c2 + square(tau) * local2)) * tau;
  }
}

data {
  int<lower=2> N;
  int<lower=1> D;
  int<lower=1> M;
  matrix[N, D] inputs;
  vector[N] y;

  int<lower=0> N_test;
  matrix[N_test, D] inputs_test;
  vector[N_test] y_test;

  vector<lower=0>[D] L;
  vector[D] x_min;
  vector[D] x_max;
  int<lower=0, upper=1> is_centered;
  int<lower=0, upper=1> prior_only;
  int<lower=1, upper=4> prior_type;
  int<lower=0, upper=1> estimate_rho;
  // 1: divide the spectral weights by the design-based normaliser, so that
  // alpha[d]^2 is the expected sample variance of component d. 
  // 0: the typical HSGP
  // parameterisation in which alpha[d] is the kernel marginal SD and is
  // therefore coupled to rho[d].
  int<lower=0, upper=1> normalize_variance;
  // Prior on the lengthscales. 
  // 1: the inverse-gamma prior. 
  // 2: lognormal(0, rho_vague_sd), which has a thicker right tail
  // 3: normal(rho_normal[1], rho_normal[2]). 
  int<lower=1, upper=3> rho_prior_type;
  real<lower=0> rho_vague_sd;
  vector<lower=0>[2] rho_normal;

  // Covariance function, see spd() above. 3 (exponential) is the default.
  int<lower=0, upper=3> covfn;

  real<lower=0, upper=1> r2_mean;
  real<lower=0> r2_precision;
  real<lower=0> dirichlet_concentration;
  vector<lower=0>[2] lengthscale_inv_gamma;
  vector<lower=0>[D] rho_fixed;
  real<lower=0> independent_alpha_scale;
  real<lower=0> sigma_scale;
  real<lower=0> intercept_scale;

  // Regularized horseshoe (prior_type 3)
  real<lower=0> hs_df;
  real<lower=0> hs_df_global;
  real<lower=0> hs_df_slab;
  real<lower=0> hs_scale_slab;
  real<lower=0> hs_p0;

  // Normal-gamma (prior_type 4)
  real<lower=0> ng_nu_scale;
}

transformed data {
  array[D] matrix[N, M] PHI;
  array[D] matrix[N_test, M] PHI_test;
  matrix[M, D] omega;
  matrix[M, D] basis_sample_variance;
  vector[D] midpoint = 0.5 * (x_min + x_max);
  int use_r2d2 = prior_type == 1;
  int use_independent = prior_type == 2;
  int use_rhs = prior_type == 3;
  int use_ng = prior_type == 4;

  for (d in 1:D) {
    vector[N] x_d = inputs[, d];
    vector[N_test] x_test_d = inputs_test[, d];
    if (is_centered) {
      x_d -= midpoint[d];
      x_test_d -= midpoint[d];
    }
    for (m in 1:M) {
      PHI[d][, m] = phi(L[d], m, x_d);
      PHI_test[d][, m] = phi(L[d], m, x_test_d);
      omega[m, d] = (m * pi()) / (2 * L[d]);
      basis_sample_variance[m, d] = variance(PHI[d][, m]);
    }
  }
}

parameters {
  array[D] vector[M] z;
  vector<lower=0>[estimate_rho ? D : 0] rho_parameter;
  real<lower=0> sigma;
  real intercept;

  // Only the block belonging to prior_type is given nonzero size, so the
  // unused arms cost neither parameters nor sampling time.
  array[use_r2d2 ? 1 : 0] real<lower=0, upper=1> R2_raw;
  simplex[use_r2d2 ? D : 1] variance_proportion;
  vector<lower=0>[use_independent ? D : 0] alpha_independent;

  // Regularized horseshoe: alpha = z * lambda_tilde * tau, z half-normal.
  vector<lower=0>[use_rhs ? D : 0] hs_z;
  vector<lower=0>[use_rhs ? D : 0] hs_local;
  array[use_rhs ? 1 : 0] real<lower=0> hs_global;
  array[use_rhs ? 1 : 0] real<lower=0> hs_slab;

  // Normal-gamma: alpha = sigma * z * sqrt(psi2), z half-normal.
  vector<lower=0>[use_ng ? D : 0] ng_z;
  vector<lower=0>[use_ng ? D : 0] ng_psi2;
  array[use_ng ? 1 : 0] real<lower=0> ng_lambda;
  array[use_ng ? 1 : 0] real<lower=0> ng_nu;
}

transformed parameters {
  real<lower=0, upper=1> R2 = use_r2d2 ? R2_raw[1] : 0.5;
  real<lower=0> global_signal_variance = square(sigma) * R2 / (1 - R2);
  vector<lower=0>[D] rho;
  vector<lower=0>[D] alpha;
  vector[N] f;

  if (estimate_rho) {
    rho = rho_parameter;
  } else {
    rho = rho_fixed;
  }

  if (use_r2d2) {
    for (d in 1:D) {
      alpha[d] = sqrt(global_signal_variance * variance_proportion[d]);
    }
  } else if (use_independent) {
    alpha = alpha_independent;
  } else if (use_rhs) {
    alpha = hs_z .* regularized_horseshoe_scale(
      hs_local, hs_global[1], square(hs_scale_slab) * hs_slab[1]
    );
  } else {
    alpha = sigma * ng_z .* sqrt(ng_psi2);
  }

  {
    matrix[N, D] f_component;
    for (d in 1:D) {
      vector[M] spectral_density;
      vector[M] spectral_scale;
      real normalizer;
      for (m in 1:M) {
        spectral_density[m] = spd(covfn, rho[d], omega[m, d]);
      }
      normalizer = normalize_variance
        ? dot_product(spectral_density, basis_sample_variance[, d])
        : 1.0;
      spectral_scale = alpha[d] * sqrt(spectral_density / normalizer);
      f_component[, d] = PHI[d] * (spectral_scale .* z[d]);
    }
    f = intercept + f_component * rep_vector(1.0, D);
  }
}

model {
  for (d in 1:D) {
    z[d] ~ std_normal();
  }
  if (estimate_rho) {
    if (rho_prior_type == 1) {
      rho_parameter ~ inv_gamma(
        lengthscale_inv_gamma[1], lengthscale_inv_gamma[2]
      );
    } else if (rho_prior_type == 2) {
      rho_parameter ~ lognormal(0, rho_vague_sd);
    } else {
      rho_parameter ~ normal(rho_normal[1], rho_normal[2]);
    }
  }
  sigma ~ normal(0, sigma_scale);
  intercept ~ normal(0, intercept_scale);

  if (use_r2d2) {
    R2_raw ~ beta(r2_mean * r2_precision, (1 - r2_mean) * r2_precision);
    variance_proportion ~ dirichlet(rep_vector(dirichlet_concentration, D));
  } else if (use_independent) {
    alpha_independent ~ normal(0, independent_alpha_scale);
  } else if (use_rhs) {
    // Piironen & Vehtari (2017): the global scale encodes a prior guess hs_p0
    // on the number of active components.
    hs_z ~ std_normal();
    hs_local ~ student_t(hs_df, 0, 1);
    hs_global ~ student_t(
      hs_df_global, 0, hs_p0 / (D - hs_p0) * sigma / sqrt(N)
    );
    hs_slab ~ inv_gamma(0.5 * hs_df_slab, 0.5 * hs_df_slab);
  } else {
    // Griffin & Brown (2010): psi2 ~ Gamma(lambda, 1/(2 gamma^2)) with
    // gamma^2 = nu / (2 lambda), so E[psi2] = nu and E[alpha^2] = sigma^2 nu.
    ng_z ~ std_normal();
    ng_lambda ~ exponential(1);
    ng_nu ~ inv_gamma(2, ng_nu_scale);
    ng_psi2 ~ gamma(ng_lambda[1], ng_lambda[1] / ng_nu[1]);
  }

  if (!prior_only) {
    y ~ normal(f, sigma);
  }
}

generated quantities {
  real implied_signal_variance = dot_self(alpha);
  real implied_R2 = implied_signal_variance
    / (implied_signal_variance + square(sigma));
  real signal_variance = variance(f);
  real R2_empirical = signal_variance / (signal_variance + square(sigma));
  vector[N] log_lik;
  vector[N_test] log_lik_test;
  vector[N_test] f_test;
  real R2_empirical_test;
  real R2_oos;

  {
    // Recomputed rather than stored: saving the D x M scales for every draw
    // dominates the output size once D is large.
    matrix[N_test, D] f_component_test;
    for (d in 1:D) {
      vector[M] spectral_density;
      vector[M] spectral_scale;
      for (m in 1:M) {
        spectral_density[m] = spd(covfn, rho[d], omega[m, d]);
      }
      spectral_scale = alpha[d] * sqrt(
        spectral_density / (normalize_variance
          ? dot_product(spectral_density, basis_sample_variance[, d])
          : 1.0)
      );
      f_component_test[, d] = PHI_test[d] * (spectral_scale .* z[d]);
    }
    f_test = intercept + f_component_test * rep_vector(1.0, D);
  }

  // Model-based Bayesian R2 on the test design. This uses the model's own
  // sigma, not the held-out residuals, so it measures the implied
  // signal-to-noise at new inputs rather than out-of-sample fit.
  R2_empirical_test = variance(f_test) / (variance(f_test) + square(sigma));

  // Residual-based Bayesian R2 (Gelman et al. 2019) on the held-out set. The
  // denominator uses the realised test residuals, so a model that predicts
  // badly is penalised. This is the genuine out-of-sample quantity.
  R2_oos = variance(f_test) / (variance(f_test) + variance(y_test - f_test));

  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y[n] | f[n], sigma);
  }
  for (n in 1:N_test) {
    log_lik_test[n] = normal_lpdf(y_test[n] | f_test[n], sigma);
  }
}
