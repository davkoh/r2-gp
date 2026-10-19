
// Hilbert space approximate GP (HSGP) multiple inputs

// Data-interface changes to be aware of when preparing inputs (compared to moHSGP):
//   inputs : matrix[N, D]      (moHSGP: vector[N])
//   y      : vector[N]         (moHSGP: matrix[N, D])
//   M      : array[D] int      (moHSGP: single int)
//   L      : vector[D]         (moHSGP: single real; L[d] = c[d] * S[d])
//   x_min, x_max, s : vector[D] (moHSGP: scalars)
//   deriv_dim : int in 1..D    (which input dim to differentiate wrt)

functions {
	// Univariate eigenfunction for one input dimension.
	vector phi(real L, int m, vector x) {
		return (1 / sqrt(L)) * sin((m * pi()) / (2 * L) * (x + L));
	}

	// Tensor-product basis matrix PHI[N, M_total].
	// Column j is the product over input dims of the univariate eigenfunctions
	// indexed by the multi-index S_idx[j, ].
	matrix build_PHI(matrix xs, array[,] int S_idx, vector L) {
		int N = rows(xs);
		int D = cols(xs);
		int Mtot = size(S_idx);
		matrix[N, Mtot] PHI;
		for (j in 1:Mtot) {
			vector[N] phi_j = rep_vector(1.0, N);
			for (d in 1:D) {
				phi_j = phi_j .* phi(L[d], S_idx[j, d], xs[, d]);
			}
			PHI[, j] = phi_j;
		}
		return PHI;
	}

	// ARD squared-exponential spectral density at frequency vector w.
	// S(w) = alpha^2 (2*pi)^{D/2} prod_d rho_d exp(-0.5 sum_d rho_d^2 w_d^2)
	real spd_se(real alpha, vector rho, vector w) {
		int D = num_elements(w);
		real quad = dot_self(rho .* w);           // sum_d rho_d^2 w_d^2
		return square(alpha) * pow(2 * pi(), 0.5 * D) * prod(rho) * exp(-0.5 * quad);
	}

	// ARD Matern spectral density at frequency vector w (nu = 1.5 or 2.5).
	// S(w) = alpha^2 * C_nu * prod_d rho_d * (2 nu + sum_d rho_d^2 w_d^2)^{-(nu + D/2)}
	// with C_nu = 2^D pi^{D/2} Gamma(nu + D/2) (2 nu)^nu / Gamma(nu).
	// For D = 1 this equals the isotropic 1D Matern SPD used in sHSGP.stan.
	real spd_matern(real alpha, vector rho, vector w, real nu) {
		int D = num_elements(w);
		real quad = dot_self(rho .* w);
		real C = pow(2.0, D) * pow(pi(), 0.5 * D) * tgamma(nu + 0.5 * D)
		         * pow(2 * nu, nu) / tgamma(nu);
		return square(alpha) * C * prod(rho) * pow(2 * nu + quad, -(nu + 0.5 * D));
	}
}

data {
	int<lower=1> N;                    // sample size
	int<lower=1> D;                    // number of input dimensions
	matrix[N, D] inputs;               // input matrix (observed; noisy if latent)
	vector[N] y;                       // single output
	array[D] int<lower=1> M;           // no. of basis functions per input dim
	vector<lower=0>[D] L;              // boundary factor per input dim (L[d] = c[d] * S[d])
	// per-dimension bounds (latent-input prior support / centering midpoint)
	vector[D] x_min;
	vector[D] x_max;
	vector<lower=0>[D] s;              // measurement-error SD per input dim (latent model)
	// switches
	int<lower=0, upper=2> covfn;       // 0: SE, 1: Matern 3/2, 2: Matern 5/2
	int<lower=0, upper=1> latent;      // 1: inputs are latent (measured with error)
	int<lower=0, upper=1> is_deriv;    // 1: use derivative-kernel spectral density
	int<lower=1, upper=D> deriv_dim;   // input dim to differentiate wrt (used if is_deriv)
	int<lower=0, upper=1> rho_prior;   // 0: Normal; 1: InvGamma (length-scales)
	int<lower=0, upper=1> is_centered; // 1: shift each input dim by its midpoint
	// prior hyperparams (two-parameter families)
	vector[2] ls_param;                // length-scale prior params (shared across dims)
	vector[2] msd_param;               // marginal SD (alpha) prior
	vector[2] esd_param;               // error SD (sigma) prior
	vector[2] intc_y;                  // intercept prior
}

transformed data {
	int M_total = prod(M);
	array[M_total, D] int S_idx;       // multi-index set of univariate indices
	matrix[M_total, D] omega;          // per-dim frequencies sqrt(lambda_{d, S_idx[j,d]})
	vector[D] mid = 0.5 * (x_min + x_max);
	matrix[N, M_total] PHI_data;       // basis for fixed (non-latent) inputs

	// Enumerate all D-tuples (mixed-radix counting; dim 1 varies fastest).
	{
		array[D] int idx = rep_array(1, D);
		for (j in 1:M_total) {
			for (d in 1:D) {
				S_idx[j, d] = idx[d];
			}
			for (d in 1:D) {
				idx[d] += 1;
				if (idx[d] <= M[d]) break;
				idx[d] = 1;
			}
		}
	}
	// Per-dimension frequencies: sqrt(lambda) = m*pi / (2 L).
	for (j in 1:M_total) {
		for (d in 1:D) {
			omega[j, d] = (S_idx[j, d] * pi()) / (2 * L[d]);
		}
	}
	// Precompute PHI for the observed inputs (used when not latent).
	{
		matrix[N, D] xs;
		if (is_centered) {
			for (d in 1:D) xs[, d] = inputs[, d] - mid[d];
		} else {
			xs = inputs;
		}
		PHI_data = build_PHI(xs, S_idx, L);
	}
}

parameters {
	vector[M_total] beta;              // basis coefficients
	vector<lower=0>[D] rho;            // ARD length-scales (one per input dim)
	real<lower=0> alpha;              // GP marginal SD
	real<lower=0> sigma;             // error SD
	real intercept;
	// latent inputs on the unit box, scaled per-dim to [x_min, x_max]
	matrix<lower=0, upper=1>[latent ? N : 0, latent ? D : 0] x01;
}

transformed parameters {
	matrix[N, D] x;                    // inputs used by the GP (latent or observed)
	matrix[N, M_total] PHI;
	vector[M_total] diagSPD;
	vector[N] f;

	if (latent) {
		matrix[N, D] x_shifted;
		for (d in 1:D) {
			x[, d] = x_min[d] + (x_max[d] - x_min[d]) * x01[, d];
		}
		if (is_centered) {
			for (d in 1:D) x_shifted[, d] = x[, d] - mid[d];
		} else {
			x_shifted = x;
		}
		PHI = build_PHI(x_shifted, S_idx, L);
	} else {
		x = inputs;
		PHI = PHI_data;
	}

	// Spectral-density weights (ARD). Derivative process (wrt deriv_dim) uses
	// the derivative-kernel spectral density, i.e. base SPD times w_deriv_dim^2.
	for (j in 1:M_total) {
		vector[D] w = omega[j, ]';
		real S;
		if (covfn == 1) {
			S = spd_matern(alpha, rho, w, 1.5);
		} else if (covfn == 2) {
			S = spd_matern(alpha, rho, w, 2.5);
		} else {
			S = spd_se(alpha, rho, w);
		}
		if (is_deriv) {
			S = S * square(w[deriv_dim]);
		}
		diagSPD[j] = sqrt(S);
	}

	f = intercept + PHI * (diagSPD .* beta);
}

model {
	// Coefficients
	beta ~ normal(0, 1);
	// Length-scale prior (shared hyperparameters across input dims)
	if (rho_prior == 0) {
		rho ~ normal(ls_param[1], ls_param[2]);
	} else {
		rho ~ inv_gamma(ls_param[1], ls_param[2]);
	}
	// Marginal SD, error SD, intercept
	alpha ~ normal(msd_param[1], msd_param[2]);
	sigma ~ normal(esd_param[1], esd_param[2]);
	intercept ~ normal(intc_y[1], intc_y[2]);
	// Latent-input measurement model.
	// x01 has an implicit uniform(0, 1) prior from its [0, 1] constraint, so x
	// is uniform on [x_min, x_max] per dimension (linear map => constant Jacobian).
	if (latent) {
		for (d in 1:D) {
			inputs[, d] ~ normal(x[, d], s[d]);
		}
	}
	// Likelihood
	y ~ normal(f, sigma);
}
