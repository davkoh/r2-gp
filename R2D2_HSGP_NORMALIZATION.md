# R2D2 additive HSGP normalization

<style>
mjx-container[display="true"] {
  display: block;
  max-width: 100%;
  overflow-x: auto;
  overflow-y: hidden;
  padding-block: 0.2rem;
}
</style>

## Purpose

`Stanmodels/aHSGP_R2D2.stan` models the predictor as a sum of $D$ univariate
exponential-kernel HSGPs,

$$
f(x_i) = \sum_{d=1}^D f_d(x_{id}),
$$

where $x_i=(x_{i1},\ldots,x_{iD})$ is the covariate vector for observation $i$
and $f_d$ is the function associated with covariate $d$. There is no
tensor-product or cross-product kernel. Under the R2D2 prior, a single global
signal variance is decomposed among the additive functions:

$$
R^2 \sim \operatorname{Beta}(\mu_R\phi_R,(1-\mu_R)\phi_R), \qquad
\alpha^2 = \sigma^2\frac{R^2}{1-R^2}, \qquad
\alpha_d^2 = \pi_d\alpha^2,
$$

where $\sigma^2$ is the observation-error variance, $\alpha^2$ is the global
signal variance, and $\pi=(\pi_1,\ldots,\pi_D)$ has a symmetric Dirichlet prior.

## Why normalization is needed

Fix an additive component $d$ and evaluate it at the observed values
$x_{1d},\ldots,x_{Nd}$, where $N$ is the number of observations. Let $M$ be
the number of HSGP basis functions and define:

- $\Phi_d\in\mathbb{R}^{N\times M}$ as the basis matrix, with entry
  $[\Phi_d]_{im}=\phi_{dm}(x_{id})$, corresponding to `PHI[d]`;
- $\phi_{dm}\in\mathbb{R}^N$ as column $m$ of $\Phi_d$;
- $\omega_{dm}$ as the frequency associated with basis function $m$,
  corresponding to `omega[m, d]`;
- $\rho_d>0$ as the kernel length scale for component $d$, corresponding to
  `rho[d]`;
- $S_{dm}=S(\omega_{dm},\rho_d)>0$ as the exponential-kernel spectral density,
  stored as `spectral_density[m]` inside the loop for component `d`;
- $S_d=(S_{d1},\ldots,S_{dM})^{\mathsf T}$ as the vector of spectral densities;
- $z_d\in\mathbb{R}^M$ as independent standard-normal coefficients, so
  $z_d\sim N(0,I_M)$, corresponding to `z[d]` in Stan.

Here $I_M$ is the $M\times M$ identity matrix and $\operatorname{diag}(S_d)$
is the diagonal matrix whose diagonal is $S_d$.

### Step 1: obtain the finite-design HSGP covariance

The unscaled HSGP values over the observed design form the vector

$$
g_d=\Phi_d\operatorname{diag}(\sqrt{S_d})z_d,
$$

where the square root is applied elementwise. Because $E_z[z_d]=0$ and
$E_z[z_dz_d^{\mathsf T}]=I_M$, the conditional mean and covariance of $g_d$
given the inputs and $\rho_d$ are

$$
E_z[g_d]=0,
\qquad
\widetilde K_d
=E_z[g_dg_d^{\mathsf T}]
=\Phi_d\operatorname{diag}(S_d)\Phi_d^{\mathsf T}.
$$

The matrix $\widetilde K_d\in\mathbb{R}^{N\times N}$ is the finite-basis HSGP
approximation to the kernel matrix. Kernel normalization states that
$K_d(x,x)=1$, so the diagonal elements of $\widetilde K_d$ are approximately
one away from boundary and basis-truncation effects. This controls variance at
each individual input. It does not imply that a realized function has sample
variance one across a finite, correlated collection of inputs.

### Step 2: write sample variance as a quadratic form

Let $\mathbf{1}_N\in\mathbb{R}^N$ be a vector of ones, let $I_N$ be the
$N\times N$ identity matrix, and define the centering matrix

$$
P_N=I_N-\frac{1}{N}\mathbf{1}_N\mathbf{1}_N^{\mathsf T}.
$$

The matrix $P_N$ subtracts a vector's sample mean and satisfies
$P_N^{\mathsf T}=P_N$ and $P_N^2=P_N$. The sample variance of $g_d$, using the
usual denominator $N-1$, is

$$
\widehat{\operatorname{Var}}_N(g_d)
=\frac{1}{N-1}\lVert P_Ng_d\rVert^2
=\frac{1}{N-1}g_d^{\mathsf T}P_Ng_d.
$$

Here $\lVert\cdot\rVert$ is the Euclidean norm, the superscript $\mathsf T$
denotes transpose, and the subscript $N$ distinguishes finite-sample variance
from pointwise GP variance.

### Step 3: take the expectation over HSGP coefficients

For a mean-zero random vector $u$ with covariance matrix $K$ and a fixed matrix
$A$, the quadratic-form identity is

$$
E[u^{\mathsf T}Au]=\operatorname{tr}(AK).
$$

Applying this identity conditional on the complete input matrix
$X\in\mathbb{R}^{N\times D}$ and on $\rho_d$ gives

$$
\begin{aligned}
Q_d(\rho_d)
&:=E_z\left[
  \widehat{\operatorname{Var}}_N(g_d)\mid X,\rho_d
\right] \\
&=\frac{1}{N-1}\operatorname{tr}(P_N\widetilde K_d) \\
&=\frac{1}{N-1}\operatorname{tr}\left(
  P_N\Phi_d\operatorname{diag}(S_d)\Phi_d^{\mathsf T}
\right).
\end{aligned}
$$

In this expression, $E_z$ means expectation only over $z_d$, and
$\operatorname{tr}(A)$ is the sum of the diagonal elements of a square matrix
$A$. The equivalent symmetric expression is
$\operatorname{tr}(P_N\widetilde K_dP_N)/(N-1)$ because $P_N^2=P_N$ and trace
is invariant under cyclic permutations.

### Step 4: reduce the trace to the quantity used in code

Because $\operatorname{diag}(S_d)$ is diagonal, the trace separates over the
$M$ basis functions:

$$
\begin{aligned}
Q_d(\rho_d)
&=\frac{1}{N-1}\sum_{m=1}^M
  S_{dm}\,\phi_{dm}^{\mathsf T}P_N\phi_{dm} \\
&=\sum_{m=1}^M S_{dm}\,
  \widehat{\operatorname{Var}}_N(\phi_{dm}).
\end{aligned}
$$

Thus, the implementation only needs the sample variance of each basis column
and its spectral density. It does not construct the $N\times N$ matrix
$\widetilde K_d$. In Stan, $Q_d(\rho_d)$ is
`hsgp_variance_normalizer[d]`, and
$\widehat{\operatorname{Var}}_N(\phi_{dm})$ is
`basis_sample_variance[m, d]`.

### Step 5: show why the unnormalized value is generally not one

Expanding $P_N$ gives

$$
Q_d(\rho_d)
=\frac{1}{N-1}\left\{
  \operatorname{tr}(\widetilde K_d)
  -\frac{1}{N}\mathbf{1}_N^{\mathsf T}
   \widetilde K_d\mathbf{1}_N
\right\}.
$$

The first term sums the $N$ pointwise variances. The second term is $N$ times
the prior variance of the sample mean of $g_d$. For a positively correlated
GP, this second term is positive: variation shared across all input locations
looks like a nearly constant shift and is removed by sample centering.
Consequently, $Q_d(\rho_d)$ can be substantially below one even when every
diagonal element of the exact kernel is one. Finite HSGP truncation and boundary
effects can also move the first term away from $N$.

### Step 6: normalize the component

The model scales the component with its R2D2 amplitude $\alpha_d$, represented
by `alpha[d]` in Stan, as

$$
f_d=\frac{\alpha_d}{\sqrt{Q_d(\rho_d)}}g_d.
$$

Since $Q_d$ depends on the design and length scale but not on the realized
coefficient vector, it follows that

$$
E_z\left[
  \widehat{\operatorname{Var}}_N(f_d)
  \mid X,\rho_d,\alpha_d
\right]
=\alpha_d^2.
$$

The coefficient vectors are independent across components, so their expected
sample cross-covariances are zero. Therefore,

$$
E_z\left[
  \widehat{\operatorname{Var}}_N\left(\sum_{d=1}^D f_d\right)
  \mathrel{\Big|}X,\rho,\alpha
\right]
=\sum_{d=1}^D\alpha_d^2
=\alpha^2,
$$

where $\rho=(\rho_1,\ldots,\rho_D)$ and
$\alpha=(\alpha_1,\ldots,\alpha_D)$. This is a design-dependent prior
normalization, not a post-sampling rescaling. Individual realized functions can
still have sample variance above or below their allocated $\alpha_d^2$.

### Mapping to the Stan coefficients

The additive Stan model uses a non-centered parameterization. Its sampled
coefficient vector is

$$
z_d\sim N(0,I_M).
$$

This vector is represented by `z[d]` in Stan. Thus, `z[d]` has the same role
that the standard-normal `beta` vector has in `sHSGPmi.stan`. The normalization
does not change the prior variance of `z[d]`; that variance remains one. Stan
instead constructs the scale vector

$$
c_d=\alpha_d\sqrt{\frac{S_d}{Q_d(\rho_d)}},
$$

represented by `spectral_scale`, where the square root and division are
elementwise. The effective HSGP basis coefficients are

$$
\beta_d^{\mathrm{eff}}=c_d\odot z_d,
$$

where $\odot$ denotes elementwise multiplication. This is the Stan expression

```stan
spectral_scale .* z[d]
```

The component evaluated at the observed inputs is

$$
f_d=\Phi_d\beta_d^{\mathrm{eff}}.
$$

This final vector is represented by `f_component[, d]`.

Equivalently, effective coefficient $m$ has the centered prior

$$
\beta_{dm}^{\mathrm{eff}}
\sim N\left(0,\frac{\alpha_d^2S_{dm}}{Q_d(\rho_d)}\right).
$$

Therefore, the normalization scales the variance of the effective basis
coefficients. Computationally, however, the model keeps the sampled `z`
coefficients standard normal and applies that scale through the deterministic
`spectral_scale` transformation.

### Effect on HSGP convergence

Let the subscript $M$ make the finite number of basis functions explicit. The
usual unnormalized HSGP covariance is

$$
\widetilde K_{d,M}
=\Phi_{d,M}\operatorname{diag}(S_{d,M})\Phi_{d,M}^{\mathsf T}.
$$

With a sufficiently large boundary and increasing $M$, the standard HSGP
approximation satisfies

$$
\widetilde K_{d,M}\longrightarrow K_d,
$$

where $K_d$ is the intended exponential-kernel matrix on the observed inputs.
The normalized covariance used by this model is

$$
K_{d,M}^{\mathrm{norm}}
=\frac{\alpha_d^2}{Q_{d,M}}\widetilde K_{d,M},
$$

where

$$
Q_{d,M}
=\frac{1}{N-1}\operatorname{tr}
  \left(P_N\widetilde K_{d,M}\right).
$$

As $M$ increases, the normalizer converges with the HSGP covariance:

$$
\begin{aligned}
Q_{d,M}&\longrightarrow Q_{d,\infty}, \\
Q_{d,\infty}&=\frac{1}{N-1}\operatorname{tr}(P_NK_d).
\end{aligned}
$$

It follows that

$$
K_{d,M}^{\mathrm{norm}}
\longrightarrow
K_{d,\infty}^{\mathrm{norm}},
\qquad
K_{d,\infty}^{\mathrm{norm}}
=\frac{\alpha_d^2}{Q_{d,\infty}}K_d.
$$

The normalization therefore does not disrupt convergence or change the
relative spectral weighting: all frequencies are multiplied by the same
scalar. The correlation shape and smoothness controlled by `rho[d]` are
preserved. What changes is the limiting amplitude convention. Instead of
interpreting `alpha[d]^2` as pointwise kernel variance, the normalized model
interprets it as expected sample variance over the supplied `inputs[, d]`.
The limiting process is a scalar-rescaled version of the conventional GP unless
$Q_{d,\infty}=1$.

## R2 prior interface

The Stan model uses mean and precision rather than direct Beta shape values:

- `r2_mean` is $\mu_R$, constrained to $(0,1)$.
- `r2_precision` is $\phi_R>0$.

For example, `r2_mean = 1.0 / 3.0` and `r2_precision = 3.0` imply
$R^2\sim\operatorname{Beta}(1,2)$. The R test script stores the same pair in
`r2_shape` for historical reasons and converts it to the two Beta shapes when
calling `rbeta`.

## Length-scale options

In Stan:

- Set `estimate_rho = 1` to estimate one length scale per function using
  `lengthscale_inv_gamma` as the inverse-gamma shape and scale.
- Set `estimate_rho = 0` to remove the length scales from the parameter space
  and use the positive values supplied in `rho_fixed`.

In the R test script, `sample_rho` selects the equivalent behavior and
`fixed_rho` supplies the common fixed value. For example:

```r
sample_rho <- FALSE
fixed_rho <- 1
```

## Interpretation and limitations

The normalization matches sample variance in expectation over HSGP coefficient
draws. It does not force every realized function to have variance exactly
$\alpha_d^2$. Consequently, empirical R2 draws still vary around the target
Beta draws, particularly with only one or two additive functions. The agreement
usually improves as independent components are added because realization-level
variance fluctuations average out.

Exact draw-by-draw matching would require rescaling each realized function by
its observed sample standard deviation. That would make the prior depend
nonlinearly on its latent coefficient draw and would no longer be the intended
Gaussian process prior.

Increasing the number of observations reduces sample noise, and increasing the
number of basis functions reduces HSGP truncation error. Neither replaces the
design normalization. The normalization should be recomputed whenever inputs,
boundaries, basis count, or length scales change; both implementations do this
automatically.

## Prior-predictive check

Run:

```sh
Rscript aHSGP_R2D2_prior_test.R
```

The script writes `aHSGP_R2D2_prior_density.png`. The dashed black curve is the
sampled target Beta prior. Colored curves are empirical R2 distributions from
realized additive HSGPs. The two facets use independent y-axis scales so the
shape within each prior family remains visible.