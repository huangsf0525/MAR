# API Reference: `fit_mar_ssvs()`

`fit_mar_ssvs()` implements the continuous spike-and-slab SSVS baseline for MAR(P).

## Function signature

```r
fit_mar_ssvs(
  Y, P = 1L,
  A_init = NULL, B_init = NULL,
  Sigma_A_init = NULL, Sigma_B_init = NULL,
  n_iter = 2250L, burnin = 750L, warm_iter = 10L,
  init_method = c("map", "gibbs", "none"),
  tau0_A = 0.05, tau1_A = 2,
  tau0_B = tau0_A, tau1_B = tau1_A,
  prior_A = list(diag = 0.4, off = 0.2),
  prior_B = list(diag = 0.4, off = 0.2),
  fixed_theta = TRUE,
  theta_A_init = NULL, theta_B_init = NULL,
  alpha_theta = 1, beta_theta = 1,
  warm_fixed_theta = TRUE,
  nu_A = NULL, nu_B = NULL,
  scale_init = 1, scale_df = NULL,
  eta_shape = 1, eta_rate = 1,
  pip_cut = 0.5, pip_rule = c(">", ">="),
  rescale = TRUE,
  scale_update = c("before_rescale", "after_rescale"),
  numerics = c("density", "log_odds"),
  coefficient_jitter = 0,
  whitening = c("original", "paper"),
  draw_initial_indicators = TRUE,
  seed = NULL, warm_seed = NULL, sampling_seed = NULL,
  keep_draws = TRUE, verbose_every = NULL,
  n_chains = 1L, workers = NULL, mcmc_clock = NULL
)
```

## Inputs

| Argument | Description |
|---|---|
| `Y` | Numeric `G x N x T` training array. |
| `P` | MAR lag order. |
| `A_init` | One A matrix or list of `P` matrices. |
| `B_init` | One B matrix or list of `P` matrices. |
| `Sigma_A_init` | Initial A covariance. |
| `Sigma_B_init` | Initial B covariance. |
| `n_iter` | Total MCMC sweeps. |
| `burnin` | Burn-in sweeps. |
| `warm_iter` | Warm-start sweeps. |
| `init_method` | `"map"`, `"gibbs"`, or `"none"`. |
| `tau0_A` | A-side spike standard deviation. |
| `tau1_A` | A-side slab standard deviation. |
| `tau0_B` | B-side spike standard deviation. |
| `tau1_B` | B-side slab standard deviation. |
| `prior_A` | A-side prior inclusion probabilities. |
| `prior_B` | B-side prior inclusion probabilities. |
| `fixed_theta` | Keep inclusion probability fixed. |
| `theta_A_init` | Initial A inclusion probability/probabilities. |
| `theta_B_init` | Initial B inclusion probability/probabilities. |
| `alpha_theta` | Beta-prior alpha for estimated theta. |
| `beta_theta` | Beta-prior beta for estimated theta. |
| `warm_fixed_theta` | Keep theta fixed during warm-up. |
| `nu_A` | Inverse-Wishart df for `Sigma_A`. |
| `nu_B` | Inverse-Wishart df for `Sigma_B`. |
| `scale_init` | Initial common scale. |
| `scale_df` | Scale-update df. |
| `eta_shape` | Gamma hyperprior shape. |
| `eta_rate` | Gamma hyperprior rate. |
| `pip_cut` | PIP cutoff. |
| `pip_rule` | `">"` or `">="`; match the threshold rule used by the run being reproduced. |
| `rescale` | Balance MAR factor scales. |
| `scale_update` | Update scale before/after rescaling. |
| `numerics` | `"density"` or `"log_odds"`. |
| `coefficient_jitter` | Optional covariance jitter. |
| `whitening` | `"original"` or `"paper"`. |
| `draw_initial_indicators` | Draw initial spike/slab indicators. |
| `seed` | Base seed. |
| `warm_seed` | Warm-up seed. |
| `sampling_seed` | Formal sampling seed. |
| `keep_draws` | Keep posterior draws. |
| `verbose_every` | Progress interval. |
| `n_chains` | Number of independent chains. |
| `workers` | PSOCK worker count. |
| `mcmc_clock` | Optional timer. |

## Prior

For coefficient $\phi$:

$$
\phi\mid\gamma
\sim
(1-\gamma)N(0,\tau_0^2)
+
\gamma N(0,\tau_1^2),
$$

with

$$
\gamma\sim\operatorname{Bernoulli}(\theta),
\qquad
\tau_1>\tau_0.
$$

The spike shrinks inactive coefficients toward zero but does not make them exactly zero.

## Outputs

| Output | Description |
|---|---|
| `method` | `"SSVS"`. |
| `P` | MAR lag order. |
| `A_hat` | Thresholded A estimate. |
| `B_hat` | Thresholded B estimate. |
| `PIP_A` | A-side PIPs. |
| `PIP_B` | B-side PIPs. |
| `A_mean_active` | A mean using only draws with `gamma = 1`. |
| `B_mean_active` | B mean using only draws with `gamma = 1`. |
| `A_mean_all` | A mean over all retained draws. |
| `B_mean_all` | B mean over all retained draws. |
| `Sigma_A_last` | Final `Sigma_A` draw. |
| `Sigma_B_last` | Final `Sigma_B` draw. |
| `scale_last` | Final scale draw. |
| `runtime_sec` | Formal MCMC-loop runtime. |
| `settings` | Stored settings. |
| `A_chain` | A posterior draws. |
| `B_chain` | B posterior draws. |
| `gamma_A_chain` | A indicator draws. |
| `gamma_B_chain` | B indicator draws. |
| `scale_chain` | Scale draws. |
| `theta_A_chain` | A theta draws. |
| `theta_B_chain` | B theta draws. |

`A_mean_active` estimates $E(A_{ij}\mid\gamma_{ij}=1,Y)$, while `A_mean_all` averages across both spike and slab states. The same distinction applies to B.
