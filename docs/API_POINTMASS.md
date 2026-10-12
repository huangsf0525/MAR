# API Reference: `fit_mar_pointmass()`

`fit_mar_pointmass()` implements the exact-zero point-mass sampler used for CGS and GGS.

## Function signature

```r
fit_mar_pointmass(
  Y, P = 1L, partition_A = NULL, partition_B = NULL,
  A_init = NULL, B_init = NULL,
  Sigma_A_init = NULL, Sigma_B_init = NULL,
  n_iter = 10000L, burnin = 5000L, warm_iter = 20L,
  tau_A = 4, tau_B = 4,
  prior_A = list(diag = 0.4, off = 0.2),
  prior_B = list(diag = 0.4, off = 0.2),
  nu_A = NULL, nu_B = NULL,
  scale_init = NULL, scale_df = NULL,
  eta_shape = 1, eta_rate = 1,
  pip_cut = 0.5, pip_rule = c(">", ">="),
  rescale = TRUE, random_lag_scan = FALSE,
  init_candidates = NULL,
  numerics = c("precision", "conditional"),
  singleton_update = c("auto", "scalar", "group_1x1"),
  whitening = c("original", "paper"),
  seed = NULL, warm_seed = NULL, sampling_seed = NULL,
  keep_draws = TRUE, verbose_every = NULL,
  n_chains = 1L, workers = NULL, mcmc_clock = NULL
)
```

## Inputs

| Argument | Description |
|---|---|
| `Y` | Numeric `G x N x T` training array. |
| `P` | MAR lag order; `1 <= P < T`. |
| `partition_A` | A-side partition; `NULL` gives scalar CGS updates. |
| `partition_B` | B-side partition; `NULL` gives scalar CGS updates. |
| `A_init` | One `G x G` matrix or a list of `P` matrices. |
| `B_init` | One `N x N` matrix or a list of `P` matrices. |
| `Sigma_A_init` | Initial positive-definite `G x G` covariance. |
| `Sigma_B_init` | Initial positive-definite `N x N` covariance. |
| `n_iter` | Total formal MCMC sweeps. |
| `burnin` | MCMC sweeps discarded as burn-in. |
| `warm_iter` | Warm-start sweeps before formal sampling. |
| `tau_A` | A-side slab scale; scalar or diagonal/off-diagonal specification. |
| `tau_B` | B-side slab scale. |
| `prior_A` | A-side prior inclusion probabilities. |
| `prior_B` | B-side prior inclusion probabilities. |
| `nu_A` | Inverse-Wishart df for `Sigma_A`; default `G + 2`. |
| `nu_B` | Inverse-Wishart df for `Sigma_B`; default `N + 2`. |
| `scale_init` | Initial common covariance-prior scale. |
| `scale_df` | Degrees of freedom used in scale update. |
| `eta_shape` | Gamma hyperprior shape. |
| `eta_rate` | Gamma hyperprior rate. |
| `pip_cut` | PIP threshold. |
| `pip_rule` | `">"` or `">="`; match the threshold rule used by the run being reproduced. |
| `rescale` | Frobenius balancing of factor pairs. |
| `random_lag_scan` | Randomize lag update order. |
| `init_candidates` | Optional candidate A/B initial states. |
| `numerics` | `"precision"` or `"conditional"`. |
| `singleton_update` | Singleton group update policy. |
| `whitening` | `"original"` or `"paper"`. |
| `seed` | Base random seed. |
| `warm_seed` | Warm-start seed. |
| `sampling_seed` | Formal sampling seed. |
| `keep_draws` | Keep posterior chains. |
| `verbose_every` | Progress-print interval. |
| `n_chains` | Number of independent chains. |
| `workers` | PSOCK worker count. |
| `mcmc_clock` | Optional external timer interface. |

## CGS

```r
fit <- fit_mar_pointmass(
  Y = Y,
  P = 1L,
  partition_A = NULL,
  partition_B = NULL
)
```

## GGS

```r
fit <- fit_mar_pointmass(
  Y = Y,
  P = 1L,
  partition_A = list(c(1, 2), c(3, 4)),
  partition_B = list(c(1, 2), c(3))
)
```

A and B dispatch independently.

## Outputs

| Output | Description |
|---|---|
| `method` | `"CGS"` if both partitions are `NULL`; otherwise `"GGS"`. |
| `P` | MAR lag order. |
| `A_hat` | Selected A estimates; list of `P` matrices. |
| `B_hat` | Selected B estimates; list of `P` matrices. |
| `PIP_A` | A-side PIPs. |
| `PIP_B` | B-side PIPs. |
| `Sigma_A_last` | Final draw of `Sigma_A`. |
| `Sigma_B_last` | Final draw of `Sigma_B`. |
| `scale_last` | Final scale draw. |
| `warm_start_index` | Selected initialization candidate index. |
| `warm_scores` | Warm-start residual scores. |
| `runtime_sec` | Formal MCMC-loop runtime in seconds. |
| `settings` | Stored estimator settings. |
| `A_chain` | A posterior draws if retained. |
| `B_chain` | B posterior draws if retained. |
| `scale_chain` | Scale posterior draws if retained. |

If `n_post = n_iter - burnin`, then:

```text
A_chain : n_post x P x G x G
B_chain : n_post x P x N x N
```

For publication selection:

```r
pip_cut = 0.5
pip_rule = ">="
```
