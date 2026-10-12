# Core Usage Guide

This document explains what the MAR(P) core estimator solves and how to operate it in a new analysis.

## 1. Model estimated by the core

$$
Y_t = \sum_{\ell=1}^{P} A_\ell Y_{t-\ell}B_\ell^\top + E_t,
$$

with

$$
E_t\sim MN(0,\Sigma_A,\Sigma_B).
$$

The training data are passed as a three-dimensional R array:

```r
dim(Y)
# G N T
```

The core estimates $A_1,\ldots,A_P$, $B_1,\ldots,B_P$, $\Sigma_A$, $\Sigma_B$, and posterior inclusion probabilities used for sparse structure selection.

## 2. Which function should be used?

Use `fit_mar_pointmass()` when inactive coefficients or groups should be exactly zero. This covers CGS and GGS.

Use `fit_mar_ssvs()` for the continuous spike-and-slab SSVS baseline, where inactive coefficients are shrunk toward zero rather than fixed at zero.

## 3. CGS versus GGS

### CGS

```r
partition_A = NULL
partition_B = NULL
```

Each coefficient is updated individually.

### GGS

Provide a partition on one or both sides:

```r
partition_A = list(...)
partition_B = list(...)
```

The A-side and B-side dispatch independently.

Example:

```r
fit <- fit_mar_pointmass(
  Y = Y,
  P = 2L,
  partition_A = list(c(1, 2), c(3, 4)),
  partition_B = NULL
)
```

This uses group-wise updates on A and scalar updates on B.

## 4. Valid partition format

For a six-dimensional side:

```r
partition_A <- list(
  group_1 = c(1, 2, 3),
  group_2 = c(4),
  group_3 = c(5, 6)
)
```

A partition must be disjoint, nonempty, and contain each index exactly once.

## 5. Data preparation

The core requires a finite numeric array with exactly three dimensions.

Example:

```r
G <- 4L
N <- 3L
T <- 100L
Y <- array(rnorm(G * N * T), dim = c(G, N, T))
```

The function requires:

```text
T > P
```

The number of conditional observations is $T-P$.

## 6. Initial values

`A_init` can be `NULL`, one `G x G` matrix reused at every lag, or a list of `P` matrices.

`B_init` works analogously with `N x N` matrices.

Example:

```r
A_init <- diag(G)
B_init <- diag(N)
```

For MAR(3):

```r
A_init <- list(A1_init, A2_init, A3_init)
B_init <- list(B1_init, B2_init, B3_init)
```

`Sigma_A_init` and `Sigma_B_init` must be symmetric positive-definite matrices.

## 7. MCMC controls

Main controls:

```r
n_iter
burnin
warm_iter
```

The retained sample size is

$$
n_{post}=n_{iter}-burnin.
$$

`warm_iter` is distinct from burn-in. It is used before the formal MCMC run for warm-start updates.

## 8. PIP selection

For manuscript reproduction:

```r
pip_cut = 0.5
pip_rule = ">="
```

so a selection unit is retained when

$$
\mathrm{PIP}\ge 0.5.
$$

## 9. Rescaling

Because

$$
A_\ell Y_{t-\ell}B_\ell^\top
=
(cA_\ell)Y_{t-\ell}(B_\ell/c)^\top,
$$

MAR factors have a scale-identification issue.

Use:

```r
rescale = TRUE
```

The implementation balances Frobenius norms while preserving the bilinear product.

## 10. Whitening

Conditional updates whiten by the covariance associated with the opposite mode.

For manuscript-style reproduction in this implementation:

```r
whitening = "paper"
```

`whitening = "original"` is retained as a legacy/diagnostic option.

## 11. Reproducible seeds

The core exposes:

```r
seed
warm_seed
sampling_seed
```

Example:

```r
fit <- fit_mar_pointmass(
  Y = Y,
  P = 1L,
  seed = 1001L,
  warm_seed = 2001L,
  sampling_seed = 3001L
)
```

## 12. Multiple chains

```r
fit <- fit_mar_pointmass(
  Y = Y,
  P = 1L,
  n_chains = 4L,
  workers = 4L,
  seed = 1000L
)
```

When `n_chains > 1`, the outer result contains `method`, `P`, `n_chains`, `workers`, `chain_seeds`, and `chains`.

Single-chain fits are stored under:

```r
fit$chains$chain_1
fit$chains$chain_2
```

## 13. Posterior draws

Use:

```r
keep_draws = TRUE
```

Point-mass fits retain `A_chain`, `B_chain`, and `scale_chain`.

SSVS additionally retains `gamma_A_chain`, `gamma_B_chain`, `theta_A_chain`, and `theta_B_chain`.

## 14. Recommended workflow

```text
1. Arrange data as G x N x T.
2. Choose P.
3. Choose CGS, GGS, or SSVS.
4. Define partitions if using GGS.
5. Choose prior inclusion probabilities and slab/spike scales.
6. Run the estimator.
7. Inspect PIP_A and PIP_B.
8. Extract A_hat and B_hat.
9. Check posterior diagnostics if draws were retained.
10. Perform forecasting/evaluation outside the core.
```
