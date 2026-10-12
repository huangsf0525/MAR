# Efficient Bayesian Structure Selection for Matrix Autoregressive Models

Documentation and local-run result summaries for Bayesian structure-selection methods for matrix autoregressive (MAR) models.

This extracted folder contains documentation with local-run result summaries. The R source and original result CSV files are stored separately on this machine. The reusable estimator core is:

```text
MAR_paper_core_standalone.R
```

It provides three Bayesian structure-selection strategies through two public R functions:

- **CGS** — componentwise Gibbs sampling with an exact point mass at zero;
- **GGS** — group-wise Gibbs sampling with exact-zero group selection;
- **SSVS** — continuous spike-and-slab stochastic search variable selection.

The simulation programs used in the paper are kept separate from the estimator interface so the core can be reused for new data sets, alternative lag orders, and other experiments.

## What problem does the core solve?

For a matrix-valued time series $Y_t\in\mathbb{R}^{G\times N}$, the MAR(P) model is

$$
Y_t = \sum_{\ell=1}^{P} A_\ell Y_{t-\ell}B_\ell^\top + E_t,
$$

where

$$
A_\ell\in\mathbb{R}^{G\times G},\qquad
B_\ell\in\mathbb{R}^{N\times N},
$$

and

$$
E_t\sim MN(0,\Sigma_A,\Sigma_B).
$$

The core estimates sparse row-side and column-side dynamics and returns posterior inclusion probabilities (PIPs) for structure selection.

It is intended for questions such as:

- Which entries of $A_\ell$ are active?
- Which entries of $B_\ell$ are active?
- Should off-diagonal coefficients be selected individually or in prespecified groups?
- How can MAR(P) structure be estimated without repeatedly sampling one very large coefficient vector?

## Public functions

Run the examples from a directory containing the core file, then source it:

```r
source("MAR_paper_core_standalone.R")
```

This defines:

```r
fit_mar_pointmass()
fit_mar_ssvs()
```

### `fit_mar_pointmass()`

Use for exact-zero point-mass selection.

- `partition_A = NULL`, `partition_B = NULL` -> **CGS**.
- Supplying one or both partitions -> **GGS** on those sides.

### `fit_mar_ssvs()`

Use for the continuous spike-and-slab **SSVS** baseline.

Inactive coefficients are shrunk toward zero through the spike distribution rather than fixed exactly at zero.

## Minimal input

Both functions require a numeric training array:

```text
G x N x T
```

and a lag order `P` with

```text
1 <= P < T
```

Example:

```r
dim(Y)
# 20 20 90

P <- 1L
```

## Quick start

### CGS

```r
fit_cgs <- fit_mar_pointmass(
  Y = Y,
  P = 1L,
  partition_A = NULL,
  partition_B = NULL,
  n_iter = 2250L,
  burnin = 750L,
  pip_cut = 0.5,
  pip_rule = ">=",
  whitening = "paper",
  rescale = TRUE,
  seed = 123L
)
```

### GGS

```r
groups_A <- list(
  group_1 = c(1, 2),
  group_2 = c(3, 4)
)

fit_ggs <- fit_mar_pointmass(
  Y = Y,
  P = 1L,
  partition_A = groups_A,
  partition_B = NULL,
  n_iter = 10000L,
  burnin = 5000L,
  pip_cut = 0.5,
  pip_rule = ">=",
  whitening = "paper",
  rescale = TRUE,
  seed = 123L
)
```

### SSVS

```r
fit_ssvs <- fit_mar_ssvs(
  Y = Y,
  P = 1L,
  tau0_A = 0.05,
  tau1_A = 2,
  tau0_B = 0.05,
  tau1_B = 2,
  n_iter = 2250L,
  burnin = 750L,
  pip_cut = 0.5,
  pip_rule = ">=",
  whitening = "paper",
  rescale = TRUE,
  seed = 123L
)
```

## Main outputs

For both estimators:

```r
fit$A_hat
fit$B_hat
fit$PIP_A
fit$PIP_B
```

Even when `P = 1`, these are lists of length `P`.

For MAR(3):

```r
fit$A_hat[[1]]  # A1
fit$A_hat[[2]]  # A2
fit$A_hat[[3]]  # A3
```

If `keep_draws = TRUE`, posterior coefficient chains have dimensions:

```text
retained_draw x lag x row x column
```

## Documented local runs

The numerical results in [Simulation 1](./docs/SIMULATION1.md) and [Simulation 2](./docs/SIMULATION2.md) come from saved local core-based runs with 10 replications per method. Simulation 2 uses the sequential run without replication-level parallelism. The original local output paths and run identifiers are recorded in the [reproducibility notes](./docs/REPRODUCIBILITY.md). These are not 100-replication manuscript results.

The local runners use a PIP cutoff of 0.5. Their exact threshold operators differ; see the simulation pages and [reproducibility notes](./docs/REPRODUCIBILITY.md). The examples above show one valid API configuration, not the full saved run configuration.

## Folder layout

```text
.
README.md
docs/
  INDEX.md
  CORE_USAGE.md
  API_POINTMASS.md
  API_SSVS.md
  METHODS.md
  INSTALLATION.md
  SIMULATION1.md
  SIMULATION2.md
  REPRODUCIBILITY.md
```

## Documentation

- [Core usage guide](./docs/CORE_USAGE.md)
- [Point-mass API: CGS/GGS](./docs/API_POINTMASS.md)
- [SSVS API](./docs/API_SSVS.md)
- [Method notes](./docs/METHODS.md)
- [Installation](./docs/INSTALLATION.md)
- [Simulation 1](./docs/SIMULATION1.md)
- [Simulation 2](./docs/SIMULATION2.md)
- [Reproducibility](./docs/REPRODUCIBILITY.md)

## Dependencies

Core estimators:

```text
R
MASS
MCMCpack
```

Local simulation runners additionally use packages such as:

```text
Matrix
mcmcse
```

## Scope of the core file

`MAR_paper_core_standalone.R` contains estimation routines only. It does not contain data-generating processes, Monte Carlo replication loops, forecast evaluation, simulation tables, file-output management, or empirical OECD preprocessing.

## Citation

If this code is used in academic work, cite the associated manuscript:

> Po-Syu Chen, Yi-Siang Ma, Ray-Bing Chen, and Shih-Feng Huang.  
> *Efficient Bayesian Structure Selection for Matrix Autoregressive Models.*

Add final journal/DOI/preprint metadata before archival release.
