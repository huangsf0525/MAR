# Supplement S2

Install `MASS`, `MCMCpack`, `mcmcse`, and `Matrix` first. Place the core file at `R/mar_estimators.R`. Start a fresh R session for this experiment from the repository root, then run both R blocks in order. The first creates 100 training arrays; the second fits each array.

## Data generation

```r
s2_cfg <- list(
  m = 24L, n = 28L, T_total = 300L, T_train = 160L,
  simulation_burn = 200L, truth_seed = 1L, factor_radius_cap = 0.98,
  data_seed_base = 0L
)

s2_group_cov <- function(d, core, rho_core, rho_rest) {
  rest <- setdiff(seq_len(d), core)
  block <- function(k, rho, scales) {
    diag(scales, k) %*% toeplitz(rho^(0:(k - 1L))) %*% diag(scales, k)
  }
  core_block <- block(
    length(core), rho_core,
    seq(1.35, 1.05, length.out = length(core))
  )
  rest_block <- block(
    length(rest), rho_rest,
    seq(.95, .75, length.out = length(rest))
  )
  out <- matrix(0, d, d)
  out[core, core] <- core_block
  out[rest, rest] <- rest_block
  v <- numeric(d)
  v[core] <- 1
  v[rest] <- .4
  v <- v / sqrt(sum(v^2))
  out <- out + (.06 * mean(diag(out))) * tcrossprod(v)
  out + 1e-10 * diag(d)
}

s2_truth <- function(cfg) {
  set.seed(cfg$truth_seed)
  A <- matrix(0, cfg$m, cfg$m)
  B <- matrix(0, cfg$n, cfg$n)
  ca <- 1:6
  ra <- 7:24
  cb <- 1:7
  rb <- 8:28
  A[ca, ca] <- 0.12
  diag(A)[ca] <- 0.42
  B[cb, cb] <- 0.10
  diag(B)[cb] <- 0.45
  qa <- c(.18, .15, .12, .10, .08, .06)
  qb <- c(.16, .14, .12, .10, .08, .06, .05)
  diag(A)[ra] <- .30
  diag(B)[rb] <- .28
  A[cbind(8:24, 7:23)] <- .05
  B[cbind(9:28, 8:27)] <- .04
  for (i in ra) A[i, ca] <- qa * (.90 + .20 * runif(length(ca)))
  raw_radius_A <- max(Mod(eigen(A, only.values = TRUE)$values))
  if (raw_radius_A >= cfg$factor_radius_cap) {
    A <- A * (cfg$factor_radius_cap / raw_radius_A)
  }
  for (i in rb) B[i, cb] <- qb * (.90 + .20 * runif(length(cb)))
  raw_radius_B <- max(Mod(eigen(B, only.values = TRUE)$values))
  if (raw_radius_B >= cfg$factor_radius_cap) {
    B <- B * (cfg$factor_radius_cap / raw_radius_B)
  }
  list(
    A = A, B = B,
    Sigma_A = s2_group_cov(cfg$m, ca, .75, .55),
    Sigma_B = s2_group_cov(cfg$n, cb, .70, .55)
  )
}

s2_simulate <- function(truth, cfg, rep_id) {
  set.seed(cfg$data_seed_base + rep_id)
  Y <- array(0, c(cfg$m, cfg$n, cfg$simulation_burn + cfg$T_total))
  LA <- t(chol(truth$Sigma_A))
  LB <- t(chol(truth$Sigma_B))
  Y[, , 1L] <- LA %*% matrix(rnorm(cfg$m * cfg$n), cfg$m, cfg$n) %*% t(LB)
  for (t in 2:dim(Y)[3L]) {
    Y[, , t] <- truth$A %*% Y[, , t - 1L] %*% t(truth$B) +
      LA %*% matrix(rnorm(cfg$m * cfg$n), cfg$m, cfg$n) %*% t(LB)
  }
  Y[, , cfg$simulation_burn + seq_len(cfg$T_total), drop = FALSE]
}

truth <- s2_truth(s2_cfg)
Y_train_list <- lapply(seq_len(100L), function(rep_id) {
  Y <- s2_simulate(truth, s2_cfg, rep_id)
  Y[, , 1:160, drop = FALSE]
})
```

## Core-function inputs

```r
source("R/mar_estimators.R")
fits_by_rep <- lapply(seq_len(100L), function(rep_id) {
  Y_train <- Y_train_list[[rep_id]]

  fit_cgs <- fit_mar_pointmass(
    Y = Y_train, P = 1L,
    partition_A = NULL, partition_B = NULL,
    A_init = diag(24), B_init = diag(28),
    Sigma_A_init = diag(24), Sigma_B_init = diag(28),
    n_iter = 10000L, burnin = 5000L, warm_iter = 20L,
    tau_A = 4, tau_B = 4,
    prior_A = list(diag = 0.4, off = 0.2),
    prior_B = list(diag = 0.4, off = 0.2),
    nu_A = 26, nu_B = 30,
    scale_init = 26, scale_df = c(26, 30),
    eta_shape = 1, eta_rate = 1,
    pip_cut = 0.5, pip_rule = ">",
    rescale = TRUE, random_lag_scan = FALSE,
    init_candidates = NULL, numerics = "precision",
    singleton_update = "auto", whitening = "original",
    seed = 100000L + rep_id, warm_seed = NULL, sampling_seed = NULL,
    keep_draws = TRUE, verbose_every = NULL,
    n_chains = 1L, workers = NULL, mcmc_clock = NULL
  )

  fit_ggs <- fit_mar_pointmass(
    Y = Y_train, P = 1L,
    partition_A = c(list(1:6), as.list(7:24)),
    partition_B = c(list(1:7), as.list(8:28)),
    A_init = diag(24), B_init = diag(28),
    Sigma_A_init = diag(24), Sigma_B_init = diag(28),
    n_iter = 10000L, burnin = 5000L, warm_iter = 20L,
    tau_A = 4, tau_B = 4,
    prior_A = list(diag = 0.4, off = 0.2),
    prior_B = list(diag = 0.4, off = 0.2),
    nu_A = 26, nu_B = 30,
    scale_init = 26, scale_df = c(26, 30),
    eta_shape = 1, eta_rate = 1,
    pip_cut = 0.5, pip_rule = ">",
    rescale = TRUE, random_lag_scan = FALSE,
    init_candidates = NULL, numerics = "precision",
    singleton_update = "scalar", whitening = "original",
    seed = 100000L + rep_id, warm_seed = NULL, sampling_seed = NULL,
    keep_draws = TRUE, verbose_every = NULL,
    n_chains = 1L, workers = NULL, mcmc_clock = NULL
  )

  list(CGS = fit_cgs, GGS = fit_ggs)
})
```
