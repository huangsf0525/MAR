# Simulation 1

Install `MASS`, `MCMCpack`, `mcmcse`, and `Matrix` first. Place the core file at `R/mar_estimators.R`. Start a fresh R session for this experiment from the repository root, then run both R blocks in order. The first creates 100 training arrays; the second fits each array.

## Data generation

```r
safe_chol <- function(S, jitter_seq = c(0, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2)) {
  p <- nrow(S)
  S <- 0.5 * (S + t(S))
  for (jj in jitter_seq) {
    out <- tryCatch(chol(S + jj * diag(p)), error = function(e) NULL)
    if (!is.null(out)) {
      return(out)
    }
  }
  stop("safe_chol() failed: matrix is not positive definite.")
}

spectral_radius <- function(M) {
  max(Mod(eigen(M, only.values = TRUE)$values))
}

rmatnorm_simple <- function(G, N, Sigma_L, Sigma_R) {
  Z <- matrix(rnorm(G * N), G, N)
  L <- t(safe_chol(Sigma_L))
  R <- t(safe_chol(Sigma_R))
  L %*% Z %*% t(R)
}

make_truth_matrices <- function(G, N,
                                target_radius = 0.85,
                                diag_A = 0.35,
                                diag_B = 0.35,
                                off_value = 0.18,
                                seed = 123) {
  set.seed(seed)

  A0 <- matrix(0, G, G)
  B0 <- matrix(0, N, N)

  diag(A0) <- seq(diag_A, diag_A * 0.45, length.out = G)
  diag(B0) <- seq(diag_B, diag_B * 0.45, length.out = N)

  if (G >= 4) {
    A0[1, 2] <- off_value
    A0[2, 1] <- -0.80 * off_value
    A0[3, 4] <- 0.70 * off_value
  }
  if (G >= 6) A0[5, 6] <- -0.65 * off_value
  if (G >= 8) A0[7, 8] <- 0.55 * off_value
  if (G >= 10) A0[9, 10] <- -0.50 * off_value
  if (G >= 12) A0[11, 12] <- 0.45 * off_value

  if (N >= 4) {
    B0[1, 3] <- 0.90 * off_value
    B0[3, 1] <- -0.75 * off_value
    B0[2, 4] <- 0.65 * off_value
  }
  if (N >= 6) {
    B0[4, 5] <- 0.60 * off_value
    B0[6, 2] <- -0.50 * off_value
  }
  if (N >= 8) B0[7, 8] <- 0.45 * off_value
  if (N >= 10) B0[10, 9] <- -0.40 * off_value
  if (N >= 12) B0[11, 12] <- 0.35 * off_value

  rho_phi <- spectral_radius(kronecker(B0, A0))
  if (is.finite(rho_phi) && rho_phi > 0) {
    sf <- sqrt(target_radius / rho_phi)
    A0 <- sf * A0
    B0 <- sf * B0
  }

  rownames(A0) <- colnames(A0) <- paste0("R", seq_len(G))
  rownames(B0) <- colnames(B0) <- paste0("C", seq_len(N))

  list(
    A0 = A0,
    B0 = B0
  )
}

simulate_mar1 <- function(A0, B0, T_total = 160, burn_sim = 100,
                          Sigma_L = NULL, Sigma_R = NULL, seed = 1234) {
  set.seed(seed)

  G <- nrow(A0)
  N <- nrow(B0)
  if (is.null(Sigma_L)) Sigma_L <- diag(G)
  if (is.null(Sigma_R)) Sigma_R <- diag(N)

  T_all <- T_total + burn_sim
  Y <- array(0, dim = c(G, N, T_all))

  for (tt in 2:T_all) {
    E_t <- rmatnorm_simple(G, N, Sigma_L, Sigma_R)
    Y[, , tt] <- A0 %*% Y[, , tt - 1] %*% t(B0) + E_t
  }

  Y <- Y[, , (burn_sim + 1):T_all, drop = FALSE]
  dimnames(Y) <- list(rownames(A0), rownames(B0), paste0("t", seq_len(T_total)))
  Y
}

truth <- make_truth_matrices(20L, 20L, target_radius = 0.85, seed = 101L)
Y_train_list <- lapply(seq_len(100L), function(rep_id) {
  Y <- simulate_mar1(truth$A0, truth$B0,
    T_total = 120L,
    burn_sim = 100L, Sigma_L = diag(20), Sigma_R = diag(20),
    seed = 101000L + rep_id
  )
  Y[, , 1:90, drop = FALSE]
})
```

## Core-function inputs

```r
source("R/mar_estimators.R")
fits_by_rep <- lapply(seq_len(100L), function(rep_id) {
  Y_train <- Y_train_list[[rep_id]]

  fit_cgs <- fit_mar_pointmass(
    Y = Y_train, P = 1L, partition_A = NULL, partition_B = NULL,
    A_init = diag(20), B_init = diag(20),
    Sigma_A_init = diag(20), Sigma_B_init = diag(20),
    n_iter = 2250L, burnin = 750L, warm_iter = 20L,
    tau_A = 2, tau_B = 2,
    prior_A = list(diag = 0.4, off = 0.2),
    prior_B = list(diag = 0.4, off = 0.2),
    nu_A = 22, nu_B = 22,
    scale_init = 20, scale_df = c(22, 22),
    eta_shape = 1, eta_rate = 1,
    pip_cut = 0.5, pip_rule = ">", rescale = TRUE,
    random_lag_scan = FALSE, init_candidates = NULL,
    numerics = "precision", singleton_update = "auto",
    whitening = "original",
    seed = 301000L + rep_id, warm_seed = NULL, sampling_seed = NULL,
    keep_draws = TRUE, verbose_every = NULL,
    n_chains = 1L, workers = NULL, mcmc_clock = NULL
  )

  fit_ssvs <- fit_mar_ssvs(
    Y = Y_train, P = 1L, A_init = diag(20), B_init = diag(20),
    Sigma_A_init = diag(20), Sigma_B_init = diag(20),
    n_iter = 2250L, burnin = 750L, warm_iter = 10L,
    init_method = "map", tau0_A = 0.05, tau1_A = 2,
    tau0_B = 0.05, tau1_B = 2,
    prior_A = list(diag = 0.4, off = 0.2),
    prior_B = list(diag = 0.4, off = 0.2),
    fixed_theta = TRUE, theta_A_init = NULL, theta_B_init = NULL,
    alpha_theta = 1, beta_theta = 1, warm_fixed_theta = TRUE,
    nu_A = 22, nu_B = 22,
    scale_init = 1, scale_df = c(22, 22),
    eta_shape = 1, eta_rate = 1,
    pip_cut = 0.5, pip_rule = ">=", rescale = TRUE,
    scale_update = "before_rescale", numerics = "density",
    coefficient_jitter = 0, whitening = "original",
    draw_initial_indicators = TRUE,
    seed = 201000L + rep_id, warm_seed = NULL, sampling_seed = NULL,
    keep_draws = TRUE, verbose_every = NULL,
    n_chains = 1L, workers = NULL, mcmc_clock = NULL
  )

  list(CGS = fit_cgs, SSVS = fit_ssvs)
})
```
