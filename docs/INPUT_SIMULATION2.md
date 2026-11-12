# Simulation 2

Install `MASS`, `MCMCpack`, `mcmcse`, and `Matrix` first. Place the core file at `R/mar_estimators.R`. Start a fresh R session for this experiment from the repository root, then run both R blocks in order. The first creates 100 training arrays; the second fits each array.

## Data generation

```r
G <- 15L
N <- 10L
L <- 3L
K_group_A <- 4L
K_group_B <- 3L
signal_scale_lag1 <- 1.20
signal_scale_lag3 <- 1.30
target_companion_rho <- 0.90
noise_scale <- 0.5

safe_chol <- function(S,
                      jitter_seq = c(
                        0, 1e-12, 1e-10, 1e-8,
                        1e-6, 1e-4, 1e-2
                      )) {
  S <- 0.5 * (S + t(S))
  p <- nrow(S)

  for (jj in jitter_seq) {
    R <- tryCatch(
      chol(S + jj * diag(p)),
      error = function(e) NULL
    )
    if (!is.null(R)) {
      return(R)
    }
  }

  stop("safe_chol() failed: matrix is not SPD after jitter.")
}

rescale_pair <- function(A, B) {
  norm_A <- sqrt(sum(A^2))
  norm_B <- sqrt(sum(B^2))

  if (!is.finite(norm_A) || !is.finite(norm_B) ||
    norm_A < .Machine$double.eps || norm_B < .Machine$double.eps) {
    return(list(A = A, B = B))
  }

  f_A <- sqrt(norm_B / norm_A)
  list(A = f_A * A, B = B / f_A)
}

make_ar1_cov <- function(p, rho, scale = NULL) {
  if (p <= 0L) {
    return(matrix(numeric(0), 0, 0))
  }
  if (is.null(scale)) scale <- rep(1, p)

  R <- toeplitz(rho^(0:(p - 1L)))
  D <- diag(scale, p)
  D %*% R %*% D
}

make_group_cov <- function(p, core_idx,
                           rho_core = 0.70,
                           rho_rest = 0.45,
                           cross_strength = 0.04,
                           jitter = 1e-10) {
  core_idx <- sort(unique(core_idx))
  rest_idx <- setdiff(seq_len(p), core_idx)

  stopifnot(length(core_idx) > 0L, all(core_idx %in% seq_len(p)))

  core_scale <- seq(1.15, 0.95, length.out = length(core_idx))
  rest_scale <- if (length(rest_idx) > 0L) {
    seq(0.90, 0.75, length.out = length(rest_idx))
  } else {
    numeric(0)
  }

  S_core <- make_ar1_cov(length(core_idx), rho_core, core_scale)
  S_rest <- make_ar1_cov(length(rest_idx), rho_rest, rest_scale)
  S_block <- as.matrix(Matrix::bdiag(S_core, S_rest))

  ordered_idx <- c(core_idx, rest_idx)
  if (!identical(ordered_idx, seq_len(p))) {
    P <- diag(p)[, ordered_idx, drop = FALSE]
    S_block <- P %*% S_block %*% t(P)
  }

  if (cross_strength > 0 && length(rest_idx) > 0L) {
    v <- numeric(p)
    v[core_idx] <- 1
    v[rest_idx] <- 0.35
    v <- v / sqrt(sum(v^2))
    S_block <- S_block +
      cross_strength * mean(diag(S_block)) * tcrossprod(v)
  }

  S_block <- 0.5 * (S_block + t(S_block)) + jitter * diag(p)
  min_eig <- min(eigen(S_block, symmetric = TRUE, only.values = TRUE)$values)

  if (!is.finite(min_eig) || min_eig <= 0) {
    S_block <- S_block + (abs(min_eig) + 1e-8) * diag(p)
  }

  S_block
}

make_true_matrix <- function(p, K,
                             diag_core,
                             off_core,
                             diag_rest,
                             core_effect,
                             lower_effect) {
  stopifnot(K >= 2L, K < p, length(core_effect) == K)

  core <- seq_len(K)
  rest <- (K + 1L):p
  M <- matrix(0, p, p)

  M[core, core] <- off_core
  diag(M)[core] <- diag_core

  for (ii in rest) {
    M[ii, core] <- core_effect * (0.90 + 0.20 * runif(K))
  }

  diag(M)[rest] <- diag_rest

  if (length(rest) >= 2L) {
    for (rr in 2:length(rest)) {
      row_id <- rest[rr]
      col_id <- rest[seq_len(rr - 1L)]
      M[row_id, col_id] <- lower_effect *
        (0.85 + 0.30 * runif(length(col_id)))
    }
  }

  M[core, rest] <- 0
  M
}

companion_matrix <- function(Phi_list) {
  L_here <- length(Phi_list)
  m <- nrow(Phi_list[[1]])

  stopifnot(
    L_here >= 1L,
    all(vapply(Phi_list, nrow, integer(1)) == m),
    all(vapply(Phi_list, ncol, integer(1)) == m)
  )

  C <- matrix(0, nrow = m * L_here, ncol = m * L_here)
  C[seq_len(m), ] <- do.call(cbind, Phi_list)

  if (L_here > 1L) {
    C[(m + 1L):(m * L_here), seq_len(m * (L_here - 1L))] <-
      diag(m * (L_here - 1L))
  }

  C
}

companion_rho <- function(Phi_list) {
  max(Mod(eigen(companion_matrix(Phi_list), only.values = TRUE)$values))
}

get_true_parameters <- function() {
  set.seed(1)

  A1 <- make_true_matrix(
    p = G,
    K = K_group_A,
    diag_core = 0.56,
    off_core = 0.22,
    diag_rest = 0.40,
    core_effect = c(0.40, 0.36, 0.32, 0.29),
    lower_effect = 0.150
  )

  B1 <- make_true_matrix(
    p = N,
    K = K_group_B,
    diag_core = 0.58,
    off_core = 0.21,
    diag_rest = 0.39,
    core_effect = c(0.38, 0.34, 0.30),
    lower_effect = 0.140
  )

  A2 <- matrix(0, G, G)
  B2 <- matrix(0, N, N)

  A3 <- make_true_matrix(
    p = G,
    K = K_group_A,
    diag_core = -0.50,
    off_core = -0.20,
    diag_rest = -0.37,
    core_effect = -c(0.38, 0.34, 0.31, 0.28),
    lower_effect = -0.140
  )

  B3 <- make_true_matrix(
    p = N,
    K = K_group_B,
    diag_core = 0.52,
    off_core = 0.19,
    diag_rest = 0.36,
    core_effect = c(0.36, 0.32, 0.29),
    lower_effect = 0.130
  )

  pair1 <- rescale_pair(A1, B1)
  pair3 <- rescale_pair(A3, B3)

  A1 <- sqrt(signal_scale_lag1) * pair1$A
  B1 <- sqrt(signal_scale_lag1) * pair1$B
  A3 <- sqrt(signal_scale_lag3) * pair3$A
  B3 <- sqrt(signal_scale_lag3) * pair3$B

  A_list <- list(A1 = A1, A2 = A2, A3 = A3)
  B_list <- list(B1 = B1, B2 = B2, B3 = B3)
  Phi_list <- Map(function(A, B) kronecker(B, A), A_list, B_list)

  rho_now <- companion_rho(Phi_list)
  shrink_step <- 0.95
  shrink_count <- 0L

  while (rho_now >= target_companion_rho && shrink_count < 200L) {
    A_list$A1 <- sqrt(shrink_step) * A_list$A1
    B_list$B1 <- sqrt(shrink_step) * B_list$B1
    A_list$A3 <- sqrt(shrink_step) * A_list$A3
    B_list$B3 <- sqrt(shrink_step) * B_list$B3
    shrink_count <- shrink_count + 1L
    Phi_list <- Map(function(A, B) kronecker(B, A), A_list, B_list)
    rho_now <- companion_rho(Phi_list)
  }

  if (!is.finite(rho_now) || rho_now >= 1) {
    stop("Unable to construct a stable MAR(3) DGP.")
  }

  Sigma_A_true <- make_group_cov(
    p = G,
    core_idx = seq_len(K_group_A),
    rho_core = 0.70,
    rho_rest = 0.45,
    cross_strength = 0.04
  )

  Sigma_B_true <- make_group_cov(
    p = N,
    core_idx = seq_len(K_group_B),
    rho_core = 0.65,
    rho_rest = 0.40,
    cross_strength = 0.04
  )

  covariance_factor_scale <- sqrt(noise_scale)
  Sigma_A_true <- covariance_factor_scale * Sigma_A_true
  Sigma_B_true <- covariance_factor_scale * Sigma_B_true

  list(
    A = A_list,
    B = B_list,
    Phi = Phi_list,
    Sigma_A_true = Sigma_A_true,
    Sigma_B_true = Sigma_B_true,
    Sigma_e = kronecker(Sigma_B_true, Sigma_A_true)
  )
}

simulate_mar3 <- function(seed, true_par, T_total, burn) {
  set.seed(seed)

  G_here <- nrow(true_par$A[[1]])
  N_here <- nrow(true_par$B[[1]])
  L_here <- length(true_par$Phi)
  m <- G_here * N_here
  total_needed <- burn + T_total

  stopifnot(total_needed > L_here)

  R_e <- safe_chol(true_par$Sigma_e)
  innovations <- matrix(
    rnorm(total_needed * m),
    nrow = total_needed,
    ncol = m
  ) %*% R_e

  Yvec <- matrix(0, nrow = total_needed, ncol = m)
  Yvec[seq_len(L_here), ] <- innovations[seq_len(L_here), ]

  for (tt in (L_here + 1L):total_needed) {
    conditional_mean <- numeric(m)

    for (ll in seq_len(L_here)) {
      conditional_mean <- conditional_mean +
        as.numeric(true_par$Phi[[ll]] %*% Yvec[tt - ll, ])
    }

    Yvec[tt, ] <- conditional_mean + innovations[tt, ]
  }

  Yvec <- Yvec[(burn + 1L):(burn + T_total), , drop = FALSE]
  array(t(Yvec), dim = c(G_here, N_here, T_total))
}

true_par <- get_true_parameters()
Y_train_list <- lapply(seq_len(100L), function(rep_id) {
  Y_arr <- simulate_mar3(
    seed = rep_id, true_par = true_par,
    T_total = 590L, burn = 400L
  )
  Y_arr[, , 1:450, drop = FALSE]
})
```

## Core-function inputs

```r
source("R/mar_estimators.R")
fits_by_rep <- lapply(seq_len(100L), function(rep_id) {
  Y_train <- Y_train_list[[rep_id]]

  fit_ggs <- fit_mar_pointmass(
    Y = Y_train, P = 3L,
    partition_A = c(
      list(core_A_1_to_4 = 1:4),
      setNames(as.list(5:15), paste0("singleton_A", 5:15))
    ),
    partition_B = c(
      list(core_B_1_to_3 = 1:3),
      setNames(as.list(4:10), paste0("singleton_B", 4:10))
    ),
    A_init = rep(list(0.30 * diag(15)), 3L),
    B_init = rep(list(0.30 * diag(10)), 3L),
    Sigma_A_init = diag(15), Sigma_B_init = diag(10),
    init_candidates = list(
      list(
        A = rep(list(0.30 * diag(15)), 3L),
        B = rep(list(0.30 * diag(10)), 3L)
      ),
      list(
        A = list(0.30 * diag(15), 0.30 * diag(15), -0.30 * diag(15)),
        B = rep(list(0.30 * diag(10)), 3L)
      )
    ),
    n_iter = 10000L, burnin = 5000L, warm_iter = 50L,
    tau_A = list(diag = 0.5, off = 0.5),
    tau_B = list(diag = 0.5, off = 0.5),
    prior_A = list(diag = 0.4, off = 0.2),
    prior_B = list(diag = 0.4, off = 0.2),
    nu_A = 17, nu_B = 12, scale_df = c(17, 12),
    eta_shape = 1, eta_rate = 1,
    scale_init = 12.5, pip_cut = 0.5, pip_rule = ">",
    rescale = TRUE, random_lag_scan = TRUE,
    numerics = "conditional", singleton_update = "scalar",
    whitening = "original", seed = NULL,
    warm_seed = 100017L + rep_id,
    sampling_seed = 101000L + rep_id,
    keep_draws = TRUE, verbose_every = NULL,
    n_chains = 1L, workers = NULL, mcmc_clock = NULL
  )

  fit_ssvs <- fit_mar_ssvs(
    Y = Y_train, P = 3L,
    A_init = 0.20 * diag(15), B_init = 0.20 * diag(10),
    Sigma_A_init = diag(15), Sigma_B_init = diag(10),
    n_iter = 10000L, burnin = 5000L, warm_iter = 50L,
    init_method = "gibbs", tau0_A = 0.1, tau1_A = 2,
    tau0_B = 0.1, tau1_B = 2,
    prior_A = list(diag = 0.4, off = 0.2),
    prior_B = list(diag = 0.4, off = 0.2),
    fixed_theta = TRUE, theta_A_init = NULL, theta_B_init = NULL,
    alpha_theta = 1, beta_theta = 1, warm_fixed_theta = TRUE,
    nu_A = 15, nu_B = 10,
    scale_init = 12.5, scale_df = c(17L, 12L),
    eta_shape = 1, eta_rate = 1,
    scale_update = "after_rescale", pip_cut = 0.5, pip_rule = ">",
    rescale = TRUE, numerics = "log_odds", coefficient_jitter = 1e-10,
    whitening = "original",
    draw_initial_indicators = FALSE, seed = 100000L + rep_id,
    warm_seed = NULL, sampling_seed = 101000L + rep_id,
    keep_draws = TRUE, verbose_every = NULL,
    n_chains = 1L, workers = NULL, mcmc_clock = NULL
  )

  list(GGS = fit_ggs, SSVS = fit_ssvs)
})
```
