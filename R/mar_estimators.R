# Required packages: MASS, MCMCpack, mcmcse, and Matrix. Install them before
# sourcing this file.

library(MASS)
library(MCMCpack)
library(mcmcse)
library(Matrix)

# Each estimator takes a numeric G x N x T_train array Y:
# Simulation 1: MAR(1), 20 x 20 x 90; see INPUT_SIMULATION1.md.
# Simulation 2: MAR(3), 15 x 10 x 450; see INPUT_SIMULATION2.md.
# Supplement S2: MAR(1), 24 x 28 x 160; see INPUT_SUPPLEMENT_S2.md.
# Each input file specifies the required training array and estimator settings.
# P: number of autoregressive lags.
# partition_A / partition_B: coefficient groups; NULL uses individual coefficients.
# A_init / B_init: starting coefficient matrices.
# Sigma_A_init / Sigma_B_init: starting innovation covariance matrices.
# n_iter: total MCMC sweeps.
# burnin: initial MCMC sweeps discarded from posterior summaries.
# warm_iter: initialization sweeps before the main MCMC run.
# tau_A / tau_B: prior scales for the A and B coefficients.
# prior_A / prior_B: prior inclusion probabilities for diagonal and off-diagonal entries.
# pip_cut / pip_rule: threshold and comparison used to select coefficients from PIPs.
# seed: seed for the estimator's main setup.
# warm_seed: seed for warm-up when supplied separately.
# sampling_seed: seed for posterior sampling when supplied separately.
# numerics: numerical form of the coefficient updates.
# whitening: choice of data whitening convention.
# singleton_update: update rule for one-coefficient groups.
# scale_update: timing of the scale update relative to coefficient rescaling.
# Exact values for each simulation are in the corresponding input Markdown file.

fit_mar_pointmass <- function(Y, P = 1L, partition_A = NULL, partition_B = NULL,
    A_init = NULL, B_init = NULL, Sigma_A_init = NULL, Sigma_B_init = NULL,
    n_iter = 10000L, burnin = 5000L, warm_iter = 20L,
    tau_A = 4, tau_B = 4, prior_A = list(diag = 0.4, off = 0.2),
    prior_B = list(diag = 0.4, off = 0.2), nu_A = NULL, nu_B = NULL,
    scale_init = NULL, scale_df = NULL, eta_shape = 1, eta_rate = 1,
    pip_cut = 0.5, pip_rule = c(">", ">="), rescale = TRUE,
    random_lag_scan = FALSE, init_candidates = NULL,
    numerics = c("precision", "conditional"),
    singleton_update = c("auto", "scalar", "group_1x1"),
    whitening = c("original", "paper"),
    seed = NULL, warm_seed = NULL, sampling_seed = NULL,
    keep_draws = TRUE, verbose_every = NULL, n_chains = 1L, workers = NULL,
    mcmc_clock = NULL) {
  # Input checks and matrix operations are private to this estimator.
  integer_arg <- function(x, name, lower = 0L) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x != floor(x) || x < lower || x > .Machine$integer.max)
      stop(name, " must be an integer >= ", lower, ".")
    as.integer(x)
  }
  if (!is.numeric(Y) || length(dim(Y)) != 3L || any(dim(Y) < 1L) || any(!is.finite(Y)))
    stop("Y must be a finite numeric G x N x T training array.")
  G <- dim(Y)[1L]; N <- dim(Y)[2L]; Tt <- dim(Y)[3L]
  P <- integer_arg(P, "P", 1L)
  n_iter <- integer_arg(n_iter, "n_iter", 1L)
  burnin <- integer_arg(burnin, "burnin")
  warm_iter <- integer_arg(warm_iter, "warm_iter")
  if (Tt <= P || n_iter <= burnin) stop("Require T > P and n_iter > burnin.")
  if (!is.numeric(pip_cut) || length(pip_cut) != 1L || !is.finite(pip_cut) ||
      pip_cut < 0 || pip_cut > 1) stop("pip_cut must be in [0, 1].")
  pip_rule <- match.arg(pip_rule)
  numerics <- match.arg(numerics)
  singleton_update <- match.arg(singleton_update)
  whitening <- match.arg(whitening)
  for (flag in c("rescale", "keep_draws")) {
    x <- get(flag)
    if (!is.logical(x) || length(x) != 1L || is.na(x)) stop(flag, " must be TRUE or FALSE.")
  }
  if (!is.null(verbose_every)) verbose_every <- integer_arg(verbose_every, "verbose_every", 1L)
  for (name in c("seed", "warm_seed", "sampling_seed"))
    if (!is.null(get(name))) assign(name, integer_arg(get(name), name))
  if (identical(n_chains, "auto")) {
    n_chains <- parallel::detectCores(logical = TRUE)
    if (is.na(n_chains)) stop("Cannot detect logical CPU cores.")
  }
  n_chains <- integer_arg(n_chains, "n_chains", 1L)
  if (n_chains > 1L) {
    if (is.null(workers)) workers <- parallel::detectCores(logical = TRUE)
    if (length(workers) == 1L && is.na(workers)) workers <- 1L
    workers <- integer_arg(workers, "workers", 1L)
    base_seed <- if (is.null(seed)) sample.int(.Machine$integer.max, 1L) else seed
    for (value in list(base_seed, warm_seed, sampling_seed)) {
      if (!is.null(value) && value > .Machine$integer.max - n_chains + 1L)
        stop("A seed plus n_chains - 1 exceeds the integer limit.")
    }
    params <- lapply(as.list(match.call(expand.dots = TRUE))[-1L],
                     eval, envir = parent.frame())
    params$n_chains <- 1L
    params$workers <- 1L
    jobs <- lapply(seq_len(n_chains), function(i) {
      args <- params
      offset <- i - 1L
      args$seed <- as.integer(as.double(base_seed) + offset)
      if (!is.null(warm_seed))
        args$warm_seed <- as.integer(as.double(warm_seed) + offset)
      if (!is.null(sampling_seed))
        args$sampling_seed <- as.integer(as.double(sampling_seed) + offset)
      args
    })
    self <- sys.function()
    cluster <- parallel::makePSOCKcluster(min(workers, n_chains))
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    chains <- parallel::parLapply(cluster, jobs, function(args, fit) do.call(fit, args),
                                  fit = self)
    names(chains) <- paste0("chain_", seq_len(n_chains))
    return(list(method = chains[[1L]]$method, P = P, n_chains = n_chains,
                workers = length(cluster), chain_seeds = vapply(jobs, function(x) x$seed, integer(1)),
                chains = chains))
  }
  positive <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0)
      stop(name, " must be a positive finite number.")
    x
  }
  spec <- function(x, name, probability = FALSE) {
    if (is.numeric(x) && length(x) == 1L) x <- c(diag = x, off = x)
    if (!(is.numeric(x) || is.list(x)) || length(x) != 2L ||
        !all(c("diag", "off") %in% names(x)))
      stop(name, " must be a scalar or list(diag=..., off=...).")
    out <- c(diag = positive(x[["diag"]], name), off = positive(x[["off"]], name))
    if (probability && any(out >= 1)) stop(name, " probabilities must be less than 1.")
    out
  }
  matrix_arg <- function(x, d, name, covariance = FALSE) {
    if (is.null(x)) x <- diag(d)
    if (!is.matrix(x) || !is.numeric(x) || !identical(dim(x), c(d, d)) || any(!is.finite(x)))
      stop(name, " must be a finite ", d, " x ", d, " matrix.")
    if (covariance) {
      if (!isSymmetric(x, tol = 1e-10) || inherits(try(chol(x), silent = TRUE), "try-error"))
        stop(name, " must be symmetric positive definite.")
    }
    x
  }
  lag_init <- function(x, d, name) {
    if (is.null(x) || is.matrix(x)) x <- rep(list(matrix_arg(x, d, name)), P)
    if (!is.list(x) || length(x) != P) stop(name, " must be one matrix or a list of P matrices.")
    lapply(x, matrix_arg, d = d, name = name)
  }
  initial_A <- lag_init(A_init, G, "A_init")
  initial_B <- lag_init(B_init, N, "B_init")
  initial_SA <- matrix_arg(Sigma_A_init, G, "Sigma_A_init", TRUE)
  initial_SB <- matrix_arg(Sigma_B_init, N, "Sigma_B_init", TRUE)
  if (is.null(nu_A)) nu_A <- G + 2
  if (is.null(nu_B)) nu_B <- N + 2
  if (positive(nu_A, "nu_A") <= G - 1 || positive(nu_B, "nu_B") <= N - 1)
    stop("Inverse-Wishart degrees of freedom must exceed dimension - 1.")
  eta_shape <- positive(eta_shape, "eta_shape"); eta_rate <- positive(eta_rate, "eta_rate")
  if (is.null(scale_df)) scale_df <- c(nu_A, nu_B)
  if (!is.numeric(scale_df) || length(scale_df) != 2L || any(!is.finite(scale_df)) || any(scale_df <= 0))
    stop("scale_df must contain two positive degrees of freedom.")
  if (is.null(scale_init)) scale_init <- (G + N) / 2
  scale_init <- positive(scale_init, "scale_init")
  priorA <- spec(prior_A, "prior_A", TRUE); priorB <- spec(prior_B, "prior_B", TRUE)
  Q <- Tt - P
  target <- Y[, (seq_len(N)), (P + 1L):Tt, drop = FALSE]
  lags <- lapply(seq_len(P), function(lag) Y[, , (P + 1L - lag):(Tt - lag), drop = FALSE])
  slice <- function(x, q) matrix(x[, , q], nrow = dim(x)[1L], ncol = dim(x)[2L])
  sym <- function(x) 0.5 * (x + t(x))
  chol_safe <- function(x, jitter = chol_jitter, scaled = chol_scaled) {
    x <- sym(x)
    factor <- if (scaled) max(abs(diag(x)), 1) else 1
    if (any(!is.finite(x))) stop("Non-finite covariance or precision matrix.")
    for (eps in jitter) {
      out <- tryCatch(chol(x / factor + eps * diag(nrow(x))), error = function(e) NULL)
      if (!is.null(out)) return(if (scaled) sqrt(factor) * out else out)
    }
    stop("Cholesky decomposition failed after jitter.")
  }
  inv <- function(x) chol2inv(chol_safe(x))
  balance_pair <- function(a, b, symmetric = FALSE) {
    na <- sqrt(sum(a^2)); nb <- sqrt(sum(b^2))
    if (!is.finite(na) || !is.finite(nb) || na < .Machine$double.eps || nb < .Machine$double.eps)
      return(list(a, b))
    factor <- sqrt(nb / na)
    a <- factor * a; b <- b / factor
    if (symmetric) { a <- sym(a); b <- sym(b) }
    list(a, b)
  }
  balance <- function(state) {
    if (!rescale) return(state)
    for (lag in seq_len(P)) {
      pair <- balance_pair(state$A[[lag]], state$B[[lag]])
      state$A[[lag]] <- pair[[1L]]; state$B[[lag]] <- pair[[2L]]
    }
    pair <- balance_pair(state$SA, state$SB, sym_cov)
    state$SA <- pair[[1L]]; state$SB <- pair[[2L]]
    state
  }
  update_scale <- function(state, map = FALSE) {
    shape <- 0.5 * (G * scale_df[1L] + N * scale_df[2L]) + eta_shape
    if (scale_matrix_trace) {
      rate <- 0.5 * (sum(diag((diag(G) / G) %*% inv(state$SA))) +
                      sum(diag((diag(N) / N) %*% inv(state$SB)))) + eta_rate
    } else {
      rate <- 0.5 * sum(diag(inv(state$SA))) / G +
        0.5 * sum(diag(inv(state$SB))) / N + eta_rate
    }
    if (map) return(if (shape > 1) (shape - 1) / rate else shape / rate)
    rgamma(1L, shape = shape, rate = rate)
  }
  selected <- function(pip) if (pip_rule == ">") pip > pip_cut else pip >= pip_cut
  # A single, stable output shape is used even when P=1.
  summarize <- function(chain, indicators, side) {
    estimates <- pips <- active_means <- all_means <- vector("list", P)
    d <- dim(chain)[3L]; n <- dim(chain)[1L]
    for (lag in seq_len(P)) {
      draws <- array(chain[, lag, , , drop = FALSE], c(n, d, d))
      on <- array(indicators[, lag, , , drop = FALSE], c(n, d, d))
      pip <- apply(on, c(2L, 3L), mean)
      active <- matrix(0, d, d)
      for (i in seq_len(d)) for (j in seq_len(d))
        if (any(on[, i, j])) active[i, j] <- mean(draws[on[, i, j], i, j])
      estimates[[lag]] <- active * selected(pip)
      pips[[lag]] <- pip; active_means[[lag]] <- active
      all_means[[lag]] <- apply(draws, c(2L, 3L), mean)
    }
    labels <- paste0(side, seq_len(P))
    names(estimates) <- names(pips) <- names(active_means) <- names(all_means) <- labels
    list(hat = estimates, PIP = pips, active = active_means, mean = all_means)
  }

  tauA <- spec(tau_A, "tau_A"); tauB <- spec(tau_B, "tau_B")
  partition <- function(x, d, name) {
    if (is.null(x)) return(lapply(seq_len(d), function(i) i))
    if (!is.list(x) || !length(x) || any(lengths(x) == 0L) ||
        !all(vapply(x, is.numeric, logical(1)))) stop(name, " must be a nonempty list of index groups.")
    idx <- unlist(x, use.names = FALSE)
    if (any(!is.finite(idx)) || any(idx != floor(idx)) || anyDuplicated(idx) ||
        !setequal(idx, seq_len(d))) stop(name, " must contain each index 1:", d, " exactly once.")
    lapply(x, function(g) sort(as.integer(g)))
  }
  partA <- partition(partition_A, G, "partition_A")
  partB <- partition(partition_B, N, "partition_B")
  modeA <- if (singleton_update == "auto")
    if (is.null(partition_A)) "scalar" else "group_1x1" else singleton_update
  modeB <- if (singleton_update == "auto")
    if (is.null(partition_B)) "scalar" else "group_1x1" else singleton_update
  if (!is.logical(random_lag_scan) || length(random_lag_scan) != 1L || is.na(random_lag_scan))
    stop("random_lag_scan must be TRUE or FALSE.")
  # Numerical policies preserve existing implementations without dispatching on P.
  precision <- numerics == "precision"
  chol_scaled <- precision; sym_cov <- precision; scale_matrix_trace <- FALSE
  chol_jitter <- c(0, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2)
  logdet <- function(x) {
    r <- if (precision) chol_safe(x) else chol_safe(x, jitter = c(0, 1e-12, 1e-10, 1e-8, 1e-6))
    2 * sum(log(diag(r)))
  }
  quadratic <- function(x, v) {
    if (precision) {
      z <- forwardsolve(t(chol_safe(x)), as.numeric(v))
      return(as.numeric(crossprod(z)))
    }
    r <- chol_safe(x, jitter = c(0, 1e-12, 1e-10, 1e-8, 1e-6))
    z <- backsolve(r, forwardsolve(t(r), v))
    as.numeric(crossprod(v, z))
  }
  logit <- function(p) {
    if (!precision) p <- pmin(pmax(p, 1e-10), 1 - 1e-10)
    log(p / (1 - p))
  }
  probability <- function(x) plogis(if (precision) pmax(pmin(x, 700), -700) else pmin(x, 700))
  update_mode <- function(response, design, Sigma, Phi, scale, groups, tau, prior, nu,
                          singleton_mode) {
    d <- nrow(response); nobs <- ncol(response)
    residual_all <- response - Phi %*% design
    all_singletons <- all(lengths(groups) == 1L)
    if (precision) {
      S <- sym(Sigma) + 1e-10 * diag(d)
      omega <- inv(S)
      gram <- crossprod(t(design))
    }
    for (i in seq_len(d)) {
      other <- setdiff(seq_len(d), i)
      if (precision) {
        w <- if (length(other)) -as.numeric(omega[i, other]) / omega[i, i] else numeric(0)
        variance <- if (length(other)) max(1 / omega[i, i], 1e-12) else max(S[i, i], 1e-12)
      } else if (length(other)) {
        Soo <- sym(Sigma[other, other, drop = FALSE]) + 1e-10 * diag(length(other))
        soi <- Sigma[other, i, drop = FALSE]
        w <- as.numeric(solve(Soo, soi))
        variance <- as.numeric(Sigma[i, i] - crossprod(soi, solve(Soo, soi)))
        if (!is.finite(variance) || variance <= 1e-12) variance <- 1e-12
      } else { w <- numeric(0); variance <- max(as.numeric(Sigma[i, i]), 1e-12) }
      adjustment <- if (length(other)) as.numeric(t(w) %*% residual_all[other, , drop = FALSE]) else rep(0, nobs)
      y <- as.numeric(response[i, ] - adjustment)
      row_residual <- y - as.numeric(Phi[i, ] %*% design)
      # MAR(1) visits singleton coefficients in row order. Original MAR(3)
      # visits the diagonal first, then each off-diagonal group.
      updates <- if (singleton_mode == "scalar" && all_singletons)
        as.list(seq_len(d)) else
        c(list(i), lapply(groups, function(block) setdiff(block, i)))
      for (idx in updates) {
        if (!length(idx)) next
        if (length(idx) == 1L &&
            (idx == i || singleton_mode == "scalar")) {
          j <- idx[[1L]]
          x <- as.numeric(design[j, ])
          partial <- row_residual + Phi[i, j] * x
          slab <- if (j == i) tau[["diag"]] else tau[["off"]]
          inclusion <- if (j == i) prior[["diag"]] else prior[["off"]]
          xx <- if (precision) gram[j, j] else sum(x^2)
          denom <- variance + slab^2 * xx
          v <- max(slab^2 * variance / denom, 1e-12)
          m <- slab^2 * sum(x * partial) / denom
          odds <- logit(inclusion) + 0.5 * (log(v) - log(slab^2)) + m^2 / (2 * v)
          Phi[i, j] <- if (rbinom(1L, 1L, probability(odds))) rnorm(1L, m, sqrt(v)) else 0
          row_residual <- partial - Phi[i, j] * x
        } else {
          Xg <- design[idx, , drop = FALSE]
          partial <- row_residual + as.numeric(Phi[i, idx, drop = FALSE] %*% Xg)
          post_precision <- Xg %*% t(Xg) / variance + diag(1 / tau[["off"]]^2, length(idx))
          V <- if (precision) chol2inv(chol_safe(post_precision)) else
            tryCatch(solve(post_precision), error = function(e) solve(post_precision + 1e-10 * diag(length(idx))))
          V <- sym(V)
          mu <- V %*% (Xg %*% matrix(partial, ncol = 1L) / variance)
          prior_logdet <- if (precision) length(idx) * log(tau[["off"]]^2) else
            logdet(diag(tau[["off"]]^2, length(idx)))
          odds <- logit(prior[["off"]]) + 0.5 * (logdet(V) - prior_logdet) + 0.5 * quadratic(V, mu)
          on <- rbinom(1L, 1L, probability(odds))
          if (!on) Phi[i, idx] <- 0 else if (length(idx) == 1L)
            Phi[i, idx] <- rnorm(1L, as.numeric(mu), sqrt(as.numeric(V))) else
            Phi[i, idx] <- as.numeric(MASS::mvrnorm(1L, as.numeric(mu), V + 1e-10 * diag(length(idx))))
          row_residual <- partial - as.numeric(Phi[i, idx, drop = FALSE] %*% Xg)
        }
      }
      residual_all[i, ] <- as.numeric(response[i, ] - Phi[i, ] %*% design)
    }
    residual <- response - Phi %*% design
    list(Phi = Phi, Sigma = MCMCpack::riwish(nu + nobs, scale * (diag(d) / d) + residual %*% t(residual)))
  }
  contribution <- function(a, b, x) {
    out <- array(0, dim(x))
    for (q in seq_len(Q)) out[, , q] <- a %*% slice(x, q) %*% t(b)
    out
  }
  residual <- function(state, exclude = integer(0)) {
    out <- target
    for (lag in setdiff(seq_len(P), exclude))
      out <- out - contribution(state$A[[lag]], state$B[[lag]], lags[[lag]])
    out
  }
  work <- function(state, lag, side) {
    response <- residual(state, lag)
    is_A <- side == "A"
    d <- if (is_A) G else N; other <- if (is_A) N else G
    chol_other <- chol_safe(if (is_A) state$SB else state$SA)
    whiten <- if (whitening == "original") t(solve(chol_other)) else solve(chol_other)
    coef <- if (is_A) state$B[[lag]] else state$A[[lag]]
    yw <- xw <- array(0, c(d, other, Q))
    for (q in seq_len(Q)) {
      y <- slice(response, q); x <- slice(lags[[lag]], q)
      if (!is_A) { y <- t(y); x <- t(x) }
      yw[, , q] <- y %*% whiten
      xw[, , q] <- x %*% t(coef) %*% whiten
    }
    list(Y = matrix(yw, nrow = d), X = matrix(xw, nrow = d))
  }
  sweep <- function(state) {
    order_A <- if (random_lag_scan) sample(seq_len(P)) else seq_len(P)
    for (lag in order_A) {
      w <- work(state, lag, "A")
      up <- update_mode(w$Y, w$X, state$SA, state$A[[lag]], state$scale,
                        partA, tauA, priorA, nu_A, modeA)
      state$A[[lag]] <- up$Phi; state$SA <- up$Sigma
    }
    order_B <- if (random_lag_scan) sample(seq_len(P)) else seq_len(P)
    for (lag in order_B) {
      w <- work(state, lag, "B")
      up <- update_mode(w$Y, w$X, state$SB, state$B[[lag]], state$scale,
                        partB, tauB, priorB, nu_B, modeB)
      state$B[[lag]] <- up$Phi; state$SB <- up$Sigma
    }
    state$scale <- update_scale(state)
    balance(state)
  }
  base <- list(A = initial_A, B = initial_B, SA = initial_SA, SB = initial_SB, scale = scale_init)
  candidates <- list(base)
  if (!is.null(init_candidates)) {
    if (!is.list(init_candidates) || !length(init_candidates)) stop("init_candidates must be a nonempty list.")
    candidates <- lapply(init_candidates, function(x) {
      if (!is.list(x) || !all(c("A", "B") %in% names(x))) stop("Each candidate needs A and B initial values.")
      state <- base
      state$A <- lag_init(x$A, G, "candidate A"); state$B <- lag_init(x$B, N, "candidate B")
      state
    })
  }
  if (!is.null(seed)) set.seed(seed)
  warmed <- lapply(candidates, function(state) {
    if (!is.null(warm_seed)) set.seed(warm_seed)
    for (i in seq_len(warm_iter)) state <- sweep(state)
    state
  })
  scores <- vapply(warmed, function(state) sum(residual(state)^2), numeric(1))
  chosen <- which.min(scores); state <- warmed[[chosen]]
  if (!is.null(sampling_seed)) set.seed(sampling_seed)
  n_post <- n_iter - burnin
  A_chain <- array(NA_real_, c(n_post, P, G, G)); B_chain <- array(NA_real_, c(n_post, P, N, N))
  scale_chain <- numeric(n_post)
  if (!is.null(mcmc_clock)) {
    mcmc_clock$start()
    on.exit(mcmc_clock$cancel(), add = TRUE)
  }
  started <- Sys.time()
  for (iter in seq_len(n_iter)) {
    state <- sweep(state)
    if (iter > burnin) {
      s <- iter - burnin
      for (lag in seq_len(P)) { A_chain[s, lag, , ] <- state$A[[lag]]; B_chain[s, lag, , ] <- state$B[[lag]] }
      scale_chain[s] <- state$scale
    }
    if (!is.null(verbose_every) && iter %% verbose_every == 0L) cat("Point-mass MAR(", P, ") ", iter, "/", n_iter, "\n", sep = "")
  }
  runtime <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  if (!is.null(mcmc_clock)) mcmc_clock$stop()
  a <- summarize(A_chain, A_chain != 0, "A"); b <- summarize(B_chain, B_chain != 0, "B")
  method <- if (is.null(partition_A) && is.null(partition_B)) "CGS" else "GGS"
  out <- list(method = method, P = P, A_hat = a$hat, B_hat = b$hat, PIP_A = a$PIP, PIP_B = b$PIP,
    Sigma_A_last = state$SA, Sigma_B_last = state$SB, scale_last = state$scale,
    warm_start_index = chosen, warm_scores = scores, runtime_sec = runtime,
    settings = list(P = P, n_iter = n_iter, burnin = burnin, warm_iter = warm_iter,
      partition_A = partA, partition_B = partB, tau_A = tauA, tau_B = tauB,
      prior_A = priorA, prior_B = priorB, pip_cut = pip_cut, pip_rule = pip_rule,
      random_lag_scan = random_lag_scan, numerics = numerics,
      singleton_update = c(A = modeA, B = modeB),
      singleton_update_policy = singleton_update, whitening = whitening,
      rescale = rescale))
  if (keep_draws) { out$A_chain <- A_chain; out$B_chain <- B_chain; out$scale_chain <- scale_chain }
  out
}

fit_mar_ssvs <- function(Y, P = 1L,
    A_init = NULL, B_init = NULL, Sigma_A_init = NULL, Sigma_B_init = NULL,
    n_iter = 2250L, burnin = 750L, warm_iter = 10L,
    init_method = c("map", "gibbs", "none"),
    tau0_A = 0.05, tau1_A = 2, tau0_B = tau0_A, tau1_B = tau1_A,
    prior_A = list(diag = 0.4, off = 0.2), prior_B = list(diag = 0.4, off = 0.2),
    fixed_theta = TRUE, theta_A_init = NULL, theta_B_init = NULL,
    alpha_theta = 1, beta_theta = 1, warm_fixed_theta = TRUE,
    nu_A = NULL, nu_B = NULL, scale_init = 1, scale_df = NULL, eta_shape = 1, eta_rate = 1,
    pip_cut = 0.5, pip_rule = c(">", ">="), rescale = TRUE,
    scale_update = c("before_rescale", "after_rescale"),
    numerics = c("density", "log_odds"), coefficient_jitter = 0,
    whitening = c("original", "paper"),
    draw_initial_indicators = TRUE, seed = NULL, warm_seed = NULL, sampling_seed = NULL,
    keep_draws = TRUE, verbose_every = NULL, n_chains = 1L, workers = NULL,
    mcmc_clock = NULL) {
  # Input checks and matrix operations are private to this estimator.
  integer_arg <- function(x, name, lower = 0L) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x != floor(x) || x < lower || x > .Machine$integer.max)
      stop(name, " must be an integer >= ", lower, ".")
    as.integer(x)
  }
  if (!is.numeric(Y) || length(dim(Y)) != 3L || any(dim(Y) < 1L) || any(!is.finite(Y)))
    stop("Y must be a finite numeric G x N x T training array.")
  G <- dim(Y)[1L]; N <- dim(Y)[2L]; Tt <- dim(Y)[3L]
  P <- integer_arg(P, "P", 1L)
  n_iter <- integer_arg(n_iter, "n_iter", 1L)
  burnin <- integer_arg(burnin, "burnin")
  warm_iter <- integer_arg(warm_iter, "warm_iter")
  if (Tt <= P || n_iter <= burnin) stop("Require T > P and n_iter > burnin.")
  if (!is.numeric(pip_cut) || length(pip_cut) != 1L || !is.finite(pip_cut) ||
      pip_cut < 0 || pip_cut > 1) stop("pip_cut must be in [0, 1].")
  pip_rule <- match.arg(pip_rule)
  numerics <- match.arg(numerics)
  whitening <- match.arg(whitening)
  for (flag in c("rescale", "keep_draws")) {
    x <- get(flag)
    if (!is.logical(x) || length(x) != 1L || is.na(x)) stop(flag, " must be TRUE or FALSE.")
  }
  if (!is.null(verbose_every)) verbose_every <- integer_arg(verbose_every, "verbose_every", 1L)
  for (name in c("seed", "warm_seed", "sampling_seed"))
    if (!is.null(get(name))) assign(name, integer_arg(get(name), name))
  if (identical(n_chains, "auto")) {
    n_chains <- parallel::detectCores(logical = TRUE)
    if (is.na(n_chains)) stop("Cannot detect logical CPU cores.")
  }
  n_chains <- integer_arg(n_chains, "n_chains", 1L)
  if (n_chains > 1L) {
    if (is.null(workers)) workers <- parallel::detectCores(logical = TRUE)
    if (length(workers) == 1L && is.na(workers)) workers <- 1L
    workers <- integer_arg(workers, "workers", 1L)
    base_seed <- if (is.null(seed)) sample.int(.Machine$integer.max, 1L) else seed
    for (value in list(base_seed, warm_seed, sampling_seed)) {
      if (!is.null(value) && value > .Machine$integer.max - n_chains + 1L)
        stop("A seed plus n_chains - 1 exceeds the integer limit.")
    }
    params <- lapply(as.list(match.call(expand.dots = TRUE))[-1L],
                     eval, envir = parent.frame())
    params$n_chains <- 1L
    params$workers <- 1L
    jobs <- lapply(seq_len(n_chains), function(i) {
      args <- params
      offset <- i - 1L
      args$seed <- as.integer(as.double(base_seed) + offset)
      if (!is.null(warm_seed))
        args$warm_seed <- as.integer(as.double(warm_seed) + offset)
      if (!is.null(sampling_seed))
        args$sampling_seed <- as.integer(as.double(sampling_seed) + offset)
      args
    })
    self <- sys.function()
    cluster <- parallel::makePSOCKcluster(min(workers, n_chains))
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    chains <- parallel::parLapply(cluster, jobs, function(args, fit) do.call(fit, args),
                                  fit = self)
    names(chains) <- paste0("chain_", seq_len(n_chains))
    return(list(method = chains[[1L]]$method, P = P, n_chains = n_chains,
                workers = length(cluster), chain_seeds = vapply(jobs, function(x) x$seed, integer(1)),
                chains = chains))
  }
  positive <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0)
      stop(name, " must be a positive finite number.")
    x
  }
  spec <- function(x, name, probability = FALSE) {
    if (is.numeric(x) && length(x) == 1L) x <- c(diag = x, off = x)
    if (!(is.numeric(x) || is.list(x)) || length(x) != 2L ||
        !all(c("diag", "off") %in% names(x)))
      stop(name, " must be a scalar or list(diag=..., off=...).")
    out <- c(diag = positive(x[["diag"]], name), off = positive(x[["off"]], name))
    if (probability && any(out >= 1)) stop(name, " probabilities must be less than 1.")
    out
  }
  matrix_arg <- function(x, d, name, covariance = FALSE) {
    if (is.null(x)) x <- diag(d)
    if (!is.matrix(x) || !is.numeric(x) || !identical(dim(x), c(d, d)) || any(!is.finite(x)))
      stop(name, " must be a finite ", d, " x ", d, " matrix.")
    if (covariance) {
      if (!isSymmetric(x, tol = 1e-10) || inherits(try(chol(x), silent = TRUE), "try-error"))
        stop(name, " must be symmetric positive definite.")
    }
    x
  }
  lag_init <- function(x, d, name) {
    if (is.null(x) || is.matrix(x)) x <- rep(list(matrix_arg(x, d, name)), P)
    if (!is.list(x) || length(x) != P) stop(name, " must be one matrix or a list of P matrices.")
    lapply(x, matrix_arg, d = d, name = name)
  }
  initial_A <- lag_init(A_init, G, "A_init")
  initial_B <- lag_init(B_init, N, "B_init")
  initial_SA <- matrix_arg(Sigma_A_init, G, "Sigma_A_init", TRUE)
  initial_SB <- matrix_arg(Sigma_B_init, N, "Sigma_B_init", TRUE)
  if (is.null(nu_A)) nu_A <- G + 2
  if (is.null(nu_B)) nu_B <- N + 2
  if (positive(nu_A, "nu_A") <= G - 1 || positive(nu_B, "nu_B") <= N - 1)
    stop("Inverse-Wishart degrees of freedom must exceed dimension - 1.")
  eta_shape <- positive(eta_shape, "eta_shape"); eta_rate <- positive(eta_rate, "eta_rate")
  if (is.null(scale_df)) scale_df <- c(nu_A, nu_B)
  if (!is.numeric(scale_df) || length(scale_df) != 2L || any(!is.finite(scale_df)) || any(scale_df <= 0))
    stop("scale_df must contain two positive degrees of freedom.")
  if (is.null(scale_init)) scale_init <- (G + N) / 2
  scale_init <- positive(scale_init, "scale_init")
  priorA <- spec(prior_A, "prior_A", TRUE); priorB <- spec(prior_B, "prior_B", TRUE)
  Q <- Tt - P
  target <- Y[, (seq_len(N)), (P + 1L):Tt, drop = FALSE]
  lags <- lapply(seq_len(P), function(lag) Y[, , (P + 1L - lag):(Tt - lag), drop = FALSE])
  slice <- function(x, q) matrix(x[, , q], nrow = dim(x)[1L], ncol = dim(x)[2L])
  sym <- function(x) 0.5 * (x + t(x))
  chol_safe <- function(x, jitter = chol_jitter, scaled = chol_scaled) {
    x <- sym(x)
    factor <- if (scaled) max(abs(diag(x)), 1) else 1
    if (any(!is.finite(x))) stop("Non-finite covariance or precision matrix.")
    for (eps in jitter) {
      out <- tryCatch(chol(x / factor + eps * diag(nrow(x))), error = function(e) NULL)
      if (!is.null(out)) return(if (scaled) sqrt(factor) * out else out)
    }
    stop("Cholesky decomposition failed after jitter.")
  }
  inv <- function(x) chol2inv(chol_safe(x))
  balance_pair <- function(a, b, symmetric = FALSE) {
    na <- sqrt(sum(a^2)); nb <- sqrt(sum(b^2))
    if (!is.finite(na) || !is.finite(nb) || na < .Machine$double.eps || nb < .Machine$double.eps)
      return(list(a, b))
    factor <- sqrt(nb / na)
    a <- factor * a; b <- b / factor
    if (symmetric) { a <- sym(a); b <- sym(b) }
    list(a, b)
  }
  balance <- function(state) {
    if (!rescale) return(state)
    for (lag in seq_len(P)) {
      pair <- balance_pair(state$A[[lag]], state$B[[lag]])
      state$A[[lag]] <- pair[[1L]]; state$B[[lag]] <- pair[[2L]]
    }
    pair <- balance_pair(state$SA, state$SB, sym_cov)
    state$SA <- pair[[1L]]; state$SB <- pair[[2L]]
    state
  }
  update_scale <- function(state, map = FALSE) {
    shape <- 0.5 * (G * scale_df[1L] + N * scale_df[2L]) + eta_shape
    if (scale_matrix_trace) {
      rate <- 0.5 * (sum(diag((diag(G) / G) %*% inv(state$SA))) +
                      sum(diag((diag(N) / N) %*% inv(state$SB)))) + eta_rate
    } else {
      rate <- 0.5 * sum(diag(inv(state$SA))) / G +
        0.5 * sum(diag(inv(state$SB))) / N + eta_rate
    }
    if (map) return(if (shape > 1) (shape - 1) / rate else shape / rate)
    rgamma(1L, shape = shape, rate = rate)
  }
  selected <- function(pip) if (pip_rule == ">") pip > pip_cut else pip >= pip_cut
  # A single, stable output shape is used even when P=1.
  summarize <- function(chain, indicators, side) {
    estimates <- pips <- active_means <- all_means <- vector("list", P)
    d <- dim(chain)[3L]; n <- dim(chain)[1L]
    for (lag in seq_len(P)) {
      draws <- array(chain[, lag, , , drop = FALSE], c(n, d, d))
      on <- array(indicators[, lag, , , drop = FALSE], c(n, d, d))
      pip <- apply(on, c(2L, 3L), mean)
      active <- matrix(0, d, d)
      for (i in seq_len(d)) for (j in seq_len(d))
        if (any(on[, i, j])) active[i, j] <- mean(draws[on[, i, j], i, j])
      estimates[[lag]] <- active * selected(pip)
      pips[[lag]] <- pip; active_means[[lag]] <- active
      all_means[[lag]] <- apply(draws, c(2L, 3L), mean)
    }
    labels <- paste0(side, seq_len(P))
    names(estimates) <- names(pips) <- names(active_means) <- names(all_means) <- labels
    list(hat = estimates, PIP = pips, active = active_means, mean = all_means)
  }

  init_method <- match.arg(init_method)
  scale_update <- match.arg(scale_update)
  for (flag in c("fixed_theta", "draw_initial_indicators")) {
    x <- get(flag)
    if (!is.logical(x) || length(x) != 1L || is.na(x)) stop(flag, " must be TRUE or FALSE.")
  }
  tau0_A <- positive(tau0_A, "tau0_A"); tau1_A <- positive(tau1_A, "tau1_A")
  tau0_B <- positive(tau0_B, "tau0_B"); tau1_B <- positive(tau1_B, "tau1_B")
  if (tau1_A <= tau0_A || tau1_B <= tau0_B) stop("Each slab scale must exceed its spike scale.")
  alpha_theta <- positive(alpha_theta, "alpha_theta"); beta_theta <- positive(beta_theta, "beta_theta")
  if (!is.numeric(coefficient_jitter) || length(coefficient_jitter) != 1L ||
      !is.finite(coefficient_jitter) || coefficient_jitter < 0)
    stop("coefficient_jitter must be a finite nonnegative number.")
  density <- numerics == "density"
  chol_scaled <- FALSE; sym_cov <- density; scale_matrix_trace <- density
  chol_jitter <- if (density) c(0, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2) else
    c(0, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2)
  prior_matrix <- function(d, prior) {
    one <- matrix(prior[["off"]], d, d); diag(one) <- prior[["diag"]]
    do.call(cbind, rep(list(one), P))
  }
  fixed_A <- prior_matrix(G, priorA); fixed_B <- prior_matrix(N, priorB)
  # Optional explicit initial probabilities also support estimating theta.
  theta_initial <- function(x, d, fixed) {
    if (is.null(x)) return(if (fixed_theta) fixed else alpha_theta / (alpha_theta + beta_theta))
    if (!is.numeric(x) || !(length(x) %in% c(1L, d * d, d * d * P)) ||
        any(!is.finite(x)) || any(x <= 0 | x >= 1)) stop("theta initial values must be probabilities of compatible size.")
    if (length(x) == 1L) return(x)
    matrix(rep_len(as.numeric(x), d * d * P), d, d * P)
  }
  theta_A <- theta_initial(theta_A_init, G, fixed_A)
  theta_B <- theta_initial(theta_B_init, N, fixed_B)
  if (fixed_theta) { fixed_A <- theta_A; fixed_B <- theta_B }
  inclusion <- function(phi, theta, tau0, tau1) {
    theta <- rep_len(as.numeric(theta), length(phi))
    if (density) {
      theta <- pmin(pmax(theta, 1e-8), 1 - 1e-8)
      log1 <- log(theta) + dnorm(phi, mean = 0, sd = tau1, log = TRUE)
      log0 <- log1p(-theta) + dnorm(phi, mean = 0, sd = tau0, log = TRUE)
      return(plogis(log1 - log0))
    }
    theta <- pmin(pmax(theta, 1e-10), 1 - 1e-10)
    odds <- log(theta / (1 - theta)) - 0.5 * (log(tau1^2) - log(tau0^2)) -
      0.5 * phi^2 * (1 / tau1^2 - 1 / tau0^2)
    plogis(pmax(pmin(odds, 700), -700))
  }
  split_lags <- function(x, d) lapply(seq_len(P), function(lag) x[, (lag - 1L) * d + seq_len(d), drop = FALSE])
  # The block diagonal design is constructed from P; no order-specific code paths.
  X_block <- array(0, c(P * G, P * N, Q))
  for (q in seq_len(Q)) for (lag in seq_len(P))
    X_block[(lag - 1L) * G + seq_len(G), (lag - 1L) * N + seq_len(N), q] <- slice(lags[[lag]], q)
  work <- function(state, side) {
    is_A <- side == "A"
    d <- if (is_A) G else N; other <- if (is_A) N else G
    coef <- do.call(cbind, if (is_A) state$B else state$A)
    chol_other <- chol_safe(if (is_A) state$SB else state$SA)
    whiten <- if (whitening == "original") t(solve(chol_other)) else solve(chol_other)
    yw <- array(0, c(d, other, Q)); xw <- array(0, c(P * d, other, Q))
    for (q in seq_len(Q)) {
      y <- slice(target, q); x <- slice(X_block, q)
      if (!is_A) { y <- t(y); x <- t(x) }
      yw[, , q] <- y %*% whiten
      xw[, , q] <- x %*% t(coef) %*% whiten
    }
    list(Y = matrix(yw, nrow = d), X = matrix(xw, nrow = P * d))
  }
  update_mode <- function(w, Sigma, Phi, theta, scale, tau0, tau1, nu, fixed, map = FALSE) {
    d <- nrow(Phi); k <- ncol(Phi); nobs <- ncol(w$Y)
    prob <- inclusion(as.vector(Phi), theta, tau0, tau1)
    gamma <- if (map) prob else rbinom(length(prob), size = 1L, prob = prob)
    prior_precision <- if (map) (1 - gamma) / (tau0^2) + gamma / (tau1^2) else
      1 / ifelse(gamma == 1L, tau1^2, tau0^2)
    gram <- tcrossprod(w$X)
    like <- kronecker(gram, inv(Sigma))
    if (!density) like <- sym(like)
    prec <- like + diag(prior_precision, d * k)
    if (density) prec <- sym(prec)
    covariance <- inv(prec)
    if (!density) covariance <- sym(covariance)
    ols <- w$Y %*% t(w$X) %*% solve(gram + 1e-8 * diag(k))
    mu <- as.numeric(covariance %*% (like %*% as.vector(ols)))
    if (map) values <- mu else {
      covariance <- sym(covariance)
      # Avoid adding even a zero matrix when retaining the original density path.
      if (coefficient_jitter > 0) covariance <- covariance + coefficient_jitter * diag(d * k)
      values <- as.numeric(MASS::mvrnorm(1L, mu = mu, Sigma = covariance))
    }
    coef <- matrix(values, d, k)
    residual <- w$Y - coef %*% w$X
    posterior_scale <- residual %*% t(residual) + scale * (diag(d) / d)
    sigma <- if (map) sym(posterior_scale / (nu + nobs + d + 1)) else
      MCMCpack::riwish(nu + nobs, if (density) sym(posterior_scale) else posterior_scale)
    next_theta <- if (!is.null(fixed)) fixed else if (map)
      (sum(gamma) + alpha_theta) / (length(gamma) + alpha_theta + beta_theta) else
      rbeta(1L, sum(gamma) + alpha_theta, length(gamma) - sum(gamma) + beta_theta)
    list(Phi = coef, Sigma = sigma, gamma = matrix(if (map) as.integer(gamma >= 0.5) else gamma, d, k), theta = next_theta)
  }
  sweep <- function(state, map = FALSE, warm = FALSE) {
    a <- update_mode(work(state, "A"), state$SA, do.call(cbind, state$A), state$theta_A,
      state$scale, tau0_A, tau1_A, nu_A, if (fixed_theta) fixed_A else if (warm && warm_fixed_theta) state$theta_A else NULL, map)
    state$A <- split_lags(a$Phi, G); state$SA <- a$Sigma; state$theta_A <- a$theta
    b <- update_mode(work(state, "B"), state$SB, do.call(cbind, state$B), state$theta_B,
      state$scale, tau0_B, tau1_B, nu_B, if (fixed_theta) fixed_B else if (warm && warm_fixed_theta) state$theta_B else NULL, map)
    state$B <- split_lags(b$Phi, N); state$SB <- b$Sigma; state$theta_B <- b$theta
    state$gamma_A <- split_lags(a$gamma, G); state$gamma_B <- split_lags(b$gamma, N)
    if (scale_update == "before_rescale") state$scale <- update_scale(state, map)
    state <- balance(state)
    if (scale_update == "after_rescale") state$scale <- update_scale(state, map)
    state
  }
  if (!is.logical(warm_fixed_theta) || length(warm_fixed_theta) != 1L || is.na(warm_fixed_theta))
    stop("warm_fixed_theta must be TRUE or FALSE.")
  state <- list(A = initial_A, B = initial_B, SA = initial_SA, SB = initial_SB,
                scale = scale_init, theta_A = theta_A, theta_B = theta_B)
  if (!is.null(seed)) set.seed(seed)
  if (draw_initial_indicators) {
    state$gamma_A <- split_lags(matrix(rbinom(G * G * P, 1, as.numeric(theta_A)), G, G * P), G)
    state$gamma_B <- split_lags(matrix(rbinom(N * N * P, 1, as.numeric(theta_B)), N, N * P), N)
  }
  if (!is.null(warm_seed)) set.seed(warm_seed)
  if (init_method != "none") for (i in seq_len(warm_iter))
    state <- sweep(state, map = init_method == "map", warm = TRUE)
  if (!is.null(sampling_seed)) set.seed(sampling_seed)
  n_post <- n_iter - burnin
  A_chain <- array(NA_real_, c(n_post, P, G, G)); B_chain <- array(NA_real_, c(n_post, P, N, N))
  gamma_A_chain <- array(FALSE, dim(A_chain)); gamma_B_chain <- array(FALSE, dim(B_chain))
  scale_chain <- theta_A_chain <- theta_B_chain <- numeric(n_post)
  if (!is.null(mcmc_clock)) {
    mcmc_clock$start()
    on.exit(mcmc_clock$cancel(), add = TRUE)
  }
  started <- Sys.time()
  for (iter in seq_len(n_iter)) {
    state <- sweep(state)
    if (iter > burnin) {
      s <- iter - burnin
      for (lag in seq_len(P)) {
        A_chain[s, lag, , ] <- state$A[[lag]]; B_chain[s, lag, , ] <- state$B[[lag]]
        gamma_A_chain[s, lag, , ] <- state$gamma_A[[lag]] == 1L
        gamma_B_chain[s, lag, , ] <- state$gamma_B[[lag]] == 1L
      }
      scale_chain[s] <- state$scale
      theta_A_chain[s] <- mean(state$theta_A); theta_B_chain[s] <- mean(state$theta_B)
    }
    if (!is.null(verbose_every) && iter %% verbose_every == 0L) cat("SSVS MAR(", P, ") ", iter, "/", n_iter, "\n", sep = "")
  }
  runtime <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  if (!is.null(mcmc_clock)) mcmc_clock$stop()
  a <- summarize(A_chain, gamma_A_chain, "A"); b <- summarize(B_chain, gamma_B_chain, "B")
  out <- list(method = "SSVS", P = P, A_hat = a$hat, B_hat = b$hat, PIP_A = a$PIP, PIP_B = b$PIP,
    A_mean_active = a$active, B_mean_active = b$active, A_mean_all = a$mean, B_mean_all = b$mean,
    Sigma_A_last = state$SA, Sigma_B_last = state$SB, scale_last = state$scale, runtime_sec = runtime,
    settings = list(P = P, n_iter = n_iter, burnin = burnin, warm_iter = warm_iter,
      init_method = init_method, tau0_A = tau0_A, tau1_A = tau1_A, tau0_B = tau0_B, tau1_B = tau1_B,
      prior_A = priorA, prior_B = priorB, fixed_theta = fixed_theta, pip_cut = pip_cut,
      pip_rule = pip_rule, numerics = numerics, whitening = whitening,
      scale_update = scale_update, rescale = rescale))
  if (keep_draws) {
    out$A_chain <- A_chain; out$B_chain <- B_chain
    out$gamma_A_chain <- gamma_A_chain; out$gamma_B_chain <- gamma_B_chain
    out$scale_chain <- scale_chain; out$theta_A_chain <- theta_A_chain; out$theta_B_chain <- theta_B_chain
  }
  out
}
