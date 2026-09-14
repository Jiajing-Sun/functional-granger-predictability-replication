# ==============================================================
# Simulation + estimation functions (lambda=20 grid) for
#   Relevant Functional Granger Causality with SN
#   Effect size: d_{uX} = || C_{uX} ||_HS^2  (p=1)
#
# Includes:
#   - Relevant SN tests (quadratic + adjusted-range)
#   - Plug-in point-null test for Delta=0 (classical "no causality")
#     via diagonal long-run variance (Newey-West, Bartlett) and a
#     Satterthwaite chi-square approximation.
#
# Base R only. Works on Windows/Mac/Linux.
# ==============================================================

# ---------- utilities ----------
safe_seed <- function(seed) {
  if (!is.null(seed)) set.seed(seed)
}

# Generate equation disturbances with matched marginal variance sigma^2.
# The fMA1 case is serially correlated and is used as a disturbance stress
# design; it is not an MDS-innovation design.
sim_innovations <- function(T, D, sigma_vec,
                            type = c("iid","fMA1","t5"),
                            psi = 0.7,
                            seed = NULL) {
  type <- match.arg(type)
  safe_seed(seed)
  if (length(sigma_vec) != D) stop("sigma_vec must have length D.")
  if (type == "iid") {
    mat <- matrix(rnorm(T * D), nrow = T, ncol = D)
    mat <- sweep(mat, 2, sigma_vec, "*")
    return(mat)
  }
  if (type == "fMA1") {
    # e_t = z_t + psi z_{t-1}, with Var(e_t)=sigma^2 by scaling z_t
    sd_z <- sigma_vec / sqrt(1 + psi^2)
    Z <- matrix(rnorm((T + 1L) * D), nrow = (T + 1L), ncol = D)
    Z <- sweep(Z, 2, sd_z, "*")
    E <- Z[2:(T+1L), , drop = FALSE] + psi * Z[1:T, , drop = FALSE]
    return(E)
  }
  if (type == "t5") {
    df <- 5
    # rt has variance df/(df-2); scale to variance 1 then multiply sigma
    scale_t <- sqrt((df - 2) / df)
    mat <- matrix(rt(T * D, df = df), nrow = T, ncol = D) * scale_t
    mat <- sweep(mat, 2, sigma_vec, "*")
    return(mat)
  }
  stop("Unknown type.")
}

# FARX(1) in coefficient space (diagonal operators)
sim_FARX1_diag <- function(n,
                           D = 21L,
                           theta = 0,
                           a_y = 0.4,
                           c_x = 0.5,
                           b_vec = 1/(1:D),
                           sigma_vec = 1/(1:D),
                           error_type = c("iid","fMA1","t5"),
                           psi = 0.7,
                           burn = 200L,
                           seed = NULL) {
  error_type <- match.arg(error_type)
  if (length(b_vec) != D) stop("b_vec must have length D.")
  if (length(sigma_vec) != D) stop("sigma_vec must have length D.")
  total <- as.integer(n + burn)

  # Equation disturbances for X and Y (serially correlated under fMA1)
  eX <- sim_innovations(total, D, sigma_vec, type = error_type, psi = psi,
                        seed = if (is.null(seed)) NULL else seed + 11L)
  eY <- sim_innovations(total, D, sigma_vec, type = error_type, psi = psi,
                        seed = if (is.null(seed)) NULL else seed + 29L)

  x <- matrix(0, nrow = total + 1L, ncol = D)
  y <- matrix(0, nrow = total + 1L, ncol = D)

  # recursion
  for (t in 2:(total + 1L)) {
    x[t, ] <- c_x * x[t-1L, ] + eX[t-1L, ]
    y[t, ] <- a_y * y[t-1L, ] + theta * (b_vec * x[t-1L, ]) + eY[t-1L, ]
  }

  X <- x[(burn + 2L):(total + 1L), , drop = FALSE]
  Y <- y[(burn + 2L):(total + 1L), , drop = FALSE]
  list(X = X, Y = Y)
}

# ---------- design-specific population relevance boundaries ----------

# Solve P = A P A' + W by vectorization. This base-R implementation is
# sufficient for the four-dimensional augmented state used below.
solve_discrete_lyapunov <- function(A, W, tol = 1e-10) {
  if (!is.matrix(A) || nrow(A) != ncol(A)) stop("A must be square.")
  if (!is.matrix(W) || any(dim(W) != dim(A))) {
    stop("W must have the same dimensions as A.")
  }
  d <- nrow(A)
  spectral_radius <- max(Mod(eigen(A, only.values = TRUE)$values))
  if (!is.finite(spectral_radius) || spectral_radius >= 1) {
    stop(sprintf("State matrix is not stable: spectral radius=%.8g.",
                 spectral_radius))
  }
  lhs <- diag(d * d) - kronecker(A, A)
  P <- matrix(solve(lhs, as.vector(W)), nrow = d, ncol = d)
  P <- (P + t(P)) / 2
  residual <- max(abs(P - (A %*% P %*% t(A) + W)))
  scale <- max(1, max(abs(P)))
  if (!is.finite(residual) || residual > tol * scale) {
    stop(sprintf("Discrete Lyapunov residual %.3e exceeds tolerance.", residual))
  }
  min_eigenvalue <- min(eigen(P, symmetric = TRUE, only.values = TRUE)$values)
  if (!is.finite(min_eigenvalue) || min_eigenvalue < -tol * scale) {
    stop(sprintf("Stationary state covariance is not positive semidefinite: min eigenvalue=%.3e.",
                 min_eigenvalue))
  }
  attr(P, "lyapunov_residual") <- residual
  attr(P, "spectral_radius") <- spectral_radius
  attr(P, "min_eigenvalue") <- min_eigenvalue
  P
}

# Existing white-disturbance closed form. It applies to both iid Gaussian and
# iid t5 designs because those designs have the same matched second moments.
duX_white_noise_closed_form <- function(theta,
                                        a_y,
                                        c_x,
                                        b_vec,
                                        sigma_vec) {
  if (length(b_vec) != length(sigma_vec)) {
    stop("b_vec and sigma_vec must have the same length.")
  }
  var_x <- sigma_vec^2 / (1 - c_x^2)
  cov_xy <- (c_x * theta * b_vec * var_x) / (1 - a_y * c_x)
  var_y <- (theta^2 * b_vec^2 * var_x * (1 + a_y * c_x) /
              (1 - a_y * c_x) + sigma_vec^2) / (1 - a_y^2)
  var_x_given_y <- var_x - cov_xy^2 / var_y
  partial_cov <- theta * b_vec * var_x_given_y
  list(
    d_uX = sum(partial_cov^2),
    components = data.frame(
      coordinate = seq_along(b_vec),
      var_x = var_x,
      var_y = var_y,
      cov_xy = cov_xy,
      cov_y_t_x_lag = a_y * cov_xy + theta * b_vec * var_x,
      cov_y_t_y_lag = a_y * var_y + theta * b_vec * cov_xy,
      partial_cov = partial_cov,
      stringsAsFactors = FALSE
    )
  )
}

# Exact second moments for the serially correlated fMA(1) disturbance design.
# For each diagonal coordinate, the augmented state is
#   s_t = (x_t, y_t, z^X_t, z^Y_t)',
# where e^j_t = z^j_t + psi z^j_{t-1}. The state covariance is obtained from
# a discrete Lyapunov equation, and C_{uX} is then the population partial
# covariance of y_t and x_{t-1} after projecting both on y_{t-1}.
duX_fma1_augmented_state <- function(theta,
                                     a_y,
                                     c_x,
                                     b_vec,
                                     sigma_vec,
                                     psi = 0.7) {
  if (length(b_vec) != length(sigma_vec)) {
    stop("b_vec and sigma_vec must have the same length.")
  }
  B <- rbind(c(1, 0), c(0, 1), c(1, 0), c(0, 1))
  components <- vector("list", length(b_vec))

  for (k in seq_along(b_vec)) {
    g <- theta * b_vec[k]
    A <- rbind(
      c(c_x, 0, psi, 0),
      c(g, a_y, 0, psi),
      c(0, 0, 0, 0),
      c(0, 0, 0, 0)
    )
    var_z <- sigma_vec[k]^2 / (1 + psi^2)
    W <- B %*% (diag(var_z, 2L)) %*% t(B)
    P <- solve_discrete_lyapunov(A, W)
    lag_cov <- A %*% P

    var_x <- P[1, 1]
    var_y <- P[2, 2]
    cov_xy <- P[1, 2]
    cov_y_t_x_lag <- lag_cov[2, 1]
    cov_y_t_y_lag <- lag_cov[2, 2]
    partial_cov <- cov_y_t_x_lag - cov_y_t_y_lag * cov_xy / var_y
    # Independent algebraic check using Cov(e^Y_t, y_{t-1}) = psi Var(z^Y_t).
    partial_cov_equivalent <-
      g * (var_x - cov_xy^2 / var_y) - psi * var_z * cov_xy / var_y
    equivalence_error <- partial_cov - partial_cov_equivalent

    components[[k]] <- data.frame(
      coordinate = k,
      var_x = var_x,
      var_y = var_y,
      cov_xy = cov_xy,
      cov_y_t_x_lag = cov_y_t_x_lag,
      cov_y_t_y_lag = cov_y_t_y_lag,
      partial_cov = partial_cov,
      partial_cov_equivalent = partial_cov_equivalent,
      equivalence_error = equivalence_error,
      stringsAsFactors = FALSE
    )
  }

  components <- do.call(rbind, components)
  if (max(abs(components$equivalence_error)) > 1e-10) {
    stop(sprintf("Partial-covariance identity check failed: max error=%.3e.",
                 max(abs(components$equivalence_error))))
  }
  list(
    d_uX = sum(components$partial_cov^2),
    d_uX_equivalent = sum(components$partial_cov_equivalent^2),
    components = components
  )
}

# Construct one relevance boundary per disturbance design and run two
# deterministic checks:
#   (i) the augmented-state formula at psi=0 equals the white-noise formula;
#  (ii) an independent partial-covariance identity reproduces each
#       design-specific d_uX(theta0) boundary;
# (iii) theta=0 gives a zero population effect for white and fMA1 designs.
design_specific_duX_thresholds <- function(theta0,
                                           a_y,
                                           c_x,
                                           b_vec,
                                           sigma_vec,
                                           psi = 0.7,
                                           tol = 1e-10) {
  white <- duX_white_noise_closed_form(
    theta = theta0, a_y = a_y, c_x = c_x,
    b_vec = b_vec, sigma_vec = sigma_vec
  )
  state_psi0 <- duX_fma1_augmented_state(
    theta = theta0, a_y = a_y, c_x = c_x,
    b_vec = b_vec, sigma_vec = sigma_vec, psi = 0
  )
  fma1 <- duX_fma1_augmented_state(
    theta = theta0, a_y = a_y, c_x = c_x,
    b_vec = b_vec, sigma_vec = sigma_vec, psi = psi
  )
  white_zero <- duX_white_noise_closed_form(
    theta = 0, a_y = a_y, c_x = c_x,
    b_vec = b_vec, sigma_vec = sigma_vec
  )
  fma1_zero <- duX_fma1_augmented_state(
    theta = 0, a_y = a_y, c_x = c_x,
    b_vec = b_vec, sigma_vec = sigma_vec, psi = psi
  )

  psi0_error <- state_psi0$d_uX - white$d_uX
  if (!is.finite(psi0_error) || abs(psi0_error) > tol) {
    stop(sprintf(
      "psi=0 augmented-state check failed: error=%.3e (tol=%.3e).",
      psi0_error, tol
    ))
  }
  zero_signal_error <- max(abs(c(white_zero$d_uX, fma1_zero$d_uX)))
  if (!is.finite(zero_signal_error) || zero_signal_error > tol) {
    stop(sprintf(
      "theta=0 population-effect check failed: max effect=%.3e (tol=%.3e).",
      zero_signal_error, tol
    ))
  }

  out <- data.frame(
    error = c("iid", "fMA1", "t5"),
    theta0 = rep(theta0, 3L),
    disturbance_psi = c(0, psi, 0),
    Delta = c(white$d_uX, fma1$d_uX, white$d_uX),
    d_uX_independent_check_at_theta0 = c(
      state_psi0$d_uX_equivalent,
      fma1$d_uX_equivalent,
      state_psi0$d_uX_equivalent
    ),
    second_moment_method = c(
      "white-noise closed form (matched Gaussian variance)",
      "augmented-state discrete Lyapunov (serial fMA1 disturbance)",
      "white-noise closed form (matched t5 variance)"
    ),
    stringsAsFactors = FALSE
  )
  out$boundary_error <- out$d_uX_independent_check_at_theta0 - out$Delta
  if (any(!is.finite(out$boundary_error)) ||
      max(abs(out$boundary_error)) > tol) {
    stop(sprintf(
      "Design-specific boundary check failed: max error=%.3e (tol=%.3e).",
      max(abs(out$boundary_error)), tol
    ))
  }
  attr(out, "psi0_white_noise_error") <- psi0_error
  attr(out, "zero_signal_error") <- zero_signal_error
  out
}

# Ridge (Tikhonov) multi-response regression: Y ≈ X %*% Beta
fit_ridge_multi <- function(X, Y, lambda) {
  n <- nrow(X)
  p <- ncol(X)
  XtX <- crossprod(X) / n
  XtY <- crossprod(X, Y) / n
  Beta <- tryCatch(
    solve(XtX + lambda * diag(p), XtY),
    error = function(e) qr.solve(XtX + lambda * diag(p), XtY)
  )
  Beta
}

# Restricted model residuals (only Y lag)
restricted_residuals <- function(Y, lambda_fun) {
  n <- nrow(Y)
  D <- ncol(Y)
  if (n < 2) stop("Need n>=2.")
  Ylag <- Y[1:(n-1L), , drop = FALSE]
  Ycur <- Y[2:n, , drop = FALSE]
  lambda <- lambda_fun(n - 1L)
  Beta <- fit_ridge_multi(Ylag, Ycur, lambda)  # D x D
  U <- Ycur - Ylag %*% Beta
  list(U = U, Ylag = Ylag, Beta = Beta, lambda = lambda)
}

# Restricted residuals for both the response and lagged predictor.
restricted_residuals_uv <- function(X, Y, lambda_fun) {
  n <- nrow(Y)
  if (nrow(X) != n) stop("X and Y must have the same number of rows.")
  if (n < 2) stop("Need n>=2.")
  Ylag <- Y[1:(n-1L), , drop = FALSE]
  Ycur <- Y[2:n, , drop = FALSE]
  Xlag <- X[1:(n-1L), , drop = FALSE]
  lambda <- lambda_fun(n - 1L)
  BetaY <- fit_ridge_multi(Ylag, Ycur, lambda)
  BetaX <- fit_ridge_multi(Ylag, Xlag, lambda)
  U <- Ycur - Ylag %*% BetaY
  V <- Xlag - Ylag %*% BetaX
  list(U = U, V = V, Ylag = Ylag, Xlag = Xlag,
       BetaY = BetaY, BetaX = BetaX, lambda = lambda)
}

# Sequential estimator for d_{uX} with RELEVANT scaling:
#   C_hat(k) = (1/n) sum_{t=2}^k u_t \otimes v_{t-1}, denom fixed at n
#   d_hat(k) = ||C_hat(k)||_F^2  (Frobenius = HS under orthonormal basis)
estimate_duX_relevant <- function(X, Y, k, denom_n, lambda_fun) {
  if (k < 2) return(NA_real_)
  Xk <- X[1:k, , drop = FALSE]
  Yk <- Y[1:k, , drop = FALSE]
  res <- restricted_residuals_uv(Xk, Yk, lambda_fun)
  U <- res$U  # (k-1) x D
  V <- res$V  # (k-1) x D, residualized lagged predictor
  N <- nrow(U)
  Chat <- crossprod(U, V) / N
  diagonal <- mean(rowSums(U^2) * rowSums(V^2))
  ((N * sum(Chat^2) - diagonal) / (N - 1)) * (N / denom_n)^2
}

# Relevant SN test for d_{uX}
relevant_SN_test_duX <- function(X, Y,
                                Delta,
                                q_quad, q_range,
                                m = 20L,
                                lambda_fun = function(n) n^(-1),
                                return_path = FALSE) {
  n <- nrow(Y)
  if (m < 3) stop("m must be >= 3.")
  lambdas <- 0.25 + 0.75 * (1:(m-1L))/m

  d_full <- estimate_duX_relevant(X, Y, k = n, denom_n = n, lambda_fun = lambda_fun)
  d_sub <- numeric(length(lambdas))

  for (i in seq_along(lambdas)) {
    k <- floor(n * lambdas[i])
    if (k < 2) {
      stop(sprintf("Prefix lambda=%g has fewer than two observations.", lambdas[i]))
    }
    d_sub[i] <- estimate_duX_relevant(X, Y, k = k, denom_n = n, lambda_fun = lambda_fun)
  }

  if (!is.finite(d_full) || any(!is.finite(d_sub))) {
    stop("All full-sample and prefix effect estimates must be finite.")
  }

  # centered subsample path (relevant SN structure)
  G <- d_sub - (lambdas^2) * d_full

  # quadratic SN: discrete-uniform nu on lambdas -> mean
  V <- sqrt(mean(G^2))

  # adjusted-range SN
  H <- max(G) - min(G)

  if (any(!is.finite(G)) || !is.finite(V) || !is.finite(H) || V <= 0 || H <= 0) {
    stop("The complete fixed-grid self-normalization path must be finite and nondegenerate.")
  }

  rej_quad  <- as.integer(d_full > Delta + q_quad  * V)
  rej_range <- as.integer(d_full > Delta + q_range * H)

  out <- list(
    reject_quad = rej_quad,
    reject_range = rej_range,
    d_hat = d_full,
    V = V,
    H = H
  )
  if (return_path) out$G <- G
  out
}

# Plug-in point-null test: Delta=0 ("no causality") for d_{uX}
# We test H0: C_{uX}=0 by the quadratic form:
#   T = N * ||C_hat||_F^2,  C_hat = (1/N) sum u_t \otimes v_{t-1}
# Under H0, vec(sqrt(N) C_hat) approx N(0, Omega). With diagonal LRV,
# T approx sum_j sigma_j^2 Z_j^2. We use Satterthwaite approximation.
plugin_test_duX_Delta0 <- function(X, Y,
                                  alpha = 0.05,
                                  lambda_fun = function(n) n^(-1),
                                  L_bw = NULL) {
  n <- nrow(Y)
  D <- ncol(Y)
  if (n < 3) stop("Need n>=3.")
  res <- restricted_residuals_uv(X, Y, lambda_fun)
  U <- res$U               # (n-1) x Dy
  V <- res$V               # (n-1) x Dx
  N <- nrow(U)             # N = n-1
  Dy <- ncol(U)
  Dx <- ncol(V)

  # build W_t = vec(u_t v_{t-1}^T) as blocks: each block is u_t * v_{t-1,j}
  W <- matrix(0, nrow = N, ncol = Dy * Dx)
  col_start <- 1L
  for (j in 1:Dx) {
    W[, col_start:(col_start + Dy - 1L)] <- U * V[, j]
    col_start <- col_start + Dy
  }

  w_bar <- colMeans(W)
  wc <- sweep(W, 2, w_bar, "-")

  # bandwidth for scalar Newey-West on each component
  if (is.null(L_bw)) {
    L <- max(1L, floor(N^(1/4)))
  } else {
    L <- as.integer(L_bw)
    L <- max(0L, min(L, N-1L))
  }

  gamma0 <- colSums(wc^2) / N
  lrv <- gamma0

  if (L >= 1L) {
    for (h in 1:L) {
      w_h <- 1 - h/(L + 1)
      gam_h <- colSums(wc[(h+1L):N, , drop = FALSE] * wc[1:(N-h), , drop = FALSE]) / N
      lrv <- lrv + 2 * w_h * gam_h
    }
  }
  lrv <- pmax(lrv, 0)

  # Satterthwaite: T ≈ scale * chisq_df
  mu1 <- sum(lrv)
  mu2 <- sum(lrv^2)
  # Guard against numerical issues
  if (mu1 <= 0 || mu2 <= 0) {
    return(list(reject = NA_integer_, stat = NA_real_, cv = NA_real_,
                df = NA_real_, scale = NA_real_, L = L))
  }

  df <- (mu1^2) / mu2
  scale <- mu2 / mu1
  cv <- scale * qchisq(1 - alpha, df = df)

  stat <- N * sum(w_bar^2)   # = || sqrt(N) * mean(W) ||^2
  reject <- as.integer(stat > cv)

  list(reject = reject, stat = stat, cv = cv, df = df, scale = scale, L = L)
}
