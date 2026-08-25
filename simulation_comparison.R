#!/usr/bin/env Rscript

## Reproducible comparison study for anomalous slopes in fixed-effects panels.
## Usage:
##   Rscript simulation_comparison.R smoke
##   Rscript simulation_comparison.R pilot
##   Rscript simulation_comparison.R final

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) tolower(args[[1]]) else "smoke"
if (!mode %in% c("smoke", "pilot", "final")) {
  stop("mode must be one of: smoke, pilot, final")
}

options(stringsAsFactors = FALSE, scipen = 99)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/",
                             mustWork = TRUE)
## Use the script directory as the project directory in the submission bundle.
study_dir <- dirname(script_path)
output_dir <- file.path(study_dir, mode)

safe_solve <- function(A, b = NULL, ridge = 1e-10) {
  A <- (A + t(A)) / 2
  ans <- tryCatch(
    if (is.null(b)) solve(A) else solve(A, b),
    error = function(e) NULL
  )
  if (!is.null(ans)) return(ans)
  A <- A + diag(ridge * max(1, mean(abs(diag(A)))), nrow(A))
  if (is.null(b)) qr.solve(A) else qr.solve(A, b)
}

toeplitz_cov <- function(p, rho = 0.5) {
  outer(seq_len(p), seq_len(p), function(j, k) rho^abs(j - k))
}

rmvn_rows <- function(nr, Sigma) {
  matrix(rnorm(nr * nrow(Sigma)), nr, nrow(Sigma)) %*% chol(Sigma)
}

generate_series <- function(T, p, design, Sigma) {
  if (design == "iid") return(rmvn_rows(T, Sigma))
  if (design == "ma") {
    e <- rmvn_rows(T + 2L, Sigma)
    return(e[3:(T + 2L), , drop = FALSE] +
             0.5 * e[2:(T + 1L), , drop = FALSE] +
             0.3 * e[1:T, , drop = FALSE])
  }
  if (design == "ar") {
    burn <- 50L
    total <- T + burn
    e <- rmvn_rows(total, Sigma)
    w <- matrix(0, total, p)
    w[1, ] <- e[1, ]
    w[2, ] <- 0.4 * w[1, ] + e[2, ]
    for (tt in 3:total) {
      w[tt, ] <- 0.4 * w[tt - 1L, ] - 0.6 * w[tt - 2L, ] + e[tt, ]
    }
    return(w[(burn + 1L):total, , drop = FALSE])
  }
  stop("Unknown predictor design: ", design)
}

make_directions <- function(n_out, K, pattern) {
  if (n_out == 0L) return(matrix(0, 0L, K))
  if (pattern == "spherical") {
    U <- matrix(rnorm(n_out * K), n_out, K)
    return(U / sqrt(rowSums(U^2)))
  }
  if (pattern == "aligned") {
    u <- rep(1 / sqrt(K), K)
    return(matrix(rep(u, each = n_out), n_out, K))
  }
  if (pattern == "sparse") {
    U <- matrix(0, n_out, K)
    coords <- rep(seq_len(K), length.out = n_out)
    signs <- rep(c(-1, 1), length.out = n_out)
    U[cbind(seq_len(n_out), coords)] <- signs
    return(U)
  }
  stop("Unknown direction pattern: ", pattern)
}

generate_panel <- function(scenario, zeta, seed) {
  set.seed(seed)
  n <- scenario$n
  T <- scenario$T
  K <- scenario$K
  m <- scenario$m
  beta <- c(4, 3, 2, 1)[seq_len(K)]
  theta <- rep(1, m)
  n_out <- if (zeta > 0) floor(n * scenario$rho) else 0L
  ## Match simulationcode.r: the first floor(n*rho) units are heterogeneous.
  out_idx <- if (n_out > 0L) seq_len(n_out) else integer(0)
  is_out <- seq_len(n) %in% out_idx

  gamma <- matrix(0, n, K)
  if (n_out > 0L) {
    U <- make_directions(n_out, K, scenario$direction)
    strength <- runif(n_out, scenario$signal_lower * zeta, zeta)
    gamma[out_idx, ] <- U * (strength * sqrt(sum(beta^2)))
  }

  Sigma_x <- toeplitz_cov(K)
  Sigma_z <- toeplitz_cov(m)
  X <- vector("list", n)
  Z <- vector("list", n)
  Y <- matrix(0, n, T)
  c_i <- runif(n, -1, 1)
  c_i[n] <- -sum(c_i[-n])

  for (i in seq_len(n)) {
    Xi <- generate_series(T, K, scenario$predictor, Sigma_x)
    Zi <- generate_series(T, m, scenario$predictor, Sigma_z)
    if (scenario$error == "normal") {
      eps <- rnorm(T, 0, 0.5)
    } else if (scenario$error == "mixnorm") {
      eta <- rnorm(T, ifelse(rbinom(T, 1, 0.5) == 1, 0.5, -0.5), 1)
      eps <- 0.5 * eta / sqrt(1.25)
    } else {
      stop("Unknown error design: ", scenario$error)
    }
    X[[i]] <- Xi
    Z[[i]] <- Zi
    Y[i, ] <- c_i[i] + as.vector(Xi %*% (beta + gamma[i, ])) +
      as.vector(Zi %*% theta) + eps
  }

  list(
    X = X, Z = Z, Y = Y, beta = beta, theta = theta,
    gamma = gamma, out_idx = out_idx, is_out = is_out,
    n = n, T = T, K = K, m = m
  )
}

prepare_panel <- function(dat) {
  n <- dat$n
  K <- dat$K
  m <- dat$m
  Xw <- vector("list", n)
  Zw <- vector("list", n)
  yw <- vector("list", n)
  Sxx <- vector("list", n)
  beta_projection_Z <- vector("list", n)
  beta_projection_y <- vector("list", n)
  sum_ztz <- matrix(0, m, m)
  sum_zty <- numeric(m)

  for (i in seq_len(n)) {
    Xi <- scale(dat$X[[i]], center = TRUE, scale = FALSE)
    Zi <- scale(dat$Z[[i]], center = TRUE, scale = FALSE)
    yi <- dat$Y[i, ] - mean(dat$Y[i, ])
    A <- crossprod(Xi)
    Ainv <- safe_solve(A)
    rz <- Zi - Xi %*% (Ainv %*% crossprod(Xi, Zi))
    ry <- yi - as.vector(Xi %*% (Ainv %*% crossprod(Xi, yi)))
    sum_ztz <- sum_ztz + crossprod(rz)
    sum_zty <- sum_zty + as.vector(crossprod(rz, ry))
    Xw[[i]] <- Xi
    Zw[[i]] <- Zi
    yw[[i]] <- yi
    Sxx[[i]] <- A
  }

  theta_hat <- as.vector(safe_solve(sum_ztz, sum_zty))
  beta_ind <- matrix(0, n, K)
  sxy <- vector("list", n)
  syy <- numeric(n)
  rss_ind <- numeric(n)
  yadj <- vector("list", n)
  for (i in seq_len(n)) {
    ai <- yw[[i]] - as.vector(Zw[[i]] %*% theta_hat)
    bi <- as.vector(safe_solve(Sxx[[i]], crossprod(Xw[[i]], ai)))
    xy <- as.vector(crossprod(Xw[[i]], ai))
    yy <- sum(ai^2)
    beta_ind[i, ] <- bi
    sxy[[i]] <- xy
    syy[i] <- yy
    rss_ind[i] <- max(0, yy - sum(xy * bi))
    yadj[[i]] <- ai
  }

  list(
    X = Xw, Z = Zw, y = yw, yadj = yadj, theta = theta_hat,
    Sxx = Sxx, sxy = sxy, syy = syy, rss_ind = rss_ind,
    beta_ind = beta_ind, n = dat$n, T = dat$T, K = K, m = m
  )
}

sum_list_matrices <- function(x, idx) {
  Reduce(`+`, x[idx])
}

sum_list_vectors <- function(x, idx) {
  Reduce(`+`, x[idx])
}

pooled_beta <- function(pre, idx) {
  if (!length(idx)) return(rep(NA_real_, pre$K))
  A <- sum_list_matrices(pre$Sxx, idx)
  b <- sum_list_vectors(pre$sxy, idx)
  as.vector(safe_solve(A, b))
}

pooled_joint_beta <- function(pre, idx) {
  if (!length(idx)) return(rep(NA_real_, pre$K))
  p <- pre$K + pre$m
  A <- matrix(0, p, p)
  b <- numeric(p)
  for (i in idx) {
    W <- cbind(pre$X[[i]], pre$Z[[i]])
    A <- A + crossprod(W)
    b <- b + as.vector(crossprod(W, pre$y[[i]]))
  }
  as.vector(safe_solve(A, b))[seq_len(pre$K)]
}

rss_at_beta <- function(pre, idx, beta) {
  sum(vapply(idx, function(i) {
    max(0, pre$syy[i] - 2 * sum(beta * pre$sxy[[i]]) +
          as.numeric(crossprod(beta, pre$Sxx[[i]] %*% beta)))
  }, numeric(1)))
}

bh_select <- function(p, alpha = 0.05) {
  which(p.adjust(p, method = "BH") <= alpha)
}

lm_method <- function(pre, ref_idx, alpha = 0.05) {
  beta_hat <- pooled_beta(pre, ref_idx)
  rss_ref <- rss_at_beta(pre, ref_idx, beta_hat)
  sigma2 <- rss_ref / max(1, length(ref_idx) * (pre$T - 1))
  stats <- numeric(pre$n)
  for (i in seq_len(pre$n)) {
    resid <- pre$yadj[[i]] - as.vector(pre$X[[i]] %*% beta_hat)
    score <- as.vector(crossprod(pre$X[[i]], resid))
    stats[i] <- max(0, as.numeric(crossprod(score, safe_solve(pre$Sxx[[i]], score))) /
                      max(sigma2, 1e-12))
  }
  p <- pchisq(stats, df = pre$K, lower.tail = FALSE)
  list(selected = bh_select(p, alpha), beta_test = beta_hat,
       p = p, stat = stats, sigma2 = sigma2)
}

proposed_method <- function(pre, tau_c = 0.50, k_c = 0.60, alpha = 0.05) {
  clean_info <- screen_clean_set(pre, tau_c = tau_c, k_c = k_c)
  clean_pre <- clean_info$indices
  ans <- lm_method(pre, clean_pre, alpha)
  ans$clean_pre <- clean_pre
  ans$c0 <- clean_info$c0
  ans$neighborhood_fraction <- clean_info$neighborhood_fraction
  ans
}

screen_clean_set <- function(pre, tau_c = 0.50, k_c = 0.60) {
  dm <- as.matrix(dist(pre$beta_ind))
  c0 <- as.numeric(quantile(dm[lower.tri(dm)], probs = tau_c, names = FALSE))
  neigh <- (rowSums(dm <= c0) - 1) / max(pre$n - 1L, 1L)
  ord <- order(-neigh, seq_len(pre$n))
  list(
    indices = ord[seq_len(floor(k_c * pre$n))],
    c0 = c0,
    neighborhood_fraction = neigh
  )
}

joint_fit <- function(pre, idx) {
  p <- pre$K + pre$m
  A <- matrix(0, p, p)
  b <- numeric(p)
  for (i in idx) {
    W <- cbind(pre$X[[i]], pre$Z[[i]])
    A <- A + crossprod(W)
    b <- b + as.vector(crossprod(W, pre$y[[i]]))
  }
  Ainv <- safe_solve(A)
  coef <- as.vector(Ainv %*% b)
  rss <- sum(vapply(idx, function(i) {
    W <- cbind(pre$X[[i]], pre$Z[[i]])
    sum((pre$y[[i]] - as.vector(W %*% coef))^2)
  }, numeric(1)))
  df <- max(1, length(idx) * (pre$T - 1L) - p)
  list(coef = coef, beta = coef[seq_len(pre$K)],
       theta = coef[pre$K + seq_len(pre$m)], Ainv = Ainv,
       sigma2 = rss / df, df = df)
}

make_slope_score_model <- function(pre, train_idx) {
  fit <- joint_fit(pre, train_idx)
  g <- do.call(rbind, lapply(train_idx, function(i) {
    W <- cbind(pre$X[[i]], pre$Z[[i]])
    resid <- pre$y[[i]] - as.vector(W %*% fit$coef)
    pre$X[[i]] * resid
  }))
  omega <- stats::cov(g)
  ridge <- 1e-6 * max(mean(diag(omega)), 1e-8)
  omega <- omega + diag(ridge, pre$K)
  list(fit = fit, omega = omega, omega_inv = safe_solve(omega))
}

slope_nonconformity_scores <- function(pre, score_model, idx = seq_len(pre$n)) {
  lapply(idx, function(i) {
    W <- cbind(pre$X[[i]], pre$Z[[i]])
    resid <- pre$y[[i]] - as.vector(W %*% score_model$fit$coef)
    g <- pre$X[[i]] * resid
    pmax(0, rowSums((g %*% score_model$omega_inv) * g))
  })
}

batch_conformal_weights <- function(n_ref, batch_size, eta) {
  r <- seq_len(n_ref)
  log_denom <- lchoose(n_ref + batch_size, batch_size)
  log_w <- lchoose(r + eta - 2L, eta - 1L) +
    lchoose(n_ref + batch_size - r - eta + 1L, batch_size - eta) -
    log_denom
  w <- exp(log_w)
  base <- exp(lchoose(n_ref + eta - 1L, eta - 1L) - log_denom)
  ## Numerical normalization enforces the exact identity sum(w) + base = 1.
  total <- sum(w) + base
  list(tail = rev(cumsum(rev(w / total))), base = base / total)
}

batch_conformal_pvalues <- function(reference_scores, test_scores, eta) {
  reference_scores <- sort(as.numeric(reference_scores))
  n_ref <- length(reference_scores)
  batch_size <- length(test_scores[[1]])
  weights <- batch_conformal_weights(n_ref, batch_size, eta)
  vapply(test_scores, function(s) {
    s_eta <- sort(s, partial = eta)[eta]
    first_ge <- findInterval(s_eta, reference_scores, left.open = TRUE) + 1L
    if (first_ge > n_ref) weights$base else
      min(1, weights$tail[first_ge] + weights$base)
  }, numeric(1))
}

batch_conformal_method <- function(pre, pre_ref, dat_ref,
                                   reference = c("oracle", "screened"),
                                   alpha = 0.05, tau_c = 0.50, k_c = 0.60,
                                   seed = 1L) {
  reference <- match.arg(reference)
  h <- floor(k_c * pre_ref$n)
  if (reference == "oracle") {
    candidates <- which(!dat_ref$is_out)
    if (length(candidates) < h) stop("Too few clean oracle-reference units")
    set.seed(seed)
    selected_ref <- sample(candidates, h, replace = FALSE)
  } else {
    selected_ref <- screen_clean_set(pre_ref, tau_c, k_c)$indices
  }
  set.seed(seed + 1L)
  selected_ref <- sample(selected_ref, length(selected_ref), replace = FALSE)
  n_train <- floor(length(selected_ref) / 2L)
  train_idx <- selected_ref[seq_len(n_train)]
  cal_idx <- selected_ref[-seq_len(n_train)]

  score_model <- make_slope_score_model(pre_ref, train_idx)
  calibration <- unlist(
    slope_nonconformity_scores(pre_ref, score_model, cal_idx),
    use.names = FALSE
  )
  test_scores <- slope_nonconformity_scores(pre, score_model)
  eta <- ceiling(pre$T / 2)
  p <- batch_conformal_pvalues(calibration, test_scores, eta)
  selected <- bh_select(p, alpha)
  purity <- mean(!dat_ref$is_out[selected_ref])
  clean_recall <- sum(!dat_ref$is_out[selected_ref]) / sum(!dat_ref$is_out)
  list(selected = selected, p = p, stat = vapply(test_scores, median, numeric(1)),
       eta = eta, reference_size = length(selected_ref),
       training_size = length(train_idx), calibration_size = length(cal_idx),
       reference_purity = purity, reference_clean_recall = clean_recall)
}

vsom_bootstrap_cutoff <- function(pre, R = 1000L, alpha = 0.05, seed = 1L) {
  W <- do.call(rbind, lapply(seq_len(pre$n), function(i) {
    cbind(pre$X[[i]], pre$Z[[i]])
  }))
  Ainv <- safe_solve(crossprod(W))
  leverage <- pmin(1 - 1e-10, pmax(0, rowSums((W %*% Ainv) * W)))
  denom <- pmax(1 - leverage, 1e-10)
  df <- max(1, pre$n * (pre$T - 1L) - ncol(W))
  percentile <- numeric(R)
  set.seed(seed)
  for (b in seq_len(R)) {
    eps <- matrix(rnorm(pre$n * pre$T), nrow = pre$T, ncol = pre$n)
    eps <- sweep(eps, 2L, colMeans(eps), FUN = "-")
    ew <- as.vector(eps)
    resid <- ew - as.vector(W %*% (Ainv %*% crossprod(W, ew)))
    sigma2 <- sum(resid^2) / df
    q <- resid^2 / (max(sigma2, 1e-12) * denom)
    percentile[b] <- as.numeric(quantile(q, 1 - alpha, names = FALSE, type = 8))
  }
  median(percentile)
}

vsom_method <- function(pre, cutoff) {
  fit <- joint_fit(pre, seq_len(pre$n))
  q_list <- vector("list", pre$n)
  for (i in seq_len(pre$n)) {
    W <- cbind(pre$X[[i]], pre$Z[[i]])
    resid <- pre$y[[i]] - as.vector(W %*% fit$coef)
    leverage <- pmin(1 - 1e-10,
                     pmax(0, rowSums((W %*% fit$Ainv) * W)))
    q_list[[i]] <- resid^2 /
      (max(fit$sigma2, 1e-12) * pmax(1 - leverage, 1e-10))
  }
  ## The original VSOM rule is observation-specific. A panel unit is flagged
  ## whenever at least one of its T squared standardized residuals exceeds the
  ## bootstrap cutoff; no BH adjustment is applied.
  stat <- vapply(q_list, max, numeric(1))
  list(selected = which(stat > cutoff), stat = stat, cutoff = cutoff,
       observation_flags = lapply(q_list, function(q) which(q > cutoff)))
}

influence_methods <- function(pre) {
  all_idx <- seq_len(pre$n)
  A0 <- sum_list_matrices(pre$Sxx, all_idx)
  b0 <- sum_list_vectors(pre$sxy, all_idx)
  beta0 <- as.vector(safe_solve(A0, b0))
  syy0 <- sum(pre$syy)
  rss0 <- max(0, syy0 - sum(b0 * beta0))
  df0 <- max(1, pre$n * (pre$T - 1) - pre$m - pre$K)
  sigma0 <- rss0 / df0
  cook <- numeric(pre$n)
  dffits_sq <- numeric(pre$n)
  for (i in all_idx) {
    Am <- A0 - pre$Sxx[[i]]
    bm <- b0 - pre$sxy[[i]]
    beta_m <- as.vector(safe_solve(Am, bm))
    delta <- beta0 - beta_m
    rss_m <- max(0, (syy0 - pre$syy[i]) - sum(bm * beta_m))
    df_m <- max(1, (pre$n - 1) * (pre$T - 1) - pre$m - pre$K)
    sigma_m <- rss_m / df_m
    cook[i] <- as.numeric(crossprod(delta, A0 %*% delta)) /
      (pre$K * max(sigma0, 1e-12))
    ## This is the squared, K-normalized DFFITS statistic. The cutoff 4/n is
    ## equivalent to the conventional unsquared cutoff 2*sqrt(K/n).
    dffits_sq[i] <- as.numeric(crossprod(delta, Am %*% delta)) /
      (pre$K * max(sigma_m, 1e-12))
  }
  cutoff <- 4 / pre$n
  list(
    Cook = list(selected = which(cook > cutoff), beta_test = beta0,
                stat = cook, cutoff = cutoff),
    DFFITS = list(selected = which(dffits_sq > cutoff), beta_test = beta0,
                  stat = dffits_sq, cutoff = cutoff)
  )
}

corrected_beta <- function(pre, selected) {
  ref <- setdiff(seq_len(pre$n), selected)
  if (!length(ref)) ref <- seq_len(pre$n)
  pooled_joint_beta(pre, ref)
}

evaluate_method <- function(name, result, pre, dat) {
  selected <- sort(unique(result$selected))
  truth <- dat$is_out
  tp <- sum(truth[selected])
  fp <- sum(!truth[selected])
  n_out <- sum(truth)
  n_clean <- sum(!truth)
  tpr <- if (n_out > 0) tp / n_out else NA_real_
  fpr <- fp / n_clean
  fdp <- fp / max(length(selected), 1L)
  beta_est <- corrected_beta(pre, selected)
  err <- beta_est - dat$beta
  data.frame(
    method = name, selected = length(selected), retained = pre$n - length(selected),
    fallback = as.integer(length(selected) == pre$n),
    TP = tp, FP = fp,
    TPR = tpr, FPR = fpr, FDR = fdp,
    RMSE = sqrt(mean(err^2)), MAE = mean(abs(err)),
    reference_purity = if (!is.null(result$reference_purity))
      result$reference_purity else NA_real_,
    reference_clean_recall = if (!is.null(result$reference_clean_recall))
      result$reference_clean_recall else NA_real_,
    stringsAsFactors = FALSE
  )
}

run_replication <- function(task) {
  dat <- generate_panel(task$scenario, task$zeta, task$seed)
  pre <- prepare_panel(dat)
  reference_seed <- if (!is.null(task$reference_seed)) task$reference_seed else
    task$seed + 600000000L
  dat_ref <- generate_panel(task$scenario, task$zeta, reference_seed)
  pre_ref <- prepare_panel(dat_ref)
  all_idx <- seq_len(dat$n)
  clean_idx <- which(!dat$is_out)

  prop <- proposed_method(pre)
  oracle <- lm_method(pre, clean_idx)
  pooled <- lm_method(pre, all_idx)
  bc_oracle <- batch_conformal_method(
    pre, pre_ref, dat_ref, reference = "oracle", seed = reference_seed + 1000L
  )
  bc_screened <- batch_conformal_method(
    pre, pre_ref, dat_ref, reference = "screened", seed = reference_seed + 2000L
  )
  cutoff <- if (!is.null(task$vsom_cutoff)) task$vsom_cutoff else
    qchisq(0.95, df = 1)
  vsom <- vsom_method(pre, cutoff)
  infl <- influence_methods(pre)
  rows <- rbind(
    evaluate_method("Proposed", prop, pre, dat),
    evaluate_method("Oracle LM--BH", oracle, pre, dat),
    evaluate_method("Pooled LM--BH", pooled, pre, dat),
    evaluate_method("Batch Conformal--BH (Oracle)", bc_oracle, pre, dat),
    evaluate_method("Batch Conformal--BH (Screened)", bc_screened, pre, dat),
    evaluate_method("VSOM", vsom, pre, dat),
    evaluate_method("Cook's distance", infl$Cook, pre, dat),
    evaluate_method("DFFITS", infl$DFFITS, pre, dat)
  )
  rows$scenario <- task$scenario$name
  rows$zeta <- task$zeta
  rows$replication <- task$replication
  rows$seed <- task$seed
  rows
}

scenario_catalog <- function() {
  base <- list(n = 200L, T = 100L, K = 4L, m = 4L, rho = 0.10,
               predictor = "iid", error = "normal", direction = "spherical",
               signal_lower = 0.10)
  list(
    iid_spherical = modifyList(base, list(
      name = "IID--spherical", predictor = "iid", error = "normal",
      direction = "spherical", rho = 0.10, T = 100L)),
    ar_spherical = modifyList(base, list(
      name = "AR--spherical", predictor = "ar", error = "mixnorm",
      direction = "spherical", rho = 0.10, T = 100L)),
    ar_aligned = modifyList(base, list(
      name = "AR--aligned", predictor = "ar", error = "mixnorm",
      direction = "aligned", rho = 0.10, T = 100L))
  )
}

make_tasks <- function(mode) {
  sc <- scenario_catalog()
  if (mode == "smoke") {
    design <- data.frame(key = "iid_spherical", zeta = 0.10, reps = 2L)
  } else if (mode == "pilot") {
    design <- expand.grid(
      key = names(sc), zeta = c(0, 0.10, 0.30, 0.50),
      stringsAsFactors = FALSE
    )
    design$reps <- 5L
  } else {
    ## This design is declared before inspecting the pilot results. The first
    ## two settings retain the spherical contamination in simulationcode.r;
    ## the third examines masking when heterogeneous slopes share a direction.
    design <- rbind(
      expand.grid(
        key = c("iid_spherical", "ar_spherical"),
        zeta = c(0, 0.10, 0.20, 0.30, 0.50),
        stringsAsFactors = FALSE
      ),
      expand.grid(
        key = "ar_aligned", zeta = c(0.10, 0.20, 0.30, 0.50),
        stringsAsFactors = FALSE
      )
    )
    design$reps <- 100L
  }
  tasks <- list()
  pos <- 1L
  for (j in seq_len(nrow(design))) {
    for (r in seq_len(design$reps[j])) {
      code <- match(design$key[j], names(sc))
      ## Positive signal levels use common random numbers within each scenario.
      ## The global-null seed is separated because no anomaly directions are drawn.
      null_shift <- if (design$zeta[j] == 0) 50000L else 0L
      tasks[[pos]] <- list(
        scenario = sc[[design$key[j]]], zeta = design$zeta[j], replication = r,
        seed = 202608190L + code * 100000L + null_shift + r
      )
      pos <- pos + 1L
    }
  }
  tasks
}

summarize_results <- function(raw) {
  keys <- unique(raw[c("scenario", "zeta", "method")])
  metrics <- c("TPR", "FPR", "FDR", "RMSE", "MAE", "selected")
  out <- vector("list", nrow(keys))
  for (j in seq_len(nrow(keys))) {
    take <- raw$scenario == keys$scenario[j] & raw$zeta == keys$zeta[j] &
      raw$method == keys$method[j]
    d <- raw[take, , drop = FALSE]
    row <- keys[j, , drop = FALSE]
    for (v in metrics) {
      x <- d[[v]]
      row[[paste0(v, "_mean")]] <- if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
      row[[paste0(v, "_sd")]] <- if (sum(!is.na(x)) <= 1L) NA_real_ else sd(x, na.rm = TRUE)
    }
    row$replications <- length(unique(d$replication))
    row$failures <- sum(!is.finite(d$RMSE))
    out[[j]] <- row
  }
  ans <- do.call(rbind, out)
  method_order <- c(
    "Proposed", "Oracle LM--BH", "Pooled LM--BH",
    "Batch Conformal--BH (Oracle)", "Batch Conformal--BH (Screened)",
    "VSOM", "Cook's distance", "DFFITS"
  )
  ans <- ans[order(ans$scenario, ans$zeta, match(ans$method, method_order)), ]
  rownames(ans) <- NULL
  ans
}

fmt_cell <- function(mu, s, digits = 3L) {
  if (!is.finite(mu)) return("--")
  if (!is.finite(s)) return(formatC(mu, format = "f", digits = digits))
  paste0(formatC(mu, format = "f", digits = digits), " (",
         formatC(s, format = "f", digits = digits), ")")
}

latex_escape_method <- function(x) {
  x <- gsub("--", "-", x, fixed = TRUE)
  x
}

write_latex_table <- function(d, metrics, caption, label, path, digits) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  cat("\\begin{table}[!htbp]\n\\centering\n", file = con)
  cat("\\caption{", caption, "}\n", sep = "", file = con)
  cat("\\label{", label, "}\n", sep = "", file = con)
  cat("\\begin{tabular}{ll", paste(rep("c", length(metrics)), collapse = ""), "}\n",
      sep = "", file = con)
  cat("\\toprule\n$\\zeta$ & Method & ", paste(metrics, collapse = " & "), " \\\\\n",
      sep = "", file = con)
  cat("\\midrule\n", file = con)
  last_z <- NA_real_
  for (i in seq_len(nrow(d))) {
    ztxt <- if (!isTRUE(all.equal(d$zeta[i], last_z))) format(d$zeta[i], trim = TRUE) else ""
    vals <- vapply(seq_along(metrics), function(j) {
      v <- metrics[j]
      fmt_cell(d[[paste0(v, "_mean")]][i], d[[paste0(v, "_sd")]][i], digits[j])
    }, character(1))
    cat(ztxt, " & ", latex_escape_method(d$method[i]), " & ",
        paste(vals, collapse = " & "), " \\\\\n", sep = "", file = con)
    last_z <- d$zeta[i]
    next_z <- if (i < nrow(d)) d$zeta[i + 1L] else NA_real_
    if (i < nrow(d) && !isTRUE(all.equal(last_z, next_z))) cat("\\midrule\n", file = con)
  }
  cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n", file = con)
}

tasks <- make_tasks(mode)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
workers <- suppressWarnings(as.integer(Sys.getenv("SIM_WORKERS", "1")))
if (!is.finite(workers) || workers < 1L) workers <- 1L
workers <- min(workers, length(tasks))
message("Mode: ", mode, "; tasks: ", length(tasks), "; workers: ", workers)

if (workers > 1L) {
  cl <- parallel::makeCluster(workers)
  parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(MASS)))
  export_names <- setdiff(ls(envir = .GlobalEnv), c("cl", "tasks", "workers"))
  parallel::clusterExport(cl, export_names, envir = .GlobalEnv)
  pieces <- parallel::parLapply(cl, tasks, run_replication)
  parallel::stopCluster(cl)
} else {
  pieces <- lapply(seq_along(tasks), function(j) {
    if (j %% 25L == 0L || j == 1L || j == length(tasks)) {
      message("Completed task ", j, " / ", length(tasks))
    }
    run_replication(tasks[[j]])
  })
}

raw <- do.call(rbind, pieces)
summary <- summarize_results(raw)
prefix <- paste0("simulation_comparison_", mode)
write.csv(raw, file.path(output_dir, paste0(prefix, "_raw.csv")), row.names = FALSE)
write.csv(summary, file.path(output_dir, paste0(prefix, "_summary.csv")), row.names = FALSE)
saveRDS(list(raw = raw, summary = summary, tasks = tasks),
        file.path(output_dir, paste0(prefix, ".rds")))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

if (mode == "final") {
  main <- summary[summary$scenario == "AR--aligned", ]
  write_latex_table(
    main, c("TPR", "FPR", "FDR"),
    "Detection performance in the aligned design. Entries are Monte Carlo means with standard deviations in parentheses over 100 replications.",
    "tab:comparison-detection", file.path(output_dir, "comparison_detection_table.tex"), c(3L, 3L, 3L)
  )
  write_latex_table(
    main, c("RMSE", "MAE"),
    "Common-slope estimation in the aligned design. Entries are Monte Carlo means with standard deviations in parentheses over 100 replications.",
    "tab:comparison-estimation", file.path(output_dir, "comparison_estimation_table.tex"), c(4L, 4L)
  )
}

message("Saved output under: ", output_dir)
