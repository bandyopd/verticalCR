VerticalCR_check_matrix <- function(x, name) {
  x <- as.matrix(x)
  if (!all(is.finite(x))) stop(name, " must contain only finite numeric values.")
  x
}

VerticalCR_check_vars <- function(data, vars, label) {
  if (is.null(vars) || length(vars) == 0L) stop(label, " must be provided.")
  miss <- setdiff(vars, names(data))
  if (length(miss) > 0L) {
    stop(label, " contains variables not found in data: ", paste(miss, collapse = ", "))
  }
  vars
}

VerticalCR_prepare_inputs <- function(data = NULL,
                                       time = NULL,
                                       status = NULL,
                                       X = NULL,
                                       Z = NULL,
                                       W = NULL,
                                       cluster = NULL,
                                       total_vars = NULL,
                                       relative_vars = NULL,
                                       time_var = "time",
                                       status_var = "status",
                                       threshold_var = "W",
                                       cluster_var = "cluster") {
  if (!is.null(data)) {
    if (!is.data.frame(data)) stop("data must be a data frame.")
    required <- c(time_var, status_var, threshold_var, cluster_var)
    miss <- setdiff(required, names(data))
    if (length(miss) > 0L) {
      stop("data is missing required columns: ", paste(miss, collapse = ", "))
    }
    total_vars <- VerticalCR_check_vars(data, total_vars, "total_vars")
    relative_vars <- VerticalCR_check_vars(data, relative_vars, "relative_vars")
    time <- data[[time_var]]
    status <- data[[status_var]]
    W <- data[[threshold_var]]
    cluster <- data[[cluster_var]]
    X <- as.matrix(data[, total_vars, drop = FALSE])
    Z <- as.matrix(data[, relative_vars, drop = FALSE])
  } else {
    if (is.null(X)) stop("X must be provided when data is NULL.")
    X <- VerticalCR_check_matrix(X, "X")
    if (is.null(Z)) Z <- X
    Z <- VerticalCR_check_matrix(Z, "Z")
    if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(ncol(X)))
    if (is.null(colnames(Z))) colnames(Z) <- paste0("z", seq_len(ncol(Z)))
    total_vars <- colnames(X)
    relative_vars <- colnames(Z)
  }
  list(
    time = time,
    status = status,
    X = X,
    Z = Z,
    W = W,
    cluster = cluster,
    total_vars = total_vars,
    relative_vars = relative_vars
  )
}

make_sum_to_zero_contrast <- function(levels) {
  levels <- as.character(levels)
  m <- length(levels)
  if (m < 2L) stop("cluster must contain at least two distinct centers.")
  cmat <- qr.Q(qr(stats::contr.helmert(m)))
  rownames(cmat) <- levels
  colnames(cmat) <- paste0("v", seq_len(m - 1L))
  cmat
}

make_cluster_indicator <- function(cluster, levels = sort(unique(cluster))) {
  cluster <- as.character(cluster)
  levels <- as.character(levels)
  miss <- setdiff(unique(cluster), levels)
  if (length(miss) > 0L) stop("cluster contains levels not present in levels.")
  out <- matrix(0, nrow = length(cluster), ncol = length(levels))
  out[cbind(seq_along(cluster), match(cluster, levels))] <- 1
  colnames(out) <- paste0("cluster_", levels)
  out
}

make_event_time_basis <- function(time, df = 5L, degree = 3L) {
  time <- as.numeric(time)
  if (!all(is.finite(time))) stop("time must contain only finite values.")
  basis <- splines::bs(time, df = df, degree = degree, intercept = TRUE)
  list(
    basis = as.matrix(basis),
    df = df,
    degree = degree,
    knots = attr(basis, "knots"),
    Boundary.knots = attr(basis, "Boundary.knots")
  )
}

predict_event_time_basis <- function(time, basis_info) {
  splines::bs(as.numeric(time),
              knots = basis_info$knots,
              Boundary.knots = basis_info$Boundary.knots,
              degree = basis_info$degree,
              intercept = TRUE)
}

VerticalCR_step_eval0 <- function(x, y, xout) {
  if (length(x) == 0L) return(rep(0, length(xout)))
  stats::stepfun(x, c(0, y), right = TRUE)(xout)
}

VerticalCR_smooth_threshold_design <- function(A, W, zeta, h) {
  A <- as.matrix(A)
  g <- stats::pnorm((as.numeric(W) - as.numeric(zeta)[1L]) / as.numeric(h)[1L])
  cbind(A, g, A * g)
}

fit_twocause_vertical <- function(time = NULL,
                                  status = NULL,
                                  X = NULL,
                                  Z = NULL,
                                  W = NULL,
                                  cluster = NULL,
                                  data = NULL,
                                  total_vars = NULL,
                                  relative_vars = NULL,
                                  time_var = "time",
                                  status_var = "status",
                                  threshold_var = "W",
                                  cluster_var = "cluster",
                                  spline_df = 5L,
                                  init_beta = NULL,
                                  init_gamma = NULL,
                                  init_kappa = NULL,
                                  init_U = NULL,
                                  init_zetah = 0.5,
                                  init_zetal = 0.5,
                                  init_rho = 0.5,
                                  init_theta = 0.5,
                                  h1 = NA_real_,
                                  h2 = NA_real_,
                                  lr = 1,
                                  mite = 50L,
                                  miter = 50L,
                                  eps = 1e-5,
                                  rho_bound = 1.0,
                                  zeta_lr_pi = 1.0) {
  prepared <- VerticalCR_prepare_inputs(
    data = data, time = time, status = status, X = X, Z = Z, W = W,
    cluster = cluster, total_vars = total_vars, relative_vars = relative_vars,
    time_var = time_var, status_var = status_var,
    threshold_var = threshold_var, cluster_var = cluster_var
  )
  time <- prepared$time
  status <- prepared$status
  X <- prepared$X
  Z <- prepared$Z
  W <- prepared$W
  cluster <- prepared$cluster

  time <- as.numeric(time)
  status <- as.integer(status)
  W <- as.numeric(W)
  cluster <- as.character(cluster)
  n <- length(time)
  if (length(status) != n || length(W) != n || length(cluster) != n) {
    stop("time, status, W, and cluster must have the same length.")
  }
  if (any(!status %in% c(0L, 1L, 2L))) {
    stop("status must be coded as 0=censored, 1=cause 1, 2=reference cause.")
  }
  if (any(time < 0) || !all(is.finite(time))) stop("time must be finite and nonnegative.")

  X <- VerticalCR_check_matrix(X, "X")
  Z <- VerticalCR_check_matrix(Z, "Z")
  if (nrow(X) != n || nrow(Z) != n) stop("X and Z must have length(time) rows.")

  basis_info <- make_event_time_basis(time, df = spline_df)
  Bt <- basis_info$basis
  Delta <- as.numeric(status > 0L)
  D <- as.numeric(status == 1L)

  levels <- sort(unique(cluster))
  Cmat <- make_sum_to_zero_contrast(levels)
  Q_indicator <- make_cluster_indicator(cluster, levels)
  Q <- Q_indicator %*% Cmat

  p <- ncol(X)
  q <- ncol(Z)
  r <- ncol(Bt)
  m <- ncol(Q)
  ibeta <- if (is.null(init_beta)) rep(0.1, 2L * p + 1L) else as.numeric(init_beta)
  igamma <- if (is.null(init_gamma)) rep(0.1, 2L * q + 1L) else as.numeric(init_gamma)
  ikappa <- if (is.null(init_kappa)) rep(0.1, r) else as.numeric(init_kappa)
  if (length(ibeta) != 2L * p + 1L) stop("init_beta has incompatible length.")
  if (length(igamma) != 2L * q + 1L) stop("init_gamma has incompatible length.")
  if (length(ikappa) != r) stop("init_kappa has incompatible length.")

  if (is.null(init_U)) {
    mean_time <- max(mean(time), 1e-8)
    u_all <- vapply(levels, function(cl) {
      idx <- cluster == cl
      log((sum(Delta[idx]) + 0.5) / (sum(time[idx]) + 0.5 * mean_time))
    }, numeric(1))
    u_all <- u_all - mean(u_all)
    iU <- drop(crossprod(Cmat, pmin(pmax(u_all, -1.5), 1.5)))
  } else {
    init_U <- as.numeric(init_U)
    if (length(init_U) == length(levels)) {
      init_U <- init_U - mean(init_U)
      iU <- drop(crossprod(Cmat, init_U))
    } else if (length(init_U) == m) {
      iU <- init_U
    } else {
      stop("init_U must have length N clusters or N - 1 contrast coefficients.")
    }
  }

  if (!is.finite(h1)) h1 <- log(n) / sqrt(n) * stats::sd(W)
  nevent <- sum(Delta)
  if (nevent <= 0L) stop("At least one observed event is required.")
  if (!is.finite(h2)) h2 <- log(nevent) / sqrt(nevent) * stats::sd(W)

  fit <- tw_vertical(
    xtime = time, Delta = Delta, D = D, X = X, Z = Z, Bt = Bt, W = W, Q = Q,
    ibeta = ibeta, igamma = igamma, ikappa = ikappa, iU = iU,
    izetah = c(init_zetah), izetal = c(init_zetal), irho = c(init_rho)[1L],
    h1 = h1, h2 = h2, itheta = init_theta, lr = lr, mite = as.integer(mite),
    miter = as.integer(miter), eps = eps, rho_bound = rho_bound,
    zeta_lr_pi = zeta_lr_pi
  )

  fit$gamma1 <- fit$gamma
  fit$se_gamma1 <- fit$se_gamma
  fit$kappa1 <- fit$kappa
  fit$se_kappa1 <- fit$se_kappa
  fit$rho1 <- as.numeric(fit$rho)[1L]
  v_hat <- as.numeric(fit$v)
  if (length(v_hat) == ncol(Cmat)) {
    fit$v <- v_hat
    fit$U <- drop(Cmat %*% v_hat)
  }
  fit$converged <- isTRUE(fit$converged) &&
    all(is.finite(unlist(fit[c("beta", "gamma", "zetah", "zetal", "rho", "theta")])))
  fit$status_coding <- "0=censored, 1=cause 1, 2=reference cause"
  fit$cluster_levels <- levels
  fit$contrast_matrix <- Cmat
  fit$frailty_constraint <- "sum-to-zero center frailties via U = C v"
  fit$Q <- Q
  fit$Q_indicator <- Q_indicator
  fit$basis_info <- basis_info
  fit$training_data <- list(time = time, status = status, X = X, Z = Z, W = W,
                            cluster = cluster)
  fit$total_vars <- prepared$total_vars
  fit$relative_vars <- prepared$relative_vars
  fit$h1 <- h1
  fit$h2 <- h2
  fit$zeta_lr_pi <- zeta_lr_pi
  class(fit) <- c("twocause_vertical_fit", class(fit))
  fit
}

fit_twocause_vertical_chunked <- function(time = NULL,
                                          status = NULL,
                                          X = NULL,
                                          Z = NULL,
                                          W = NULL,
                                          cluster = NULL,
                                          data = NULL,
                                          total_vars = NULL,
                                          relative_vars = NULL,
                                          time_var = "time",
                                          status_var = "status",
                                          threshold_var = "W",
                                          cluster_var = "cluster",
                                          chunk_size = 2L,
                                          max_chunks = 50L,
                                          checkpoint_dir = NULL,
                                          ...) {
  args <- list(...)
  args$miter <- as.integer(chunk_size)
  fit <- NULL
  for (chunk in seq_len(max_chunks)) {
    if (!is.null(fit)) {
      args$init_beta <- fit$beta
      args$init_gamma <- fit$gamma
      args$init_kappa <- fit$kappa
      args$init_U <- fit$v
      args$init_zetah <- fit$zetah
      args$init_zetal <- fit$zetal
      args$init_rho <- fit$rho
      args$init_theta <- fit$theta
    }
    fit <- do.call(fit_twocause_vertical,
                   c(list(data = data, time = time, status = status,
                          X = X, Z = Z, W = W, cluster = cluster,
                          total_vars = total_vars,
                          relative_vars = relative_vars,
                          time_var = time_var,
                          status_var = status_var,
                          threshold_var = threshold_var,
                          cluster_var = cluster_var), args))
    fit$chunk <- chunk
    if (!is.null(checkpoint_dir)) {
      if (!dir.exists(checkpoint_dir)) dir.create(checkpoint_dir, recursive = TRUE)
      saveRDS(fit, file.path(checkpoint_dir, sprintf("twocause_chunk_%03d.rds", chunk)))
    }
    if (isTRUE(fit$converged)) break
  }
  fit
}

extract_twocause_parameters <- function(fit, covariate_names = NULL,
                                        total_names = NULL,
                                        relative_names = NULL) {
  if (!inherits(fit, "twocause_vertical_fit")) stop("fit must be from fit_twocause_vertical().")
  p <- (length(fit$beta) - 1L) / 2L
  q <- (length(fit$gamma) - 1L) / 2L
  if (!is.null(covariate_names)) {
    total_names <- covariate_names
    relative_names <- covariate_names
  }
  if (is.null(total_names)) total_names <- fit$total_vars
  if (is.null(relative_names)) relative_names <- fit$relative_vars
  if (is.null(total_names)) total_names <- paste0("x", seq_len(p))
  if (is.null(relative_names)) relative_names <- paste0("z", seq_len(q))
  beta_names <- c(paste0("beta1_", total_names[seq_len(p)]),
                  "beta2",
                  paste0("beta3_", total_names[seq_len(p)]))
  gamma_names <- c(paste0("gamma1_", relative_names[seq_len(q)]),
                   "gamma2",
                   paste0("gamma3_", relative_names[seq_len(q)]))
  par <- c(fit$beta, fit$zetah, fit$gamma, fit$kappa, fit$zetal, fit$rho, fit$theta)
  se <- c(fit$se_beta, fit$se_zetah, fit$se_gamma, fit$se_kappa,
          fit$se_zetal, fit$se_rho1, fit$se_theta)
  nm <- c(beta_names, "zeta_lambda", gamma_names,
          paste0("kappa", seq_along(fit$kappa)),
          "zeta_pi", "rho", "theta")
  z <- par / se
  data.frame(
    parameter = nm,
    estimate = as.numeric(par),
    se = as.numeric(se),
    z = as.numeric(z),
    p_value = 2 * stats::pnorm(-abs(z)),
    row.names = NULL
  )
}

predict_twocause_vertical <- function(fit, newdata = NULL, times,
                                      use_EB_frailty = TRUE) {
  if (!inherits(fit, "twocause_vertical_fit")) stop("fit must be from fit_twocause_vertical().")
  if (is.null(newdata)) {
    dat <- fit$training_data
  } else {
    dat <- newdata
  }
  X <- VerticalCR_check_matrix(dat$X, "newdata$X")
  Z <- if (is.null(dat$Z)) X else VerticalCR_check_matrix(dat$Z, "newdata$Z")
  W <- as.numeric(dat$W)
  cluster <- as.character(dat$cluster)
  times <- sort(unique(as.numeric(times)))
  train <- fit$training_data

  Xt_train <- VerticalCR_smooth_threshold_design(train$X, train$W, fit$zetah, fit$h1)
  u_train <- setNames(as.numeric(fit$U), fit$cluster_levels)[train$cluster]
  u_train[is.na(u_train)] <- 0
  eta_train <- drop(Xt_train %*% fit$beta) + u_train

  event_times <- sort(unique(train$time[train$status > 0L]))
  dH0 <- numeric(length(event_times))
  for (j in seq_along(event_times)) {
    tj <- event_times[j]
    d <- sum(train$time == tj & train$status > 0L)
    risk <- train$time >= tj
    dH0[j] <- d / sum(exp(pmin(pmax(eta_train[risk], -35), 35)))
  }
  keep <- event_times <= max(times)
  event_times <- event_times[keep]
  dH0 <- dH0[keep]

  Xt_new <- VerticalCR_smooth_threshold_design(X, W, fit$zetah, fit$h1)
  Zt_new <- VerticalCR_smooth_threshold_design(Z, W, fit$zetal, fit$h2)
  u_new <- rep(0, nrow(X))
  if (isTRUE(use_EB_frailty)) {
    u_by_cluster <- setNames(as.numeric(fit$U), fit$cluster_levels)
    u_new <- u_by_cluster[cluster]
    u_new[is.na(u_new)] <- 0
  }

  eta <- drop(Xt_new %*% fit$beta) + u_new
  Bt_grid <- as.matrix(predict_event_time_basis(event_times, fit$basis_info))
  F1 <- F2 <- matrix(0, nrow = nrow(X), ncol = length(times))
  for (i in seq_len(nrow(X))) {
    xi1 <- drop(Zt_new[i, , drop = FALSE] %*% fit$gamma1) +
      drop(Bt_grid %*% fit$kappa1) + fit$rho1 * u_new[i]
    pi1 <- 1 / (1 + exp(-pmin(pmax(xi1, -35), 35)))
    dHtot <- exp(pmin(pmax(eta[i], -35), 35)) * dH0
    Sprev <- exp(-c(0, cumsum(dHtot)[-length(dHtot)]))
    Fstep1 <- cumsum(Sprev * pi1 * dHtot)
    Fstep2 <- cumsum(Sprev * (1 - pi1) * dHtot)
    F1[i, ] <- VerticalCR_step_eval0(event_times, Fstep1, times)
    F2[i, ] <- VerticalCR_step_eval0(event_times, Fstep2, times)
  }
  list(times = times, F1 = F1, F2 = F2)
}

loglik_twocause_vertical <- function(fit) {
  if (!inherits(fit, "twocause_vertical_fit")) stop("fit must be from fit_twocause_vertical().")
  dat <- fit$training_data
  X <- dat$X
  Z <- dat$Z
  time <- dat$time
  status <- dat$status
  Delta <- as.numeric(status > 0L)
  D <- as.numeric(status == 1L)
  u_by_cluster <- setNames(as.numeric(fit$U), fit$cluster_levels)
  u <- u_by_cluster[dat$cluster]
  u[is.na(u)] <- 0

  Xt <- VerticalCR_smooth_threshold_design(X, dat$W, fit$zetah, fit$h1)
  eta <- drop(Xt %*% fit$beta) + u
  eta <- pmin(pmax(eta, -35), 35)

  cox_pl <- 0
  event_times <- sort(unique(time[Delta > 0]))
  for (tj in event_times) {
    fail <- time == tj & Delta > 0
    risk <- time >= tj
    d <- sum(fail)
    cox_pl <- cox_pl + sum(eta[fail]) - d * log(sum(exp(eta[risk])))
  }

  Zt <- VerticalCR_smooth_threshold_design(Z, dat$W, fit$zetal, fit$h2)
  Bt <- fit$basis_info$basis
  xi <- drop(Zt %*% fit$gamma1) + drop(Bt %*% fit$kappa1) + fit$rho1 * u
  xi <- pmin(pmax(xi, -35), 35)
  rel_ll <- sum(Delta * (D * xi - log1p(exp(xi))))

  theta <- max(as.numeric(fit$theta)[1L], 1e-8)
  frailty_ll <- -0.5 * (length(fit$U) * log(2 * pi * theta) + sum(fit$U^2) / theta)
  out <- cox_pl + rel_ll + frailty_ll
  attr(out, "components") <- c(cox_partial = cox_pl,
                               relative = rel_ll,
                               frailty = frailty_ll)
  out
}

select_spline_df_twocause <- function(time = NULL,
                                      status = NULL,
                                      X = NULL,
                                      Z = NULL,
                                      W = NULL,
                                      cluster = NULL,
                                      data = NULL,
                                      total_vars = NULL,
                                      relative_vars = NULL,
                                      time_var = "time",
                                      status_var = "status",
                                      threshold_var = "W",
                                      cluster_var = "cluster",
                                      candidate_df = 3:8,
                                      criterion = c("BIC", "AIC"),
                                      keep_fits = FALSE,
                                      ...) {
  criterion <- match.arg(criterion)
  candidate_df <- sort(unique(as.integer(candidate_df)))
  if (any(candidate_df < 2L)) stop("candidate_df must contain integers >= 2.")
  n <- if (!is.null(data)) nrow(data) else length(time)
  if (!is.finite(n) || n <= 0L) stop("Cannot determine sample size.")

  rows <- vector("list", length(candidate_df))
  fits <- if (keep_fits) vector("list", length(candidate_df)) else NULL
  for (i in seq_along(candidate_df)) {
    df <- candidate_df[i]
    fit <- tryCatch(
      fit_twocause_vertical(
        time = time, status = status, X = X, Z = Z, W = W, cluster = cluster,
        data = data, total_vars = total_vars, relative_vars = relative_vars,
        time_var = time_var, status_var = status_var,
        threshold_var = threshold_var, cluster_var = cluster_var,
        spline_df = df, ...
      ),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      rows[[i]] <- data.frame(
        spline_df = df,
        converged = FALSE,
        loglik = NA_real_,
        n_parameters = NA_integer_,
        AIC = NA_real_,
        BIC = NA_real_,
        error = conditionMessage(fit)
      )
      next
    }
    ll <- loglik_twocause_vertical(fit)
    k <- length(fit$beta) + length(fit$gamma) + length(fit$kappa) +
      length(fit$v) + 4L
    rows[[i]] <- data.frame(
      spline_df = df,
      converged = isTRUE(fit$converged),
      loglik = as.numeric(ll),
      n_parameters = k,
      AIC = -2 * as.numeric(ll) + 2 * k,
      BIC = -2 * as.numeric(ll) + log(n) * k,
      error = NA_character_
    )
    if (keep_fits) fits[[i]] <- fit
  }

  table <- do.call(rbind, rows)
  score <- table[[criterion]]
  ok <- is.finite(score)
  if (!any(ok)) stop("All spline_df candidates failed.")
  best_row <- which(ok)[which.min(score[ok])]
  out <- list(
    criterion = criterion,
    selected_df = table$spline_df[best_row],
    summary = table,
    best_fit = if (keep_fits) fits[[best_row]] else NULL,
    fits = fits
  )
  class(out) <- "twocause_spline_df_selection"
  out
}
