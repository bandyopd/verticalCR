manuscript_or_default <- function(x, y) {
  if (is.null(x)) y else x
}

manuscript_simulate_twocause_data <- function(
    N = 40L,
    ni = c(rep(40L, 20L), rep(20L, 20L)),
    theta = 0.5,
    zeta_lambda = 0.35,
    zeta_pi = 0.35,
    censoring_rate_target = 0.15,
    weibull_a = 0.08,
    weibull_b = 1.20,
    beta1 = c(0.40, -0.30),
    beta2 = 0.50,
    beta3 = c(0.50, -0.50),
    gamma1 = c(0.30, -0.20),
    gamma2 = 0.50,
    gamma3 = c(0.50, -0.50),
    rho = 0.50,
    tau0 = 6,
    U_cluster = NULL,
    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (length(ni) == 1L) ni <- rep(as.integer(ni), N)
  if (length(ni) != N) stop("ni must have length 1 or N.")
  n <- sum(ni)
  cluster <- rep(seq_len(N), times = ni)
  if (is.null(U_cluster)) U_cluster <- stats::rnorm(N, 0, sqrt(theta))
  U_true <- U_cluster[cluster]

  X <- matrix(stats::rnorm(n * 2L), ncol = 2L)
  colnames(X) <- c("X1", "X2")
  W <- stats::runif(n)
  I_lambda <- as.numeric(W > zeta_lambda)
  eta <- drop(X %*% beta1) + (beta2 + drop(X %*% beta3)) * I_lambda + U_true
  Tstar <- (-log(pmax(stats::runif(n), .Machine$double.eps)) /
              (weibull_a * exp(pmin(pmax(eta, -35), 35))))^(1 / weibull_b)

  I_pi <- as.numeric(W > zeta_pi)
  f1 <- 0.3 * sin(pi * Tstar / tau0)
  xi <- f1 + drop(X %*% gamma1) +
    (gamma2 + drop(X %*% gamma3)) * I_pi + rho * U_true
  p1 <- 1 / (1 + exp(-pmin(pmax(xi, -35), 35)))
  D <- ifelse(stats::runif(n) <= p1, 1L, 2L)

  censor_u <- stats::runif(n)
  censor_rate_at <- function(cmax) mean(censor_u * cmax < Tstar)
  lo <- max(min(Tstar[Tstar > 0], na.rm = TRUE) * 1e-6, 1e-8)
  hi <- max(max(Tstar, na.rm = TRUE) / min(1 - censoring_rate_target, 0.95), 1)
  while (censor_rate_at(hi) > censoring_rate_target && hi < 1e8) hi <- hi * 2
  for (iter in seq_len(60L)) {
    mid <- (lo + hi) / 2
    if (censor_rate_at(mid) > censoring_rate_target) lo <- mid else hi <- mid
  }
  cmax <- (lo + hi) / 2
  C <- censor_u * cmax
  event <- as.integer(Tstar <= C)
  time <- pmin(Tstar, C)
  cause <- ifelse(event == 1L, D, 0L)

  data <- data.frame(
    id = seq_len(n),
    cluster = cluster,
    time = time,
    status = event,
    cause = cause,
    X1 = X[, 1],
    X2 = X[, 2],
    W = W,
    U_true = U_true,
    Tstar = Tstar,
    C = C
  )
  true_parameters <- list(
    N = N, ni = ni, n = n, theta = theta,
    zeta_lambda = zeta_lambda, zeta_pi = zeta_pi,
    a = weibull_a, b = weibull_b,
    beta1 = beta1, beta2 = beta2, beta3 = beta3,
    gamma1 = gamma1, gamma2 = gamma2, gamma3 = gamma3,
    rho = rho, U_cluster = U_cluster, tau0 = tau0,
    censoring_rate_target = censoring_rate_target,
    censoring_rate_observed = mean(event == 0)
  )
  list(data = data, true_parameters = true_parameters)
}

manuscript_true_cif_twocause <- function(newdata, times, true_par) {
  times <- sort(unique(as.numeric(times)))
  grid <- seq(0, max(times), length.out = max(800L, 200L * length(times)))
  dt <- c(0, diff(grid))
  out1 <- out2 <- matrix(0, nrow(newdata), length(times))
  for (i in seq_len(nrow(newdata))) {
    x <- c(newdata$X1[i], newdata$X2[i])
    w <- newdata$W[i]
    u <- newdata$U_true[i]
    il <- as.numeric(w > true_par$zeta_lambda)
    eta <- sum(x * true_par$beta1) +
      (true_par$beta2 + sum(x * true_par$beta3)) * il + u
    exp_eta <- exp(pmin(pmax(eta, -35), 35))
    lambda <- true_par$a * true_par$b * grid^(true_par$b - 1) * exp_eta
    surv <- exp(-true_par$a * grid^true_par$b * exp_eta)
    ip <- as.numeric(w > true_par$zeta_pi)
    xi <- 0.3 * sin(pi * grid / true_par$tau0) +
      sum(x * true_par$gamma1) +
      (true_par$gamma2 + sum(x * true_par$gamma3)) * ip +
      true_par$rho * u
    pi1 <- 1 / (1 + exp(-pmin(pmax(xi, -35), 35)))
    f1 <- cumsum(surv * lambda * pi1 * dt)
    f2 <- cumsum(surv * lambda * (1 - pi1) * dt)
    out1[i, ] <- stats::approx(grid, f1, xout = times, rule = 2)$y
    out2[i, ] <- stats::approx(grid, f2, xout = times, rule = 2)$y
  }
  list(F1 = out1, F2 = out2)
}

manuscript_array_to_flist <- function(x) {
  if (is.list(x)) return(x)
  if (length(dim(x)) == 3L) return(list(F1 = x[, , 1], F2 = x[, , 2]))
  stop("CIF object must be a list or a 3D array.")
}

manuscript_metric_row <- function(test_data, pred, true_cif, times, t0 = 3) {
  pred <- manuscript_array_to_flist(pred)
  true_cif <- manuscript_array_to_flist(true_cif)
  j0 <- which.min(abs(times - t0))
  imse1 <- mean((pred$F1 - true_cif$F1)^2, na.rm = TRUE)
  imse2 <- mean((pred$F2 - true_cif$F2)^2, na.rm = TRUE)
  y1 <- as.numeric(test_data$time <= t0 & test_data$cause == 1L)
  y2 <- as.numeric(test_data$time <= t0 & test_data$cause == 2L)
  bs1 <- mean((y1 - pred$F1[, j0])^2, na.rm = TRUE)
  bs2 <- mean((y2 - pred$F2[, j0])^2, na.rm = TRUE)
  auc_one <- function(y, risk) {
    case <- y == 1
    control <- test_data$time > t0
    if (!any(case) || !any(control)) return(NA_real_)
    mean(outer(risk[case], risk[control], ">")) +
      0.5 * mean(outer(risk[case], risk[control], "=="))
  }
  c(BS1_t0 = bs1, BS2_t0 = bs2,
    IBS1 = mean(colMeans((outer(y1, rep(1, length(times))) - pred$F1)^2), na.rm = TRUE),
    IBS2 = mean(colMeans((outer(y2, rep(1, length(times))) - pred$F2)^2), na.rm = TRUE),
    AUC1_t0 = auc_one(y1, pred$F1[, j0]),
    AUC2_t0 = auc_one(y2, pred$F2[, j0]),
    IMSE_1 = imse1, IMSE_2 = imse2, IMSE_all = mean(c(imse1, imse2)))
}

manuscript_hard_design <- function(dat, threshold = FALSE, zeta = NULL) {
  x <- as.matrix(dat[, c("X1", "X2")])
  if (!threshold) return(as.data.frame(x))
  g <- as.numeric(dat$W > zeta)
  out <- cbind(x, Iw = g, x * g)
  colnames(out) <- c("X1", "X2", "Iw", "X1_Iw", "X2_Iw")
  as.data.frame(out)
}

manuscript_fit_reduced <- function(train, model, zeta_grid, spline_df = 5L) {
  train$status_event <- as.integer(train$cause > 0L)
  frailty <- model == "V1"
  threshold <- model == "V2"
  candidate_zeta <- if (threshold) zeta_grid else NA_real_
  best <- NULL
  for (zeta in candidate_zeta) {
    z <- if (is.na(zeta)) NULL else zeta
    xd <- manuscript_hard_design(train, threshold, z)
    df_total <- data.frame(time = train$time, event = train$status_event,
                           cluster = factor(train$cluster), xd)
    rhs <- paste(colnames(xd), collapse = " + ")
    if (frailty) rhs <- paste(rhs, "+ frailty(cluster, distribution='gaussian')")
    cox <- survival::coxph(stats::as.formula(paste("survival::Surv(time,event)~", rhs)),
                           data = df_total, x = TRUE, model = TRUE)
    ev <- train[train$status_event == 1L, , drop = FALSE]
    xd_ev <- manuscript_hard_design(ev, threshold, z)
    bt_raw <- splines::bs(ev$time, df = spline_df, intercept = FALSE)
    basis_info <- list(knots = attr(bt_raw, "knots"),
                       Boundary.knots = attr(bt_raw, "Boundary.knots"),
                       degree = attr(bt_raw, "degree"),
                       intercept = attr(bt_raw, "intercept"))
    bt <- as.matrix(bt_raw)
    colnames(bt) <- paste0("bt", seq_len(ncol(bt)))
    df_rel <- data.frame(y = as.integer(ev$cause == 1L), xd_ev, bt)
    glm <- stats::glm(stats::as.formula(paste("y~", paste(setdiff(names(df_rel), "y"), collapse = "+"))),
                      data = df_rel, family = stats::binomial(),
                      control = stats::glm.control(maxit = 100, epsilon = 1e-5))
    obj <- as.numeric(cox$loglik[length(cox$loglik)]) + as.numeric(stats::logLik(glm))
    one <- list(cox = cox, glm = glm, threshold = threshold, frailty = frailty,
                zeta = z, basis_info = basis_info,
                objective = obj)
    if (is.null(best) || one$objective > best$objective) best <- one
  }
  best$model_class <- "reduced"
  best
}

manuscript_predict_reduced <- function(fit, newdata, times) {
  bh <- survival::basehaz(fit$cox, centered = FALSE)
  event_times <- bh$time[bh$time <= max(times)]
  h0 <- stats::stepfun(bh$time, c(0, bh$hazard), right = TRUE)(event_times)
  dH0 <- c(h0[1], diff(h0))
  xd <- manuscript_hard_design(newdata, fit$threshold, fit$zeta)
  beta <- stats::coef(fit$cox)
  beta <- beta[intersect(names(beta), colnames(xd))]
  eta <- if (length(beta)) drop(as.matrix(xd[, names(beta), drop = FALSE]) %*% beta) else rep(0, nrow(newdata))
  bt <- splines::bs(event_times, knots = fit$basis_info$knots,
                    Boundary.knots = fit$basis_info$Boundary.knots,
                    degree = fit$basis_info$degree, intercept = fit$basis_info$intercept)
  f1 <- f2 <- matrix(0, nrow(newdata), length(times))
  cum1 <- cum2 <- matrix(0, nrow(newdata), length(event_times))
  for (j in seq_along(event_times)) {
    nd <- cbind(xd, as.data.frame(matrix(rep(bt[j, ], each = nrow(newdata)),
                                         nrow = nrow(newdata))))
    colnames(nd) <- c(colnames(xd), paste0("bt", seq_len(ncol(bt))))
    p1 <- stats::predict(fit$glm, newdata = nd, type = "response")
    dH <- exp(pmin(pmax(eta, -35), 35)) * dH0[j]
    if (j == 1L) {
      cs1 <- p1 * dH
      cs2 <- (1 - p1) * dH
    } else {
      prevH <- rowSums(cbind(cs1, cs2))
      surv <- exp(-prevH)
      cs1 <- cs1 + surv * p1 * dH
      cs2 <- cs2 + surv * (1 - p1) * dH
    }
    cum1[, j] <- cs1
    cum2[, j] <- cs2
  }
  for (k in seq_along(times)) {
    idx <- which(event_times <= times[k])
    if (length(idx)) {
      jj <- idx[length(idx)]
      f1[, k] <- cum1[, jj]
      f2[, k] <- cum2[, jj]
    }
  }
  list(F1 = pmin(pmax(f1, 0), 1), F2 = pmin(pmax(f2, 0), 1))
}

manuscript_extract_v3_parameters <- function(fit, true_par, rep_id) {
  beta <- as.numeric(fit$beta)
  gamma <- as.numeric(fit$gamma1)
  est <- c(beta1_X1 = beta[1], beta1_X2 = beta[2], beta2 = beta[3],
           beta3_X1 = beta[4], beta3_X2 = beta[5],
           zeta_lambda = as.numeric(fit$zetah)[1],
           gamma1_X1 = gamma[1], gamma1_X2 = gamma[2], gamma2 = gamma[3],
           gamma3_X1 = gamma[4], gamma3_X2 = gamma[5],
           zeta_pi = as.numeric(fit$zetal)[1],
           rho = as.numeric(fit$rho)[1], theta = as.numeric(fit$theta)[1])
  tru <- c(beta1_X1 = true_par$beta1[1], beta1_X2 = true_par$beta1[2],
           beta2 = true_par$beta2, beta3_X1 = true_par$beta3[1],
           beta3_X2 = true_par$beta3[2], zeta_lambda = true_par$zeta_lambda,
           gamma1_X1 = true_par$gamma1[1], gamma1_X2 = true_par$gamma1[2],
           gamma2 = true_par$gamma2, gamma3_X1 = true_par$gamma3[1],
           gamma3_X2 = true_par$gamma3[2], zeta_pi = true_par$zeta_pi,
           rho = true_par$rho, theta = true_par$theta)
  se <- c(beta1_X1 = fit$se_beta[1], beta1_X2 = fit$se_beta[2],
          beta2 = fit$se_beta[3], beta3_X1 = fit$se_beta[4],
          beta3_X2 = fit$se_beta[5], zeta_lambda = fit$se_zetah[1],
          gamma1_X1 = fit$se_gamma1[1], gamma1_X2 = fit$se_gamma1[2],
          gamma2 = fit$se_gamma1[3], gamma3_X1 = fit$se_gamma1[4],
          gamma3_X2 = fit$se_gamma1[5], zeta_pi = fit$se_zetal[1],
          rho = fit$se_rho1[1], theta = fit$se_theta[1])
  data.frame(rep = rep_id, Parameter = names(est), True = tru[names(est)],
             Estimate = est, SE = se[names(est)], row.names = NULL)
}

manuscript_summarize_table1 <- function(raw) {
  do.call(rbind, lapply(unique(raw$Parameter), function(p) {
    d <- raw[raw$Parameter == p, , drop = FALSE]
    err <- d$Estimate - d$True
    data.frame(Parameter = p, True = d$True[1],
               Bias = mean(err, na.rm = TRUE),
               `Empirical SD` = stats::sd(d$Estimate, na.rm = TRUE),
               `Mean SE` = mean(d$SE, na.rm = TRUE),
               RMSE = sqrt(mean(err^2, na.rm = TRUE)),
               Coverage = mean(abs(err) <= 1.96 * d$SE, na.rm = TRUE),
               check.names = FALSE)
  }))
}

manuscript_summarize_table2 <- function(raw) {
  models <- c("V0_basic_vertical", "V1_vertical_frailty",
              "V2_vertical_threshold", "V3_proposed")
  do.call(rbind, lapply(models, function(m) {
    d <- raw[raw$model == m, , drop = FALSE]
    data.frame(Model = m,
               Frailty = ifelse(m %in% c("V1_vertical_frailty", "V3_proposed"), "Yes", "No"),
               Threshold = ifelse(m %in% c("V2_vertical_threshold", "V3_proposed"), "Yes", "No"),
               BS1_t0_mean = mean(d$BS1_t0, na.rm = TRUE),
               BS2_t0_mean = mean(d$BS2_t0, na.rm = TRUE),
               IBS1_mean = mean(d$IBS1, na.rm = TRUE),
               IBS2_mean = mean(d$IBS2, na.rm = TRUE),
               AUC1_t0_mean = mean(d$AUC1_t0, na.rm = TRUE),
               AUC2_t0_mean = mean(d$AUC2_t0, na.rm = TRUE),
               IMSE1_mean = mean(d$IMSE_1, na.rm = TRUE),
               IMSE2_mean = mean(d$IMSE_2, na.rm = TRUE),
               IMSE_all_mean = mean(d$IMSE_all, na.rm = TRUE),
               row.names = NULL)
  }))
}

manuscript_run_table_setting <- function(R = 500L,
                                         N = 40L,
                                         ni = c(rep(40L, 20L), rep(20L, 20L)),
                                         censoring = 0.15,
                                         theta = 0.5,
                                         zeta = 0.35,
                                         workers = 1L,
                                         seed_train_base = 10000L,
                                         seed_test_base = 20000L,
                                         fit_control = list()) {
  times <- seq(0.25, 6, by = 0.25)
  zeta_grid <- seq(0.20, 0.80, by = 0.02)
  fit_control <- utils::modifyList(list(
    spline_df = 5L, mite = 100L, miter = 100L, eps = 1e-5,
    rho_bound = 1.0, zeta_lr_pi = 0.5
  ), fit_control)
  one_rep <- function(r) {
    train_obj <- manuscript_simulate_twocause_data(
      N = N, ni = ni, theta = theta, zeta_lambda = zeta, zeta_pi = zeta,
      censoring_rate_target = censoring, seed = seed_train_base + r)
    test_obj <- manuscript_simulate_twocause_data(
      N = N, ni = ni, theta = theta, zeta_lambda = zeta, zeta_pi = zeta,
      censoring_rate_target = censoring,
      U_cluster = train_obj$true_parameters$U_cluster,
      seed = seed_test_base + r)
    train <- train_obj$data
    test <- test_obj$data
    true_cif <- manuscript_true_cif_twocause(test, times, test_obj$true_parameters)
    xtrain <- as.matrix(train[, c("X1", "X2")])
    fit_v3 <- fit_twocause_vertical(
      time = train$time,
      status = ifelse(train$cause == 0L, 0L, ifelse(train$cause == 1L, 1L, 2L)),
      X = xtrain, Z = xtrain, W = train$W, cluster = train$cluster,
      spline_df = fit_control$spline_df,
      h1 = manuscript_or_default(fit_control$h1, NA_real_),
      h2 = manuscript_or_default(fit_control$h2, NA_real_),
      mite = fit_control$mite,
      miter = fit_control$miter, eps = fit_control$eps,
      rho_bound = fit_control$rho_bound, zeta_lr_pi = fit_control$zeta_lr_pi)
    params <- manuscript_extract_v3_parameters(fit_v3, train_obj$true_parameters, r)
    models <- list(
      V0_basic_vertical = manuscript_fit_reduced(train, "V0", zeta_grid, fit_control$spline_df),
      V1_vertical_frailty = manuscript_fit_reduced(train, "V1", zeta_grid, fit_control$spline_df),
      V2_vertical_threshold = manuscript_fit_reduced(train, "V2", zeta_grid, fit_control$spline_df),
      V3_proposed = fit_v3
    )
    preds <- lapply(names(models), function(m) {
      if (m == "V3_proposed") {
        pred <- predict_twocause_vertical(
          models[[m]],
          newdata = list(X = as.matrix(test[, c("X1", "X2")]),
                         Z = as.matrix(test[, c("X1", "X2")]),
                         W = test$W, cluster = test$cluster),
          times = times, use_EB_frailty = TRUE)
      } else {
        pred <- manuscript_predict_reduced(models[[m]], test, times)
      }
      data.frame(rep = r, model = m,
                 t(manuscript_metric_row(test, pred, true_cif, times, 3)),
                 row.names = NULL)
    })
    list(parameters = params, prediction = do.call(rbind, preds))
  }
  reps <- if (workers > 1L) {
    cl <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    root <- Sys.getenv("VerticalCR_PROJECT_ROOT", unset = "")
    parallel::clusterExport(cl, varlist = c("root"), envir = environment())
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages({
        library(Rcpp)
        library(RcppArmadillo)
        library(survival)
        library(splines)
      })
      if (nzchar(root)) {
        dll <- file.path(root, "VerticalCR", "src", "VerticalCR.dll")
        if (file.exists(dll) && !"VerticalCR" %in% names(getLoadedDLLs())) {
          dyn.load(dll)
        }
        source(file.path(root, "VerticalCR", "R", "RcppExports.R"),
               local = globalenv())
        source(file.path(root, "VerticalCR", "R", "model_estimation.R"),
               local = globalenv())
        source(file.path(root, "VerticalCR", "R", "manuscript_reproduction.R"),
               local = globalenv())
        assign("tw_vertical",
               function(xtime, Delta, D, X, Z, Bt, W, Q, ibeta, igamma,
                        ikappa, iU, izetah, izetal, irho, h1, h2, itheta,
                        lr = 1, mite = 20L, miter = 20L, eps = 1E-3,
                        rho_bound = 1.0, zeta_lr_pi = 1.0) {
                 .Call("_VerticalCR_tw_vertical", xtime, Delta, D, X, Z,
                       Bt, W, Q, ibeta, igamma, ikappa, iU, izetah,
                       izetal, irho, h1, h2, itheta, lr, mite, miter,
                       eps, rho_bound, zeta_lr_pi, PACKAGE = "VerticalCR")
               },
               envir = globalenv())
      } else {
        library(VerticalCR)
      }
      NULL
    })
    parallel::parLapply(cl, seq_len(R), one_rep)
  } else {
    lapply(seq_len(R), one_rep)
  }
  raw_parameters <- do.call(rbind, lapply(reps, `[[`, "parameters"))
  raw_prediction <- do.call(rbind, lapply(reps, `[[`, "prediction"))
  list(table1 = manuscript_summarize_table1(raw_parameters),
       table2 = manuscript_summarize_table2(raw_prediction),
       raw_parameters = raw_parameters,
       raw_prediction = raw_prediction)
}

manuscript_run_factorial_tables <- function(R = 500L,
                                            workers = 1L,
                                            scenarios = expand.grid(
                                              censoring = c(0.15, 0.30),
                                              theta = c(0.5, 1.0),
                                              zeta = c(0.35, 0.50),
                                              KEEP.OUT.ATTRS = FALSE),
                                            start = 1L,
                                            end = nrow(scenarios)) {
  rows1 <- rows2 <- list()
  for (s in seq.int(start, min(end, nrow(scenarios)))) {
    sc <- scenarios[s, ]
    res <- manuscript_run_table_setting(
      R = R, workers = workers, censoring = sc$censoring,
      theta = sc$theta, zeta = sc$zeta)
    t1 <- res$table1
    t2 <- res$table2
    t1$censoring <- sc$censoring
    t1$theta_true <- sc$theta
    t1$zeta_true <- sc$zeta
    t2$censoring <- sc$censoring
    t2$theta_true <- sc$theta
    t2$zeta_true <- sc$zeta
    rows1[[length(rows1) + 1L]] <- t1
    rows2[[length(rows2) + 1L]] <- t2
  }
  list(table1 = do.call(rbind, rows1), table2 = do.call(rbind, rows2))
}
