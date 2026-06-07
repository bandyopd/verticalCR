\name{fit_twocause_vertical}
\alias{fit_twocause_vertical}
\alias{fit_twocause_vertical_chunked}
\alias{extract_twocause_parameters}
\alias{predict_twocause_vertical}
\alias{loglik_twocause_vertical}
\alias{select_spline_df_twocause}
\alias{make_sum_to_zero_contrast}
\alias{make_cluster_indicator}
\alias{make_event_time_basis}
\alias{predict_event_time_basis}
\title{Fit and use a two-cause vertical threshold frailty model}
\description{
High-level utilities for fitting the two-cause vertical competing risks model
with threshold effects and shared center frailties. The frailty vector is
represented by a sum-to-zero contrast, \eqn{U = C v}, so that
\eqn{\sum_i U_i = 0}.
}
\usage{
fit_twocause_vertical(time = NULL, status = NULL, X = NULL, Z = NULL,
  W = NULL, cluster = NULL, data = NULL, total_vars = NULL,
  relative_vars = NULL, time_var = "time", status_var = "status",
  threshold_var = "W", cluster_var = "cluster",
  spline_df = 5L, init_beta = NULL, init_gamma = NULL,
  init_kappa = NULL, init_U = NULL, init_zetah = 0.5,
  init_zetal = 0.5, init_rho = 0.5, init_theta = 0.5,
  h1 = NA_real_, h2 = NA_real_, lr = 1, mite = 50L,
  miter = 50L, eps = 1e-05, rho_bound = 1,
  zeta_lr_pi = 1)

fit_twocause_vertical_chunked(time = NULL, status = NULL, X = NULL,
  Z = NULL, W = NULL, cluster = NULL, data = NULL,
  total_vars = NULL, relative_vars = NULL, time_var = "time",
  status_var = "status", threshold_var = "W",
  cluster_var = "cluster",
  chunk_size = 2L, max_chunks = 50L, checkpoint_dir = NULL, ...)

extract_twocause_parameters(fit, covariate_names = NULL,
  total_names = NULL, relative_names = NULL)

predict_twocause_vertical(fit, newdata = NULL, times,
  use_EB_frailty = TRUE)

loglik_twocause_vertical(fit)

select_spline_df_twocause(time = NULL, status = NULL, X = NULL,
  Z = NULL, W = NULL, cluster = NULL, data = NULL,
  total_vars = NULL, relative_vars = NULL, time_var = "time",
  status_var = "status", threshold_var = "W",
  cluster_var = "cluster", candidate_df = 3:8,
  criterion = c("BIC", "AIC"), keep_fits = FALSE, ...)
}
\arguments{
\item{time}{Observed follow-up times.}
\item{status}{Event status coded as 0 for censored, 1 for cause 1, and 2 for
the reference cause.}
\item{X}{Numeric matrix for the total hazard component.}
\item{Z}{Numeric matrix for the relative cause component. Defaults to
\code{X}.}
\item{W}{Threshold variable, on the same scale as the starting thresholds.}
\item{cluster}{Cluster or center identifier.}
\item{data}{Optional data frame. When supplied, use \code{time_var},
\code{status_var}, \code{threshold_var}, \code{cluster_var},
\code{total_vars}, and \code{relative_vars} to define the model variables.}
\item{total_vars}{Character vector naming covariates in the total hazard
component. Required when \code{data} is supplied.}
\item{relative_vars}{Character vector naming covariates in the relative hazard
component. Required when \code{data} is supplied.}
\item{time_var, status_var, threshold_var, cluster_var}{Column names in
\code{data}.}
\item{spline_df}{Degrees of freedom for the B-spline time basis in the relative
hazard component.}
\item{init_beta, init_gamma, init_kappa, init_U, init_zetah, init_zetal,
init_rho, init_theta}{Initial values. If \code{init_U} has length equal to the
number of clusters it is centered and projected to contrast coordinates; if it
has length \eqn{N-1}, it is treated as contrast coefficients.}
\item{h1, h2}{Bandwidths for the total-hazard and relative-hazard threshold
smooths. If omitted, data-adaptive defaults are used.}
\item{lr}{Newton update learning rate.}
\item{mite}{Maximum inner iterations.}
\item{miter}{Maximum outer iterations.}
\item{eps}{Convergence tolerance.}
\item{rho_bound}{Absolute bound applied to the frailty loading parameter.}
\item{zeta_lr_pi}{Additional learning-rate multiplier for the relative-hazard
threshold update.}
\item{chunk_size}{Number of outer iterations per chunk for chunked fitting.}
\item{max_chunks}{Maximum number of chunks.}
\item{checkpoint_dir}{Optional directory for chunk checkpoint RDS files.}
\item{...}{Additional arguments passed to \code{fit_twocause_vertical}.}
\item{fit}{A fitted object returned by \code{fit_twocause_vertical}.}
\item{covariate_names}{Optional names for the columns of \code{X} and
\code{Z}.}
\item{total_names, relative_names}{Optional names for the total-hazard and
relative-hazard covariates. These override names stored in the fitted object.}
\item{newdata}{Optional list with elements \code{X}, \code{Z}, \code{W}, and
\code{cluster}. If omitted, predictions are returned for the training data.}
\item{times}{Prediction times.}
\item{use_EB_frailty}{Whether to use empirical Bayes frailty estimates for
clusters observed during training.}
\item{candidate_df}{Integer vector of candidate spline degrees of freedom.}
\item{criterion}{Model-selection criterion used for selecting
\code{spline_df}; either \code{"BIC"} or \code{"AIC"}.}
\item{keep_fits}{Whether to keep all candidate fitted objects and return the
best fitted object.}
}
\value{
\code{fit_twocause_vertical} returns a list containing parameter estimates,
standard errors, frailty estimates, bandwidths, contrast matrices, basis
information, and convergence diagnostics. \code{extract_twocause_parameters}
returns a data frame. \code{predict_twocause_vertical} returns a list with
\code{times}, \code{F1}, and \code{F2}. \code{select_spline_df_twocause}
returns a list containing the selected spline dimension and a comparison table.
}
\examples{
\dontrun{
fit <- fit_twocause_vertical(
  data = dat,
  time_var = "time",
  status_var = "status",
  threshold_var = "W",
  cluster_var = "center",
  total_vars = c("x1", "x2"),
  relative_vars = c("z1", "z2"),
  lr = 0.7,
  zeta_lr_pi = 0.5
)

extract_twocause_parameters(fit)
predict_twocause_vertical(fit, times = seq(0.25, 6, by = 0.25))
}
}
