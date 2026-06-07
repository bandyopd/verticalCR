# verticaltry 1.1.0

## New features

- Added high-level two-cause vertical competing risks model fitting through
  `fit_twocause_vertical()`.
- Added chunked warm-start fitting for large real-data analyses through
  `fit_twocause_vertical_chunked()`.
- Added named covariate interfaces separating total-hazard and relative-cause
  design matrices.
- Added spline-dimension selection with AIC/BIC-style criteria through
  `select_spline_df_twocause()`.
- Added cumulative-incidence prediction through `predict_twocause_vertical()`.
- Added parameter extraction and manuscript simulation helper functions.

## Implementation notes

- Center frailties are represented by the sum-to-zero contrast
  `U = C v`, so that `sum(U) = 0`.
- The package uses Rcpp and RcppArmadillo for the core optimization routines.
- Low-level legacy Rcpp functions remain exported for reproducibility, but the
  recommended user-facing functions are the high-level wrappers documented in
  `README.md`.
