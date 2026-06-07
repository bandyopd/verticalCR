# verticaltry

`verticaltry` contains the model-estimation code used for the manuscript
analyses of frailty-assisted vertical competing risks models with threshold
effects.

The current high-level interface focuses on the two-cause vertical model used
in the final simulation and real-data analyses. The frailty vector is fitted
with the sum-to-zero constraint

```text
U = C v,    sum_i U_i = 0,
```

where `C` is a Helmert contrast matrix generated from the cluster labels.

## Installation

From GitHub:

```r
install.packages("remotes")
remotes::install_github("BingWang-AHU/verticaltry")
```

From a source tarball:

```r
install.packages(c("Rcpp", "RcppArmadillo", "survival", "splines", "nnet"))
install.packages("verticaltry_1.1.0.tar.gz", repos = NULL, type = "source")
```

During development, install directly from the package directory:

```r
install.packages(".", repos = NULL, type = "source")
```

If R is not on the system `PATH`, run the commands with the full Rscript path,
for example:

```powershell
& "D:\App\R-4.6.0\bin\x64\Rscript.exe" -e "install.packages('.', repos=NULL, type='source')"
```

## Data format

The main fitting function expects:

- `time`: observed follow-up time.
- `status`: `0` for censoring, `1` for cause 1, and `2` for the reference
  cause.
- `X`: numeric matrix for the total hazard component.
- `Z`: numeric matrix for the relative cause probability component. If omitted,
  `Z = X`.
- `W`: threshold variable. The threshold grid and starting values must be on the
  same scale as `W`; the manuscript analyses use a normalized `W` in `[0, 1]`.
- `cluster`: transplant-center or cluster identifier.

## Basic fitting example with named covariates

The recommended interface uses a data frame and names the total-hazard and
relative-hazard covariates separately. This avoids accidental reuse of the same
design matrix in both model components.

```r
library(verticaltry)

fit <- fit_twocause_vertical(
  data = dat,
  time_var = "time",
  status_var = "status",
  threshold_var = "W",
  cluster_var = "center",
  total_vars = c("x_total1", "x_total2"),
  relative_vars = c("z_rel1", "z_rel2"),
  spline_df = 5,
  init_zetah = 0.5,
  init_zetal = 0.5,
  init_theta = 0.5,
  lr = 0.7,
  mite = 80,
  miter = 80,
  eps = 1e-5,
  rho_bound = 1,
  zeta_lr_pi = 0.5
)

fit$converged
fit$zetah
fit$zetal
fit$theta
```

The older matrix interface is still supported:

```r
fit <- fit_twocause_vertical(
  time = dat$time,
  status = dat$status,
  X = as.matrix(dat[, c("x_total1", "x_total2")]),
  Z = as.matrix(dat[, c("z_rel1", "z_rel2")]),
  W = dat$W,
  cluster = dat$center
)
```

## Choosing the spline dimension

The relative hazard component uses a B-spline approximation for the unknown
time function. Candidate degrees of freedom can be compared by an AIC/BIC-style
criterion:

```r
df_sel <- select_spline_df_twocause(
  data = dat,
  time_var = "time",
  status_var = "status",
  threshold_var = "W",
  cluster_var = "center",
  total_vars = c("x_total1", "x_total2"),
  relative_vars = c("z_rel1", "z_rel2"),
  candidate_df = 3:8,
  criterion = "BIC",
  lr = 0.7,
  mite = 50,
  miter = 50,
  zeta_lr_pi = 0.5
)

df_sel$selected_df
df_sel$summary
```

Set `keep_fits = TRUE` if you want the best fitted object returned as
`df_sel$best_fit`.

## Parameter table

```r
tab <- extract_twocause_parameters(
  fit
)
print(tab)
```

The output contains the estimate, standard error, Wald statistic, and Wald
p-value for the regression parameters, thresholds, spline coefficients, frailty
loading, and frailty variance.

## Prediction

For prediction on centers observed in the training data, empirical Bayes
frailty estimates can be used:

```r
pred <- predict_twocause_vertical(
  fit,
  newdata = list(
    X = as.matrix(test[, c("x1", "x2")]),
    Z = as.matrix(test[, c("x1", "x2")]),
    W = test$W,
    cluster = test$center
  ),
  times = seq(0.25, 6, by = 0.25),
  use_EB_frailty = TRUE
)

dim(pred$F1)
dim(pred$F2)
```

For new centers not present in the training data, the prediction automatically
uses frailty zero for those centers.

## Chunked optimization for full-data analyses

For large data sets, use chunked optimization. This repeatedly calls the C++
optimizer for a small number of outer iterations and warm-starts the next chunk.
It can save checkpoints after each chunk.

```r
fit <- fit_twocause_vertical_chunked(
  time = dat$time,
  status = dat$status,
  X = as.matrix(dat[, c("x1", "x2")]),
  Z = as.matrix(dat[, c("x1", "x2")]),
  W = dat$W,
  cluster = dat$center,
  chunk_size = 2,
  max_chunks = 60,
  checkpoint_dir = "outputs/checkpoints_twocause",
  spline_df = 5,
  lr = 0.7,
  mite = 20,
  eps = 1e-5,
  rho_bound = 1,
  zeta_lr_pi = 0.5
)
```

The final returned object has the same structure as `fit_twocause_vertical()`.

## Low-level C++ interface

The package also exports the raw Rcpp functions:

- `tw_vertical()`
- `tw_verticalite2()`
- `tw_verticalite3()`
- `tw_updateb()`
- `tw_updateg()`
- `tw_updateu()`

These are useful for reproducing legacy scripts, but routine analyses should
use `fit_twocause_vertical()` or `fit_twocause_vertical_chunked()`.
