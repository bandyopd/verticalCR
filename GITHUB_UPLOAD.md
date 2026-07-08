# GitHub upload checklist for `VerticalCR`

This package is intended to accompany the manuscript on frailty-assisted
vertical competing risks models with threshold effects.

## Repository contents

Keep the following files/directories under version control:

- `DESCRIPTION`
- `NAMESPACE`
- `README.md`
- `NEWS.md`
- `INSTALL.md`
- `R/`
- `src/*.cpp`
- `src/Makevars`
- `src/Makevars.win`
- `man/`
- `.Rbuildignore`
- `.gitignore`

Do not commit generated build artifacts:

- `*.tar.gz`
- `*.zip`
- `*.Rcheck/`
- `..Rcheck/`
- `src/*.o`
- `src/*.so`
- `src/*.dll`
- `src/symbols.rds`
- `.Rhistory`
- `.Rproj.user/`

These files are excluded by `.gitignore` and/or `.Rbuildignore`.

## Recommended repository description

Frailty-assisted vertical competing risks models with threshold effects.

## Suggested GitHub topics

`r-package`, `survival-analysis`, `competing-risks`, `frailty-model`,
`threshold-model`, `rcpp`, `transplantation`

## Build and check locally

From the parent directory of the package, run:

```r
install.packages(c("Rcpp", "RcppArmadillo", "survival", "splines", "nnet"))
```

Then from a shell:

```sh
R CMD build VerticalCR
R CMD check VerticalCR_1.1.0.tar.gz
```

On Windows, ensure that Rtools is fully available on the build shell path. The
build needs both the compiler (`g++`) and MSYS utilities such as `sed`,
`basename`, and `rm`.

For the local setup used in this project, the build command may require paths
similar to:

```powershell
$env:BINPREF = "D:/App/rtools45/x86_64-w64-mingw32.static.posix/bin/"
$env:COMPILER_PATH = "D:/App/rtools45/x86_64-w64-mingw32.static.posix/bin"
$env:PATH = "D:/App/rtools45/usr/bin;D:/App/rtools45/x86_64-w64-mingw32.static.posix/bin;" + $env:PATH
& "D:/App/R-4.6.0/bin/x64/R.exe" CMD build .
```

If compilation reaches the C++ files but fails with messages such as
`sed: command not found`, `basename: No such file or directory`, or
`rm: command not found`, the Rtools MSYS directory is not visible to the shell
used by `R CMD build`.

## Minimal smoke test after installation

```r
library(VerticalCR)

set.seed(1)
dat <- manuscript_simulate_twocause_data(
  N = 6,
  ni = rep(20, 6),
  theta = 0.5,
  censoring_rate = 0.30
)

fit <- fit_twocause_vertical(
  data = dat,
  time_var = "time",
  status_var = "status",
  threshold_var = "W",
  cluster_var = "cluster",
  total_vars = c("X1", "X2"),
  relative_vars = c("Z1", "Z2"),
  spline_df = 4,
  mite = 10,
  miter = 10,
  eps = 1e-4
)

extract_twocause_parameters(fit)
```

## Upload steps

1. Create an empty GitHub repository named `VerticalCR`.
2. In the local package directory, set the remote:

   ```sh
   git remote add origin https://github.com/<user-or-org>/VerticalCR.git
   ```

   If an origin already exists, update it:

   ```sh
   git remote set-url origin https://github.com/<user-or-org>/VerticalCR.git
   ```

3. Commit the package files:

   ```sh
   git add DESCRIPTION NAMESPACE README.md NEWS.md INSTALL.md GITHUB_UPLOAD.md \
       .Rbuildignore .gitignore R src man
   git commit -m "Prepare VerticalCR package for manuscript release"
   ```

4. Push:

   ```sh
   git branch -M main
   git push -u origin main
   ```

