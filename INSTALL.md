# Installation notes

This package contains Rcpp/RcppArmadillo source code. A working C++ toolchain
is required when installing from source.

## Recommended Windows setup

Use the Rscript path configured for this project:

```powershell
$env:PATH = "D:\App\rtools45\usr\bin;D:\App\rtools45\x86_64-w64-mingw32.static.posix\bin;" + $env:PATH
$env:BINPREF = "D:/App/rtools45/x86_64-w64-mingw32.static.posix/bin/"
$env:COMPILER_PATH = "D:/App/rtools45/x86_64-w64-mingw32.static.posix/bin"

& "D:\App\R-4.6.0\bin\x64\Rscript.exe" -e "Rcpp::compileAttributes('VerticalCR')"
& "D:\App\R-4.6.0\bin\x64\R.exe" CMD build VerticalCR
& "D:\App\R-4.6.0\bin\x64\R.exe" CMD INSTALL VerticalCR_1.1.0.tar.gz
```

If `R CMD build` fails at the linking stage with messages such as
`sed: command not found`, `rm: command not found`, or `basename: No such file or
directory`, the Rtools MSYS utilities are not visible to the `make` subprocess.
Check that the following files exist and that their directory is at the front of
`PATH`:

```text
D:\App\rtools45\usr\bin\sed.exe
D:\App\rtools45\usr\bin\rm.exe
D:\App\rtools45\usr\bin\basename.exe
D:\App\rtools45\x86_64-w64-mingw32.static.posix\bin\g++.exe
```

The package source itself is independent of those local paths; they are only
needed for compiling the C++ shared library.

## Development workflow

After changing files in `src/`, regenerate the Rcpp wrappers:

```r
Rcpp::compileAttributes("VerticalCR")
```

The high-level fitting interface is in `R/model_estimation.R`. The low-level
C++ implementation used by the two-cause model is in `src/twocase-blockc.cpp`.
