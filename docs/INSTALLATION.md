# Installation

The repository is written in **R**.

## Core dependencies

```text
MASS
MCMCpack
```

The local simulation runners additionally use packages such as:

```text
Matrix
mcmcse
```

## Minimal installation

```r
install.packages(c("MASS", "MCMCpack"))
```

For the local simulation runners:

```r
install.packages(c("MASS", "Matrix", "MCMCpack", "mcmcse"))
```

## Non-interactive installation

A repository helper may install missing packages and redirect installation to `R_LIBS_USER` when the default site library is not writable.

Run:

```bash
Rscript install_packages.R
```

The `bash` label only indicates a terminal command block; the executed program is R.

## Load the estimator

```r
source("MAR_paper_core_standalone.R")
```

Verify:

```r
exists("fit_mar_pointmass")
exists("fit_mar_ssvs")
```

Both should return `TRUE`.

## Publication environment

For archival reproducibility, record:

```r
R.version.string
sessionInfo()
```

An `renv.lock` file may also be added to pin package versions.
