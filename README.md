# glmSTARMA

## Spatio-temporal Models based on (Double) Generalized Linear Models
<!-- badges: start -->
[![R-CMD-check](https://github.com/stmaletz/glmSTARMA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/stmaletz/glmSTARMA/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/stmaletz/glmSTARMA/graph/badge.svg)](https://app.codecov.io/gh/stmaletz/glmSTARMA)
[![CRAN status](https://www.r-pkg.org/badges/version/glmSTARMA)](https://CRAN.R-project.org/package=glmSTARMA)
<!-- badges: end -->

## Installation

The official version of glmSTARMA can be installed using the R Console:

``` r
install.packages("glmSTARMA")
```

You can install the development version of glmSTARMA from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("stmaletz/glmSTARMA")
```
## Structure

This package implements spatio-temporal autoregressive moving average (STARMA)–type models that account for both spatial and temporal dependence. Spatial dependence is introduced via spatial weight matrices, while temporal dependence is modeled through lagged observations and lagged values of the linear predictor.

The core model-fitting functions are `glmstarma()` and `dglmstarma()`:
`glmstarma()` fits a model for the conditional mean of a spatio-temporal process, whereas `dglmstarma()` extends this framework by jointly modeling the conditional mean and the conditional dispersion.
The mean model generalizes spatio-temporal Poisson autoregressive models and supports a range of distributions from the exponential dispersion family. The dispersion model can be interpreted as a spatio-temporal extension of GARCH or log-GARCH–type models.
For simulation studies, synthetic data can be generated using `glmstarma.sim()` and `dglmstarma.sim()`.

Das Paket enthält außerdem drei vorverarbeitete Datensätze. Die Rohdaten, sowie die Skripte die zur Vorverarbeitung verwendet wurden befinden sich in den Unterordnern des Verzeichnisses `data-raw` von diesem Repository:

* `rota`: Rota Virus Infections in Germany
* `chickenpox`: Chickenpox Infections in Hungary
* `SST`: Sea Surface Temperature Anomalies in the Pacific


## Example








