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

This package also includes preprocessed datasets. The raw data and the scripts used for preprocessing can be found in the subdirectories of the `data-raw` directory of this repository.
The following datasets are included:

### `rota`: Rotavirus Infections in Germany

This dataset contains weekly counts of reported rotavirus infections in the 411 German counties from January 2001 to December 2024, along with population data and additional covariates.

The infection data were retrieved from the Robert Koch Institute (RKI) via the SurvStat@RKI 2.0 application  
(https://survstat.rki.de). Population data at the county level were obtained from the Federal Statistical Office of Germany (Destatis)  
(https://www.destatis.de). To disaggregate population data for Berlin into its districts, census data from 2022 were used  
(https://ergebnisse.zensus2022.de). Since population figures are only available on an annual basis, they were linearly interpolated to obtain weekly values.

The directory `data-raw/rota/` additionally contains shapefiles for the German counties, which were copied from  
https://github.com/ostojanovic/BSTIM. The original shapefiles are provided by the German Federal Agency for Cartography and Geodesy  
(Bundesamt für Kartographie und Geodäsie; © GeoBasis-DE / BKG 2018) and are licensed under the  
[Data licence Germany – attribution – version 2.0](https://www.govdata.de/dl-de/by-2-0).  
The row-normalized adjacency matrices were constructed based on these shapefiles.

The population and census data (retrieved on 2025-12-08) are licensed under the  
[Data licence Germany – attribution – version 2.0](https://www.govdata.de/dl-de/by-2-0).

The infection data are subject to the terms of use of SurvStat@RKI  
(https://survstat.rki.de/Content/Instruction/DataUsage.aspx).  
The data included in this package are aggregated and preprocessed for methodological research purposes.  

Please cite the infection data as follows:
- Robert Koch-Institut: SurvStat@RKI 2.0, https://survstat.rki.de, retrieved on 2025-02-03.

### `chickenpox`: Chickenpox Infections in Hungary

This dataset contains weekly counts of reported chickenpox infections in the 20 regions of Hungary from January 2005 to December 2015.

The infection data were obtained from the UCI Machine Learning Repository  
(https://archive.ics.uci.edu/dataset/580/hungarian+chickenpox+cases) and complemented with population data from the Hungarian Central Statistical Office  
(https://www.ksh.hu/stadat_files/nep/en/nep0034.html). Since population data are only available on an annual basis, they were linearly interpolated to obtain weekly values.

The data preprocessing scripts are available in the `data-raw/chickenpox/` directory of this repository.

The data are licensed under the Creative Commons Attribution 4.0 International  
[(CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/legalcode).


### `SST`: Sea Surface Temperature Anomalies in the Pacific

This dataset is a transformed subset of the `SST_df` dataset included in the (archived) R package `STRbook`, which accompanies the books listed below.  
The package remains available on GitHub at https://github.com/andrewzm/STRbook.

The data included in this package have been subsetted and reformatted for use in spatio-temporal modeling examples. The original `STRbook` package is distributed under the GNU General Public License (GPL).

For further details on the original dataset and its scientific background, see:
- Wikle, C. K., Zammit-Mangion, A., and Cressie, N. (2019). *Spatio-Temporal Statistics with R*. Chapman & Hall/CRC, Boca Raton, FL.
- Cressie, N., and Wikle, C. K. (2011). *Statistics for Spatio-Temporal Data*. John Wiley & Sons.



## Example
We illustrate the usage of the `glmSTARMA` package using the `rota` dataset, which contains weekly counts of reported rotavirus infections in the 411 German counties from January 2001 to December 2024.

We first fit a log-linear model assuming marginal conditional Poisson distributions for the counts, as proposed by Maletz et al. (2024). Specifically, we assume

$$
Y_{i,t} \mid \mathcal{F}_{t-1} \sim \text{Poisson}(\mu_{i,t}), \quad i = 1, \ldots, 411, \quad t = 1, \ldots, 1252,
$$

and define $\boldsymbol{\psi}_t := \log(\boldsymbol{\mu}_t)$ elementwise. The linear predictor is given by

$$
\boldsymbol{\psi}_t
= \beta_0 \mathbf{1}
+ \sum_{i = 1}^4 \sum_{\ell = 0}^2 \beta_{i,\ell} \mathbf{W}^{(\ell)} \log(\mathbf{Y}_{t - i} + \mathbf{1})
+ \sum_{k = 1}^{6} \gamma_k \mathbf{X}_{k, t},
$$

where $\mathbf{X}_{k,t}$, $k = 1, \ldots, 6$, are covariates describing population size, indicators for former German Democratic Republic regions, yearly seasonal effects, and vaccination effects for eastern and western Germany.  
The matrices $\mathbf{W}^{(0)}$, $\mathbf{W}^{(1)}$, and $\mathbf{W}^{(2)}$ denote the identity matrix, a row-normalized adjacency matrix, and a row-normalized second-order adjacency matrix, respectively.

The model can be fitted as follows:

```r
library(glmSTARMA)
data("rota")

covariates <- list(
  population   = population_germany,
  gdr          = TimeConstant(gdr_feature),
  season_cos   = SpatialConstant(cos(2 * pi / 52 * 1:1252)),
  season_sin   = SpatialConstant(sin(2 * pi / 52 * 1:1252)),
  vaccine_west = (gdr_feature == 0) %*% t(seq(ncol(rota)) >= 654),
  vaccine_east = (gdr_feature > 0) %*% t(seq(ncol(rota)) >= 654)
)

fit <- glmstarma(
  rota,
  list(past_obs = rep(2, 4)),
  wlist = W_germany,
  covariates = covariates,
  family = vpoisson("log")
)
```
A related model is described by Jahn et al. (2023), who propose a similar specification but use a different link function and do not transform past observations. This model can be fitted by adapting the code as follows:
```r
fit <- glmstarma(
  rota,
  list(past_obs = rep(2, 4)),
  wlist = W_germany,
  covariates = covariates,
  family = vpoisson("softplus")
)
```
### References
-   Jahn, M., Weiß, C.H., Kim, H.Y. (2023), Approximately linear INGARCH models for spatio-temporal counts, Journal of the Royal Statistical Society Series C: Applied Statistics, 72(2), 476-497, [DOI: 10.1093/jrsssc/qlad018](doi.org/10.1093/jrsssc/qlad018)
-   Maletz, S., Fokianos, K., & Fried, R. (2024). Spatio-Temporal Count Autoregression. Data Science in Science, 3(1). [DOI: 10.1080/26941899.2024.2425171](doi.org/10.1080/26941899.2024.2425171)
-   Smyth, G.K. (1989), Generalized Linear Models with Varying Dispersion. Journal of the Royal Statistical Society: Series B (Methodological), 51(1), 47-60. [DOI: 10.1111/j.2517-6161.1989.tb01747.x](doi.org/10.1111/j.2517-6161.1989.tb01747.x)






