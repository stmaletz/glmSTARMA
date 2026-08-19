# glmsSTARMA 1.2.0

## New features
- Added the `plot()` method for spatially constant covariates.
- Output of `generateW()` is now assigned the class `generateW` for which `print()` and `plot()` methods are provided.

## Improvements
- Updated documentation and help-page titles.
- Updated the names of the distributions in the `stfamily` class used internally. 
- Improved output of the `print()` method for spatially constant and time constant covariates.


# glmsSTARMA 1.1.0

## New Features
- Added the argument `const` to the `vgamma()` family function. It is used in the transformation of past observations when the log-link is used: past observations are now transformed as `log(y + const)` instead of `log(y)`.
- Added the argument `dispersion_constant` to the `dglmstarma.control()` and `glmstarma_sim.control()` functions. Its value is passed to the `const` argument of `vgamma()` when the log-link is used for the dispersion model in `dglmstarma()`.

## Minor Improvements
- The log-likelihood of the `QuasiBinomial` family is now more robust against numerical issues.
- Removed redundant checks for `NA` and infinite values in the `ts` argument of the `glmstarma()` and `dglmstarma()` functions.
- Removed unused C++ methods that were left over from a preliminary version of the package.
- Added `withr` to the `Suggests` field in the `DESCRIPTION` file; it is now used in the `test-data.R` test.

## Bug Fixes
- Fixed a bug that prevented a dispersion model with an identity link from being fitted.
- Fixed a bug in the `variance()` and `dev.resids()` functions of the `vbinomial()` and `vquasibinomial()` family functions.
- Fixed a bug in the `dev.resids()` function of the `vbinomial()` and `vquasibinomial()` family functions that caused problems when observed values were on the boundary of the support.
- Fixed a bug in the `variance()` function of the `vnormal()` family function, which now returns the correct dimensions when `ignore_dispersion = TRUE` and `mu` is not a matrix.
- The default dispersion parameter for `vquasibinomial()` is now correctly set to `1` instead of `0` (not allowed) if the user does not specify a value and it must be estimated.
- Fixed a bug in the `vnegative.binomial()` family that did not allow zero values in the dispersion argument. For zero-values in the dispersion, the Negative Binomial distribution reduces to a Poisson distribution.
- Fixed a bug in the `variance_fun` of the `InversGauss` family in `C++`.
- Changed the observation transformation for `vpoisson("sqrt")`, `vquasipoisson("sqrt")`, and `vnegative.binomial("sqrt")` to `sqrt(x)`, as `sqrt(x + 3/8) * 2` was causing problems.
- Fixed a bug in the `link_trafo()` and `derivative_link_trafo()` functions of the `SoftClippingBinomial` and `SoftClippingQuasiBinomial` family functions in `C++`, which caused problems when different values of `n` were used across locations.
- Fixed a bug in `dglmstarma.control()` where `parameter_init_dispersion` was not checked.
- Fixed a bug in `predict.dglmstarma()` that ignored the specified copula for prediction type `"sample"`.
- Fixed bugs in `glmstarma_predict()` and `dglmstarma_predict()` that always caused warnings when time-varying covariates were included in the model.
- Fixed a bug in `predict.glmstarma()` that prevented predictions for models with a homogeneous intercept.
- Fixed a bug in `fitted.dglmstarma()` that prevented initial time points from being removed when `drop_init = TRUE`.
- Fixed a bug that prevented residuals of a model from being scaled by the variance.
- `print.summary.glmstarma()` and `print.summary.dglmstarma()` can now handle empty models.
- `delete_glmSTARMA_data()` now gives correct messages when datasets do not exist in the cache.
- The documentation of the `SST` dataset now lists the correct number of locations/grid points.
- Removed an incorrect reference to a `vgarch` family from the `stfamily` help page.
