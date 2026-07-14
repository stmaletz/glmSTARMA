# glmsSTARMA 1.0.2

## New Features
- Added argument `const` to the `vgamma()` family function, which is used in the transformation of past observations when the log-link is used. Past observations are now transformed as `log(y + const)` instead of `log(y)`.
- Added argument `dispersion_constant` to the `dglmstarma.control()` and `glmstarma_sim.control()` function. The value of this argument is passed to the `const` argument of the `vgamma()` family function when the log-link is used for the dispersion model in the `dglmstarma`-function.
- Added `withr` to the `Suggests` field in the DESCRIPTION file. The package is now used in the `test-data.R` test.
- Removed redundant checks for NA and infinite values in the `ts` argument of the `glmstarma()` and `dglmstarma()` functions.



## Bug Fixes
- Fixed bug, that prevented a dispersion model with identity link to be fitted
- The documentation of the `SST` dataset now lists the correct number of locations/grid points
- `delete_glmSTARMA_data()` now gives correct messages when datasets do not exist in the cache.
- Fixed bug in the `variance()` and `dev.resids()` functions of the `vbinomial()` and `vquasibinomial()` family functions
- Fixed bug in the `variance()` function of the `vnormal()` family function, which now returns the correct dimensions when the `ignore_dispersion` argument is set to `TRUE` and `mu` is not a matrix.
- Fixed bug bug in `fitted.dglmstarma()` which prevented to remove initial time points when `drop_init = TRUE`.
- `print.summary.glmstarma()` and `print.summary.dglmstarma()` can now handle empty models
- Fixed bug in `predict.glmstarma()` that prevented predictions for models with homogeneous intercept.
- Fixed bugs in `glmstarma_predict` and `dglmstarma_predict` that always caused warnings, when time-varying covariates were included in the model
- Fixed bug in `predict.dglmstarma()` that ignored specified copulas for prediction type `"sample"`
- Fixed bug in `dglmstarma.control()`, where `parameter_init_dispersion` was not checked
- Removed the incorrect reference to a `vgarch` family from the `stfamily` help page
- Changed observation transformation in case of the `vpoisson("sqrt")`, `vquasipoisson("sqrt")` and `vnegative.binomial("sqrt")` to `sqrt(x)` as `sqrt(x + 3 / 8) * 2` was causing problems.

# glmsSTARMA 1.0.1

## Bug Fixes
- Fixed bug, that residuals of a model could not be scaled by the variance