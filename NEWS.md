# glmsSTARMA 1.0.2

## New Features
- Added argument `const` to the `vgamma()` family function, which is used in the transformation of past observations when the log-link is used. Past observations are now transformed as `log(y + const)` instead of `log(y)`.
- Added argument `dispersion_constant` to the `dglmstarma.control()` and `glmstarma_sim.control()` function. The value of this argument is passed to the `const` argument of the `vgamma()` family function when the log-link is used for the dispersion model in the `dglmstarma`-function.


## Bug Fixes
- Fixed bug, that prevented a dispersion model with identity link to be fitted
- The documentation of the `SST` dataset now lists the correct number of locations/grid points
- `delete_glmSTARMA_data()` now gives correct messages when datasets do not exist in the cache.
- Fixed bug in the `variance()` and `dev.resids()` functions of the `vbinomial()` and `vquasibinomial()` family functions
- Fixed bug in the `variance()` function of the `vnormal()` family function, which now returns the correct dimensions when the `ignore_dispersion` argument is set to `TRUE` and `mu` is not a matrix.
- Fixed bug bug in `fitted.dglmstarma()` which prevented to remove initial time points when `drop_init = TRUE`.


# glmsSTARMA 1.0.1

## Bug Fixes
- Fixed bug, that residuals of a model could not be scaled by the variance