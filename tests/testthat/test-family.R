# Test family functions

# test that vbinomial and vquasibinomial size is valid
test_that("size parameter is valid (vbinomial)", {
    expect_error(vbinomial(size = -1))
    expect_error(vbinomial(size = 0))
    expect_error(vbinomial(size = 1.5))
    mat <- matrix(c(1, 1, 1, 1), nrow = 2)
    expect_error(vbinomial(size = mat))
    # correct cases
    expect_s3_class(vbinomial(size = 1), "stfamily")
    expect_s3_class(vbinomial(size = c(1, 2, 3)), "stfamily")
    expect_s3_class(vbinomial(size = 5), "stfamily")
})

test_that("size parameter is valid (vbinomial)", {
    expect_error(vbinomial(size = -1))
    expect_error(vbinomial(size = 0))
    expect_error(vbinomial(size = 1.5))
    mat <- matrix(c(1, 1, 1, 1), nrow = 2)
    expect_error(vbinomial(size = mat))
    # correct cases
    expect_s3_class(vbinomial(size = 1), "stfamily")
    expect_s3_class(vbinomial(size = c(1, 2, 3)), "stfamily")
    expect_s3_class(vbinomial(size = 5), "stfamily")
})

# test that link is valid for the family
test_that("link parameter is valid (vpoisson)", {
  expect_error(vpoisson(link = "invalid_link"))
  expect_error(vpoisson(link = c("log", "identity")))
  # valid case
  expect_s3_class(vpoisson(link = "log"), "stfamily")
  expect_s3_class(vpoisson(link = "identity"), "stfamily")
  expect_s3_class(vpoisson(link = "sqrt"), "stfamily")
  expect_s3_class(vpoisson(link = "softplus"), "stfamily")
})

test_that("link parameter is valid (vquasipoisson)", {
  expect_error(vquasipoisson(link = "invalid_link"))
  expect_error(vquasipoisson(link = c("log", "identity")))
  # valid case
  expect_s3_class(vquasipoisson(link = "log"), "stfamily")
  expect_s3_class(vquasipoisson(link = "identity"), "stfamily")
  expect_s3_class(vquasipoisson(link = "sqrt"), "stfamily")
  expect_s3_class(vquasipoisson(link = "softplus"), "stfamily")
})

test_that("link parameter is valid (vbinomial)", {
  expect_error(vbinomial(link = "invalid_link"))
  expect_error(vbinomial(link = c("softclipping", "identity")))
  # valid case
  expect_s3_class(vbinomial(link = "softclipping"), "stfamily")
  expect_s3_class(vbinomial(link = "identity"), "stfamily")
  expect_s3_class(vbinomial(link = "logit"), "stfamily")
  expect_s3_class(vbinomial(link = "probit"), "stfamily")
})

test_that("link parameter is valid (vquasibinomial)", {
  expect_error(vquasibinomial(link = "invalid_link"))
  expect_error(vquasibinomial(link = c("softclipping", "identity")))
  # valid case
  expect_s3_class(vquasibinomial(link = "softclipping"), "stfamily")
  expect_s3_class(vquasibinomial(link = "identity"), "stfamily")
  expect_s3_class(vquasibinomial(link = "logit"), "stfamily")
  expect_s3_class(vquasibinomial(link = "probit"), "stfamily")
})

test_that("link parameter is valid (vgamma)", {
  expect_error(vgamma(link = "invalid_link"))
  expect_error(vgamma(link = c("inverse", "log")))
  # valid case
  expect_s3_class(vgamma(link = "inverse"), "stfamily")
  expect_s3_class(vgamma(link = "log"), "stfamily")
  expect_s3_class(vgamma(link = "identity"), "stfamily")
})

test_that("link parameter is valid (vinverse.gaussian)", {
  expect_error(vinverse.gaussian(link = "invalid_link"))
  expect_error(vinverse.gaussian(link = c("1/mu^2", "inverse")))
  # valid case
  expect_s3_class(vinverse.gaussian(link = "1/mu^2"), "stfamily")
  expect_s3_class(vinverse.gaussian(link = "inverse"), "stfamily")
  expect_s3_class(vinverse.gaussian(link = "identity"), "stfamily")
  expect_s3_class(vinverse.gaussian(link = "log"), "stfamily")
})

test_that("link parameter is valid (vnormal)", {
  expect_error(vnormal(link = "invalid_link"))
  expect_error(vnormal(link = c("identity", "log")))
  # valid case
  expect_s3_class(vnormal(link = "identity"), "stfamily")
  expect_s3_class(vnormal(link = "log"), "stfamily")
  expect_s3_class(vnormal(link = "inverse"), "stfamily")
})


# test that dispersion is NULL or a (vector or matrix) of positive number
test_that("dispersion parameter is valid (vquasipoisson)", {
  expect_error(vquasipoisson(dispersion = -1))
  expect_error(vquasipoisson(dispersion = 0))
  mat <- matrix(c(1, -0.5, 0.2, 1), nrow = 2)
  expect_error(vquasipoisson(dispersion = mat))
  mat[2, 1] <- NA
  expect_error(vquasipoisson(dispersion = mat))
  # valid case
  expect_s3_class(vquasipoisson(dispersion = c(1, 2, 3)), "stfamily")
  expect_s3_class(vquasipoisson(dispersion = matrix(c(1, 0.5, 0.5, 1), nrow = 2)), "stfamily")
  expect_s3_class(vquasipoisson(dispersion = 2), "stfamily")
  expect_s3_class(vquasipoisson(), "stfamily")
})

test_that("dispersion parameter is valid (vnegative.binomial)", {
  expect_error(vnegative.binomial(dispersion = -1))
  expect_error(vnegative.binomial(dispersion = 0))
  mat <- matrix(c(1, -0.5, 0.2, 1), nrow = 2)
  expect_error(vnegative.binomial(dispersion = mat))
  mat[2, 1] <- NA
  expect_error(vnegative.binomial(dispersion = mat))
  # valid case
  expect_s3_class(vnegative.binomial(dispersion = c(1, 2, 3)), "stfamily")
  expect_s3_class(vnegative.binomial(dispersion = matrix(c(1, 0.5, 0.5, 1), nrow = 2)), "stfamily")
  expect_s3_class(vnegative.binomial(dispersion = 2), "stfamily")
  expect_s3_class(vnegative.binomial(), "stfamily")
})

# warning if dispersion > size
test_that("dispersion parameter is valid (vquasibinomial)", {
  expect_error(vquasibinomial(dispersion = -1))
  expect_error(vquasibinomial(dispersion = 0))
  mat <- matrix(c(1, -0.5, 0.2, 1), nrow = 2)
  expect_error(vquasibinomial(dispersion = mat))
  mat[2, 1] <- NA
  expect_error(vquasibinomial(dispersion = mat))
  # valid case
  expect_s3_class(vquasibinomial(dispersion = c(1, 2, 3), size = 5), "stfamily")
  expect_s3_class(vquasibinomial(dispersion = matrix(c(1, 0.5, 0.5, 1), nrow = 2), size = 5), "stfamily")
  expect_s3_class(vquasibinomial(dispersion = 2, size = 5), "stfamily")
  expect_s3_class(vquasibinomial(), "stfamily")
  expect_s3_class(vquasibinomial(dispersion = c(1,2,3), size = c(3,4,5)), "stfamily")
  expect_s3_class(vquasibinomial(dispersion = matrix(c(1.5,0.5,0.5,1), nrow=2), size = c(3,4)), "stfamily")
  # warning case
  expect_warning(out <- vquasibinomial(dispersion = 5, size = 3))
  expect_s3_class(out, "stfamily")
  expect_warning(out <- vquasibinomial(dispersion = c(2,3,4), size = c(3,2,1)))
  expect_s3_class(out, "stfamily")
  expect_warning(out <- vquasibinomial(dispersion = matrix(c(2,0.5,0.5,2), nrow=2), size = c(3, 1)))
  expect_s3_class(out, "stfamily")
  #error when dismatch between size and dispersion
  expect_error(vquasibinomial(dispersion = c(1,2), size = 3:5))
  expect_error(vquasibinomial(dispersion = matrix(c(1,0.5,0.5,1), nrow=2), size = c(3,4,5)))  
})

test_that("dispersion parameter is valid (vgamma)", {
  expect_error(vgamma(dispersion = -1))
  expect_error(vgamma(dispersion = 0))
  mat <- matrix(c(1, -0.5, 0.2, 1), nrow = 2)
  expect_error(vgamma(dispersion = mat))
  mat[2, 1] <- NA
  expect_error(vgamma(dispersion = mat))
  # valid case
  expect_s3_class(vgamma(dispersion = c(1, 2, 3)), "stfamily")
  expect_s3_class(vgamma(dispersion = matrix(c(1, 0.5, 0.5, 1), nrow = 2)), "stfamily")
  expect_s3_class(vgamma(dispersion = 2), "stfamily")
  expect_s3_class(vgamma(), "stfamily")
})


test_that("dispersion parameter is valid (vinverse.gaussian)", {
  expect_error(vinverse.gaussian(dispersion = -1))
  expect_error(vinverse.gaussian(dispersion = 0))
  mat <- matrix(c(1, -0.5, 0.2, 1), nrow = 2)
  expect_error(vinverse.gaussian(dispersion = mat))
  mat[2, 1] <- NA
  expect_error(vinverse.gaussian(dispersion = mat))
  # valid case
  expect_s3_class(vinverse.gaussian(dispersion = c(1, 2, 3)), "stfamily")
  expect_s3_class(vinverse.gaussian(dispersion = matrix(c(1, 0.5, 0.5, 1), nrow = 2)), "stfamily")
  expect_s3_class(vinverse.gaussian(dispersion = 2), "stfamily")
  expect_s3_class(vinverse.gaussian(), "stfamily")
})

test_that("dispersion parameter is valid (vnormal)", {
  expect_error(vnormal(dispersion = -1))
  expect_error(vnormal(dispersion = 0))
  mat <- matrix(c(1, -0.5, 0.2, 1), nrow = 2)
  expect_error(vnormal(dispersion = mat))
  mat[2, 1] <- NA
  expect_error(vnormal(dispersion = mat))
  # valid case
  expect_s3_class(vnormal(dispersion = c(1, 2, 3)), "stfamily")
  expect_s3_class(vnormal(dispersion = matrix(c(1, 0.5, 0.5, 1), nrow = 2)), "stfamily")
  expect_s3_class(vnormal(dispersion = 2), "stfamily")
  expect_s3_class(vnormal(), "stfamily")
})

# test that const is valid for the family
test_that("const parameter is valid (quasi)-poisson", {
  expect_error(vpoisson(const = -1, link = "log"))
  expect_error(vquasipoisson(const = -1, link = "log"))
  expect_error(vpoisson(const = -1, link = "softplus"))
  expect_error(vquasipoisson(const = -1, link = "softplus"))
  expect_error(vpoisson(const = 0, link = "log"))
  expect_error(vquasipoisson(const = 0, link = "log"))
  expect_error(vpoisson(const = 0, link = "softplus"))
  expect_error(vquasipoisson(const = 0, link = "softplus"))
  expect_error(vpoisson(const = NA, link = "log"))
  expect_error(vquasipoisson(const = NA, link = "log"))
  expect_error(vpoisson(const = NA, link = "softplus"))
  expect_error(vquasipoisson(const = NA, link = "softplus"))
  expect_error(vpoisson(const = 1:2, link = "log"))
  expect_error(vquasipoisson(const = 1:2, link = "log"))
  expect_error(vpoisson(const = 1:2, link = "softplus"))
  expect_error(vquasipoisson(const = 1:2, link = "softplus"))
  expect_error(vpoisson(const = "a", link = "log"))
  expect_error(vquasipoisson(const = "a", link = "log"))
  expect_error(vpoisson(const = "a", link = "softplus"))
  expect_error(vquasipoisson(const = "a", link = "softplus"))

  # valid case
  expect_s3_class(vpoisson(const = 1, link = "log"), "stfamily")
  expect_s3_class(vquasipoisson(const = 1, link = "log"), "stfamily")
  expect_s3_class(vpoisson(const = 1, link = "softplus"), "stfamily")
  expect_s3_class(vquasipoisson(const = 1, link = "softplus"), "stfamily")
  expect_s3_class(vpoisson(const = 0.1, link = "log"), "stfamily")
  expect_s3_class(vquasipoisson(const = 0.1, link = "log"), "stfamily")
  expect_s3_class(vpoisson(const = 0.1, link = "softplus"), "stfamily")
  expect_s3_class(vquasipoisson(const = 0.1, link = "softplus"), "stfamily")
  expect_s3_class(vpoisson(link = "log"), "stfamily")
  expect_s3_class(vquasipoisson(link = "log"), "stfamily")
  expect_s3_class(vpoisson(link = "softplus"), "stfamily")
  expect_s3_class(vquasipoisson(link = "softplus"), "stfamily")

  expect_s3_class(vpoisson(const = 0, link = "identity"), "stfamily")
  expect_s3_class(vquasipoisson(const = 0, link = "identity"), "stfamily")
  expect_s3_class(vpoisson(const = 0, link = "sqrt"), "stfamily")
  expect_s3_class(vquasipoisson(const = 0, link = "sqrt"), "stfamily")
})

test_that("const parameter is valid (quasi)-binomial", {
  expect_error(vbinomial(const = -1, link = "softclipping"))
  expect_error(vbinomial(const = -1, link = "softclipping"))
  expect_error(vquasibinomial(const = -1, link = "softclipping"))
  expect_error(vbinomial(const = 0, link = "softclipping"))
  expect_error(vquasibinomial(const = 0, link = "softclipping"))
  expect_error(vbinomial(const = NA, link = "softclipping"))
  expect_error(vquasibinomial(const = NA, link = "softclipping"))
  expect_error(vbinomial(const = 1:2, link = "softclipping"))
  expect_error(vquasibinomial(const = 1:2, link = "softclipping"))
  expect_error(vbinomial(const = "a", link = "softclipping"))
  expect_error(vquasibinomial(const = "a", link = "softclipping"))
  # valid case
  expect_s3_class(vbinomial(const = 1, link = "softclipping"), "stfamily")
  expect_s3_class(vquasibinomial(const = 1, link = "softclipping"), "stfamily")
  expect_s3_class(vbinomial(const = 1, link = "softclipping"), "stfamily")
  expect_s3_class(vquasibinomial(const = 1, link = "softclipping"), "stfamily")
  expect_s3_class(vbinomial(const = 0.1, link = "softclipping"), "stfamily")
  expect_s3_class(vquasibinomial(const = 0.1, link = "softclipping"), "stfamily")
  expect_s3_class(vbinomial(const = 0.1, link = "softclipping"), "stfamily")
  expect_s3_class(vquasibinomial(const = 0.1, link = "softclipping"), "stfamily")
  expect_s3_class(vbinomial(link = "softclipping"), "stfamily")
  expect_s3_class(vquasibinomial(link = "softclipping"), "stfamily")
  expect_s3_class(vbinomial(link = "softclipping"), "stfamily")
  expect_s3_class(vquasibinomial(link = "softclipping"), "stfamily")
  expect_s3_class(vbinomial(const = 0, link = "identity"), "stfamily")
  expect_s3_class(vquasibinomial(const = 0, link = "identity"), "stfamily")
  expect_s3_class(vbinomial(const = 0, link = "logit"), "stfamily")
  expect_s3_class(vquasibinomial(const = 0, link = "logit"), "stfamily")
})

# test that copula object can be created
test_that("copula specification is valid", {
  expect_error(vpoisson(copula = "invalid_copula"))
  expect_error(vquasipoisson(copula = "invalid_copula"))
  expect_error(vnegative.binomial(copula = "invalid_copula"))
  expect_error(vbinomial(copula = "invalid_copula"))
  expect_error(vquasibinomial(copula = "invalid_copula"))
  expect_error(vgamma(copula = "invalid_copula"))
  expect_error(vinverse.gaussian(copula = "invalid_copula"))
  expect_error(vnormal(copula = "invalid_copula"))

  expect_error(vpoisson(copula = "normal"))
  expect_error(vquasipoisson(copula = "normal"))
  expect_error(vnegative.binomial(copula = "normal"))
  expect_error(vbinomial(copula = "normal"))
  expect_error(vquasibinomial(copula = "normal"))
  expect_error(vgamma(copula = "normal"))
  expect_error(vinverse.gaussian(copula = "normal"))
  expect_error(vnormal(copula = "normal"))

  expect_error(vpoisson(copula = "normal", copula_param = NA))
  expect_error(vquasipoisson(copula = "normal", copula_param = NA))
  expect_error(vnegative.binomial(copula = "normal", copula_param = NA))
  expect_error(vbinomial(copula = "normal", copula_param = NA))
  expect_error(vquasibinomial(copula = "normal", copula_param = NA))
  expect_error(vgamma(copula = "normal", copula_param = NA))
  expect_error(vinverse.gaussian(copula = "normal", copula_param = NA))
  expect_error(vnormal(copula = "normal", copula_param = NA))
  
  expect_error(vpoisson(copula = "normal", copula_param = -1))
  expect_error(vquasipoisson(copula = "normal", copula_param = -1))
  expect_error(vnegative.binomial(copula = "normal", copula_param = -1))
  expect_error(vbinomial(copula = "normal", copula_param = -1))
  expect_error(vquasibinomial(copula = "normal", copula_param = -1))
  expect_error(vgamma(copula = "normal", copula_param = -1))
  expect_error(vinverse.gaussian(copula = "normal", copula_param = -1))
  expect_error(vnormal(copula = "normal", copula_param = -1))

  expect_error(vpoisson(copula = "normal", copula_param = 1:2))
  expect_error(vquasipoisson(copula = "normal", copula_param = 1:2))
  expect_error(vnegative.binomial(copula = "normal", copula_param = 1:2))
  expect_error(vbinomial(copula = "normal", copula_param = 1:2))
  expect_error(vquasibinomial(copula = "normal", copula_param = 1:2))
  expect_error(vgamma(copula = "normal", copula_param = 1:2))
  expect_error(vinverse.gaussian(copula = "normal", copula_param = 1:2))
  expect_error(vnormal(copula = "normal", copula_param = 1:2))


  # valid case
  expect_s3_class(vpoisson(copula = "normal", copula_param = 1), "stfamily")
  expect_s3_class(vquasipoisson(copula = "normal", copula_param = 1), "stfamily")
  expect_s3_class(vnegative.binomial(copula = "normal", copula_param = 1), "stfamily")
  expect_s3_class(vbinomial(copula = "normal", copula_param = 1), "stfamily")
  expect_s3_class(vquasibinomial(copula = "normal", copula_param = 1), "stfamily")
  expect_s3_class(vgamma(copula = "normal", copula_param = 1), "stfamily")
  expect_s3_class(vinverse.gaussian(copula = "normal", copula_param = 1), "stfamily")
  expect_s3_class(vnormal(copula = "normal", copula_param = 1), "stfamily")
  expect_s3_class(vpoisson(copula = "t", copula_param = 1), "stfamily")
  expect_s3_class(vquasipoisson(copula = "t", copula_param = 1), "stfamily")
  expect_s3_class(vnegative.binomial(copula = "t", copula_param = 1), "stfamily")
  expect_s3_class(vbinomial(copula = "t", copula_param = 1), "stfamily")
  expect_s3_class(vquasibinomial(copula = "t", copula_param = 1), "stfamily")
  expect_s3_class(vgamma(copula = "t", copula_param = 1), "stfamily")
  expect_s3_class(vinverse.gaussian(copula = "t", copula_param = 1), "stfamily")
  expect_s3_class(vnormal(copula = "t", copula_param = 1), "stfamily")
  expect_s3_class(vpoisson(copula = "clayton", copula_param = 1), "stfamily")
  expect_s3_class(vquasipoisson(copula = "clayton", copula_param = 1), "stfamily")
  expect_s3_class(vnegative.binomial(copula = "clayton", copula_param = 1), "stfamily")
  expect_s3_class(vbinomial(copula = "clayton", copula_param = 1), "stfamily")
  expect_s3_class(vquasibinomial(copula = "clayton", copula_param = 1), "stfamily")
  expect_s3_class(vgamma(copula = "clayton", copula_param = 1), "stfamily")
  expect_s3_class(vinverse.gaussian(copula = "clayton", copula_param = 1), "stfamily")
  expect_s3_class(vnormal(copula = "clayton", copula_param = 1), "stfamily")
  expect_s3_class(vpoisson(copula = "frank", copula_param = 1), "stfamily")
  expect_s3_class(vquasipoisson(copula = "frank", copula_param = 1), "stfamily")
  expect_s3_class(vnegative.binomial(copula = "frank", copula_param = 1), "stfamily")
  expect_s3_class(vbinomial(copula = "frank", copula_param = 1), "stfamily")
  expect_s3_class(vquasibinomial(copula = "frank", copula_param = 1), "stfamily")
  expect_s3_class(vgamma(copula = "frank", copula_param = 1), "stfamily")
  expect_s3_class(vinverse.gaussian(copula = "frank", copula_param = 1), "stfamily")
  expect_s3_class(vnormal(copula = "frank", copula_param = 1), "stfamily")
  expect_s3_class(vpoisson(copula = "gumbel", copula_param = 1), "stfamily")
  expect_s3_class(vquasipoisson(copula = "gumbel", copula_param = 1), "stfamily")
  expect_s3_class(vnegative.binomial(copula = "gumbel", copula_param = 1), "stfamily")
  expect_s3_class(vbinomial(copula = "gumbel", copula_param = 1), "stfamily")
  expect_s3_class(vquasibinomial(copula = "gumbel", copula_param = 1), "stfamily")
  expect_s3_class(vgamma(copula = "gumbel", copula_param = 1), "stfamily")
  expect_s3_class(vinverse.gaussian(copula = "gumbel", copula_param = 1), "stfamily")
  expect_s3_class(vnormal(copula = "gumbel", copula_param = 1), "stfamily")
  expect_s3_class(vpoisson(copula = "joe", copula_param = 1), "stfamily")
  expect_s3_class(vquasipoisson(copula = "joe", copula_param = 1), "stfamily")
  expect_s3_class(vnegative.binomial(copula = "joe", copula_param = 1), "stfamily")
  expect_s3_class(vbinomial(copula = "joe", copula_param = 1), "stfamily")
  expect_s3_class(vquasibinomial(copula = "joe", copula_param = 1), "stfamily")
  expect_s3_class(vgamma(copula = "joe", copula_param = 1), "stfamily")
  expect_s3_class(vinverse.gaussian(copula = "joe", copula_param = 1), "stfamily")
  expect_s3_class(vnormal(copula = "joe", copula_param = 1), "stfamily")
  expect_s3_class(vpoisson(), "stfamily")
  expect_s3_class(vquasipoisson(), "stfamily")
  expect_s3_class(vnegative.binomial(), "stfamily")
  expect_s3_class(vbinomial(), "stfamily")
  expect_s3_class(vquasibinomial(), "stfamily")
  expect_s3_class(vgamma(), "stfamily")
  expect_s3_class(vinverse.gaussian(), "stfamily")
  expect_s3_class(vnormal(), "stfamily")
})


# test that all sampling methods are valid
test_that("sampling method is valid (vpoisson)", {
  expect_error(vpoisson(sampling_method = "invalid_method"))
  # valid case
  expect_s3_class(vpoisson(sampling_method = "inversion"), "stfamily")
  expect_s3_class(vpoisson(sampling_method = "poisson_process"), "stfamily")
})
test_that("sampling method is valid (vquasipoisson)", {
  expect_error(vquasipoisson(sampling_method = "invalid_method"))
  # valid case
  expect_s3_class(vquasipoisson(sampling_method = "build_up"), "stfamily")
  expect_s3_class(vquasipoisson(sampling_method = "chop_down"), "stfamily")
  expect_s3_class(vquasipoisson(sampling_method = "branching"), "stfamily")
  expect_s3_class(vquasipoisson(sampling_method = "negbin"), "stfamily")
  # warning case
  expect_warning(out <- vquasipoisson(sampling_method = "negbin", dispersion = 0.5))
  expect_s3_class(out, "stfamily")
  expect_warning(out <- vquasipoisson(sampling_method = "branching", dispersion = 0.5))
  expect_s3_class(out, "stfamily")
})



