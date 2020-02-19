###############################################################################
context("similar solutions as glm")

test_that("logitreg coefs are similar to glm solution", {
  data <- logitregsolution::sim_data(seed = 123, q_numeric = 1)

  # ~ -1 because data$x already includes an intercept.
  # equivalent:
  # coef_glm <- coef(glm(data$y ~ data$x[, -1], family=binomial()))
  coef_glm <- with(
    data,
    coef(glm(y ~ -1 + x, family = binomial()))
  )

  coef <- with(
    data,
    logitreg(response = y, design = x)
  )$coefficients

  expect_equivalent(
    signif(coef, 3),
    signif(coef_glm, 3)
  )
})

###############################################################################
context("data-trouble examples")

test_that("logitreg can deal with NAs", {
  # data needed for tests should be put in inst/testdata
  # we can load the data like this:
  load(system.file("testdata", "data-trouble.Rdata", package = "logitregsolution"))
  coef_glm <- with(
    trouble1,
    coef(glm(y ~ -1 + x, family = binomial()))
  )
  coef <- with(
    trouble1,
    logitreg(response = y, design = x)$coefficients
  )

  expect_equivalent(
    signif(coef, 3),
    signif(coef_glm, 3)
  )
})

test_that("logitreg warns about convergence", {
  load(system.file("testdata", "data-trouble.Rdata", package = "logitregsolution"))
  expect_warning(
    with(
      trouble2,
      logitreg(response = y, design = x)
    ),
    "did not converge"
  )
})
