load((system.file("testdata", "logitreg-data-trouble.Rdata",
  package = "logitreg")))

test_that("NA omitting", {
  expect_equivalent(fit_logitreg(trouble1$x, response = trouble1$y)$data,
    na.omit(cbind(trouble1$x, trouble1$y)))
})

test_that("fail if linearly separable", {
  expect_error(suppressWarnings(fit_logitreg(trouble2$x,
    response = trouble2$y)))
})

test_that("warn if not converged", {
  expect_warning(fit_logitreg(trouble2$x, response = trouble2$y,
    check_separable = FALSE),
    regexp = "Numerical optimization did not converge.")
})
