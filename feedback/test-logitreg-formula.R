###############################################################################
context("formula interface")

test_that("glm and logitreg agree", {
  data <- logitregsolution::sim_data(seed = 123, q_numeric = 1, dataframe = TRUE)
  coef_glm <- coef(glm(y ~ ., data = data, family = binomial()))
  coef <- logitreg(y ~ ., data = data)$coefficients
  expect_equivalent(signif(coef, 3), signif(coef_glm, 3))
})

test_that("works for NAs", {
  data <- logitregsolution::sim_data(seed = 123, q_numeric = 1, dataframe = TRUE)
  data[1, 2] <- NA
  coef_glm <- coef(glm(y ~ ., data = data, family = binomial()))
  coef <- logitreg(y ~ ., data = data)$coefficients
  expect_equivalent(signif(coef, 3), signif(coef_glm, 3))
})

test_that("works for factor responses", {
  data <- logitregsolution::sim_data(seed = 123, q_numeric = 1, dataframe = TRUE)
  data$y <- factor(data$y)
  coef_glm <- coef(glm(y ~ ., data = data, family = binomial()))
  coef <- logitreg(y ~ ., data = data)$coefficients
  expect_equivalent(signif(coef, 3), signif(coef_glm, 3))
})


test_that("works for formula variables in the global workspace", {
  y <- rbinom(n = 20, size = 1, p = .5)
  x <- runif(20)
  coef_glm <- coef(glm(y ~ x, family = binomial()))
  coef <- logitreg(y ~ x)$coefficients
  expect_equivalent(signif(coef, 3), signif(coef_glm, 3))
})
