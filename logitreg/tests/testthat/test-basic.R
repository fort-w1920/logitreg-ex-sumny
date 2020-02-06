test_that("basic fitting works", {
  dat_frame <- sim_data(seed = 123, dataframe = TRUE)
  dat_list <- sim_data(seed = 123)
  glm0 <- glm(response ~ X1 + X2 + X3, family = binomial, data = dat_frame)
  glmx <- fit_logitreg(dat_list$design, response = dat_list$response)
  expect_equivalent(signif(glm0$coefficients, digits = 2L),
    (signif(glmx$coefficients, digits = 2L)))
})
