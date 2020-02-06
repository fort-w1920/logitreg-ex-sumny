test_that("formula works", {
  dat_frame <- sim_data(seed = 123, dataframe = TRUE)
  dat_list <- sim_data(seed = 123)
  glm_default <- fit_logitreg(dat_list$design, response = dat_list$response)
  glm_formula <- fit_logitreg(response ~ X1 + X2 + X3, data = dat_frame)
  expect_identical(glm_default$coefficients, glm_formula$coefficients)
})
