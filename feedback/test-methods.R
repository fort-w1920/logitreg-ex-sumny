###############################################################################
context("predict & fitted")

test_that("predict options work", {
  data <- logitregsolution::sim_data(seed = 1121223, n = 100)
  m <- logitreg(response = data$y, design = data$x)

  probs_pr <- predict(m, type = "response", newdata = data$x)
  link_pr <- predict(m, type = "link", newdata = data$x)
  terms_pr <- predict(m, type = "terms", newdata = data$x)

  expect_true(min(probs_pr) >= 0)
  expect_true(max(probs_pr) <= 1)
  expect_equal(probs_pr, logistic(link_pr))
  expect_equal(rowSums(terms_pr), drop(link_pr))
})

test_that("predict fallback on fitted works (and fitted works)", {
  data <- logitregsolution::sim_data(seed = 2231, n = 100)
  m <- logitreg(response = data$y, design = data$x)

  probs_pr <- predict(m, type = "response", newdata = data$x)
  link_pr <- predict(m, type = "link", newdata = data$x)
  terms_pr <- predict(m, type = "terms", newdata = data$x)

  probs_pr2 <- predict(m, type = "response")
  link_pr2 <- predict(m, type = "link")
  terms_pr2 <- predict(m, type = "terms")

  probs_f <- fitted(m, type = "response")
  link_f <- fitted(m, type = "link")
  terms_f <- fitted(m, type = "terms")

  expect_equal(probs_pr, probs_pr2)
  expect_equal(link_pr, link_pr2)
  expect_equal(terms_pr, terms_pr2)
  expect_equal(probs_pr, probs_f)
  expect_equal(link_pr, link_f)
  expect_equal(terms_pr, terms_f)
})

###############################################################################
context("plot")

test_that("plot works", {
  data <- logitregsolution::sim_data(seed = 1121223, n = 100)
  m <- logitreg(response = data$y, design = data$x)
  plot(m, type = "boxplot")
  plot(m)
})
