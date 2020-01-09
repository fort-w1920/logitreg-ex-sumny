################################################################################

# Generate synthetic data sets for tests.  q_numeric gaussian variables (~N(0,1))
# and q_factor dummy variables.  Returns list(x, y, b): design, response, true
# coefs or data.frame with attribute('b').
sim_data <- function(n = 1500, numerics = 3, factors = 0, seed = NULL, dataframe = FALSE) {
  # set RNG seed for reproducibility:
  if (!is.null(seed))
    set.seed(seed)

  covariates <- 1 + numerics + factors

  design <- matrix(0, nrow = n, ncol = covariates)
  design[, 1] <- 1
  design[, seq_len(numerics) + 1] <- matrix(rnorm(n * numerics), nrow = n)
  if (factors) {
    # add binary factors
    dummies <- matrix(sample(c(0, 1), n * factors, replace = TRUE), nrow = n)
    design[, -seq_len(1 + numerics)] <- dummies
  }

  coefs <- runif(covariates, min = -3, max = 3)

  probabilities <- logistic(design %*% coefs)
  response <- rbinom(n, p = probabilities, size = 1)

  if (!dataframe) {
    return(list(design = design, response = response, coefs = coefs))
  }
  structure(
      data.frame(response = response, design[, -1, drop = FALSE]),
      coefs = coefs)
}



################################################################################
## Saving you from actually doing any maths:

logistic <- function(x) plogis(x)

neg_loglik <- function(coefs, design, response) {
  probabilities <- logistic(design %*% coefs)
  - sum(response * log(probabilities) + (1 - response) * log(1 - probabilities))
}

neg_loglik_deriv <- function(coefs, design, response) {
  probabilities <- logistic(design %*% coefs)
  - t(response - probabilities) %*% design
}


