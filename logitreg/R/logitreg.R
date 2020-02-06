#' Generate synthetic data sets for tests
#'
#' Generates synthetic data sets for logistic regression.
#'
#' @param n number of observations, integer
#' @param numerics number of numeric variables, integer
#' @param factors number of factor variables, integer
#' @param seed seed for the RNG, numeric
#' @param dataframe should a dataframe be returned consisting of the response
#'   and the variables (dropping the intercept column and appending the coefs
#'   as an attribute) or should everything be returned as a list (design,
#'   repsonse, coefs)? logical
#' @return list of the design, response and the coefs, or a dataframe including
#'   the response and the design (minus the intercept column) with the coefs
#'   appended as an attribute
#' @importFrom stats rbinom rnorm runif
#' @md
#' @export
sim_data <- function(n = 1500L, numerics = 3L, factors = 0L, seed = NULL,
  dataframe = FALSE) {
  # set RNG seed for reproducibility
  if (!is.null(seed))
    set.seed(seed)

  covariates <- 1L + numerics + factors

  design <- matrix(0L, nrow = n, ncol = covariates)
  design[, 1L] <- 1L
  design[, seq_len(numerics) + 1L] <- matrix(rnorm(n * numerics), nrow = n)
  if (factors) {
    # add binary factors
    dummies <- matrix(sample(c(0L, 1L), n * factors, replace = TRUE), nrow = n)
    design[, -seq_len(1L + numerics)] <- dummies
  }

  coefs <- runif(covariates, min = -3L, max = 3L)

  probabilities <- logistic(design %*% coefs)
  response <- rbinom(n, prob = probabilities, size = 1L)

  if (!dataframe) {
    return(list(design = design, response = response, coefs = coefs))
  }
  structure(
      data.frame(response = response, design[, -1, drop = FALSE]),
      coefs = coefs)
}

#' Compute the logistic transformation
#'
#' Simple wrapper around [stats::plogis].
#'
#' @param x value to transform, numeric
#' @return logistic transformed value, numeric
#' @importFrom stats plogis
#' @md
#' @export
logistic <- function(x) plogis(x)

#' Compute the negative log likelihood of some data under a logit model or its
#' gradient.
#'
#' Based on a parameter vector, design matrix and a response vector, the
#' negative log likelihood under the logit model or its gradient is computed.
#'
#' @param coefs parameter vector, numeric vector
#' @param design design matrix, numeric matrix
#' @param response response vector, numeric vector
#' @return negative log likelihood, numeric vector
#' @md
#' @export
neg_loglik <- function(coefs, design, response) {
  probabilities <- logistic(design %*% coefs)
  - sum(response * log(probabilities) +
    (1L - response) * log(1L - probabilities))
}

#' @rdname neg_loglik
#' @return gradient of the negative log likelihood, numeric vector
neg_loglik_deriv <- function(coefs, design, response) {
  probabilities <- logistic(design %*% coefs)
  - t(response - probabilities) %*% design
}

#' Basic fitting function for a logit model.
#'
#' Estimate the parameter vector under a logit model providing a design matrix
#' and a response vector. Linear separability of the data is checked solving a
#' linear program via Rglpk (if `check_separable = TRUE`, default).
#'
#' @param design design matrix, numeric matrix or model formula, formula
#' @param data data.frame
#' @param response response vector, integer vector
#' @param hessian should the hessian be computed? logical
#' @param check_separable should a test of linearly separability pe performed?
#' logical
#' @param \dots additional arguments passed to the methods or `optim`'s control
#' argument
#' @return an object of class `logitreg`, i.e., a list including the estimated
#' coefficients, the fitted values, the data used for fitting, the vcov (if
#' hessian = TRUE, otherwise a matrix of NA's), the log likelihood of the
#' fitted model, the number of parameters estimated, the convergence code of
#' optim, and the number of gradient evaluations used during optimization.
#' Several basic methods are available, namely `coef`, `vcov`, `logLik`,
#' `summary`, `plot`, `predict`, and `fitted`.
#' @importFrom checkmate test_integer assert_matrix assert_integer
#' assert_logical assert_formula assert_data_frame
#' @importFrom stats na.omit optim model.frame model.matrix
#' @importFrom utils tail
#' @examples
#' data <- sim_data(seed = 123, dataframe = TRUE)
#' glmx <- fit_logitreg(response ~ X1 + X2 + X3, data = data)
#' summary(glmx)
#' @md
#' @export
fit_logitreg <- function(design = NULL, ...) {
  UseMethod("fit_logitreg", design)
}

#' @rdname fit_logitreg
#' @export
fit_logitreg.formula <- function(design = NULL, data = NULL, hessian = TRUE,
  check_separable = TRUE, ...) {
  assert_formula(design)
  assert_data_frame(data)
  model_frame <- model.frame(design, data = data, na.action = NULL)
  model_matrix <- model.matrix(design, data = data)
  fit_logitreg(design = model_matrix, response = model_frame[[1L]],
    hessian = hessian, check_separable = check_separable, ...)
}

#' @rdname fit_logitreg
#' @export
fit_logitreg.default <- function(design = NULL, response = NULL, hessian = TRUE,
  check_separable = TRUE, ...) {
  # input checks
  if(!test_integer(response, lower = 0L, upper = 1L)) {
    warning("Response vector should be an integer vector of zeros and ones.",
      " I WILL FIX THIS BY BRUTE ROUNDING.")
    response <- as.integer(round(response, digits = 0L))
  }
  assert_matrix(design, mode = "numeric", nrows = length(response),
    min.cols = 1L)
  assert_logical(hessian, any.missing = FALSE, len = 1L)
  assert_logical(check_separable, any.missing = FALSE, len = 1L)

  # handle NA's
  drop <- c(which(is.na(design), arr.ind = TRUE)[, 1L], which(is.na(response)))
  if (length(drop) > 0L) {
    design <- design[-drop, ]
    response <- response[-drop]
  }

  # check for linear separability using linear programming via Rglpk
  if (check_separable) {
    if (is_separable(design, response = response)) {
       stop("Data is perfectly linearly separable.",
         " Try lasso or ridge logit models.")
    }
  }

  # pass dots to optim's control
  control <- list(...)

  # optimize
  coefs <- rep_len(0L, length.out = NCOL(design))
  optim_results <- optim(coefs, fn = neg_loglik, gr = neg_loglik_deriv,
    method = "BFGS", control = control, hessian = hessian, design = design,
    response = response)

  if (optim_results$convergence == 1L) {
    warning("Numerical optimization did not converge.")
  }

  # results
  coefficients <- optim_results$par
  df <- length(coefficients)
  vcov <- if (hessian) {
    chol2inv(chol(optim_results$hessian))
  } else {
    matrix(NA, nrow = length(coefs), ncol = length(coefs))
  }

  # return list
  ret <- list(coefficients = coefficients,
              fitted = logistic(design %*% coefficients),
              data = cbind(design = design, response = response),
              vcov = vcov,
              loglik = - optim_results$value,
              df = df,
              code = optim_results$convergence,
              iterations = tail(na.omit(optim_results$counts), 1L))

  structure(ret, class = "logitreg")
}

#' Test for linear separability.
#'
#' Test whether data is linear separable by solving a linear program via Rglpk.
#'
#' @inheritParams fit_logitreg
#' @return a logical whether the data is linear separable
#' @importFrom Rglpk Rglpk_solve_LP
#' @md
#' @export
is_separable <- function(design, response) {
  # some indices
  N <- NROW(design)
  P <- NCOL(design) + 1L
  # map response to 1 and -1
  response_directional <- ifelse(response == 1L, yes = 1L, no = -1L)
  # left hand side
  A <- cbind(design * response_directional, - response_directional)
  # objective function
  objective <- rep_len(0L, length.out = P)
  # right hand side
  b <- rep_len(-1L, length.out = N)
  # no boundaries
  bounds <- list(
    lower = list(ind = seq_len(P), val = rep_len(-Inf, length.out = P)),
    upper = list(ind = seq_len(P), val = rep_len(Inf, length.out = P))
  )
  solution <- Rglpk_solve_LP(objective, mat = A,
    dir = rep_len("<=", length.out = N), rhs = b, bounds = bounds)
  ifelse(solution$status == 0L, yes = TRUE, no = FALSE)
}

#' @export
coef.logitreg <- function(object, ...) {
  object$coefficients
}

#' @export
vcov.logitreg <- function(object, ...) {
  object$vcov
}

#' @export
logLik.logitreg <- function(object, ...) {
  structure(object$loglik, df = object$df, class = "logLik")
}

#' @importFrom stats vcov pnorm
#' @export
summary.logitreg <- function(object, vcov. = NULL, ...) {
  # coefficients
  cf <- object$coefficients

  # covariance matrix
  if (is.null(vcov.)) 
    vc <- vcov(object)
  else {
    if (is.function(vcov.)) vc <- vcov.(object)
    else vc <- vcov.
  }
  
  # Wald test of each coefficient
  cf <- cbind(cf, sqrt(diag(vc)), cf/sqrt(diag(vc)), 2L * pnorm(-abs(cf/sqrt(diag(vc)))))
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

  object$coefficients <- cf
  class(object) <- "summary.logitreg"
  return(object)
}

#' @importFrom stats printCoefmat
#' @export
print.summary.logitreg <- function(x,
  digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"), ...) {
  if (is.null(x$call)) {
    cat("\nLogit model\n\n")
  } else {
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  }

  cat("Coefficients:\n")
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)

  cat("\nLog-likelihood:", format(signif(x$loglik, digits)),
      "(df =", paste(x$df, ")", sep = ""), "\n")
  cat("Number of iterations in numerical optimization:", x$iterations, "\n\n")
  invisible(x)
}

#' @importFrom graphics plot
#' @importFrom ROCR prediction performance
#' @export
plot.logitreg <- function(x, ...) {
  predROCR <- prediction(x$fitted, x$data[, NCOL(x$data)])
  perfROCR <- performance(predROCR, "tpr", "fpr")
  plot(perfROCR, ...)
  
}

#' @export
predict.logitreg <- function(object, design, ...) {
  cf <- object$coefficients
  assert_matrix(design, mode = "numeric", min.rows = 1L, ncols = length(cf))
  as.vector(design %*% cf)

}

#' @export
fitted.logitreg <- function(object, ...) {
  as.vector(object$fitted)
}
