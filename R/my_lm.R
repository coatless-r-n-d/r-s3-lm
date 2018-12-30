#' Custom Linear Regression Estimation
#' 
#' Provides a custom implementation of _R_'s [stats::lm()] linear regression
#' function that uses the S3 system.
#' 
#' @param x   Either a design matrix with dimensions \eqn{n \times p}{n x p} or [`formula`] with the `data` parameter specified
#' @param ... Not used.
#' 
#' @return
#'  An object of class `my_lm` that contains:
#'  
#'  - `coefficients`: Estimated parameter values of \eqn{\hat{\beta}}{beta_hat} 
#'     with dimensions \eqn{p \times 1}{p x 1}
#'  - `cov_mat`: Covariance matrix of estimated parameter values 
#'     with dimensions \eqn{p \times p}{p x p}.
#'  - `sigma`: Standard deviation of residuals
#'  - `df`: Degrees of Freedom given by \eqn{df = N - p}
#'  - `fitted.values`: Fitted Values given by \eqn{\hat{y} = X\hat{\beta}}{y_hat = X*beta_hat}
#'  - `residuals`: Residuals given by \eqn{e = y - \hat{y}}{e = y - y_hat}
#'  - `call`: Information on how the [my_lm()] function was called.
#'
#' @details
#' 
#' Given a response vector \eqn{y} with dimensions \eqn{n \times 1}{n x 1}, 
#' a design matrix \eqn{X} with dimensions \eqn{n \times p}{n x p},
#' a vector of parameters \eqn{\beta}{beta} with dimensions \eqn{p \times 1}{p x 1},
#' and an error vector \eqn{\epsilon}{e} with dimensions \eqn{n \times 1}{n x 1} from \eqn{\epsilon \sim N\left( {0,{\sigma ^2}} \right)}{e ~ N(0, sigma^2)},
#' the standard linear regression model can be stated as:
#' 
#' \deqn{y  = {X'}\beta + \epsilon}{y = X'beta + e}
#' 
#' The least ordinary squares (OLS) solutions are then:
#' 
#' \deqn{\hat \beta  = {\left( {{X'}X} \right)^{ - 1}}{X'}y}{beta_hat = (X'X)^(-1) X' y}
#' \deqn{\operatorname{cov} \left( {\hat \beta } \right) = {\sigma ^2}{\left( {{X^T}X} \right)^{ - 1}}}{cov(beta_hat) = sigma^2 * (X' X)^(-1)}
#' 
#' @examples 
#' ## Matrix interface
#' 
#' # Create a design matrix
#' x = cbind(1, mtcars$disp)
#' 
#' # Extract response
#' y = mtcars$mpg
#' 
#' # Calculate outcome
#' my_model = my_lm(x, y)
#' 
#' ## Formula interface
#' 
#' # Calculate 
#' my_model = my_lm(mpg ~ disp, data = mtcars)
#' 
#' @export
#' @seealso [summary.my_lm()], [print.my_lm()]
my_lm = function(x, ...) { 
  UseMethod("my_lm") 
}


#' @rdname my_lm
#' @export
my_lm.default = function(x, ...) {
  
  stop("Must have either a matrix or a formula. Cannot operate on ", class(x), ".")
  
}

#' @param formula An object of class [formula()] which provides a symbolic description of the model to be fitted. 
#' @param data    An optional data frame, list or environment.
#' @rdname my_lm
#' @importFrom stats model.frame model.matrix model.response
#' @export
my_lm.formula = function(formula, data = list(), ...) {
  
  # Create a model frame
  model_info = model.frame(formula = formula, data = data)
  
  # Extract the model matrix
  x = model.matrix(formula, data = model_info)
  
  # Extract response
  y = model.response(model_info)
  
  # Compute the regression
  model = compute_my_lm(x, y)
  
  # Add function call details
  model$call = match.call()
  
  # Add formula
  model$formula = formula
  
  # Return built object
  model
}

#' @param y A `vector` of responses with dimensions \eqn{n \times 1}{n x 1}.
#' @rdname my_lm
#' @export
my_lm.matrix = function(x, y, ...) {
  
  stopifnot(is.matrix(x), is.numeric(x))
  stopifnot(is.numeric(y))
  
  model = compute_my_lm(x, y)
  model$call = match.call()
  
  model
}


compute_my_lm = function(x, y, ...) {
  ## Translation of (X'X)^-1 X' y
  # beta_hat = solve(t(x) %*% x) %*% t(x) %*% y
  
  ## More stable calculation of beta_hat
  # QR-decomposition of x
  qr_x = qr(x)
  
  # Computes a more stable variant of (X'X)^-1 X'y
  beta_hat = solve.qr(qr_x, y)
  
  # Compute the Degrees of Freedom
  df = nrow(x) - ncol(x)   # n - p 
  
  # Compute the Variance of the Residuals
  sigma2 = sum((y - x %*% beta_hat) ^ 2) / df
  
  ## More stable computation of the Covariance Matrix
  
  ## Translation of Cov(beta_hat) = sigma^2 * (X' X)^(-1)
  # cov_mat = sigma2 * solve(t(x) %*% x)
  
  # More stable computation of the covariance matrix
  cov_mat = sigma2 * chol2inv(qr_x$qr)
  
  # Make name symmetric in covariance matrix
  rownames(cov_mat) = colnames(x)
  colnames(cov_mat) = colnames(x)
  
  # Compute Fitted Values
  fitted.values = as.vector(x %*% beta_hat)
  residuals = y - fitted.values
  call = match.call()
  
  # Construct object
  output = structure(list(coefficients  = beta_hat, 
                          cov_mat = cov_mat, 
                          sigma   = sqrt(sigma2), 
                          df      = df,
                          fitted.values = fitted.values,
                          residuals = residuals,
                          call = call),
                     class = c("my_lm", "list"))
  
  # Return calculation
  output
}
