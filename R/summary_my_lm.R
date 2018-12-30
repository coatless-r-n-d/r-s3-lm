#' Compute Inference for Custom Linear Regression
#' 
#' Calculates standard error, t-test values, and p-values for each parameter
#' estimate.
#'
#' @param object An object with the class of `my_lm`.
#' @param ...    Additional options.
#' 
#' @return
#'  An object of class `summary.my_lm` that contains:
#'  
#'  - `est_mat`: A `matrix` with the following **named** columns...
#'      - `Estimate`: The parameter (beta) estimated values
#'      - `Std. Err`: The standard error values for the parameter estimates.
#'      - `t value`: The t-test value for the parameter estimate.
#'      - `Pr(>|t|)`: The p-value for the t-test value on the parameter.
#'  - `call`: Information on how the [my_lm()] function was called.
#'
#' @examples 
#' # Compute the model
#' model = my_lm(mpg ~ disp, data = mtcars)
#' 
#' # Calculate the inference values
#' model_inference = summary(model)
#' 
#' # Display inference values
#' model_inference
#' 
#' @seealso [my_lm()], [print.summary.my_lm()]
#' @rdname summary_lm
#' @importFrom stats pt
#' @export
summary.my_lm = function(object, ...) {
  # Note that summary(object, ...) instead of summary(x, ...)!
  
  # Calculate inference values for the model
  beta_hat = object$coefficients                 # Estimates of the parameters
  sterr = sqrt(diag(object$cov_mat))             # Standard error of Beta hat
  t_value = beta_hat / sterr                     # t-Test value
  p_value = 2 * pt(-abs(t_value), df = object$df) # p-value
  
  # Make output matrix
  est_mat = cbind(
    "Estimate" = beta_hat,
    "Std. Err" = sterr,
    "t value" = t_value,
    "Pr(>|t|)" = p_value
  )
  
  # Name the variables of the matrix
  rownames(est_mat) = rownames(object$cov_mat)
  
  # Create the object
  output = structure(list(est_mat = est_mat, call = object$call),
                     class = "summary.my_lm")
  
  output
}

#' @param x An object with the class of `summary.my_lm`.
#' @rdname summary_lm
#' @importFrom stats printCoefmat
#' @export
print.summary.my_lm = function(x, ...) {
  # Note that print(x,...)!!
  
  # Display function call
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  # Create summary matrix
  printCoefmat(x$est_mat,
               P.values   = TRUE,   # Format P-value in last column
               has.Pvalue = TRUE)   # Add a significance star column
  
  # Exercise: Add in the model fit information
}
