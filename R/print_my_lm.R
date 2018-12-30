#' Print Model Information for Custom Linear Regression
#' 
#' Provides improved output for displaying the fitted model using
#' the custom [my_lm()] function. 
#' 
#' @param x   An object of class `my_lm`
#' @param ... Not used.
#' 
#' @examples 
#' model = my_lm(mpg ~ disp, data = mtcars)
#' 
#' # Explicit call of the print function
#' print(model)
#' 
#' # Implicit call of the print function
#' model
#' 
#' @export
print.my_lm = function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
}