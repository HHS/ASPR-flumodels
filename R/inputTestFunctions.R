#' @title checkNonNegativeNumber
#' @description Stops execution with an error message if the input is not a non-negative number
#' @keywords internal
checkNonNegativeNumber <- function (x) {
  if ((length(x) != 1) || (x  < 0)) {
    stop(paste(deparse(substitute(x)), "must be a non-negative number."), call. = FALSE)
  }
}

#' @title checkNonNegative
#' @description Stops execution with an error message if the input is not non-negative
#' @keywords internal
checkNonNegative <- function(x) {
  if (!all(x >= 0)) {
    stop(paste(deparse(substitute(x)), "must be non-negative."), call. = FALSE)
  }
}

#' @title checkPositiveNumber
#' @description Stops execution with an error message if the input is not a positive number
#' @keywords internal
checkPositiveNumber <- function (x) {
  if ((length(x) != 1) || (x  < 0)) {
    stop(paste(deparse(substitute(x)), "must be a positive number."), call. = FALSE)
  }
}

#' @title checkPositive
#' @description Stops execution with an error message if the input is not positive
#' @keywords internal
checkPositive <- function(x) {
  if (!all(x > 0)) {
    stop(paste(deparse(substitute(x)), "must be positive."), call. = FALSE)
  }
}

#' @title checkBetween0and1
#' @description Stops execution with an error message if the input is not between 0 and 1
#' @keywords internal
checkBetween0and1 <- function (x) {
	if(!all(x >= 0) || !all(x <= 1)){
		stop(paste(deparse(substitute(x)), "must be between 0 and 1."), call. = FALSE)
	}	
}

#' @title checkDimensionsMatch
#' @description Stops execution with an error message if the length of x is not equal to either 1 or the length of y
#' @keywords internal
checkDimensionsMatch <- function (x, y) {
	if((length(x) > 1) && (length(x) != length(y))){
		stop(paste(deparse(substitute(x)), "dimensions do not match", deparse(substitute(y))), call. = FALSE)
	}	
}