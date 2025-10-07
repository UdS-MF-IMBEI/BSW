#' @title Estimating bootstrap statistics of \code{bsw()}
#' @description \code{bootbsw()} applies nonparametric bootstrapping to an object of class \code{"bsw"}
#' and computes bias-corrected accelerated confidence intervals (BCa) for the estimated Relative Risk.
#'
#' @usage bootbsw(object, ci_level = 0.95, R = 1000L, maxit = NULL, conswitch = NULL)
#'
#' @param object An object of the class \code{"bsw"}.
#' @param ci_level A value between 0 and 1 indicating the confidence interval.
#' Provides bias-corrected accelerated bootstrap confidence intervals
#' of the original estimated model parameters of \code{bsw()}.
#' @param R A positive integer greater than or equal to 1000 giving the number of bootstrap replicates.
#' @param maxit A positive integer giving the maximum number of iterations in the \code{bsw()} algorithm.
#' If \code{NULL} (the default), the value stored in the \code{bsw} object is passed internally to \code{bootbsw()}.
#' @param conswitch Specifies how the constraint matrix is constructed:
#' \describe{
#'   \item{1 (default)}{Generates all possible combinations of minimum and maximum values for the predictors (excluding the intercept), resulting in \eqn{2^{m-1}} constraints.
#'   This formulation constrains model predictions within the observed data range, making it suitable for both risk factor identification and prediction (prognosis).}
#'   \item{0}{Uses the raw design matrix \code{x} as the constraint matrix, resulting in \eqn{n} constraints.
#'   This is primarily suitable for identifying risk factors, but not for prediction tasks, as predictions are not bounded to realistic ranges.}
#' }
#' If \code{NULL} (the default), the value stored in the \code{bsw} object is passed internally to \code{bootbsw()}.
#'
#' @return An object of class \code{"bsw_boot"}, which is a list containing:
#' \describe{
#'   \item{Call_bsw}{The original call to the \code{bsw()} function used to fit the model.}
#'   \item{Successful_Bootstraps}{The number of bootstrap replicates that were completed successfully.}
#'   \item{message}{A character string with a status message indicating how many bootstrap samples succeeded.}
#'   \item{Coefficients}{A matrix with the original estimated model parameters (Orig. Est.), the mean of the bootstrap estimates (Boot. Est.),
#' the standard error of the bootstrap estimates (Boot. SE), the difference between the bootstrap mean and the original estimate (bias),
#' the Risk Difference (equal to the estimate; RD), and the bias-corrected accelerated confidence intervals at the specified level.}
#'   \item{Bootstrap_Object}{An object of class \code{"boot"} (from the \pkg{boot} package) containing the full bootstrap output,
#'     including replicates and metadata. This can be used for further analyses or plotting.}
#' }
#'
#' @importFrom boot boot boot.ci
#' @importFrom checkmate assert_class assert_integer assert_numeric
#' @importFrom stats complete.cases sd
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' x <- rnorm(100, 50, 10)
#' y <- rbinom(100, 1, exp(-4 + x * 0.04))
#' fit <- bsw(formula = y ~ x, data = data.frame(y = y, x = x))
#' result <- bootbsw(fit, ci_level = 0.90)
#' print(result)
#' }
#'
#' @author Julius Johannes Weise, Thomas Wolf, Stefan Wagenpfeil
#' @export

#Definition of the new function bootbsw
bootbsw <- function(object, ci_level = 0.95, R = 1000L, maxit = NULL, conswitch = NULL) {

  #Input validation
  checkmate::assert_class(object, "bsw")
  checkmate::assert_numeric(ci_level, lower = 0, upper = 1, any.missing = FALSE, len = 1)
  if (ci_level <= 0 || ci_level >= 1) {
    stop("ci_level must be strictly between 0 and 1 (exclusive)")
  }
  checkmate::assert_integer(R, lower = 1000, any.missing = FALSE, len = 1)

  if (!is.null(maxit)) {
    checkmate::assert_integer(maxit, lower = 200L, any.missing = FALSE, len = 1)
  }

  if (!is.null(conswitch) && !conswitch %in% c(0, 1)) {
    stop("\"conswitch\" must be 0 or 1")
  }

  #Determine number of coefficients
  coef_len <- length(object@coefficients)

  #Import bsw model information
  call_list <- as.list(object@call)

  if (is.null(maxit)) {
    maxit <- if (!is.null(call_list$maxit)) eval(call_list$maxit) else 200L
  }
  
  if (is.null(conswitch)) {
    conswitch <- if (!is.null(call_list$conswitch)) eval(call_list$conswitch) else 1
  }

  #Bootstrap help function
  f <- function(data, indices) {
    dat <- data[indices, ]
    model <- tryCatch(bsw(formula = object@formula, data = dat, maxit = maxit,
                          conswitch = conswitch), error = function(e) NULL)
    if (is.null(model)) return(rep(NA, coef_len))
    return(coef(model))
  }

  #Perform bootstrap
  fit <- boot::boot(data = object@data, statistic = f, R = R)

  #Calculation of confidence intervals for each bootstrap result
  ci <- lapply(seq_len(ncol(fit$t)), function(i) {
  ci_result <- tryCatch(boot::boot.ci(fit, conf = ci_level, type = "bca", index = i), error = function(e) NULL)
  if (is.null(ci_result) || is.null(ci_result$bca)) {
    return(c(NA, NA))
  } else {
    return(ci_result$bca[4:5]) #Extract the lower (4) and upper (5) confidence interval
  }
  })
  
  #Calculate the exponent of the confidence intervals
  ci_exp <- lapply(ci, function(interval) {
    exp(interval)
  })

  #List of confidence intervals in a matrix format
  ci_exp_matrix <- do.call(rbind, ci_exp)

  #Extract the percentage level
  ci_level_lower <- (1 - ci_level) / 2 * 100  #Lower percentage level
  ci_level_upper <- 100 - ci_level_lower      #Upper percentage level

  #Dynamic naming of the confidence interval columns
  colnames(ci_exp_matrix) <- c(
    paste(format(ci_level_lower, digits = 3, scientific = FALSE, trim = TRUE), "%", sep = ""),
    paste(format(ci_level_upper, digits = 3, scientific = FALSE, trim = TRUE), "%", sep = "")
  )

  #Calculation of the bootstrap results
  bootstrap_replicates <- fit$t
  successful_bootstraps <- sum(complete.cases(bootstrap_replicates))

  #Calculation of mean values and standard errors
  original_Est <- fit$t0
  bootstrap_means <- colMeans(bootstrap_replicates, na.rm = TRUE)
  bias <- bootstrap_means - original_Est
  bootstrap_se <- apply(bootstrap_replicates, 2, sd, na.rm = TRUE)

  #Output success message
  message_text <- if (successful_bootstraps == R) {
    sprintf("%d of %d bootstraps completed successfully.", R, R)
  } else {
    sprintf("ATTENTION: Only %d out of %d bootstraps were completed successfully.
            This may be due to a restriction in one or more variables in the model.",
            successful_bootstraps, R)
  }

  #Create table
  coef.table <- cbind(as.matrix(original_Est), as.matrix(bootstrap_means), as.matrix(bootstrap_se), as.matrix(bias), as.matrix(exp(bootstrap_means)), ci_exp_matrix)
  colnames(coef.table) <- c("Orig. Est.", "Boot. Est.", "Boot. SE", "Bias","Boot. RR", colnames(ci_exp_matrix))


  result <- list(
    Call_bsw = object@call,
    Successful_Bootstraps = successful_bootstraps,
    message = message_text,
    Coefficients = coef.table,
    Bootstrap_Object = fit
  )

  class(result) <- "bsw_boot"

  return(result)
}
