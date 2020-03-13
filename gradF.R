#' @title Deriving the first derivatives of the log likelihood function of the log-binomial model in \code{bsw()}
#' @description \code{gradF()} derives the first derivatives of the log likelihood function of the log-binomial model.
#' @usage gradF(theta, y, x)
#' @param theta A numeric vector containing the initial values of the model parameters.
#' @param y A numeric vector containing the dependent variable of the model.
#' @param x The model matrix.
#' @return A numeric vector containing the first derivatives of the log likelihood function of the log-binomial model.
#' @author Adam Bekhit, Jakob Schöpe
#' @export

gradF <- function(theta, y, x) {
  p <- exp(x %*% theta)
  p[p >= 1] <- 1 - 1e-5
  s <- (y - p) / (1 - p)
  deriv1 <- as.vector(t(x) %*% s)
  return(deriv1)
}
