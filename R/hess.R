#' @title Deriving the second partial derivatives of the log likelihood function of the log-binomial model in \code{bsw()} (Hessian matrix)
#' @description \code{hess()} derives the second partial derivatives of the log likelihood function of the log-binomial model.
#' @usage hess(theta, y, x)
#' @param theta A numeric vector containing the initial values of the model parameters.
#' @param y A numeric vector containing the dependent variable of the model.
#' @param x The model matrix.
#' @return A numeric matrix containing the second partial derivatives of the log likelihood function of the log-binomial model (Hessian matrix).
#' @author Adam Bekhit, Jakob Schöpe
#' @export

hess <- function(theta, y, x) {
  p <- exp(x %*% theta)
  p[p >= 1] <- 1 - 1e-5
  s <- p * (y - 1) / (1 - p)^2
  im <- 0
  for(i in 1:nrow(x)){
    im <- im + x[i,] %*% t(x[i,]) * s[i]
  }
  colnames(im) <- names(theta)
  rownames(im) <- names(theta)
  return(-im)
}
