#' @title Fitting a log-binomial model using the Bekhit-Schöpe-Wagenpfeil (BSW) algorithm 
#' @description \code{bsw()} fits a log-binomial model using a modified Newton-type algorithm (BSW algorithm) for solving the maximum likelihood estimation problem under linear inequality constraints.
#' @usage bsw(formula, data, maxit = 200L)
#' @param formula An object of class \code{"formula"} (or one that can be coerced to that class): a symbolic description of the model to be fitted. 
#' @param data A data frame containing the variables in the model.
#' @param maxit A positive integer giving the maximum number of iterations. 
#' @return An object of S4 class \code{"bsw"} containing the following slots:
#' \item{call}{An object of class \code{"call"}.}
#' \item{formula}{An object of class \code{"formula"}.}
#' \item{coefficients}{A numeric vector containing the estimated model parameters.}
#' \item{iter}{A positive integer indicating the number of iterations.}
#' \item{converged}{A logical constant that indicates whether the model has converged.}
#' \item{y}{A numerical vector containing the dependent variable of the model.}
#' \item{x}{The model matrix.}
#' \item{data}{A data frame containing the variables in the model.}
#' @references Wagenpfeil S (1996) Dynamische Modelle zur Ereignisanalyse. Herbert Utz Verlag Wissenschaft, Munich, Germany
#'
#' Wagenpfeil S (1991) Implementierung eines SQP-Verfahrens mit dem Algorithmus von Ritter und Best. Diplomarbeit, TUM, Munich, Germany 
#' @author Adam Bekhit, Jakob Schöpe
#' @examples
#' set.seed(123)
#' x <- rnorm(100, 50, 10)
#' y <- rbinom(100, 1, exp(-4 + x * 0.04)) 
#' fit <- bsw(formula = y ~ x, data = data.frame(y = y, x = x))
#' summary(fit)
#' @export

bsw <- function(formula, data, maxit = 200L) {
  call <- match.call()
  if (!inherits(x = formula, what = "formula")) {
    stop("\"formula\" must be of class \"formula\"")
  }
  
  else if (!is.integer(maxit)) {
    stop("\"maxit\" must be a positive integer")
  }
  
  else if (length(maxit) != 1L) {
    stop("single positive integer for \"maxit\" expected")
  }
  
  else {
    data <- stats::model.frame(formula = formula, data = data)
    y <- unname(stats::model.matrix(stats::as.formula(paste("~", all.vars(formula)[1])), data = data)[,-1])
    x <- stats::model.matrix(object = formula, data = data)
    theta <- c(log(mean(y)), rep(0, times = ncol(x) - 1))
    Amat <- constr(x)
    bvec <- rep(0, times = nrow(Amat))
    converged <- FALSE
    iter <- 0
  
    while(isFALSE(converged) & iter < maxit) {
      iter <- iter + 1
      Dmat <- Matrix::nearPD(hess(theta, y, x))$mat
      dvec <- gradF(theta, y, x) + t(theta) %*% Dmat
      fit <- quadprog::solve.QP(Dmat = Dmat, dvec = dvec, Amat = t(Amat), bvec = bvec)
      converged <- all(abs(fit$solution - theta) < 1e-4)
      theta <- fit$solution
      names(theta) <- colnames(x)
    }
    if (iter == maxit & converged == FALSE) {
      stop("Maximum number of iterations reached without convergence")
    }
    return(methods::new(Class = "bsw", call = call, formula = formula, coefficients = theta, iter = iter, converged = converged, y = y, x = x, data = data))
  }
}

