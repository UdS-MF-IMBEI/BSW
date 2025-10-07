#' @title Fitting a log-binomial model using the Bekhit-Schöpe-Wagenpfeil (BSW) algorithm
#' @description \code{bsw()} fits a log-binomial model using a modified Newton-type algorithm (BSW algorithm) for solving the maximum likelihood estimation problem under linear inequality constraints.
#' @usage bsw(formula, data, maxit = 200L, conswitch = 1)
#' @param formula An object of class \code{"formula"} (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data A data frame containing the variables in the model.
#' @param maxit A positive integer giving the maximum number of iterations.
#' @param conswitch Specifies how the constraint matrix is constructed:
#' \describe{
#'   \item{1 (default)}{Generates all possible combinations of minimum and maximum values for the predictors (excluding the intercept), resulting in \eqn{2^{m-1}} constraints.
#'   This formulation constrains model predictions within the observed data range, making it suitable for both risk factor identification and prediction (prognosis).}
#'   \item{0}{Uses the raw design matrix \code{x} as the constraint matrix, resulting in \eqn{n} constraints.
#'   This is primarily suitable for identifying risk factors, but not for prediction tasks, as predictions are not bounded to realistic ranges.}
#' }
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
#' fit <- bsw(formula = y ~ x, conswitch = 1, data = data.frame(y = y, x = x))
#' summary(fit)
#' @export

bsw <- function(formula, data, maxit = 200L, conswitch = 1) {
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

  else if (!is.null(conswitch) && !conswitch %in% c(0, 1)) {
    stop("\"conswitch\" must be 0 or 1")
  }

  else {
    data <- stats::model.frame(formula = formula, data = data)
    y <- unname(stats::model.matrix(stats::as.formula(paste("~", all.vars(formula)[1])), data = data)[,-1])
    x <- stats::model.matrix(object = formula, data = data)
    armijo_params <- list(c1 = 1e-4, rho = 0.5, min_alpha = 1e-12)
    theta <- c(log(mean(y)), rep(0, times = ncol(x) - 1))
    Amat <- constr(x, version = conswitch)
    bvec <- rep(0, times = nrow(Amat))
    converged <- FALSE
    iter <- 0

    while(isFALSE(converged) & iter < maxit) {
      iter <- iter + 1
      Dmat <- Matrix::nearPD(hess(theta, y, x))$mat
      grad <- gradF(theta, y, x)
      dvec <- grad + t(theta) %*% Dmat
      fit <- quadprog::solve.QP(Dmat = Dmat, dvec = dvec, Amat = t(Amat), bvec = bvec)
      p <- fit$solution - theta

      # Armijo line search
      alpha <- 1
      armijo_steps <- 0
      f_old <- obj_value(theta, y, x)
      gTp <- sum(grad * p)
      f_new <- obj_value(theta + alpha * p, y, x)

      while (f_new < f_old + armijo_params$c1 * alpha * gTp) {
        alpha <- alpha * armijo_params$rho
        armijo_steps <- armijo_steps + 1
        if (alpha < armijo_params$min_alpha) stop("Armijo line search failed: step size too small")
        f_new <- obj_value(theta + alpha * p, y, x)
      }

      theta_new <- theta + alpha * p
      converged <- all(abs(theta_new - theta) < 1e-4)
      theta <- theta_new
      names(theta) <- colnames(x)
    }
    if (iter == maxit & converged == FALSE) {
      stop("Maximum number of iterations reached without convergence")
    }
    return(methods::new(Class = "bsw", call = call, formula = formula, coefficients = theta, iter = iter, converged = converged, y = y, x = x, data = data))
  }
}
