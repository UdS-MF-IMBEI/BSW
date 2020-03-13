#' @title Summarizing the estimated model parameters of \code{bsw()}
#' @description For objects of class \code{"bsw"}, \code{summary()} summarizes the estimated model parameters of \code{bsw()}.
#' @param object An object of class \code{"bsw"}.
#' @return A list containing the following elements:
#' \item{coefficients}{A numeric vector containing the estimated model parameters.}
#' \item{std.err}{A numeric vector containing the estimated standard errors of the model parameters.}
#' \item{z.value}{A numeric vector containing the estimated z test statistic of the model parameters.}
#' \item{p.value}{A numeric vector containing the estimated p values of the model parameters.}
#' @aliases summary,bsw-method
#' @author Adam Bekhit, Jakob Schöpe
#' @export

setMethod(f = "summary",
          signature = "bsw",
          definition = function(object) {
            cf <- coef(object)
            ci <- confint(object)
            se <- sqrt(diag(solve(hess(cf, object@y, object@x))))
            z <- cf / se
            p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
            coef.table <- cbind(as.matrix(cf), as.matrix(se), as.matrix(z), as.matrix(p), as.matrix(exp(cf)), exp(ci))
            colnames(coef.table) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "RR", colnames(ci))
            
            cat("Call:\n")
            print(object@call)
            cat("\nConvergence:", object@converged)
            cat("\nCoefficients:\n")
            print(coef.table)
            cat("\nIterations:", object@iter, "\n")
            return(invisible(list(coefficients = cf, std.err = se, z.value = z, p.value = p)))
          }
)
