#' @title Estimating confidence intervals of the estimated model parameters of \code{bsw()}
#' @description For objects of class \code{"bsw"}, \code{confint()} estimates confidence intervals of the estimated model parameters of \code{bsw()}.
#' @param object An object of class \code{"bsw"}.
#' @param parm A specification of which model parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all model parameters are considered.
#' @param level A numeric value that indicates the level of confidence.
#' @param method A character giving the estimation method of the confidence intervals (\code{"bca"} or \code{"wald"}).
#' @param R A positive integer giving the number of bootstrap replicates.
#' @details \code{confint} provides Wald (default) and bias-corrected accelerated bootstrap confidence intervals of the estimated model parameters of \code{bsw()}.
#' @return A matrix with columns giving the lower and upper confidence limits of each estimated model parameter.
#' @aliases confint,bsw-method
#' @author Adam Bekhit, Jakob Schöpe
#' @export

setMethod(f = "confint",
          signature = "bsw",
          definition = function(object, parm, level = .95, method = "wald", R = 1000L) {
            if (!is.numeric(level)) {
              stop("\"level\" must be a numeric value")
            }
            
            else if (length(level) != 1L) {       
              stop("single numeric value for \"level\" expected")          
            }
            
            else if (!is.character(method)) {
              stop("\"method\" must be a character string")     
            }
                    
            else if (length(method) != 1L) {
              stop("single character string for \"method\" expected")
            }
            
            else if (!(method %in% c("bca", "wald"))) {
              stop("\"method\" is misspecified. Currently available confidence interval estimation procedures are: \"bca\" and \"wald\"")
            }
            
            else if (!is.integer(R)) {
              stop("\"R\" must be a positive integer")
            }
  
            else if (length(R) != 1L) {
              stop("single positive integer for \"R\" expected")
            }
  
            else if (R < 1000L) {
              stop("\"R\" must be a positive integer equal to or greater than 1000")
            }
                    
            else {
              cf <- coef(object)
              pnames <- names(cf)
              if (missing(parm)) { 
                parm <- pnames
              }
              else if (is.numeric(parm)) { 
                parm <- pnames[parm]
              }
              alpha <- (1 - level) / 2
              p <- c(alpha, 1 - alpha)                     
              ci <- array(data = NA, dim = c(length(parm), 2L), dimnames = list(parm, paste(x = format(x = 100 * p, digits = 3, scientific = FALSE, trim = TRUE), "%", sep = "")))
        
              if (method == "bca") {
                f <- function(formula, data, parm, indices) {
                  dat <- data[indices,]
                  fit <- bsw(formula = formula, data = dat)
                  return(coef(fit)[parm])
                }
                b <- boot::boot(data = object@data, statistic = f, R = R, formula = object@formula, parm = parm)
                ci[] <- matrix(unlist(lapply(1:ncol(b$t), function(i) {boot::boot.ci(b, conf = level, type = "bca", index = i)$bca[4:5]})), ncol = 2, byrow = TRUE)
              }              
              
              if (method == "wald") {
                se <- sqrt(diag(solve(hess(cf, object@y, object@x))))[parm]
                ci[] <- cf[parm] + se %o% stats::qnorm(p)
              }
              return(ci)
            }
          }
)
