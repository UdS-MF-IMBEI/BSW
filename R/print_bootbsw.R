#' @method print bsw_boot
#' @export
print.bsw_boot <- function(x, ...) {


  cat("\nBootstrap Summary:\n")
  cat(x$message, "\n")

  cat("\nBootstrap Coefficients:\n")
  print(round(x$Coefficients, 4))

  invisible(x)
}
