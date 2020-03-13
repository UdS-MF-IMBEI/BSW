#' @title Setting the linear inequality constraints for \code{bsw()}
#' @description \code{constr()} sets the linear inequality constraints for \code{bsw()}.
#' @usage constr(x)
#' @param x A model matrix.
#' @return A matrix containing the linear inequality constraints for \code{bsw()}.
#' @author Adam Bekhit, Jakob Schöpe
#' @export

constr <- function(x) {
  colMax <- matrixStats::colMaxs(x)
  colMin <- matrixStats::colMins(x)
  const <- expand.grid(lapply(2:length(colMax), function(i) {c(colMin[i], colMax[i])}))
  Amat <- unname(as.matrix(cbind(rep(1, times = nrow(const)), const)))
  return(-Amat)
}
