obj_value <- function(theta, y, x) {
  eta <- x %*% theta
  p <- exp(eta)
  p[p >= 1] <- 1 - 1e-5
  sum(y * eta + (1 - y) * log(1 - p))
}
