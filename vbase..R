vbase <- function(v, p) {
  n <- ncol(v)
  y <- lapply(1:p, function(s)
    matrix(rnorm((n - s + 1) * (n - s)), n - s + 1, n - s))
  y <- lapply(y, function(y)
    apply(y, 2, function(y)
      y - mean(y)))
  y <- lapply(1:p, function(s)
    rbind(matrix(0, s - 1, n - s), y[[s]]))
  for (s in 1:p) {
    yy <- y[[s]]
    for (t in 1:(n - s)) {
      if (t == 1) {
        yy[, t] <- yy[, t] / sqrt(sum(yy[, t] * (v %*% yy[, t])))
        next
      }
      for (r in 1:(t - 1)) {
        yy[, t] <- yy[, t] - sum(yy[, t] * (v %*% yy[, r])) * yy[, r]
      }
      yy[, t] <- yy[, t] / sqrt(sum(yy[, t] * (v %*% yy[, t])))
    }
    y[[s]] <- yy
  }
  return(y)
}
