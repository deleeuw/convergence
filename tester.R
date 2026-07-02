tester <- function (w, delta, p) {
  h <- smacof (delta,
               w,
               p = p,
               pca = TRUE,
               eps = 1e-15)
  x <- h$x
  n <- nrow (x)
  s <- eigen(tcrossprod (x))
  e <- s$vectors
  d <- sqrt (s$values[1:p])
  k <- svd(x)$u
  kperp <- e[, -(1:p)]
  g <- dPiA (x)
  z <- k %*% diag (1 / d)
  for (i in 1:p) {
    y <- matrix (0, n, p)
    y[, i] <- z[, i]
    y <- as.vector (y)
    print (max (abs (y - g %*% y)))
  }
  print("*********")
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      y <- matrix (0, n, p)
      y[, i] <- -z[, j]
      y[, j] <- z[, i]
      y <- as.vector (y)
      print (max (abs (y - g %*% y)))
    }
  }
  print("*********")
  for (i in 1:(n-p)) {
    for (j in 1:p) {
      y <- matrix (0, n, p)
      y[, j] <- kperp[, i]
      y <- as.vector (y)
      print (max (abs (y - g %*% y)))
    }
  }
  print("*********")
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      y <- matrix (0, n, p)
      y[, i] <- -x[, j]
      y[, j] <- x[, i]
      y <- as.vector (y)
      print (max (abs (g %*% y)))
    }
  }
}