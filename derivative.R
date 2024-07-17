
dGammaA <- function (x, delta, w) {
  n <- nrow (x)
  p <- ncol (x)
  h <- matrix (0, n * p, n * p)
  d <- dist (x)
  delta <- delta / sqrt (sum (w * delta ^ 2) / 2)
  dmat <- as.matrix (d)
  wmat <- as.matrix (w)
  emat <- as.matrix (delta)
  v <- -wmat
  diag(v) <- -rowSums (v)
  b <- as.matrix(-w * delta / d)
  diag(b) <- -rowSums (b)
  vinv <- solve(v + (1 / n)) - (1 / n)
  for (s in 1:p) {
    for (t in 1:p) {
      gst <- matrix (0, n, n)
      for (i in 1:n) {
        for (j in 1:n) {
          if (i == j)
            next
          gs <- x[i, s] - x[j, s]
          gt <- x[i, t] - x[j, t]
          gst[i, j] <-
            -wmat[i, j] * emat[i, j] * gs * gt / (dmat[i, j] ^ 3)
        }
      }
      diag(gst) <- -rowSums (gst)
      h[(s - 1) * n + 1:n, (t - 1) * n + 1:n] <- -vinv %*% gst
    }
  }
  for (s in 1:p) {
    kn <- (s - 1) * n + 1:n
    h[kn, kn] <- h[kn, kn] + vinv %*% b
  }
  return (h)
}

dPiA <- function (x) {
  n <- nrow (x)
  p <- ncol (x)
  sx <- svd (x)
  xu <- sx$u
  xv <- sx$v
  xd <- sx$d
  h <- matrix (0, n * p, n * p)
  for (s in 1:p) {
    for (i in 1:n) {
      ir <- (s - 1) * n + i
      e <- matrix (0, n, p)
      m <- matrix (0, p, p)
      e[i, s] <- 1
      u <- crossprod (xu, e %*% xv)
      for (k in 1:p) {
        for (l in 1:p) {
          if (k == l)
            next
          m[k, l] <-
            (xd[k] * u[k, l] + xd[l] * u[l, k]) / (xd[k] ^ 2 - xd[l] ^ 2)
        }
      }
      h[, ir] <- as.vector (e %*% xv - xu %*% diag (xd) %*% m)
    }
  }
  return (h)
}

dGammaN <- function (x, delta, w) {
  n <- nrow (x)
  p <- ncol (x)
  delta <- delta / sqrt (sum (w * delta ^ 2) / 2)
  v <- - as.matrix (w)
  diag (v) <- - rowSums(v)
  vinv <- solve (v + (1 / n)) - (1 / n)
  guttman <- function (x) {
    d <- dist (matrix (x, n, p))
    b <- - as.matrix (w * delta / d)
    diag (b) <- - rowSums (b)
    z <- vinv %*% b %*% matrix (x, n, p)
    return (as.vector (z))
  }
  return (jacobian (guttman, as.vector(x)))
}

dPiGammaN <- function (x, delta, w) {
  n <- nrow (x)
  p <- ncol (x)
  delta <- delta / sqrt (sum (w * delta ^ 2) / 2)
  v <- - as.matrix (w)
  diag (v) <- - rowSums(v)
  vinv <- solve (v + (1 / n)) - (1 / n)
  guttman <- function (x) {
    d <- dist (matrix (x, n, p))
    b <- - as.matrix (w * delta / d)
    diag (b) <- - rowSums (b)
    z <- vinv %*% b %*% matrix (x, n, p)
    z <- z %*% svd(z)$v
    return (as.vector (z))
  }
  return (jacobian (guttman, as.vector(x)))
}

dPiGammaA <- function (x, delta, w) {
  n <- nrow (x)
  guttman <- function (x, delta, w) {
    d <- dist (x)
    b <- as.matrix (- w * delta / d)
    diag (b) <- - rowSums (b)
    v <- - as.matrix (w)
    diag (v) <- - rowSums (v)
    vinv <- solve(v + (1 / n)) - (1 / n)
    return (vinv %*% b %*% x)
  }
  h <- dGammaA (x, delta, w)
  g <- dPiA (guttman (x, delta, w))
  return (g %*% h)
}

