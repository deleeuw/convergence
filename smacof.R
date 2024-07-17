
torgerson <- function (delta, p = 2) {
  h <- as.matrix(delta ^ 2)
  n <- nrow (h)
  j <- diag (n) - (1 / n)
  h <- -(j %*% h %*% j) / 2
  e <- eigen (h)
  return (e$vectors[, 1:p] %*% diag(sqrt (e$values[1:p])))
}

smacof <-
  function (delta,
            w,
            p = 2,
            xold = torgerson (delta, p),
            pca = FALSE,
            verbose = FALSE,
            eps = 1e-15,
            itmax = 10000) {
    n <- round ((1 + sqrt (1 + 8 * length (delta))) / 2)
    v <- -as.matrix(w)
    diag(v) <- -rowSums(v)
    vinv <- solve(v + (1 / n)) - (1 / n)
    delta <- delta / sqrt (sum (w * delta ^ 2) / 2)
    dold <- dist (xold)
    sold <- sum (w * (delta - dold) ^ 2) / 2
    eold <- Inf
    itel <- 1
    repeat {
      b <- as.matrix (-w * delta / dold)
      diag (b) <- -rowSums(b)
      xnew <- vinv %*% b %*% xold
      if (pca) {
        xsvd <- svd (xnew)
        xnew <- xnew %*% xsvd$v
      }
      dnew <- dist (xnew)
      snew <- sum (w * (delta - dnew) ^ 2) / 2
      enew <- sqrt (sum (v * tcrossprod(xold - xnew)))
      rnew <- enew ^ (1 / itel)
      qnew <- enew / eold
      if (verbose) {
        cat(
          "itel ",
          formatC(itel, digits = 4, format = "d"),
          "loss ",
          formatC(snew, digits = 15, format = "f"),
          "chan ",
          formatC(enew, digits = 15, format = "f"),
          "rcnf ",
          formatC(rnew, digits = 15, format = "f"),
          "qcnf ",
          formatC(qnew, digits = 15, format = "f"),
          "\n"
        )
      }
      if ((enew < eps) || (itel == itmax))
        break
      xold <- xnew
      dold <- dnew
      sold <- snew
      eold <- enew
      itel <- itel + 1
    }
    out <-
      list (
        itel = itel,
        x = xnew,
        s = snew,
        q = qnew,
        r = rnew,
        b = vinv %*% b
      )
    return (out)
  }
