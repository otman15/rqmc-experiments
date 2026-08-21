# SumUeU:
# For each point u, returns
# sum_j u_j * exp(u_j) - s.
sumueu <- function(x) {
  rowSums(x * exp(x)) - ncol(x)
}

smoothperb4 <- function(x, expw) {
  omega <- seq_len(ncol(x))^(-expw)
  
  terms <- 1 + sweep(
    30 * x^2 * (1 - x)^2 - 1,
    2,
    omega,
    `*`
  )
  
  apply(terms, 1, prod) - 1
}


mc2 <- function(x) {
  s <- ncol(x)
  row_prods <- apply((s - x) / (s - 0.5), 1, prod)
  row_prods - 1
}


polynomial <- function(x) {
  s <- ncol(x)
  a <- seq_len(s) / s
  
  terms <- 1 + sweep(x - 0.5, 2, a, `*`)
  
  apply(terms, 1, prod) - 1
}


oscillatory <- function(x) {
  # SSJ Oscillatory model, centered to match the Genz oscillatory form.
  s <- ncol(x)
  a <- seq_len(s) / s
  mean <- cos((s + 1) / 4) * prod(2 * sin(a / 2) / a)
  
  cos(rowSums(sweep(x, 2, a, `*`))) - mean
}


gaussian <- function(x) {
  # SSJ Gaussian model with Genz peak sign convention and exact centering.
  s <- ncol(x)
  mean <- (sqrt(pi) / 2 * (2 * pnorm(sqrt(2)) - 1))^s
  
  exp(-rowSums(x^2)) - mean
}


normal_scores <- function(x) {
  # Maps uniforms to normal scores, matching SSJ NormalDist.inverseF01.
  qnorm(pmin(pmax(x, .Machine$double.xmin), 1 - .Machine$double.eps))
}


smoothgauss <- function(x) {
  # SSJ SmoothGauss benchmark: smooth centered function of normal scores.
  z <- normal_scores(x)
  s <- ncol(z)
  
  pnorm(1 + rowSums(z) / sqrt(s)) - pnorm(1 / sqrt(2))
}


piecelingauss <- function(x) {
  # SSJ PieceLinGauss benchmark: centered positive part of normal sum.
  z <- normal_scores(x)
  s <- ncol(z)
  a <- pnorm(-1) - dnorm(1)
  
  pmax(rowSums(z) / sqrt(s) - 1, 0) + a
}


indsumnormal <- function(x) {
  # SSJ IndSumNormal benchmark: centered indicator of a normal sum.
  z <- normal_scores(x)
  s <- ncol(z)
  
  as.numeric(rowSums(z) / sqrt(s) >= 1) - pnorm(-1)
}
