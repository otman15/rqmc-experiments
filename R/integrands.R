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