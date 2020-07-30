random_init_probabilities <- function(k) {
  values <- runif(k, 1, 100)
  vsum <- sum(values)
  values/vsum
}

log_sum_exp <- function(x) {
  m <- max(x)
  if (is.infinite(m) && m < 0) return(-Inf)
  m + log(sum(exp(x - m)))
}