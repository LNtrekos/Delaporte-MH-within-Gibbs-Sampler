# =============================================================================
# File: component_functions.R
# Description:
#   - Core model component functions used by the Gibbs sampler.
#   - Functions:
#       * log_post_r      : Log posterior for r (up to constant).
#       * update_r_mh     : One Metropolis–Hastings update for r.
#       * sample_Yi1      : Sample Y1_i given X_i, r, p, lambda.
#       * update_Y1       : Vectorized update for Y1 across all observations.
#       * simulate_delaporte : Generate synthetic Delaporte-distributed data.
# Notes:
#   - Implements the distribution-specific pieces of the sampler.
#   - Depends on numerical_functions.R (clamp01, softmax_log).
# =============================================================================


# log posterior for r | (Y1, p) up to constant
log_post_r <- function(r, Y1, p, a_r = 0.01, b_r = 0.01) {
  p <- clamp01(p)
  n <- length(Y1)
  sum(lgamma(Y1 + r) - lgamma(r)) + n * (r * log(p)) + a_r * log(r) - b_r * r
}

# one MH update for r with log-normal RW
update_r_mh <- function(r_curr, Y1, p, v, a_r = 0.01, b_r = 0.01) {
  
  r_prop <- rlnorm(1, meanlog = log(r_curr), sdlog = v)
  
  logpost_curr <- log_post_r(r_curr, Y1, p, a_r, b_r)
  logpost_prop <- log_post_r(r_prop, Y1, p, a_r, b_r)
  log_alpha <- logpost_prop - logpost_curr
  
  if (is.na(log_alpha)) log_alpha <- -Inf
  if (log(runif(1)) < log_alpha) {
    list(r_new = r_prop, accepted = TRUE)
  } else {
    list(r_new = r_curr, accepted = FALSE)
  }
}

# sample Y1_i | xi, r, p, lambda using finite support with stable logs
sample_Yi1 <- function(xi, r, p, lambda, samples = 1, eps = 1e-12) {
  p <- clamp01(p, eps)
  lam <- max(lambda, eps)
  ys <- 0:xi
  k  <- xi - ys                       # Poisson remainder
  li <- lgamma(ys + r) - lgamma(r) - lgamma(ys + 1) +
    r * log(p) + ys * log1p(-p) +
    ifelse(k == 0, 0, k * log(lam)) - lgamma(k + 1)
  probs <- softmax_log(li)
  sample(ys, size = samples, replace = TRUE, prob = probs)
}

update_Y1 <- function(X, r, p, lambda) {
  n <- length(X)
  out <- integer(n)
  for (j in seq_len(n)) out[j] <- sample_Yi1(X[j], r, p, lambda)
  out
}

simulate_delaporte <- function(n = 100, r = 2, p = 0.5, lambda = 3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Y1 <- rnbinom(n, size = r, prob = p)
  Y2 <- rpois(n, lambda)
  Y1 + Y2
}


