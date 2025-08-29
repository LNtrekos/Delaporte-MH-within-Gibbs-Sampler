# =============================================================================
# File: numerical_functions.R
# Project: Delaporte MH-within-Gibbs Sampler
# Author: Leonidas Ntrekos
# Description:
#   - Small, dependency-free numerical utilities used across the sampler.
#   - Functions:
#       * clamp01: Clamp probabilities to (eps, 1-eps).
#       * softmax_log: Stable softmax for log-weights.
#       * (optionally others: logsumexp, clamp_pos, safe_log).
# Notes:
#   - These are low-level helpers; no model logic here.
# =============================================================================

clamp01 <- function(p, eps = 1e-12) pmin(pmax(p, eps), 1 - eps)
  

softmax_log <- function(l) {
  m <- max(l)
  if (!is.finite(m)) return(rep(1/length(l), length(l)))  # all -Inf -> uniform
  w <- exp(l - m); w / sum(w)
}



