# =============================================================================
# File: gibbs_sampler.R
# Description:
#   - Main Gibbs sampler and single sweep step.
#   - Functions:
#       * gibbs_step       : One iteration of the MH-within-Gibbs sampler.
#       * gibbs_delaporte  : Full sampler with burn-in, tuning, and optional
#                            live plotting of chains.
# Notes:
#   - This is the main entry point of the project.
#   - Depends on:
#       numerical_functions.R (clamp01, etc.)
#       component_functions.R (log_post_r, update_r_mh, update_Y1, ...)
# =============================================================================


# single Gibbs sweep step; see "Inference (MH-within-Gibbs)" in the report
gibbs_step <- function(X, r, p, lambda, Y1,
                       a_r = 0.01, b_r = 0.01,
                       a_p = 1,    b_p = 1,
                       a_l = 0.01, b_l = 0.01,
                       v_r = 0.3) {
  n <- length(X)
  
  # lambda | .  (Gamma, rate parametrization): eq. lambda-full
  shape_l <- sum(X - Y1) + a_l
  rate_l  <- n + b_l
  lambda  <- rgamma(1, shape = shape_l, rate = rate_l)
  lambda  <- max(lambda, 1e-12)   # <-- avoid 0
  
  # p | . (Beta): eq. p-full
  alpha_p <- n * r + a_p
  beta_p  <- sum(Y1) + b_p
  p       <- rbeta(1, alpha_p, beta_p)
  p       <- clamp01(p)           # <-- avoid 0/1
  
  # r | .  (MH, folded target): eq. mh-r (folded form)
  mh <- update_r_mh(r, Y1, p, v = v_r, a_r = a_r, b_r = b_r)
  r  <- mh$r_new
  
  # Y1 | .  (finite enumeration with stable softmax): eq. y1-logw
  Y1 <- update_Y1(X, r, p, lambda)
  
  list(r = r, p = p, lambda = lambda, Y1 = Y1, accepted_r = mh$accepted)
}

# main sampler with optional live plotting of chains
gibbs_delaporte <- function(
    X, iterations = 5000, burn = 1000,
    init = NULL,
    hyper = list(a_r = 0.01, b_r = 0.01,
                 a_p = 1,    b_p = 1,
                 a_l = 0.01, b_l = 0.01),
    tune = list(v_r = 0.3),
    seed = NULL,
    show.plot = FALSE,          # <— NEW
    plot_every = 500,           # <— how often to refresh the plot
    plot_subset_by = 10         # <— plot every k-th point to avoid heavy redraws
) {
  stopifnot(iterations >= 2, burn >= 0, burn < iterations)
  if (!is.null(seed)) set.seed(seed)
  
  X <- as.integer(X)
  n <- length(X)
  
  # init (simple, reproducible)
  if (is.null(init)) {
    r <- 1.0; p <- 0.5
    lambda <- max(mean(X) - r * (1 - p) / p, 0.1)
    Y1 <- pmin(X, rnbinom(n, size = r, prob = p))
  } else {
    r <- as.numeric(init$r); p <- as.numeric(init$p); lambda <- as.numeric(init$lambda)
    Y1 <- as.integer(init$Y1)
    if (length(Y1) != n) stop("init$Y1 length mismatch with X")
  }
  
  draws <- matrix(NA_real_, nrow = iterations, ncol = 3,
                  dimnames = list(NULL, c("r","p","lambda")))
  acc_r <- 0
  
  # small helper for plotting
  plot_chains <- function(iter_now) {
    
    idx <- seq(1, iter_now, by = max(1, as.integer(plot_subset_by)))
    mat <- cbind(lambda = draws[idx, "lambda"],
                 p      = draws[idx, "p"],
                 r      = draws[idx, "r"])
    
    # guard against NAs at very early iters
    ok <- complete.cases(mat)
    mat <- mat[ok, , drop = FALSE]
    
    if (nrow(mat) == 0) return(invisible())
    
    op <- par(mfrow = c(1,1))
    on.exit(par(op), add = TRUE)
    
    matplot(idx[ok], mat, type = "l",
            xlab = sprintf("iteration (every %d)", plot_subset_by),
            ylab = "value",
            col  = c("royalblue","darkgreen","purple"),
            lty  = c(1,2,3), lwd = 2)
    
    legend("right",
           legend = c(expression(lambda), expression(p), expression(r)),
           col = c("royalblue","darkgreen","purple"),
           lty = c(1,2,3), lwd = 2, bty = "n")
    
    title(sprintf("Chains up to iter %d  |  MH accept(r): %.1f%%",
                  iter_now, 100 * acc_r / iter_now), cex.main = 0.9)
  }
  
  for (it in seq_len(iterations)) {
    st <- gibbs_step(X, r, p, lambda, Y1,
                     a_r = hyper$a_r, b_r = hyper$b_r,
                     a_p = hyper$a_p, b_p = hyper$b_p,
                     a_l = hyper$a_l, b_l = hyper$b_l,
                     v_r = tune$v_r)
    r <- st$r; p <- st$p; lambda <- st$lambda; Y1 <- st$Y1
    acc_r <- acc_r + as.integer(st$accepted_r)
    draws[it, ] <- c(r, p, lambda)
    
    if (show.plot && (it %% plot_every == 0 || it == iterations)) {
      cat(sprintf("iteration: %d, MH accept(r): %.1f%%\n", it, 100 * acc_r / it))
      plot_chains(it)
    }
  }
  
  list(draws = as.data.frame(draws),
       accept_rate_r = acc_r / iterations,
       last_state = list(r = r, p = p, lambda = lambda, Y1 = Y1),
       burn = burn)
}


