# Delaporte MH-within-Gibbs 
# -----------------------------------------------------------------------------
#   - lambda | .  ~ Gamma(sum(X - Y1) + a_l,  n + b_l)
#   - p      | .  ~ Beta(n*r + a_p,          sum(Y1) + b_p)
#   - r      | .  via one MH step with log target:
#       log pi_folded(r|.) = sum[lgamma(Y1+r) - lgamma(r)] + n*r*log p + a_r*log r - b_r*r
#   - Y1_i   | .  by finite enumeration with stable log-weights
# -----------------------------------------------------------------------------

# Numerical Functions:
{
  clamp01 <- function(p, eps = 1e-12) pmin(pmax(p, eps), 1 - eps)
  
  # numerically-stable softmax for log-weights
  softmax_log <- function(l) {
    m <- max(l)
    if (!is.finite(m)) return(rep(1/length(l), length(l)))  # all -Inf -> uniform
    w <- exp(l - m); w / sum(w)
  }
}

# Components
{
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
  
  
}

# Gibbs sampler:
{
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
  
  
}

# Display functions:
{
  summarize_posterior <- function(draws_df, burn = 0) {
    D <- draws_df[(burn + 1):nrow(draws_df), , drop = FALSE]
    qfun <- function(x) as.numeric(quantile(x, c(0.025, 0.975)))
    data.frame(
      Parameter = colnames(D),
      Mean = vapply(D, mean, numeric(1)),
      Variance = vapply(D, var, numeric(1)),
      CI_2.5 = vapply(D, function(x) qfun(x)[1], numeric(1)),
      CI_97.5 = vapply(D, function(x) qfun(x)[2], numeric(1)),
      row.names = NULL
    )
  }
  
  plot_chains_final <- function(fit, subset_by = 10, burn = fit$burn, save = NULL) {
    stopifnot(is.data.frame(fit$draws), burn >= 0, burn < nrow(fit$draws))
    idx <- seq(burn + 1, nrow(fit$draws), by = max(1L, as.integer(subset_by)))
    Y   <- as.matrix(fit$draws[idx, c("lambda", "p", "r")])  # order matches legend
    
    op <- par(mfrow = c(1,1)); on.exit(par(op), add = TRUE)
    matplot(idx, Y, type = "l",
            xlab = sprintf("iteration (every %d)", subset_by),
            ylab = "value",
            col  = c("royalblue","darkgreen","purple"),
            lty  = c(1,2,3), lwd = 2)
    legend("topleft",
           legend = c(expression(lambda), expression(p), expression(r)),
           col = c("royalblue","darkgreen","purple"),
           lty = c(1,2,3), lwd = 2, bty = "n")
    title(sprintf("Chains after burn-in (start = %d)", burn), cex.main = 0.9)
    
    if (!is.null(save)) {  # optional: save to PDF
      dev.copy(pdf, save, width = 7, height = 4.5); dev.off()
    }
  }
  
  
}


# ---------------------- example usage ----------------------
if (sys.nframe() == 0) {
  X <- simulate_delaporte(n = 100, r = 2, p = 0.5, lambda = 3, seed = 69)
  hist(X, freq = FALSE, col = "royalblue", 
       main = "Histogram of Dellaporte Disturbution")
  lines(density(X), col = "blue", lwd = 2)
  
  fit <- gibbs_delaporte(
    X, iterations = 100000, burn = 10000,
    hyper = list(a_r = 0.01, b_r = 0.01, a_p = 1, b_p = 1, a_l = 0.01, b_l = 0.01),
    tune  = list(v_r = 0.2), seed = NULL,
    show.plot = TRUE,      # turn on live plotting
    plot_every = 1000,      # refresh every 500 iterations
    plot_subset_by = 50    # plot every 10th point to keep it light
  )
  cat(sprintf("MH accept (r): %.1f%%\n", 100 * fit$accept_rate_r))
  posterior_summaries <- summarize_posterior(fit$draws, burn = fit$burn)
  print(posterior_summaries %>%
          mutate(across(where(is.numeric), ~round(., 3))))
  fit$burn = 20000
  plot_chains_final(fit, subset_by = 1)
  
  
  # Plots: 
  # Set up 3x2 plot layout
  par(mfrow = c(3, 2), mar = c(4, 4, 3, 2))  # Adjust margins
  
  # ====================== r ======================
  # Trace plot of r
  plot(fit$draws$r, type = "l", main = bquote("Trace of " * r),  
       col = "royalblue", xlab = "Iteration (t)", ylab = bquote(r^{(t)}))
  
  # Histogram + Density + Posterior summaries of r
  hist(fit$draws$r, col = "lightblue", freq = FALSE,  
       ylab = "Density", xlab = bquote(r^{(t)}),  
       main = bquote("Posterior of " * r))
  lines(density(fit$draws$r), lwd = 2, col = "blue")
  abline(v = posterior_summaries$Mean[1], col = "red", lwd = 2)
  abline(v = posterior_summaries$CI_2.5[1], col = "darkred", lty = 2)
  abline(v = posterior_summaries$CI_97.5[1], col = "darkred", lty = 2)
  legend("topright", 
         legend = c(paste("Mean =", round(posterior_summaries$Mean[1], 3)),
                    paste("95% CI = (", round(posterior_summaries$CI_2.5[1], 3), 
                          ",", round(posterior_summaries$CI_97.5[1], 3), ")")),
         bty = "n", cex = 0.8)
  
  # ====================== p ======================
  # Trace plot of p
  plot(fit$draws$p, type = "l", main = bquote("Trace of " * p),  
       col = "darkgreen", xlab = "Iteration (t)", ylab = bquote(p^{(t)}))
  
  # Histogram + Density + Posterior summaries of p
  hist(fit$draws$p, col = "lightgreen", freq = FALSE,  
       ylab = "Density", xlab = bquote(p^{(t)}),  
       main = bquote("Posterior of " * p))
  lines(density(fit$draws$p), lwd = 2, col = "green")
  abline(v = posterior_summaries$Mean[2], col = "red", lwd = 2)
  abline(v = posterior_summaries$CI_2.5[2], col = "darkred", lty = 2)
  abline(v = posterior_summaries$CI_97.5[2], col = "darkred", lty = 2)
  legend("topright", 
         legend = c(paste("Mean =", round(posterior_summaries$Mean[2], 3)),
                    paste("95% CI = (", round(posterior_summaries$CI_2.5[2], 3), 
                          ",", round(posterior_summaries$CI_97.5[2], 3), ")")),
         bty = "n", cex = 0.8)
  
  # ====================== lambda ======================
  # Trace plot of lambda
  plot(fit$draws$lambda, type = "l", main = bquote("Trace of " * lambda),  
       col = "purple", xlab = "Iteration (t)", ylab = bquote(lambda^{(t)}))
  
  # Histogram + Density + Posterior summaries of lambda
  hist(fit$draws$lambda, col = "lavender", freq = FALSE,  
       ylab = "Density", xlab = bquote(lambda^{(t)}),  
       main = bquote("Posterior of " * lambda))
  lines(fit$draws$lambda, lwd = 2, col = "purple")
  abline(v = posterior_summaries$Mean[3], col = "red", lwd = 2)
  abline(v = posterior_summaries$CI_2.5[3], col = "darkred", lty = 2)
  abline(v = posterior_summaries$CI_97.5[3], col = "darkred", lty = 2)
  legend("topright", 
         legend = c(paste("Mean =", round(posterior_summaries$Mean[3], 3)),
                    paste("95% CI = (", round(posterior_summaries$CI_2.5[3], 3), 
                          ",", round(posterior_summaries$CI_97.5[3], 3), ")")),
         bty = "n", cex = 0.8)
  
  # Reset plot layout
  par(mfrow = c(1, 1))
  }

