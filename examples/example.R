# =============================================================================
# File: example.R
# Project: Delaporte MH-within-Gibbs Sampler
# Author: Leonidas Ntrekos
# Description:
#   - Reproducible example script demonstrating sampler usage.
#   - Workflow:
#       1. Source numerical, component, sampler, and display files.
#       2. Simulate synthetic Delaporte data.
#       3. Run the Gibbs sampler (gibbs_delaporte).
#       4. Summarize posterior draws with summarize_posterior.
#       5. Plot chains with plot_chains_final (optional).
# Notes:
#   - Designed to be run as a standalone script: `Rscript example.R`.
#   - Produces console output + diagnostic plots.
# =============================================================================

# ---- Load project functions ----
source("R/numerical_functions.R")
source("R/component_functions.R")
source("R/gibbs_sampler.R")
source("R/display_functions.R")

# Simulate synthetic data
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