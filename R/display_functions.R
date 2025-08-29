# =============================================================================
# File: display_functions.R
# Project: Delaporte MH-within-Gibbs Sampler
# Author: Leonidas Ntrekos
# Description:
#   - Posterior summarization and visualization helpers.
#   - Functions:
#       * summarize_posterior  : Compute posterior means, variances, CIs.
#       * plot_chains_final    : Plot chains after burn-in and optionally save.
# Notes:
#   - These functions are for diagnostics and reporting only.
#   - They should not alter sampler state or generate random draws.
# =============================================================================


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
