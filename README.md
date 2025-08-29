# Delaporte MH-within-Gibbs Sampler

A Bayesian inference implementation for the **Delaporte distribution** using a **Metropolis-Hastings within Gibbs sampling** approach in R.


## Project Structure

```
├── R/
│   ├── numerical_functions.R    # Low-level numerical utilities
│   ├── component_functions.R    # Distribution-specific functions
│   ├── gibbs_sampler.R         # Main sampler implementation
│   └── display_functions.R     # Posterior analysis and plotting
├── examples/
│   ├── example.R               # Standalone demonstration script
│   └── Code_all_together.R     # Self-contained version
├── LICENSE                     # MIT License
├── .gitignore                 # R-specific ignore patterns
└── README.md                  # This file
```

## Installation & Usage

### Prerequisites

```r
# Base R (≥ 4.0.0) with standard libraries
# No additional packages required
```

### Quick Start

```r
# Source the required functions
source("R/numerical_functions.R")
source("R/component_functions.R") 
source("R/gibbs_sampler.R")
source("R/display_functions.R")

# Or run the complete example
source("examples/example.R")
```

### Basic Example

```r
# Simulate synthetic data
X <- simulate_delaporte(n = 100, r = 2.0, p = 0.6, lambda = 1.5)

# Run the Gibbs sampler
fit <- gibbs_delaporte(
  X, 
  iterations = 10000, 
  burn = 2000,
  seed = 42
)

# Examine results
cat("Acceptance rate (r):", round(100 * fit$accept_rate_r, 1), "%\n")
posterior_summary <- summarize_posterior(fit$draws, burn = fit$burn)
print(posterior_summary)

# Plot diagnostics
plot_chains_final(fit)
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
