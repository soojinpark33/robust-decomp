# ============================================================================
# Causal Decomposition Analysis with Sensitivity Analysis
# ============================================================================
# 
# Description: This script performs causal decomposition analysis comparing
#              multiple methods (Difference in Coefficients, Oaxaca-Blinder,
#              and Causal Decomposition) with sensitivity analysis and covariate benchmarking.
#
# License: This code is provided for research purposes.
# ============================================================================

# Load required packages
library(haven)
library(causal.decomp)
library(ggplot2)
library(ggrepel)

# ============================================================================
# CONFIGURATION SECTION
# ============================================================================
# Modify these settings to match your data and analysis needs

# Data file path (relative to working directory or absolute path)
data_file <- "synthetic_data.dta"

# Variable names (modify to match your dataset)
group <- "black"                 # Group/exposure variable (should be factor)
mediator_var <- "educy"          # Mediator variable
outcome_var <- "cog27"           # Outcome variable
base_cov <- "age"                # Baseline covariate (will be centered)

# Intermediate covariates for adjustment
# Note: base_cov_centered and base_cov_centered_sq will be created automatically
intermediate_covs <- c("chdSES", "RTHLTHCH", "Pdivorce16")

# Benchmark covariates for sensitivity analysis
benchmark_covariates <- c("age_centered", "chdSES", "RTHLTHCH", "Pdivorce16")

# Bootstrap settings
n_boot <- 1000                   # Number of bootstrap iterations
conf_level <- 0.95               # Confidence level

# Sensitivity analysis settings
max_rsq <- 0.3                   # Maximum R² for sensitivity analysis grid

# SMI settings
smi_sims <- 1000                 # Number of simulations for SMI
smi_seed <- 32                    # Seed for SMI
smi_conditional <- TRUE           # Use conditional SMI

# Plot settings
plot_width <- 10                 # Plot width in inches
plot_height <- 8                 # Plot height in inches
plot_dpi <- 300                  # Plot resolution

# Output file names
output_plot_reduction <- "sensitivity_reduction.png"

# ============================================================================
# DATA LOADING AND PREPARATION
# ============================================================================

cat("\n========================================\n")
cat("Causal Decomposition Analysis\n")
cat("========================================\n\n")

# Read the Stata file
cat("Loading data from:", data_file, "\n")
if(!file.exists(data_file)) {
  stop(paste("Data file not found:", data_file, 
             "\nPlease check the file path in the configuration section."))
}
data <- read_dta(data_file)
cat("Data loaded successfully. Dimensions:", nrow(data), "rows,", ncol(data), "columns\n")

# Convert haven_labelled objects to standard R types
# This is necessary because read_dta() creates haven_labelled objects that some packages don't handle
cat("\n--- Converting haven_labelled objects to standard R types ---\n")
for(var in names(data)) {
  if(inherits(data[[var]], "haven_labelled")) {
    # Convert to numeric (preserves the numeric values, drops labels)
    data[[var]] <- as.numeric(data[[var]])
    cat("Converted", var, "from haven_labelled to numeric\n")
  }
}

# Display data summary
cat("\n--- Data Summary ---\n")
cat("First few rows:\n")
print(head(data))
cat("\nData structure:\n")
str(data, give.attr = FALSE)

# Check for missing values (should be none for synthetic data)
cat("\n--- Missing Data Check ---\n")
missing_summary <- sapply(data, function(x) sum(is.na(x)))
if(sum(missing_summary) > 0) {
  warning("Missing values detected in data. Consider using imputation.")
  print(missing_summary[missing_summary > 0])
} else {
  cat("No missing values detected. Imputation not needed.\n")
}

# Verify required variables exist
key_vars <- c(group, base_cov, mediator_var, outcome_var, intermediate_covs)
key_vars <- key_vars[key_vars %in% names(data)]
if(length(key_vars) < length(c(group, base_cov, mediator_var, outcome_var, intermediate_covs))) {
  missing_vars <- setdiff(c(group, base_cov, mediator_var, outcome_var, intermediate_covs), key_vars)
  stop(paste("The following required variables are missing from the data:", 
                paste(missing_vars, collapse = ", ")))
}

# Center baseline covariate and create squared term
base_col <- data[[base_cov]]
if(is.null(base_col)) {
  stop(paste("Baseline covariate", base_cov, "not found in data"))
}
data$age_centered <- base_col - mean(base_col, na.rm = TRUE)
data$age_centered_sq <- data$age_centered^2

# Convert group variable to factor
if(!group %in% names(data)) {
  stop(paste("Group variable", group, "not found in data"))
}
data[[group]] <- as.factor(data[[group]])



# ============================================================================
# VARIABLE DEFINITIONS
# ============================================================================
# - Group/Exposure: [group] (factor)
# - Mediator: [mediator_var] (continuous)
# - Outcome: [outcome_var] (continuous)
# - Baseline covariate: age (centered)
# - Additional covariates: [paste(intermediate_covs, collapse = ", ")]
# ============================================================================


# ============================================================================
# Bootstrap Confidence Intervals and Comparison Table
# ============================================================================

cat("\n\n========================================\n")
cat("Computing 95% Confidence Intervals and Comparison Table\n")
cat("========================================\n")

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# Bootstrap function for confidence intervals only (not estimates)
bootstrap_ci <- function(data, func, n_boot = 1000, conf_level = 0.95) {
  n <- nrow(data)
  boot_estimates <- replicate(n_boot, {
    boot_idx <- sample(1:n, n, replace = TRUE)
    boot_data <- data[boot_idx, ]
    tryCatch(func(boot_data), error = function(e) NA)
  })
  boot_estimates <- boot_estimates[!is.na(boot_estimates)]
  alpha <- 1 - conf_level
  ci_lower <- quantile(boot_estimates, alpha/2, na.rm = TRUE)
  ci_upper <- quantile(boot_estimates, 1 - alpha/2, na.rm = TRUE)
  return(list(ci_lower = ci_lower,
              ci_upper = ci_upper))
}

# Function for Difference in Coefficients method
# Returns: c(total_effect, direct_effect, indirect_effect, prop_explained)
func_dc <- function(d) {
  # Build formula strings dynamically
  covs_str <- paste(c("age_centered", "age_centered_sq", intermediate_covs), collapse = " + ")
  formula_t <- as.formula(paste(outcome_var, "~", group, "+", covs_str))
  formula_d <- as.formula(paste(outcome_var, "~", group, "+", mediator_var, "+", covs_str))
  
  model_t <- lm(formula_t, data = d)
  model_d <- lm(formula_d, data = d)
  # Extract coefficient for group variable (factor variable - coefficient name will be "group1" or similar)
  treat_coef_t <- coef(model_t)[grep(paste0("^", group), names(coef(model_t)))]
  treat_coef_d <- coef(model_d)[grep(paste0("^", group), names(coef(model_d)))]
  te <- if(length(treat_coef_t) > 0) treat_coef_t[1] else NA
  de <- if(length(treat_coef_d) > 0) treat_coef_d[1] else NA
  ie <- te - de
  prop_exp <- if(!is.na(te) && te != 0) (ie / te) * 100 else NA
  return(c(te, de, ie, prop_exp))
}

# Function for Oaxaca-Blinder (non-reference group as reference)
# Returns: c(total_difference, explained, unexplained, prop_explained)
func_kob_nonblack <- function(d) {
  # Get group levels (assuming factor with levels 0 and 1, or numeric 0/1)
  treat_levels <- levels(d[[group]])
  if(is.null(treat_levels)) {
    # If not a factor, assume 0 and 1
    d_black <- d[d[[group]] == 1, ]
    d_nonblack <- d[d[[group]] == 0, ]
  } else {
    d_black <- d[d[[group]] == treat_levels[2], ]
    d_nonblack <- d[d[[group]] == treat_levels[1], ]
  }
  
  if(nrow(d_black) < 2 || nrow(d_nonblack) < 2) return(c(NA, NA, NA, NA))
  
  mean_diff <- mean(d_black[[outcome_var]], na.rm = TRUE) - mean(d_nonblack[[outcome_var]], na.rm = TRUE)
  covs_str <- paste(c("age_centered", "age_centered_sq", intermediate_covs, mediator_var), collapse = " + ")
  formula_kob <- as.formula(paste(outcome_var, "~", covs_str))
  model_b <- lm(formula_kob, data = d_black)
  model_nb <- lm(formula_kob, data = d_nonblack)
  
  X_b <- colMeans(model.matrix(model_b), na.rm = TRUE)[mediator_var]
  X_nb <- colMeans(model.matrix(model_nb), na.rm = TRUE)[mediator_var]
  # Use non-reference group coefficients (beta_nb) for decomposition
  beta_nb <- coef(model_nb)[mediator_var]
  
  explained <- sum((X_b - X_nb) * beta_nb)
  unexplained <- mean_diff - explained
  prop_exp <- (explained / mean_diff) * 100
  
  return(c(mean_diff, explained, unexplained, prop_exp))
}

# Function for Oaxaca-Blinder (reference group as reference)
# Returns: c(total_difference, explained, unexplained, prop_explained)
func_kob_black <- function(d) {
  # Get group levels
  treat_levels <- levels(d[[group]])
  if(is.null(treat_levels)) {
    d_black <- d[d[[group]] == 1, ]
    d_nonblack <- d[d[[group]] == 0, ]
  } else {
    d_black <- d[d[[group]] == treat_levels[2], ]
    d_nonblack <- d[d[[group]] == treat_levels[1], ]
  }
  
  if(nrow(d_black) < 2 || nrow(d_nonblack) < 2) return(c(NA, NA, NA, NA))
  
  mean_diff <- mean(d_black[[outcome_var]], na.rm = TRUE) - mean(d_nonblack[[outcome_var]], na.rm = TRUE)
  covs_str <- paste(c("age_centered", "age_centered_sq", intermediate_covs, mediator_var), collapse = " + ")
  formula_kob <- as.formula(paste(outcome_var, "~", covs_str))
  model_b <- lm(formula_kob, data = d_black)
  model_nb <- lm(formula_kob, data = d_nonblack)
  
  X_b <- colMeans(model.matrix(model_b), na.rm = TRUE)[mediator_var]
  X_nb <- colMeans(model.matrix(model_nb), na.rm = TRUE)[mediator_var]
  # Use reference group coefficients (beta_b) for decomposition
  beta_b <- coef(model_b)[mediator_var]
  
  explained <- sum((X_b - X_nb) * beta_b)
  unexplained <- mean_diff - explained
  prop_exp <- (explained / mean_diff) * 100
  
  return(c(mean_diff, explained, unexplained, prop_exp))
}

# Compute estimates from original data first
cat("\n--- Computing Estimates from Original Data ---\n")

# Difference in Coefficients
cat("Computing: Difference in Coefficients...\n")
dc_orig <- func_dc(data)
dc_total_orig <- dc_orig[1]
dc_direct_orig <- dc_orig[2]
dc_indirect_orig <- dc_orig[3]
dc_prop_orig <- dc_orig[4]

# Oaxaca-Blinder (non-black reference)
cat("Computing: Oaxaca-Blinder (non-black reference)...\n")
kob_nb_orig <- func_kob_nonblack(data)
kob_nb_total_orig <- kob_nb_orig[1]
kob_nb_exp_orig <- kob_nb_orig[2]
kob_nb_unexp_orig <- kob_nb_orig[3]
kob_nb_prop_orig <- kob_nb_orig[4]

# Oaxaca-Blinder (black reference)
cat("Computing: Oaxaca-Blinder (black reference)...\n")
kob_b_orig <- func_kob_black(data)
kob_b_total_orig <- kob_b_orig[1]
kob_b_exp_orig <- kob_b_orig[2]
kob_b_unexp_orig <- kob_b_orig[3]
kob_b_prop_orig <- kob_b_orig[4]

# Run bootstrap for confidence intervals only
cat("\n--- Running Bootstrap for Confidence Intervals (", n_boot, " iterations) ---\n", sep = "")
set.seed(123)

# Difference in Coefficients
cat("Bootstrap: Difference in Coefficients...\n")
dc_boot <- bootstrap_ci(data, function(d) func_dc(d)[1], n_boot = n_boot, conf_level = conf_level)  # Total
dc_direct_boot <- bootstrap_ci(data, function(d) func_dc(d)[2], n_boot = n_boot, conf_level = conf_level)  # Direct
dc_indirect_boot <- bootstrap_ci(data, function(d) func_dc(d)[3], n_boot = n_boot, conf_level = conf_level)  # Indirect
dc_prop_boot <- bootstrap_ci(data, function(d) func_dc(d)[4], n_boot = n_boot, conf_level = conf_level)  # % Explained

# Oaxaca-Blinder (non-reference group as reference)
cat("Bootstrap: Oaxaca-Blinder (non-reference group as reference)...\n")
kob_nb_boot <- bootstrap_ci(data, function(d) func_kob_nonblack(d)[1], n_boot = n_boot, conf_level = conf_level)  # Total
kob_nb_exp_boot <- bootstrap_ci(data, function(d) func_kob_nonblack(d)[2], n_boot = n_boot, conf_level = conf_level)  # Explained
kob_nb_unexp_boot <- bootstrap_ci(data, function(d) func_kob_nonblack(d)[3], n_boot = n_boot, conf_level = conf_level)  # Unexplained
kob_nb_prop_boot <- bootstrap_ci(data, function(d) func_kob_nonblack(d)[4], n_boot = n_boot, conf_level = conf_level)  # % Explained

# Oaxaca-Blinder (reference group as reference)
cat("Bootstrap: Oaxaca-Blinder (reference group as reference)...\n")
kob_b_boot <- bootstrap_ci(data, function(d) func_kob_black(d)[1], n_boot = n_boot, conf_level = conf_level)  # Total
kob_b_exp_boot <- bootstrap_ci(data, function(d) func_kob_black(d)[2], n_boot = n_boot, conf_level = conf_level)  # Explained
kob_b_unexp_boot <- bootstrap_ci(data, function(d) func_kob_black(d)[3], n_boot = n_boot, conf_level = conf_level)  # Unexplained
kob_b_prop_boot <- bootstrap_ci(data, function(d) func_kob_black(d)[4], n_boot = n_boot, conf_level = conf_level)  # % Explained

# Causal Decomposition Analysis using Sequential Mediation Imputation (SMI)
cat("\n--- Causal Decomposition Analysis (SMI) ---\n")
covs_str_m <- paste(c("age_centered", "age_centered_sq"), collapse = " + ")
formula_m <- as.formula(paste(mediator_var, "~", group, "+", covs_str_m))
fit.m1 <- lm(formula_m, data = data)

covs_str_y <- paste(c("age_centered", "age_centered_sq", intermediate_covs), collapse = " + ")
formula_y <- as.formula(paste(outcome_var, "~", group, "*", mediator_var, "+", covs_str_y))
fit.y1 <- lm(formula_y, data = data)

res.1a1 <- smi(fit.m = fit.m1, fit.y = fit.y1, sims = smi_sims, conditional = smi_conditional, 
               covariates = c("age_centered"), group = group, seed = smi_seed)

# Display results
print(res.1a1)

cat("Extracting smi results...\n")
smi_total <- res.1a1$result[1,1]
smi_direct <- res.1a1$result[2,1]
smi_indirect <- res.1a1$result[3,1]
smi_prop <- (smi_indirect / smi_total) * 100

# Get CIs from smi if available, otherwise use bootstrap
  smi_total_ci <- res.1a1$result[1,2:3]
  smi_direct_ci <- res.1a1$result[2,2:3]
  smi_indirect_ci <- res.1a1$result[3,2:3]


# Create comparison table using original estimates
comparison_table <- data.frame(
  Component = c("Total Disparity", 
                "Explained Component (Disparity Reduction)",
                "Unexplained Component (Disparity Remaining)",
                "% Explained"),
  
  # Difference in Coefficients (using original estimates)
  DC_Estimate = c(dc_total_orig, dc_indirect_orig, dc_direct_orig, dc_prop_orig),
  DC_CI_Lower = c(dc_boot$ci_lower, dc_indirect_boot$ci_lower, dc_direct_boot$ci_lower, dc_prop_boot$ci_lower),
  DC_CI_Upper = c(dc_boot$ci_upper, dc_indirect_boot$ci_upper, dc_direct_boot$ci_upper, dc_prop_boot$ci_upper),
  
  # Oaxaca-Blinder (non-black reference) - using original estimates
  KOB_NB_Estimate = c(kob_nb_total_orig, kob_nb_exp_orig, kob_nb_unexp_orig, kob_nb_prop_orig),
  KOB_NB_CI_Lower = c(kob_nb_boot$ci_lower, kob_nb_exp_boot$ci_lower, kob_nb_unexp_boot$ci_lower, kob_nb_prop_boot$ci_lower),
  KOB_NB_CI_Upper = c(kob_nb_boot$ci_upper, kob_nb_exp_boot$ci_upper, kob_nb_unexp_boot$ci_upper, kob_nb_prop_boot$ci_upper),
  
  # Oaxaca-Blinder (black reference) - using original estimates
  KOB_B_Estimate = c(kob_b_total_orig, kob_b_exp_orig, kob_b_unexp_orig, kob_b_prop_orig),
  KOB_B_CI_Lower = c(kob_b_boot$ci_lower, kob_b_exp_boot$ci_lower, kob_b_unexp_boot$ci_lower, kob_b_prop_boot$ci_lower),
  KOB_B_CI_Upper = c(kob_b_boot$ci_upper, kob_b_exp_boot$ci_upper, kob_b_unexp_boot$ci_upper, kob_b_prop_boot$ci_upper),
  
  # Causal Decomposition (smi) - using original estimates from smi
  SMI_Estimate = c(smi_total, smi_indirect, smi_direct, smi_prop),
  SMI_CI_Lower = c(smi_total_ci[1], smi_indirect_ci[1], smi_direct_ci[1], NA),
  SMI_CI_Upper = c(smi_total_ci[2], smi_indirect_ci[2], smi_direct_ci[2], NA)
)

# Format the table for better display
cat("\n\n========================================\n")
cat("COMPARISON TABLE: Three Decomposition Methods\n")
cat("========================================\n\n")

# Print formatted table
for(i in seq_len(nrow(comparison_table))) {
  cat(sprintf("%-40s\n", comparison_table$Component[i]))
  cat(sprintf("  Difference in Coefficients: %.3f (95%% CI: %.3f, %.3f)\n",
              comparison_table$DC_Estimate[i], 
              comparison_table$DC_CI_Lower[i], 
              comparison_table$DC_CI_Upper[i]))
  cat(sprintf("  KOB (non-black ref):        %.3f (95%% CI: %.3f, %.3f)\n",
              comparison_table$KOB_NB_Estimate[i], 
              comparison_table$KOB_NB_CI_Lower[i], 
              comparison_table$KOB_NB_CI_Upper[i]))
  cat(sprintf("  KOB (black ref):            %.3f (95%% CI: %.3f, %.3f)\n",
              comparison_table$KOB_B_Estimate[i], 
              comparison_table$KOB_B_CI_Lower[i], 
              comparison_table$KOB_B_CI_Upper[i]))
  cat(sprintf("  Causal Decomp (smi):        %.3f (95%% CI: %.3f, %.3f)\n",
              comparison_table$SMI_Estimate[i], 
              comparison_table$SMI_CI_Lower[i], 
              comparison_table$SMI_CI_Upper[i]))
  cat("\n")
}

# Also print as data frame
print(comparison_table)


# ============================================================================
# Sensitivity Analysis with Cinelli's Covariate Benchmarks
# ============================================================================

cat("\n\n========================================\n")
cat("Sensitivity Analysis with Covariate Benchmarks\n")
cat("========================================\n")

# Run sensitivity analysis
treat_levels <- levels(data[[group]])
if(is.null(treat_levels) || length(treat_levels) < 2) {
  stop("Group variable must be a factor with at least 2 levels")
}
sel.lev.group <- treat_levels[2]  # Select non-reference level

sens_results <- sensitivity(
  boot.res = res.1a1,
  fit.y = fit.y1,
  fit.m = fit.m1,
  mediator = mediator_var,
  covariates = c("age_centered"),
  group = group,
  sel.lev.group = sel.lev.group,
  conf.level = conf_level,
  max.rsq = max_rsq
)

# Extract sensitivity analysis results
# Based on actual structure: sens_results has rm, ry, red, rem matrices
rm_values <- sens_results$rm  # R-squared with mediator (0 to 0.3)
ry_values <- sens_results$ry  # R-squared with outcome (0 to 0.3)
red_matrix <- -sens_results$red  # Disparity reduction matrix (31x31)
rem_matrix <- sens_results$rem  # Disparity remaining matrix (31x31)
lower_red <- -sens_results$lower_red  # Lower CI for reduction
lower_rem <- -sens_results$lower_rem  # Lower CI for remaining


# Calculate Cinelli's covariate benchmarks
# Partial R-squared of each covariate with exposure (A) and outcome (Y)
cat("\n--- Calculating Covariate Benchmarks ---\n")

# Observed covariates to benchmark
# Filter to only include covariates that exist in the data
benchmark_covariates <- unique(benchmark_covariates)
benchmark_covariates <- benchmark_covariates[benchmark_covariates %in% names(data)]
cat("Benchmark covariates:", paste(benchmark_covariates, collapse = ", "), "\n")
if(length(benchmark_covariates) == 0) {
  warning("No benchmark covariates found in data. Please check benchmark_covariates in configuration.")
}

# Calculate Cinelli's benchmarks for each covariate
# Need: r2.y (partial R2 with outcome) and r2.m (partial R2 with mediator)
benchmarks <- data.frame(
  covariate = character(),
  r2_y = numeric(),  # Partial R-squared with outcome
  r2_m = numeric(),  # Partial R-squared with mediator
  stringsAsFactors = FALSE
)

for(cov in benchmark_covariates) {
  if(cov %in% names(data)) {
    # Cinelli & Hazlett benchmark: Equations 54 and 56 from supplementary material
    # R = treatment variable, M = mediator variable, Y = outcome variable
    # X, C = other covariates (age_centered, age_centered_sq, and user-specified covariates)
    # X_j = cov (the benchmark covariate)
    
    # Define other covariates (excluding the benchmark covariate)
    all_covs <- c("age_centered", "age_centered_sq", intermediate_covs)
    other_covs <- setdiff(all_covs, cov)
    
    # ========================================================================
    # Equation 54: Partial R² for mediator
    # R²_{M~U|R,X,C} = R²_{M~X_j} / (1 - R²_{M~R+X+C})
    # ========================================================================
    
    # R²_{M~X_j}: R² from regressing mediator on just X_j
    model_m_xj <- lm(as.formula(paste(mediator_var, "~", cov)), data = data)
    r2_m_xj <- summary(model_m_xj)$r.squared
    
    # R²_{M~R+X+C}: R² from regressing mediator on R + X + C (excluding X_j)
    if(length(other_covs) > 0) {
      formula_m_rxc <- as.formula(paste(mediator_var, "~", group, "+", paste(other_covs, collapse = " + ")))
    } else {
      formula_m_rxc <- as.formula(paste(mediator_var, "~", group))
    }
    model_m_rxc <- lm(formula_m_rxc, data = data)
    r2_m_rxc <- summary(model_m_rxc)$r.squared
    
    # Calculate partial R² for mediator
    if(r2_m_rxc < 1) {
      r2_m <- r2_m_xj / (1 - r2_m_rxc)
    } else {
      r2_m <- 0
    }
    
    # ========================================================================
    # Equation 56: Partial R² for outcome
    # R²_{Y~U|R,X,M,C} = (|R_{Y~U|R,X,C}| - |R_{Y~M|R,X,C} * R_{M~U|R,X,C}|) / sqrt(1 - R²_{Y~M|R,X,C})
    # where R²_{Y~U|R,X,C} = R²_{Y~X_j} / (1 - R²_{Y~R+X+C})
    # ========================================================================
    
    # R²_{Y~X_j}: R² from regressing outcome on just X_j
    model_y_xj <- lm(as.formula(paste(outcome_var, "~", cov)), data = data)
    r2_y_xj <- summary(model_y_xj)$r.squared
    
    # R²_{Y~R+X+C}: R² from regressing outcome on R + X + C (excluding X_j and M)
    if(length(other_covs) > 0) {
      formula_y_rxc <- as.formula(paste(outcome_var, "~", group, "+", paste(other_covs, collapse = " + ")))
    } else {
      formula_y_rxc <- as.formula(paste(outcome_var, "~", group))
    }
    model_y_rxc <- lm(formula_y_rxc, data = data)
    r2_y_rxc <- summary(model_y_rxc)$r.squared
    
    # R²_{Y~U|R,X,C}
    if(r2_y_rxc < 1) {
      r2_y_u_rxc <- r2_y_xj / (1 - r2_y_rxc)
    } else {
      r2_y_u_rxc <- 0
    }
    
    # R²_{Y~M|R,X,C}: R² from regressing outcome on mediator, conditional on R + X + C
    if(length(other_covs) > 0) {
      formula_y_m_rxc <- as.formula(paste(outcome_var, "~", group, "+", mediator_var, "+", paste(other_covs, collapse = " + ")))
    } else {
      formula_y_m_rxc <- as.formula(paste(outcome_var, "~", group, "+", mediator_var))
    }
    model_y_m_rxc <- lm(formula_y_m_rxc, data = data)
    r2_y_m_rxc <- summary(model_y_m_rxc)$r.squared
    
    # To compute correlations, we need residuals from partial regressions
    # R_{Y~U|R,X,C}: correlation between residuals of Y~R+X+C and U~R+X+C
    # Since U is unobserved, we use X_j as proxy: correlation between residuals of Y~R+X+C and X_j~R+X+C
    residuals_y_rxc <- residuals(model_y_rxc)
    
    # Model X_j on R+X+C to get residuals
    if(length(other_covs) > 0) {
      formula_xj_rxc <- as.formula(paste(cov, "~", group, "+", paste(other_covs, collapse = " + ")))
    } else {
      formula_xj_rxc <- as.formula(paste(cov, "~", group))
    }
    model_xj_rxc <- lm(formula_xj_rxc, data = data)
    residuals_xj_rxc <- residuals(model_xj_rxc)
    
    # R_{Y~U|R,X,C} (using X_j as proxy for U)
    r_y_u_rxc <- cor(residuals_y_rxc, residuals_xj_rxc, use = "complete.obs")
    
    # R_{Y~M|R,X,C}: correlation between residuals of Y~R+X+C and M~R+X+C
    residuals_m_rxc <- residuals(model_m_rxc)
    r_y_m_rxc <- cor(residuals_y_rxc, residuals_m_rxc, use = "complete.obs")
    
    # R_{M~U|R,X,C}: correlation between residuals of M~R+X+C and U~R+X+C (using X_j as proxy)
    r_m_u_rxc <- cor(residuals_m_rxc, residuals_xj_rxc, use = "complete.obs")
    
    # Calculate partial R² for outcome using Equation 56
    numerator <- abs(r_y_u_rxc) - abs(r_y_m_rxc * r_m_u_rxc)
    denominator <- sqrt(1 - r2_y_m_rxc)
    
    if(denominator > 0 && numerator >= 0) {
      r2_y <- (numerator / denominator)^2
    } else {
      r2_y <- 0
    }
    
    # Ensure non-negative and handle edge cases
    r2_y <- max(0, min(r2_y, 1))
    r2_m <- max(0, min(r2_m, 1))
    
    benchmarks <- rbind(benchmarks, data.frame(
      covariate = cov,
      r2_y = r2_y,
      r2_m = r2_m
    ))
    
    cat(sprintf("%s: Partial R2 with outcome = %.4f, Partial R2 with mediator = %.4f\n", 
                cov, r2_y, r2_m))
  }
}

# Remove any duplicate rows from benchmarks (in case code was run multiple times)
benchmarks <- unique(benchmarks)

# Rename covariates for display in plots
benchmarks$covariate_display <- benchmarks$covariate
benchmarks$covariate_display[benchmarks$covariate == "age_centered"] <- "age"
benchmarks$covariate_display[benchmarks$covariate == "Pdivorce16"] <- "pdivorced"
benchmarks$covariate_display[benchmarks$covariate == "RTHLTHCH"] <- "maltreat"
# Keep chdSES as is

cat("\n--- Final Benchmark Summary ---\n")
cat("Number of unique benchmarks:", nrow(benchmarks), "\n")
print(benchmarks)

# Create custom sensitivity plots with Cinelli's benchmarks
library(ggplot2)
library(ggrepel)

# Prepare data for contour plots
# Create data frames from the matrices
red_df <- expand.grid(
  r2_y = ry_values,
  r2_m = rm_values
)
red_df$reduction <- as.vector(red_matrix)
red_df$lower_ci <- as.vector(lower_red)

rem_df <- expand.grid(
  r2_y = ry_values,
  r2_m = rm_values
)
rem_df$remaining <- as.vector(rem_matrix)
rem_df$lower_ci <- as.vector(lower_rem)

# Calculate reduction values at benchmark points
# Find closest match in red_df for each benchmark point
benchmarks$reduction_value <- NA
for(i in seq_len(nrow(benchmarks))) {
  # Find closest r2_y and r2_m values
  idx_y <- which.min(abs(red_df$r2_y - benchmarks$r2_y[i]))
  idx_m <- which.min(abs(red_df$r2_m - benchmarks$r2_m[i]))
  # Find row with both closest values
  distances <- sqrt((red_df$r2_y - benchmarks$r2_y[i])^2 + (red_df$r2_m - benchmarks$r2_m[i])^2)
  closest_idx <- which.min(distances)
  benchmarks$reduction_value[i] <- red_df$reduction[closest_idx]
}

# Plot 1: Disparity Reduction (contour plot)
p1 <- ggplot(red_df, aes(x = r2_y, y = r2_m, z = reduction)) +
  geom_contour_filled(alpha = 0.7, bins = 15) +
  scale_fill_grey(start = 0.9, end = 0.3,
                  guide = guide_legend(reverse = TRUE),
                  labels = function(x) {
                    # Extract all numbers from labels (handles intervals like "(0.1,0.2]" or "(-0.1,0.2]")
                    sapply(x, function(val) {
                      val_str <- as.character(val)
                      # Extract numbers with signs: pattern matches -?[0-9.]+ (including negative sign)
                      num_matches <- regmatches(val_str, gregexpr("-?[0-9]+\\.[0-9]+|-?[0-9]+", val_str))
                      nums <- as.numeric(unlist(num_matches))
                      
                      if(length(nums) >= 2) {
                        # For intervals, use the number furthest from zero to preserve sign
                        # This correctly handles both negative and positive intervals
                        furthest <- nums[which.max(abs(nums))]
                        sprintf("%.2f", furthest)
                      } else if(length(nums) == 1) {
                        # Single number (preserve sign)
                        sprintf("%.2f", nums[1])
                      } else {
                        # Fallback - try to extract as numeric directly
                        num_val <- suppressWarnings(as.numeric(val_str))
                        if(!is.na(num_val)) {
                          sprintf("%.2f", num_val)
                        } else {
                          as.character(val)
                        }
                      }
                    })
                  }) +
  geom_contour(color = "white", linewidth = 0.5, bins = 15) +
  
  # Add solid line where estimate (reduction) equals zero
  geom_contour(data = red_df, aes(x = r2_y, y = r2_m, z = reduction), 
               breaks = 0, color = "black", linetype = "solid", 
               linewidth = 1.2, inherit.aes = FALSE) +
  
  # Add dotted line where confidence interval covers zero (lower_ci = 0)
  geom_contour(data = red_df, aes(x = r2_y, y = r2_m, z = lower_ci), 
               breaks = 0, color = "black", linetype = "dashed", 
               linewidth = 1, inherit.aes = FALSE) +
  
  # Add benchmark points (inherit.aes = FALSE to avoid inheriting z aesthetic)
  geom_point(data = benchmarks, 
             aes(x = r2_y, y = r2_m),
             inherit.aes = FALSE,
             color = "black", fill = "black",
             size = 5, shape = 17, stroke = 2) +
  
  # Add covariate labels
  geom_text_repel(data = benchmarks,
                  aes(x = r2_y, y = r2_m, label = covariate_display),
                  inherit.aes = FALSE,
                  color = "black",
                  size = 5, fontface = "bold",
                  box.padding = 0.5,
                  point.padding = 0.3,
                  min.segment.length = 0.2,
                  segment.color = "black",
                  segment.size = 0.5,
                  direction = "y") +
  
  labs(
    x = "Partial R² with Outcome (r².y)",
    y = "Partial R² with Mediator (r².m)",
    fill = "Disparity\nReduction"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 13, color = "gray40"),
    axis.title = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(size = 14),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = "right"
  ) +
  scale_x_continuous(limits = c(0, 0.3), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.3), expand = c(0, 0))

print(p1)



# Save plots
ggsave(output_plot_reduction, p1, width = plot_width, height = plot_height, dpi = plot_dpi)
cat("\nPlot saved as '", output_plot_reduction, "'\n", sep = "")

# Print benchmark summary
cat("\n--- Covariate Benchmark Summary ---\n")
print(benchmarks)
