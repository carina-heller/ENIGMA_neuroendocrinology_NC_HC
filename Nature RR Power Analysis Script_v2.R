#-------------------------------------------------------------------------------
              ## Nature Registered Report Power Analysis Script ##
#-------------------------------------------------------------------------------

# Clear environment  
rm(list = ls())

# Set working directory  
setwd("~/Documents/R Working Directory/ENIGMA")

################################################################################  
## Power analysis for Menstrual Phase LMM  
## Model: outcome ~ phase + age + age² + ICV + (1|site) + (1|participant)  
## Phase: 3 levels (early follicular, late follicular, luteal)  
## FDR correction: Benjamini-Hochberg across 34 Desikan atlas regions  
## 3 separate families: surface area, volume, cortical thickness  
## 5 true effect regions, 29 null regions per family  
## Power definitions: any, all, average  
################################################################################

## ---- packages ----  
library(lme4)  
library(lmerTest)  
library(MASS)  # for mvrnorm (correlated regions)  
library(ggplot2)  
set.seed(123)

################################################################################  
## 1. SIMULATION PARAMETERS  
################################################################################

# Effect sizes (Cohen's d for phase effect)  
effect_sizes <- c(0.20, 0.30, 0.40, 0.48)

# Sample sizes to test (number of unique participants)  
sample_sizes <- c(660, 1000, 1500, 2000, 2500, 3000, 3500)

# Number of sites  
n_sites <- 66

# Longitudinal data structure  
prop_longitudinal <- 0.40  # ~40% have multiple time points  
mean_timepoints_if_longitudinal <- 2.5  # Average 2-3 time points for longitudinal participants

# Age distributions to test  
age_distributions <- list(  
  uniform = function(n) runif(n, 18, 45),
  
  normal_young = function(n) {  
    age <- rnorm(n, mean = 28, sd = 8)  
    pmin(pmax(age, 18), 45)  
  },
  
  normal_middle = function(n) {  
    age <- rnorm(n, mean = 32, sd = 8)  
    pmin(pmax(age, 18), 45)  
  }  
)

# Random effects - typical neuroimaging values  
site_icc <- 0.10        # Sites explain 10% of variance  
participant_icc <- 0.40 # Participants explain 40% of variance  
total_var <- 1          # Standardized outcome

# Calculate variance components  
var_participant <- participant_icc * total_var  
var_site <- site_icc * total_var  
var_residual <- total_var - var_participant - var_site

sigma_participant <- sqrt(var_participant)  
sigma_site <- sqrt(var_site)  
sigma_eps <- sqrt(var_residual)

# Fixed effects (standardized)  
beta_intercept <- 0      # intercept  
beta_age <- 0.2          # small age effect  
beta_age2 <- -0.1        # small quadratic age effect  
beta_ICV <- 0.3          # moderate ICV effect

# ICV correlation with outcome  
icv_correlation <- 0.4

# Significance level (FDR threshold)  
alpha <- 0.05

# Number of simulations  
n_sims <- 100

# Desikan atlas parameters  
n_regions <- 34          # regions per family  
n_true <- 5              # true effect regions per family  
n_null <- n_regions - n_true  # null regions per family  
region_correlation <- 0.50   # inter-region correlation

# Brain measure families  
brain_measures <- c("surface_area", "volume", "cortical_thickness",   
                    "sulcal_depth", "gyrification")

# Calculate expected observations per participant  
expected_obs_per_participant <- (1 - prop_longitudinal) * 1 +  
  prop_longitudinal * mean_timepoints_if_longitudinal

cat("Expected observations per participant:", round(expected_obs_per_participant, 2), "\n")

################################################################################  
## 2. FUNCTION: Generate participants per site  
################################################################################  
generate_site_allocation <- function(n_participants, n_sites) {  
  # Generate random proportions  
  props <- runif(n_sites)  
  props <- props / sum(props)
  
  # Allocate participants  
  n_per_site <- floor(props * n_participants)
  
  # Ensure minimum of 4 per site  
  n_per_site <- pmax(n_per_site, 10)
  
  # Adjust to match exactly n_participants  
  current_total <- sum(n_per_site)
  
  if (current_total < n_participants) {  
    diff <- n_participants - current_total  
    sites_to_add <- sample(1:n_sites, diff, replace = TRUE)  
    for (s in sites_to_add) {  
      n_per_site[s] <- n_per_site[s] + 1  
    }  
  } else if (current_total > n_participants) {  
    diff <- current_total - n_participants  
    while (diff > 0) {  
      eligible_sites <- which(n_per_site > 10)  
      if (length(eligible_sites) == 0) {  
        eligible_sites <- 1:n_sites  
      }  
      largest_site <- eligible_sites[which.max(n_per_site[eligible_sites])]  
      n_per_site[largest_site] <- n_per_site[largest_site] - 1  
      diff <- diff - 1  
    }  
  }
  
  if (sum(n_per_site) != n_participants) {  
    stop("Site allocation failed: ", sum(n_per_site), " != ", n_participants)  
  }
  
  return(n_per_site)  
}

################################################################################  
## 3. FUNCTION: Generate time points per participant  
################################################################################  
generate_time_points <- function(n_participants) {  
  is_longitudinal <- rbinom(n_participants, 1, prop_longitudinal)  
  n_time_points <- ifelse(is_longitudinal,  
                          sample(2:3, n_participants, replace = TRUE),  
                          1)  
  return(n_time_points)  
}

################################################################################  
## 4. FUNCTION: Build correlation matrix for regions  
################################################################################  
# Compound symmetry correlation matrix:  
# All regions have correlation = region_correlation with each other  
# Diagonal = 1  
build_correlation_matrix <- function(n_regions, region_correlation) {  
  cor_mat <- matrix(region_correlation, nrow = n_regions, ncol = n_regions)  
  diag(cor_mat) <- 1  
  return(cor_mat)  
}

################################################################################  
## 5. FUNCTION: Simulate one dataset & test  
## Returns 3 power indicators (any, all, average) for one brain measure family  
################################################################################

sim_one <- function(n_participants, n_sites, effect_size, age_fun) {
  
  tryCatch({
    
    # -------------------------------------------------------------------  
    # STEP 1: Build participant structure  
    # -------------------------------------------------------------------
    
    # Allocate participants to sites  
    n_per_site <- generate_site_allocation(n_participants, n_sites)
    
    # Generate time points per participant  
    n_time_points_per_participant <- generate_time_points(n_participants)
    
    # Create participant-level data  
    participant_data <- data.frame(  
      participant_id = 1:n_participants,  
      site = rep(1:n_sites, times = n_per_site),  
      n_time_points = n_time_points_per_participant,  
      baseline_age = age_fun(n_participants),  
      u_participant = rnorm(n_participants, 0, sigma_participant)  
    )
    
    # Expand to observation level  
    obs_data <- do.call(rbind, lapply(1:n_participants, function(i) {  
      n_obs <- participant_data$n_time_points[i]  
      data.frame(  
        participant_id = rep(i, n_obs),  
        time_point = 1:n_obs,  
        time_since_baseline = (1:n_obs - 1) * (1/12)  
      )  
    }))
    
    # Merge with participant data  
    data <- merge(obs_data, participant_data, by = "participant_id")  
    n_total_obs <- nrow(data)
    
    # -------------------------------------------------------------------  
    # STEP 2: Generate covariates  
    # -------------------------------------------------------------------
    
    # Site random effects  
    u_site <- rnorm(n_sites, 0, sigma_site)  
    data$random_site <- u_site[data$site]
    
    # Age at each time point  
    data$age <- data$baseline_age + data$time_since_baseline  
    data$age_scaled <- scale(data$age)[,1]  
    data$age2_scaled <- data$age_scaled^2
    
    # ICV (constant within participant, correlated with outcome)  
    participant_baseline <- with(participant_data,  
                                 beta_intercept +  
                                   beta_age * scale(baseline_age)[,1] +  
                                   beta_age2 * scale(baseline_age)[,1]^2 +  
                                   u_participant)
    
    participant_icv <- participant_baseline * icv_correlation +  
      rnorm(n_participants, 0, sqrt(1 - icv_correlation^2))
    
    data$icv <- participant_icv[data$participant_id]  
    data$icv_scaled <- scale(data$icv)[,1]
    
    # Phase: 3 levels (early follicular = reference, late follicular, luteal)  
    data$phase <- sample(c("early_follicular", "late_follicular", "luteal"),  
                         n_total_obs, replace = TRUE)  
    data$phase <- factor(data$phase,  
                         levels = c("early_follicular", "late_follicular", "luteal"))
    
    # Dummy code phase contrasts  
    # beta_phase applies to BOTH late follicular and luteal (same Cohen's d)  
    beta_phase <- effect_size * sigma_eps  
    data$phase_late_fol <- as.numeric(data$phase == "late_follicular")  
    data$phase_luteal <- as.numeric(data$phase == "luteal")
    
    # -------------------------------------------------------------------  
    # STEP 3: Generate correlated outcomes across 34 regions  
    # -------------------------------------------------------------------
    
    # Build correlation matrix for regions  
    cor_mat <- build_correlation_matrix(n_regions, region_correlation)
    
    # Generate correlated residuals across regions for all observations  
    # mvrnorm produces n_total_obs x n_regions matrix of correlated errors  
    correlated_errors <- mvrnorm(n = n_total_obs,  
                                 mu = rep(0, n_regions),  
                                 Sigma = cor_mat * sigma_eps^2)
    
    # Linear predictor (shared across all regions as baseline)  
    # This is the systematic part WITHOUT phase effect  
    mu_base <- beta_intercept +  
      beta_age * data$age_scaled +  
      beta_age2 * data$age2_scaled +  
      beta_ICV * data$icv_scaled +  
      data$random_site +  
      data$u_participant
    
    # Phase contribution:  
    # - True regions (1:n_true): phase effect = beta_phase  
    # - Null regions (n_true+1:n_regions): phase effect = 0  
    phase_effect <- beta_phase * (data$phase_late_fol + data$phase_luteal)
    
    # Build outcome matrix: n_total_obs x n_regions  
    Y <- matrix(NA, nrow = n_total_obs, ncol = n_regions)
    
    for (r in 1:n_regions) {  
      if (r <= n_true) {  
        # True effect region  
        Y[, r] <- mu_base + phase_effect + correlated_errors[, r]  
      } else {  
        # Null region (no phase effect)  
        Y[, r] <- mu_base + correlated_errors[, r]  
      }  
    }
    
    # -------------------------------------------------------------------  
    # STEP 4: Fit LMM for each region & extract F-test p-value  
    # -------------------------------------------------------------------
    
    p_values <- numeric(n_regions)
    
    for (r in 1:n_regions) {
      
      dat_r <- data.frame(  
        y = Y[, r],  
        phase = data$phase,  
        age = data$age_scaled,  
        age2 = data$age2_scaled,  
        ICV = data$icv_scaled,  
        site = factor(data$site),  
        participant_id = factor(data$participant_id)  
      )
      
      fit <- suppressWarnings(  
        lmer(y ~ phase + age + age2 + ICV + (1|site) + (1|participant_id),  
             data = dat_r, REML = TRUE)  
      )
      
      # Overall F-test for phase factor (2 df: late_fol and luteal contrasts)  
      anova_fit <- anova(fit, type = "III")  
      p_values[r] <- anova_fit["phase", "Pr(>F)"]  
    }
    
    # -------------------------------------------------------------------  
    # STEP 5: Apply BH FDR correction across 34 regions  
    # -------------------------------------------------------------------
    
    p_adjusted <- p.adjust(p_values, method = "BH")
    
    # Which regions survive FDR correction?  
    significant <- p_adjusted < alpha
    
    # True regions are 1:n_true  
    true_region_significant <- significant[1:n_true]
    
    # -------------------------------------------------------------------  
    # STEP 6: Calculate 3 power definitions  
    # -------------------------------------------------------------------
    
    # Any: at least 1 true region survives  
    power_any <- as.numeric(any(true_region_significant))
    
    # All: all 5 true regions survive  
    power_all <- as.numeric(all(true_region_significant))
    
    # Average: proportion of true regions that survive  
    power_avg <- mean(true_region_significant)
    
    return(c(power_any = power_any,  
             power_all = power_all,  
             power_avg = power_avg))
    
  }, error = function(e) {  
    cat("\nError in simulation:", e$message, "\n")  
    return(c(power_any = NA, power_all = NA, power_avg = NA))  
  })  
}

################################################################################  
## 6. POWER ANALYSIS LOOP  
################################################################################

# Create results data frame  
# One row per combination of age distribution, effect size, sample size  
results <- expand.grid(  
  age_dist = names(age_distributions),  
  effect_size = effect_sizes,  
  sample_size = sample_sizes,  
  power_any = NA_real_,  
  power_all = NA_real_,  
  power_avg = NA_real_,  
  avg_total_obs = NA_real_,  
  stringsAsFactors = FALSE  
)

cat("\n=== STARTING POWER ANALYSIS ===\n")  
cat("Phase levels: early follicular (ref), late follicular, luteal\n")  
cat("Regions per family:", n_regions, "| True regions:", n_true, "| Null regions:", n_null, "\n")  
cat("Brain measure families:", paste(brain_measures, collapse = ", "), "\n")  
cat("FDR method: Benjamini-Hochberg | Alpha:", alpha, "\n")  
cat("Inter-region correlation:", region_correlation, "\n")  
cat("Proportion with longitudinal data:", round(prop_longitudinal * 100), "%\n")  
cat("This will run", nrow(results) * n_sims, "total simulations\n\n")

for (i in 1:nrow(results)) {
  
  age_dist_name <- results$age_dist[i]  
  d <- results$effect_size[i]  
  n <- results$sample_size[i]
  
  cat("Age dist:", age_dist_name,  
      "| Effect size:", d,  
      "| N participants:", n,  
      "| Expected total obs:", round(n * expected_obs_per_participant), "...\n")
  
  age_fun <- age_distributions[[age_dist_name]]
  
  # Run simulations  
  sim_results <- replicate(n_sims, sim_one(n, n_sites, d, age_fun))  
  # sim_results is 3 x n_sims matrix (power_any, power_all, power_avg)
  
  # Handle NAs  
  valid_cols <- !is.na(sim_results[1,])  
  n_valid <- sum(valid_cols)
  
  if (n_valid < n_sims) {  
    cat("  Warning:", n_sims - n_valid, "simulations failed\n")  
  }
  
  if (n_valid > 0) {  
    results$power_any[i] <- mean(sim_results["power_any", valid_cols])  
    results$power_all[i] <- mean(sim_results["power_all", valid_cols])  
    results$power_avg[i] <- mean(sim_results["power_avg", valid_cols])  
  }
  
  results$avg_total_obs[i] <- round(n * expected_obs_per_participant)  
}

################################################################################  
## 7. RESULTS & SAVE  
################################################################################

print(results)  
write.csv(results, "power_analysis_results_FDR_3phase.csv", row.names = FALSE)

################################################################################  
## 8. VISUALIZATION  
################################################################################

# Reshape results to long format for plotting  
results_long <- reshape(results,  
                        varying = c("power_any", "power_all", "power_avg"),  
                        v.names = "power",  
                        timevar = "power_type",  
                        times = c("Any true region", "All true regions", "Average proportion"),  
                        direction = "long")

# Power type labels  
results_long$power_type <- factor(results_long$power_type,  
                                  levels = c("Any true region",  
                                             "All true regions",  
                                             "Average proportion"))

# Plot for each age distribution  
for (age_dist_name in names(age_distributions)) {
  
  dat_plot <- subset(results_long, age_dist == age_dist_name)
  
  p <- ggplot(dat_plot, aes(x = sample_size, y = power,  
                            color = factor(effect_size),  
                            group = factor(effect_size))) +  
    geom_line(size = 1) +  
    geom_point(size = 2) +  
    geom_hline(yintercept = 0.80, linetype = "dashed", color = "gray50") +  
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +  
    facet_wrap(~power_type, ncol = 3) +  
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +  
    scale_color_brewer(palette = "RdYlBu", direction = -1,  
                       name = "Cohen's d") +  
    labs(title = paste("Power Analysis -", age_dist_name, "age distribution"),  
         subtitle = paste("3-phase model | BH FDR correction | 5/34 true regions",  
                          "| Inter-region r =", region_correlation),  
         x = "Number of Unique Participants",  
         y = "Statistical Power") +  
    theme_minimal() +  
    theme(legend.position = "right",  
          plot.title = element_text(hjust = 0.5),  
          plot.subtitle = element_text(hjust = 0.5))
  
  print(p)
  
  ggsave(paste0("power_plot_FDR_3phase_", age_dist_name, ".png"),  
         p, width = 12, height = 5, dpi = 300)  
}

# Combined plot across all age distributions (faceted by age dist and power type)  
p_combined <- ggplot(results_long, aes(x = sample_size, y = power,  
                                       color = factor(effect_size),  
                                       group = factor(effect_size))) +  
  geom_line(size = 0.8) +  
  geom_point(size = 1.5) +  
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "gray50") +  
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +  
  facet_grid(power_type ~ age_dist) +  
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +  
  scale_color_brewer(palette = "RdYlBu", direction = -1,  
                     name = "Cohen's d") +  
  labs(title = "Power Analysis: 3-Phase Menstrual Cycle Model with BH FDR Correction",  
       subtitle = paste("60% cross-sectional, 40% longitudinal | 5/34 true regions",  
                        "| Inter-region r =", region_correlation,  
                        "| Separate families per brain measure"),  
       x = "Number of Unique Participants",  
       y = "Statistical Power") +  
  theme_minimal() +  
  theme(legend.position = "bottom",  
        plot.title = element_text(hjust = 0.5),  
        plot.subtitle = element_text(hjust = 0.5))

print(p_combined)  
ggsave("power_plot_FDR_3phase_combined.png", p_combined, width = 14, height = 10, dpi = 300)

################################################################################  
## 9. SUMMARY TABLE at N = 3500  
################################################################################

summary_3500 <- subset(results, sample_size == 3500)  
summary_3500 <- summary_3500[order(summary_3500$age_dist,  
                                   summary_3500$effect_size), ]

cat("\n=== POWER AT N = 3500 PARTICIPANTS ===\n")  
cat("Expected total observations:", unique(summary_3500$avg_total_obs), "\n")  
cat("Participant ICC:", participant_icc, "| Site ICC:", site_icc, "\n")  
cat("FDR method: BH | Alpha:", alpha, "\n")  
cat("True regions:", n_true, "/ Total regions:", n_regions, "\n\n")  
print(summary_3500)

################################################################################  
## 10. MINIMUM DETECTABLE EFFECT SIZE  
################################################################################

cat("\n=== MINIMUM DETECTABLE EFFECT SIZE AT N = 3500 ===\n")  
for (age_dist_name in names(age_distributions)) {  
  cat("\nAge distribution:", age_dist_name, "\n")  
  dat_sub <- subset(summary_3500, age_dist == age_dist_name)
  
  for (power_type in c("power_any", "power_all", "power_avg")) {  
    power_vals <- dat_sub[[power_type]]
    
    if (any(power_vals >= 0.80, na.rm = TRUE)) {  
      min_d_80 <- min(dat_sub$effect_size[power_vals >= 0.80], na.rm = TRUE)  
      cat(" ", power_type, ": Cohen's d >=", min_d_80, "for 80% power\n")  
    } else {  
      cat(" ", power_type, ": 80% power not achieved at any effect size tested\n")  
    }
    
    if (any(power_vals >= 0.85, na.rm = TRUE)) {  
      min_d_80 <- min(dat_sub$effect_size[power_vals >= 0.85], na.rm = TRUE)  
      cat(" ", power_type, ": Cohen's d >=", min_d_80, "for 85% power\n")  
    } else {  
      cat(" ", power_type, ": 85% power not achieved at any effect size tested\n")  
    }
    
    if (any(power_vals >= 0.95, na.rm = TRUE)) {  
      min_d_95 <- min(dat_sub$effect_size[power_vals >= 0.95], na.rm = TRUE)  
      cat(" ", power_type, ": Cohen's d >=", min_d_95, "for 95% power\n")  
    } else {  
      cat(" ", power_type, ": 95% power not achieved at any effect size tested\n")  
    }  
  }  
}

################################################################################  
## 11. DESIGN & VARIANCE SUMMARY  
################################################################################

cat("\n=== VARIANCE DECOMPOSITION ===\n")  
cat("Participant ICC:", participant_icc, "\n")  
cat("Site ICC:", site_icc, "\n")  
cat("Residual variance proportion:", 1 - participant_icc - site_icc, "\n")

cat("\n=== DESIGN SUMMARY ===\n")  
cat("Phase levels: early follicular (ref), late follicular, luteal\n")  
cat("Cross-sectional participants:", round((1 - prop_longitudinal) * 100), "%\n")  
cat("Longitudinal participants:", round(prop_longitudinal * 100), "%\n")  
cat("Time points for longitudinal: 2-3\n")  
cat("Average observations per participant:", round(expected_obs_per_participant, 2), "\n")  
cat("Brain measure families:", paste(brain_measures, collapse = ", "), "\n")  
cat("Regions per family:", n_regions, "\n")  
cat("True effect regions per family:", n_true, "\n")  
cat("Null regions per family:", n_null, "\n")  
cat("Inter-region correlation:", region_correlation, "\n")  
cat("FDR method: Benjamini-Hochberg\n")  
cat("FDR threshold:", alpha, "\n")  

