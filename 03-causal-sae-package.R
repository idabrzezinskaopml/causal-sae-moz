# Packages
library(dplyr)
library(tidyverse)
library(remotes)
library(mquantreg)
library(CausalSAE)
library(sf)
library(haven)
library(grf)
library(tidyr)
library(ggplot2)

# Load population and sample data
sample_data <- readRDS("./output/sample_data_cov.rds")
population_data <- readRDS("./output/pop_ea_cov.rds")

# Add District identifier. First 4 digits = District
population_data$ID_DIST <- substr(population_data$BA_ID, 1, 4)

# Check identifiers
pop_dist <- unique(population_data$ID_DIST)
sample_dist <- unique(sample_data$ID_DIST)

# 41 vs 42?? not sure why 
setdiff(pop_dist, sample_dist)
setdiff(sample_dist, pop_dist)

# "0414" district ID is in population for some reason, but not sample
# Dropping it for now - to be investigated 
population_data <- subset(
  population_data,
  ID_DIST != "0414"
)

# Get sample weights for IOF 2022
iof_2022_raw <-  read_dta("Z:/01 Mozambique/01 Data/01 Poverty/04 Household data MZ/02 WB data/iof2022_v3.dta")

# Leave only hhid and sample weights - aggregate at hh level 
hh_weights <- iof_2022_raw %>%
  select(hhid, weight_hh) %>%
  group_by(hhid) %>%
  summarise(
    weight_hh = max(weight_hh, na.rm = TRUE),
    .groups = "drop"
  )

# Remove large IOF file
rm(iof_2022_raw)

# Attach weights to the sample data
sample_data <- left_join(sample_data, hh_weights)

# Move to the front
sample_data <- sample_data %>%
  relocate(weight_hh, .before = hhid)

# Fit logistic model for propensity score
# Average household size from census at district level (divide population)

# Create non-sampled data object
sample_eas <- unique(sample_data$AES)

non_sampled_data <- population_data %>%
  filter(!AES %in% sample_eas)

## TO CHECK: does it really make sense to remove sampled EAs??

obj_p_score_EBLUP <- list(data_p_score = sample_data)
class(obj_p_score_EBLUP) <- "EBLUP"

ps_hat_EBLUP <-  p_score(obj_p_score = obj_p_score_EBLUP,
                         model_formula = flood_any_5km_no_eloise ~ mean.buildings_count_BCB_gl_100m_v1_1 +
                           mean.MOZ_level0_100m_2015_2030 + mean.highway_dist_osm_2023_100m_v1 + 
                           mean.ppt_2022_yravg_tc_100m_v1 + mean.tavg_2022_tlst_100m_v1 +
                           mean.viirs_nvf_2022_100m_v1 + flow_acc +
                           (1|ID_DIST))

# Attach propensity score to sample data
sample_data$ps_hat_EBLUP <- ps_hat_EBLUP

# treatment: (consumption x treatment ) / propensity score - mean of the district
# control: consumption x (1 - treatment) / 1 - propensity - mean of the districts 
# T - C for all districts 

district_ate <- sample_data %>%
  group_by(ID_DIST) %>%
  summarise(
    # IPW estimator
    T_mean_ipw = mean(
      (log_cons * flood_any_5km_no_eloise) / ps_hat_EBLUP,
      na.rm = TRUE
    ),
    C_mean_ipw = mean(
      (log_cons * (1 - flood_any_5km_no_eloise)) / (1 - ps_hat_EBLUP),
      na.rm = TRUE
    ),
    IPW_ATE = T_mean_ipw - C_mean_ipw,
    
    # Survey-weighted IPW estimator
    T_mean_wipw = sum(
      weight_hh * log_cons * flood_any_5km_no_eloise / ps_hat_EBLUP,
      na.rm = TRUE
    ) /
      sum(
        weight_hh * flood_any_5km_no_eloise / ps_hat_EBLUP,
        na.rm = TRUE
      ),
    
    C_mean_wipw = sum(
      weight_hh * log_cons * (1 - flood_any_5km_no_eloise) / (1 - ps_hat_EBLUP),
      na.rm = TRUE
    ) /
      sum(
        weight_hh * (1 - flood_any_5km_no_eloise) / (1 - ps_hat_EBLUP),
        na.rm = TRUE
      ),
    
    WIPW_ATE = T_mean_wipw - C_mean_wipw,
    
    .groups = "drop"
  )

# Map the naive IPW estimator
ipw_plot_data <- district_ate %>%
  select(IPW_ATE, WIPW_ATE) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Method",
    values_to = "ATE"
  ) %>%
  mutate(
    Method = recode(
      Method,
      IPW_ATE = "Naive IPW",
      WIPW_ATE = "Survey-weighted naive IPW"
    ),
    District = rep(seq_len(nrow(district_ate)), 2)
  )

naive_plot <- ggplot(ipw_plot_data,
       aes(x = District, y = ATE, colour = Method)) +
  geom_point(
    position = position_dodge(width = 0.4),
    size = 3
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    colour = "black"
  ) +
  labs(
    x = "District",
    y = "ATE",
    colour = "Method"
  ) +
  theme_classic()

ggsave(naive_plot,
       filename = "./plots/naive_ATE.png")


## Example ###

m = 50
ni = rep(10, m)
Ni = rep(200, m)
N = sum(Ni)
n = sum(ni)

X <- generate_X(
  n = N,
  p = 1,
  covariance_norm = NULL,
  cov_type = "unif",
  seed = 1
)

X_outcome <- generate_X(
  n = N,
  p = 1,
  covariance_norm = NULL,
  cov_type = "lognorm",
  seed = 1
)

populations <- generate_pop(X, X_outcome,
                            coeffs = get_default_coeffs(),
                            errors_outcome = get_default_errors_outcome(),
                            rand_eff_outcome = get_default_rand_eff_outcome(),
                            rand_eff_p_score = get_default_rand_eff_p_score(),
                            regression_type = "continuous",
                            Ni_size  = 200,
                            m = 50,
                            no_sim = 1,
                            seed = 10)

samples <- generate_sample(populations, ni_size = 10,
                           sample_part = "sampled",
                           get_index = TRUE)

data_sample <- data.frame(samples[[1]]$samp_data)
index_sample <- samples[[1]]$index_s
data_out_of_sample <- populations[-index_sample, ]

hte_OR <- hte(type_hte = "OR",
              data_sample,
              data_out_of_sample,
              params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                               method = "EBLUP",
                               type_model = "gaussian"))

hte_NIPW <- hte(type_hte = "NIPW",
                data_sample,
                data_out_of_sample,
                params_impute_y = list(model_formula = y ~ X1 + Xo1 + (1 + A||group),
                                       method = "EBLUP",
                                       type_model = "gaussian"),
                params_p_score =  list(model_formula = A ~ X1 + Xo1 + (1|group),
                                       method = "EBLUP"))

hte_AIPW <- hte(type_hte = "AIPW",
                data_sample,
                data_out_of_sample,
                params_impute_y = list(model_formula = y ~ X1 + Xo1 + (1 + A||group),
                                       method = "EBLUP",
                                       type_model = "gaussian"),
                params_p_score =  list(model_formula = A ~ X1 + Xo1 + (1|group),
                                       method = "EBLUP"),
                params_OR = list(model_formula = y ~ X1 + Xo1 + (1 + A||group),
                                 method = "MQ",
                                 type_model = "continuous"))

## Try to recreate the exact data structure as in example 
## and run functions from the CausalSAE package

# non_sampled_data
# sample_data

# A - treatment variable
# y - outcome 

# Keep only common variables and transform to data frame (not spatial) +
# rename all columns 

predictors <- c(
  "mean.buildings_count_BCB_gl_100m_v1_1",
  "mean.MOZ_level0_100m_2015_2030",
  "mean.highway_dist_osm_2023_100m_v1",
  "mean.ppt_2022_yravg_tc_100m_v1",
  "mean.tavg_2022_tlst_100m_v1",
  "mean.viirs_nvf_2022_100m_v1",
  "flow_acc"
)

non_sampled_data_est <- population_data %>%
  st_drop_geometry() %>%
  as.data.frame() %>%
  select(all_of(predictors), ID_DIST, flood_any_5km_no_eloise) %>%
  rename(group = ID_DIST,
         A = flood_any_5km_no_eloise)

sample_data_est <- sample_data %>%
  st_drop_geometry() %>%
  as.data.frame() %>%
  select(all_of(predictors), ID_DIST, flood_any_5km_no_eloise, log_cons) %>%
  rename(group = ID_DIST,
         A = flood_any_5km_no_eloise,
         y = log_cons)

# Recode A and group as integer
non_sampled_data_est$A <- as.integer(non_sampled_data_est$A)

# This seems to be a glitch in the package - where y is hard-coded into the package
# Adding a null value for now
non_sampled_data_est$y <- NA_real_



hte_OR_ida <- hte(type_hte = "OR",
              sample_data_est,
              non_sampled_data_est,
              params_OR = list(model_formula = y ~ mean.buildings_count_BCB_gl_100m_v1_1 +
                                 mean.MOZ_level0_100m_2015_2030 + mean.highway_dist_osm_2023_100m_v1 + 
                                 mean.ppt_2022_yravg_tc_100m_v1 + mean.tavg_2022_tlst_100m_v1 +
                                 mean.viirs_nvf_2022_100m_v1 + flow_acc + (1|group),
                               method = "EBLUP",
                               type_model = "gaussian"))

hte_NIPW_ida <- hte(type_hte = "NIPW",
                    sample_data_est,
                    non_sampled_data_est,
                params_impute_y = list(model_formula = y ~ mean.buildings_count_BCB_gl_100m_v1_1 +
                                         mean.MOZ_level0_100m_2015_2030 + mean.highway_dist_osm_2023_100m_v1 + 
                                         mean.ppt_2022_yravg_tc_100m_v1 + mean.tavg_2022_tlst_100m_v1 +
                                         mean.viirs_nvf_2022_100m_v1 + flow_acc + (1 + A||group),
                                       method = "EBLUP",
                                       type_model = "gaussian"),
                params_p_score =  list(model_formula = A ~ mean.buildings_count_BCB_gl_100m_v1_1 +
                                         mean.MOZ_level0_100m_2015_2030 + mean.highway_dist_osm_2023_100m_v1 + 
                                         mean.ppt_2022_yravg_tc_100m_v1 + mean.tavg_2022_tlst_100m_v1 +
                                         mean.viirs_nvf_2022_100m_v1 + flow_acc + (1|group),
                                       method = "EBLUP"))


hte_AIPW_ida <- hte(type_hte = "AIPW",
                    sample_data_est,
                    non_sampled_data_est,
                params_impute_y = list(model_formula = y ~ mean.buildings_count_BCB_gl_100m_v1_1 +
                                         mean.MOZ_level0_100m_2015_2030 + mean.highway_dist_osm_2023_100m_v1 + 
                                         mean.ppt_2022_yravg_tc_100m_v1 + mean.tavg_2022_tlst_100m_v1 +
                                         mean.viirs_nvf_2022_100m_v1 + flow_acc + (1 + A||group),
                                       method = "EBLUP",
                                       type_model = "gaussian"),
                params_p_score =  list(model_formula = A ~ mean.buildings_count_BCB_gl_100m_v1_1 +
                                         mean.MOZ_level0_100m_2015_2030 + mean.highway_dist_osm_2023_100m_v1 + 
                                         mean.ppt_2022_yravg_tc_100m_v1 + mean.tavg_2022_tlst_100m_v1 +
                                         mean.viirs_nvf_2022_100m_v1 + flow_acc + (1|group),
                                       method = "EBLUP"),
                params_OR = list(model_formula = y ~ mean.buildings_count_BCB_gl_100m_v1_1 +
                                   mean.MOZ_level0_100m_2015_2030 + mean.highway_dist_osm_2023_100m_v1 + 
                                   mean.ppt_2022_yravg_tc_100m_v1 + mean.tavg_2022_tlst_100m_v1 +
                                   mean.viirs_nvf_2022_100m_v1 + flow_acc + (1 + A||group),
                                 method = "RF",
                                 tune_RF = F,
                                 clust_RF = F,
                                 type_model = "continuous"))


write.csv(hte_OR_ida, "./output/hte_OR.csv", row.names = FALSE)
write.csv(hte_NIPW_ida, "./output/hte_NIPW.csv", row.names = FALSE)
write.csv(hte_AIPW_ida, "./output/hte_AIPW.csv", row.names = FALSE)

## Plot results
library(dplyr)
library(ggplot2)

plot_data <- bind_rows(
  `CSAE-OR`   = hte_OR_ida,
  `CSAE-NIPW` = hte_NIPW_ida,
  `CSAE-AIPW` = hte_AIPW_ida,
  .id = "Method"
) %>%
  group_by(Method) %>%
  mutate(estimate = row_number()) %>%
  ungroup()

# Remove outliers
plot_data <- plot_data %>%
  filter(tau > -2)

p <- ggplot(plot_data,
       aes(x = estimate, y = tau, colour = Method)) +
  geom_point(
    position = position_dodge(width = 0.25),
    size = 3
  ) +
  scale_x_continuous(
    breaks = 1:n_distinct(plot_data$estimate),
    labels = 1:n_distinct(plot_data$estimate)
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    colour = "black"
  ) +
  labs(
    x = NULL,
    y = "ATE",
    colour = "Method"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11)
  )

ggsave(p, 
       filename = "./plots/causal_ATE_plot_without_outliers.png")


