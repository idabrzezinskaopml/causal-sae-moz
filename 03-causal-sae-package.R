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
library(glmnet)
library(MASS)

# Load population and sample data
sample_data <- readRDS("./output/sample_data_cov_std.rds")
population_data <- readRDS("./output/pop_ea_cov_std.rds")

# Note that sample data is still at the household level
# While population data is at EA level

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

# Aggregate
ea_sample_data <- sample_data %>%
  st_drop_geometry() %>%
  group_by(AES, BA_ID, ID_DIST, PA_ID, urban, province_name, Cod_provincia) %>%
  summarise(
    log_cons = mean(log_cons, na.rm = TRUE),
    weight_hh = mean(weight_hh, na.rm = T),
    
    # treatment
    flood_any_5km_no_eloise = first(flood_any_5km_no_eloise),
    
    # population
    population = first(population),
    
    # flow accumulation
    #flow_acc = first(flow_acc),
    
    # all geocovariates
    across(
      starts_with("mean."),
      ~ mean(.x, na.rm = TRUE)
    ),
    
    n_hh = n(),
    .groups = "drop"
  )


# Fit logistic model for propensity score
# Average household size from census at district level (divide population)

# Create non-sampled data object
sample_eas <- unique(ea_sample_data$AES)

non_sampled_data <- population_data %>%
  filter(!AES %in% sample_eas)

### COVARIATE SELECTION ###
lasso_df <- ea_sample_data %>%
  dplyr::select(log_cons, starts_with("mean"))

set.seed(123)
lasso_df <- na.omit(lasso_df)

lambdas_to_try <- 10 ^ seq(-3, 5, length.out = 100)

# Setting alpha = 1 implements lasso regression
lasso_cv <- cv.glmnet(lasso_df %>% dplyr::select(-log_cons) %>% as.matrix(),
                      lasso_df$log_cons,
                      alpha = 1, lambda = lambdas_to_try,
                      standardize = TRUE, nfolds = 10)

# Plot cross-validation results
plot(lasso_cv)

# Best cross-validated lambda
lambda_cv <- lasso_cv$lambda.min
lambda_cv <- 0.03

# Fit final model, get its sum of squared residuals and multiple R-squared
model_cv <- glmnet(lasso_df %>% dplyr::select(-log_cons) %>% as.matrix(),
                   lasso_df$log_cons,
                   alpha = 1, lambda = lambda_cv, standardize = TRUE)


lasso_predictors <- coef(model_cv, s = lambda_cv) %>%
  as.matrix() %>%
  as.data.frame() %>%
  rename(coef = 1) %>%
  mutate(var = rownames(.)) %>%
  as_tibble() %>%
  filter(coef != 0) %>%
  filter(var != "(Intercept)") %>%
  pull(var)


design_matrix <- model.matrix(as.formula(paste("~", paste(lasso_predictors, collapse = " + "), collapse = " ")),
                              data = lasso_df)

eigenvalues <- eigen(t(design_matrix) %*% design_matrix)$values

# Take a look
lasso_predictors


## Manual selection of predictors

xvars <- c(
  "mean.buildings_count_BCB_gl_100m_v1_1",
 "mean.MOZ_level0_100m_2015_2030",
 "mean.highway_dist_osm_2023_100m_v1",
  "mean.ppt_2022_yravg_tc_100m_v1",
  "mean.tavg_2022_tlst_100m_v1",
  "mean.viirs_nvf_2022_100m_v1"
)

## TO CHECK: does it really make sense to remove sampled EAs??

obj_p_score_EBLUP <- list(data_p_score = ea_sample_data)
class(obj_p_score_EBLUP) <- "EBLUP"

form <- as.formula(
  paste(
    "flood_any_5km_no_eloise ~",
    paste(c(lasso_predictors, "(1|ID_DIST)"), collapse = " + ")
  )
)


ps_hat_EBLUP <-  p_score(obj_p_score = obj_p_score_EBLUP,
                         model_formula = form)

# Attach propensity score to sample data
ea_sample_data$ps_hat_EBLUP <- ps_hat_EBLUP

# treatment: (consumption x treatment ) / propensity score - mean of the district
# control: consumption x (1 - treatment) / 1 - propensity - mean of the districts 
# T - C for all districts 

district_ate <- ea_sample_data %>%
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
  dplyr::select(IPW_ATE, WIPW_ATE) %>%
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
  "mean.viirs_nvf_2022_100m_v1"
)

non_sampled_data_est <- population_data %>%
  st_drop_geometry() %>%
  as.data.frame() %>%
  dplyr::select(all_of(lasso_predictors), ID_DIST, flood_any_5km_no_eloise) %>%
  rename(group = ID_DIST,
         A = flood_any_5km_no_eloise)

sample_data_est <- ea_sample_data %>%
  st_drop_geometry() %>%
  as.data.frame() %>%
  dplyr::select(all_of(lasso_predictors), ID_DIST, flood_any_5km_no_eloise, log_cons) %>%
  rename(group = ID_DIST,
         A = flood_any_5km_no_eloise,
         y = log_cons)

# Recode A and group as integer
non_sampled_data_est$A <- as.integer(non_sampled_data_est$A)

# This seems to be a glitch in the package - where y is hard-coded into the package
# Adding a null value for now
non_sampled_data_est$y <- NA_real_

form_OR <- as.formula(
  paste(
    "y ~",
    paste(predictors, collapse = " + "),
    "+ (1|group)"
  )
)

form_OR


hte_OR_ida <- hte(type_hte = "OR",
              sample_data_est,
              non_sampled_data_est,
              params_OR = list(model_formula = form_OR,
                               method = "EBLUP",
                               type_model = "gaussian"))

# Fixed effects from LASSO
x_string <- paste(predictors, collapse = " + ")

# Outcome/imputation model
form_y <- as.formula(
  paste(
    "y ~",
    x_string,
    "+ (1 + A || group)"
  )
)

# Propensity score model
form_ps <- as.formula(
  paste(
    "A ~",
    x_string,
    "+ (1 | group)"
  )
)


hte_NIPW_ida <- hte(
  type_hte = "NIPW",
  sample_data_est,
  non_sampled_data_est,
  params_impute_y = list(
    model_formula = form_y,
    method = "EBLUP",
    type_model = "gaussian"
  ),
  params_p_score = list(
    model_formula = form_ps,
    method = "EBLUP"
  )
)

fixed_effects <- paste(predictors, collapse = " + ")

form_impute_y <- as.formula(
  paste(
    "y ~",
    fixed_effects,
    "+ (1 + A||group)"
  )
)

form_p_score <- as.formula(
  paste(
    "A ~",
    fixed_effects,
    "+ (1|group)"
  )
)

form_OR <- as.formula(
  paste(
    "y ~",
    fixed_effects,
    "+ (1 + A||group)"
  )
)


fixed_effects <- paste(predictors, collapse = " + ")

form_impute_y <- as.formula(
  paste("y ~", fixed_effects, "+ (1 + A||group)")
)

form_p_score <- as.formula(
  paste("A ~", fixed_effects, "+ (1|group)")
)

form_OR <- as.formula(
  paste("y ~", fixed_effects, "+ (1 + A||group)")
)

hte_AIPW_ida <- hte(
  type_hte = "AIPW",
  sample_data_est,
  non_sampled_data_est,
  
  params_impute_y = list(
    model_formula = form_impute_y,
    method = "EBLUP",
    type_model = "gaussian"
  ),
  
  params_p_score = list(
    model_formula = form_p_score,
    method = "EBLUP"
  ),
  
  params_OR = list(
    model_formula = form_OR,
    method = "RF",
    tune_RF = FALSE,
    clust_RF = FALSE,
    type_model = "continuous"
  )
)


write.csv(hte_OR_ida, "./output/hte_OR.csv", row.names = FALSE)
write.csv(hte_NIPW_ida, "./output/hte_NIPW.csv", row.names = FALSE)
write.csv(hte_AIPW_ida, "./output/hte_AIPW.csv", row.names = FALSE)

#### MANUAL CODE FROM KASIA ####


csae_estimators <- function(data_sample, data_out_of_sample,
                            y = "y",
                            A = "corp",
                            group = "group") {
  
  # Covariates X
  #xvars <- c(
  #  "mean.buildings_count_BCB_gl_100m_v1_1",
  #  "mean.MOZ_level0_100m_2015_2030",
  #  "mean.highway_dist_osm_2023_100m_v1",
  #  "mean.ppt_2022_yravg_tc_100m_v1",
  #  "mean.tavg_2022_tlst_100m_v1",
  #  "mean.viirs_nvf_2022_100m_v1"
  #)
  
  xvars <- lasso_predictors
  
  # Mark sampled / nonsampled observations
  data_sample$.S <- 1
  data_out_of_sample$.S <- 0
  
  # Entire population N
  pop <- bind_rows(data_sample, data_out_of_sample)
  
  # ----------------------------------------------------------
  # 1. Outcome nuisance model: mu_n(X,A)
  #    Estimated ONLY in sampled data (n)
  # ----------------------------------------------------------
  
  # Allows treatment effect to vary with X
  f_y <- reformulate(
    paste0(A, " * (", paste(xvars, collapse = " + "), ")"),
    response = y
  )
  
  fit_y <- lm(f_y, data = data_sample)
  
  # mu_n(X,A): prediction under observed treatment
  pop$muA_n <- predict(fit_y, newdata = pop)
  
  # mu_1,n(X)
  dat1 <- pop
  dat1[[A]] <- 1
  pop$mu1_n <- predict(fit_y, newdata = dat1)
  
  # mu_0,n(X)
  dat0 <- pop
  dat0[[A]] <- 0
  pop$mu0_n <- predict(fit_y, newdata = dat0)
  
  # Y-tilde:
  # observed Y for sampled units;
  # imputed Y for nonsampled units
  pop$Ytilde <- ifelse(
    pop$.S == 1,
    pop[[y]],
    pop$muA_n
  )
  
  # ----------------------------------------------------------
  # 2. Propensity score: e_a,N(X)
  #    Estimated using ENTIRE population N
  # ----------------------------------------------------------
  
  f_e <- reformulate(xvars, response = A)
  
  fit_e <- glm(
    f_e,
    data = pop,
    family = binomial()
  )
  
  pop$e1_N <- predict(fit_e, newdata = pop, type = "response")
  pop$e0_N <- 1 - pop$e1_N
  
  # Optional numerical protection
  eps <- 1e-6
  pop$e1_N <- pmin(pmax(pop$e1_N, eps), 1 - eps)
  pop$e0_N <- pmin(pmax(pop$e0_N, eps), 1 - eps)
  
  # ----------------------------------------------------------
  # 3. CSAE estimators within each area/group j
  # ----------------------------------------------------------
  
  result <- pop %>%
    group_by(.data[[group]]) %>%
    summarise(
      
      N_j = n(),
      
      # ---------------------------------
      # CSAE-OR
      # P_Nj[mu_1,n - mu_0,n]
      # ---------------------------------
      OR = mean(mu1_n - mu0_n, na.rm = TRUE),
      
      # ---------------------------------
      # CSAE-IPW
      # ---------------------------------
      IPW =
        mean(.data[[A]] * Ytilde / e1_N, na.rm = TRUE) -
        mean((1 - .data[[A]]) * Ytilde / e0_N, na.rm = TRUE),
      
      # ---------------------------------
      # CSAE-NIPW
      # normalized / Hájek IPW
      # ---------------------------------
      NIPW =
        sum(.data[[A]] * Ytilde / e1_N, na.rm = TRUE) /
        sum(.data[[A]] / e1_N, na.rm = TRUE) -
        sum((1 - .data[[A]]) * Ytilde / e0_N, na.rm = TRUE) /
        sum((1 - .data[[A]]) / e0_N, na.rm = TRUE),
      
      # ---------------------------------
      # CSAE-AIPW
      # ---------------------------------
      AIPW =
        mean(
          .data[[A]] / e1_N * (Ytilde - mu1_n) + mu1_n,
          na.rm = TRUE
        ) -
        mean(
          (1 - .data[[A]]) / e0_N * (Ytilde - mu0_n) + mu0_n,
          na.rm = TRUE
        ),
      
      .groups = "drop"
    )
  
  return(list(
    estimates = result,
    population_data = pop,
    outcome_model = fit_y,
    propensity_model = fit_e
  ))
}

test_model <- csae_estimators(data_sample = sample_data_est,
                data_out_of_sample = non_sampled_data_est,
                A = "A")

test_result <- test_model$estimates

p_model <- test_model$propensity_model


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
       filename = "./plots/causal_ATE_plot_package.png")


library(tidyverse)

## MAP THE MANUAL CODE VERSION 

# Convert to long format
plot_data <- test_result %>%
  pivot_longer(
    cols = c(OR, NIPW, AIPW), # IPW, OR, NIPW, AIPW
    names_to = "Method",
    values_to = "Estimate"
  ) %>%
  mutate(
    group = factor(group, levels = unique(group))
  )# %>%
#filter(!(Method == "IPW" & Estimate < -15)) %>%
#filter(!(Method == "NIPW" & Estimate < -1)) %>%
#filter(!(Method == "OR" & Estimate < -1)) %>%
#filter(!(Method == "AIPW" & Estimate < -1))


# Plot
p <- ggplot(
  plot_data,
  aes(x = group, y = Estimate, colour = Method)
) +
  geom_point(
    position = position_dodge(width = 0.5),
    size = 3
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    colour = "black"
  ) +
  labs(
    x = "Area",
    y = "ATE",
    colour = "Method"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(
      size = 11,
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(size = 11)
  )

ggsave(
  filename = "./plots/ATE_comparison_plot.png",
  plot = p,
  width = 8,
  height = 5
)

p

