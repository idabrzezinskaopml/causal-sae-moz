library(dplyr)

csae_estimators <- function(data_sample, data_out_of_sample,
                            y = "y",
                            A = "corp",
                            group = "group") {
  
  # Covariates X
  xvars <- c(
    "mean.buildings_count_BCB_gl_100m_v1_1",
    "mean.MOZ_level0_100m_2015_2030",
    "mean.highway_dist_osm_2023_100m_v1",
    "mean.ppt_2022_yravg_tc_100m_v1",
    "mean.tavg_2022_t1st_100m_v1",
    "mean.viirs_nvf_2022_100m_v1",
    "flow_acc"
  )
  
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