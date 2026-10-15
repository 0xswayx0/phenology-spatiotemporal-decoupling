# ==============================================================================
# A STATISTICAL FRAMEWORK FOR DECOUPLING SPATIAL AND ENVIRONMENTAL DRIVERS OF
# PHENOLOGY FROM PATCHY HISTORICAL DATASETS
# AUTHORS: Swantje Grabener, Claudia Drees, Michael Staab, Wiebke Schuett, 
# Pascale Zumstein, Thorsten Assmann
# ==============================================================================
# Target Journal: Methods in Ecology and Evolution (Workflow Article)
#
# Overview of the Pipeline:
#   This script executes the complete, reproducible statistical workflow across
#   three ground beetle species (Coleoptera: Carabidae) with contrasting life
#   histories and monitoring data constraints:
#     1. Pterostichus madidus   (Internal baseline template; Dual-tensor GAMM)
#     2. Notiophilus biguttatus (External Wijster template;  Dual-tensor GAMM)
#     3. Calathus micropterus   (Internal baseline template; Single-tensor GAMM)
#
# Required Input File:
#   - "Phenology_Pipeline_Data.RData" (contains pre-standardized dynamic thermal
#     drivers 'therm_daily', ECN site coordinates 'site_coords', ECN species 
#     count datasets 'ecn_species_data', and Wijster reference data 'wijster_nb')
#
# Outputs Generated:
#   - "Table_1_Model_Selection_Summary.csv" (AIC, Deviance Explained, Concurvity)
#   - Fitted GAMM objects (.rds) and publication figures (Figures 2, 4, and 5)
# ==============================================================================

# Load required packages
library(dplyr)
library(lubridate)
library(ggplot2)
library(mgcv)
library(patchwork)

# ==============================================================================
# STAGE 0: LOAD REPOSITORY DATA & DEFINE HELPER FUNCTIONS
# ==============================================================================
cat("Loading 'Phenology_Pipeline_Data.RData'...\n")
if (!file.exists("Phenology_Pipeline_Data.RData")) {
  stop("Input file 'Phenology_Pipeline_Data.RData' not found in working directory.")
}
load("Phenology_Pipeline_Data.RData")

#' Scan Trapping Time Series for Continuous Year-Round Sequences (>= 365 Days)
#'
#' Identifies contiguous trapping runs where internal hiatuses between 
#' consecutive trap deployments do not exceed 30 days (Manuscript Section 2.2).
find_continuous_runs <- function(df, start_col = "OpenDate", end_col = "CloseDate", 
                                 subsite_col = "SubSite", n_col = "N") {
  subsites <- unique(df[[subsite_col]])
  candidate_intervals <- list()
  
  for (sub in subsites) {
    df_s <- df[df[[subsite_col]] == sub, ]
    df_s <- df_s[order(df_s[[start_col]]), ]
    if (nrow(df_s) < 2) next
    
    run_start <- df_s[[start_col]][1]
    run_end   <- df_s[[end_col]][1]
    run_rows  <- c(1)
    
    for (i in 2:nrow(df_s)) {
      gap <- as.numeric(df_s[[start_col]][i] - run_end)
      if (gap <= 30) {
        run_end  <- max(run_end, df_s[[end_col]][i])
        run_rows <- c(run_rows, i)
      } else {
        dur <- as.numeric(run_end - run_start)
        if (dur >= 365) {
          candidate_intervals[[length(candidate_intervals) + 1]] <- list(
            SubSite = sub, Start = run_start, End = run_end,
            DurationDays = dur, Data = df_s[run_rows, ], 
            TotalN = sum(df_s[[n_col]][run_rows])
          )
        }
        run_start <- df_s[[start_col]][i]
        run_end   <- df_s[[end_col]][i]
        run_rows  <- c(i)
      }
    }
    dur <- as.numeric(run_end - run_start)
    if (dur >= 365) {
      candidate_intervals[[length(candidate_intervals) + 1]] <- list(
        SubSite = sub, Start = run_start, End = run_end,
        DurationDays = dur, Data = df_s[run_rows, ], 
        TotalN = sum(df_s[[n_col]][run_rows])
      )
    }
  }
  return(candidate_intervals)
}

# ==============================================================================
# MASTER WORKFLOW FUNCTION (MANUSCRIPT SECTIONS 2.1 - 2.5 & 3.1 - 3.3)
# ==============================================================================
run_phenology_workflow <- function(sp_name, sp_short, raw_ecn_df, 
                                   external_df = NULL, model_type = "dual", 
                                   out_fig_name) {
  
  cat(sprintf("\n==================================================================\n"))
  cat(sprintf("Executing Workflow for: %s (%s-tensor configuration)\n", sp_name, model_type))
  cat(sprintf("==================================================================\n"))
  
  # ----------------------------------------------------------------------------
  # STAGE 1: Effort Standardization & Abundance Filtering (Section 2.1 & 3.1)
  # ----------------------------------------------------------------------------
  # Calculate trapping interval duration, midpoint Day of Year (DOY), and 
  # log-effort offset: ln(Effort) = ln(NTraps * Duration).
  dat <- raw_ecn_df %>%
    mutate(OpenDate  = ymd(OpenDate), 
           CloseDate = ymd(CloseDate)) %>%
    mutate(Duration  = as.numeric(CloseDate - OpenDate),
           DOY       = pmin(365, yday(OpenDate + floor(Duration / 2))),
           Year      = year(OpenDate + floor(Duration / 2))) %>%
    filter(!is.na(Duration) & Duration > 0) %>%
    mutate(Effort    = NTraps * Duration,
           LogEffort = log(Effort),
           CPUE      = N / Effort) %>%
    # Join causally aligned, DOY-standardized dynamic thermal drivers (Section 2.3)
    left_join(therm_daily, by = c("Site", "Year", "DOY")) %>%
    left_join(site_coords, by = "Site") %>%
    filter(!is.na(Latitude) & !is.na(T_roll30) & !is.na(GDD5_roll60) & !is.na(GDD5_YTD))
  
  # Apply empirical abundance filter: retain subsite-years with >= 10 individuals
  # to prevent spline overfitting to stochastic single-catch noise (Appendix S1).
  sy_grid <- dat %>%
    group_by(Site, SubSite, Year) %>%
    summarize(TotalN = sum(N), .groups = "drop") %>%
    mutate(FilterPass = TotalN >= 10)
  
  cat(sprintf("1. Retained %d subsite-years passing abundance filter (N >= 10).\n", 
              sum(sy_grid$FilterPass)))
  
  # ----------------------------------------------------------------------------
  # STAGE 2: Baseline Whole-Year Phenological Template (Section 2.2, Step 1)
  # ----------------------------------------------------------------------------
  # Prioritize internal year-round trapping sequences (if TotalN >= 10); 
  # otherwise fall back to external continuous reference data (Wijster).
  cands_internal <- find_continuous_runs(dat)
  max_n_internal <- if (length(cands_internal) > 0) max(sapply(cands_internal, function(x) x$TotalN)) else 0
  
  if (max_n_internal >= 10) {
    best_idx <- which.max(sapply(cands_internal, function(x) x$TotalN))
    df_base  <- cands_internal[[best_idx]]$Data
    cat(sprintf("2. Built internal baseline template from SubSite %s (N = %d).\n", 
                cands_internal[[best_idx]]$SubSite, cands_internal[[best_idx]]$TotalN))
  } else {
    wij_prep <- external_df %>%
      mutate(StartDate = ymd(StartDate), 
             EndDate   = ymd(EndDate)) %>%
      mutate(TrapDays  = as.numeric(EndDate - StartDate),
             DOY       = pmin(365, yday(StartDate + floor(TrapDays / 2))),
             Effort    = TrapDays,
             LogEffort = log(Effort),
             SubSite   = "Wijster_Station")
    cands_ext <- find_continuous_runs(wij_prep, start_col = "StartDate", 
                                      end_col = "EndDate", subsite_col = "SubSite", n_col = "N")
    best_idx  <- which.max(sapply(cands_ext, function(x) x$TotalN))
    df_base   <- cands_ext[[best_idx]]$Data
    cat(sprintf("2. Built external baseline template from Wijster Station (N = %d).\n", 
                cands_ext[[best_idx]]$TotalN))
  }
  
  # Fit cyclic baseline GAM with Tweedie distribution and REML penalty inflation (gamma = 1.4)
  baseline_gam <- gam(N ~ s(DOY, bs = "cc", k = 15) + offset(LogEffort), 
                      data = df_base, family = tw(), method = "REML", 
                      gamma = 1.4, knots = list(DOY = c(0.5, 365.5)))
  
  pred_grid     <- data.frame(DOY = 1:365, LogEffort = 0)
  pred_grid$fit <- as.numeric(exp(predict(baseline_gam, newdata = pred_grid)))
  
  # Derive empirical detection threshold (tau_det): minimum predicted rate at N >= 1
  X_catch        <- predict(baseline_gam, 
                            newdata = data.frame(DOY = df_base$DOY[df_base$N > 0], LogEffort = 0), 
                            type = "lpmatrix")
  phys_threshold <- min(as.numeric(exp(X_catch %*% coef(baseline_gam))))
  
  # Truncate days below tau_det to structural zeros and normalize to daily probability density
  pred_grid$fit[pred_grid$fit < phys_threshold] <- 0
  template_density <- pred_grid$fit / sum(pred_grid$fit)
  
  # ----------------------------------------------------------------------------
  # STAGE 3: Automated Gap Detection & Synthetic Augmentation (Section 2.2, Step 2)
  # ----------------------------------------------------------------------------
  df_filtered <- dat %>% inner_join(sy_grid %>% filter(FilterPass), by = c("Site", "SubSite", "Year"))
  sy_unique   <- unique(df_filtered[, c("SubSite", "Year")])
  n_sy        <- nrow(sy_unique)
  
  # Construct binary daily surveillance coverage matrix (365 days x n_subsite_years)
  active_matrix <- matrix(FALSE, nrow = 365, ncol = n_sy)
  for (j in 1:n_sy) {
    df_j <- df_filtered %>% filter(SubSite == sy_unique$SubSite[j], Year == sy_unique$Year[j])
    for (r in 1:nrow(df_j)) {
      d_start <- min(365, max(1, yday(df_j$OpenDate[r])))
      d_end   <- min(365, max(1, yday(df_j$CloseDate[r])))
      if (year(df_j$OpenDate[r]) == year(df_j$CloseDate[r])) {
        if (d_start <= d_end) active_matrix[d_start:d_end, j] <- TRUE
      } else {
        active_matrix[d_start:365, j] <- TRUE
        active_matrix[1:d_end, j]     <- TRUE
      }
    }
  }
  
  # Flag contiguous calendar windows with < 5% network surveillance coverage
  gap_doys       <- which(rowMeans(active_matrix) < 0.05)
  breaks         <- c(0, which(diff(gap_doys) > 1), length(gap_doys))
  inferred_gaps  <- lapply(1:(length(breaks) - 1), function(k) {
    chunk <- gap_doys[(breaks[k] + 1):breaks[k + 1]]
    data.frame(Start = min(chunk), End = max(chunk))
  })
  global_gaps <- bind_rows(inferred_gaps)
  
  # Project observed rates onto baseline template via OLS through the origin (kappa_{l,y})
  # and inject down-weighted (w = 0.5) continuous Tweedie pseudo-observations
  synthetic_list <- list()
  for (i in 1:nrow(sy_grid)) {
    if (!sy_grid$FilterPass[i]) next
    sub   <- sy_grid$SubSite[i]
    yr    <- sy_grid$Year[i]
    df_sy <- df_filtered %>% filter(SubSite == sub, Year == yr)
    
    p_temp <- template_density[df_sy$DOY]
    if (sum(p_temp^2) == 0) next
    k_scale         <- sum(df_sy$CPUE * p_temp) / sum(p_temp^2)
    scaled_template <- template_density * k_scale
    
    global_gaps_sy <- global_gaps %>% filter(!(Start >= min(df_sy$DOY) & End <= max(df_sy$DOY)))
    if (nrow(global_gaps_sy) == 0) next
    
    for (g in 1:nrow(global_gaps_sy)) {
      g_start <- global_gaps_sy$Start[g]
      g_end   <- global_gaps_sy$End[g]
      if (g_end < g_start + 7) next
      
      syn_doys <- seq(g_start + 7, g_end, by = 14)
      ntraps   <- mean(df_sy$NTraps, na.rm = TRUE)
      syn_eff  <- ntraps * 14
      
      synthetic_list[[length(synthetic_list) + 1]] <- data.frame(
        Site      = df_sy$Site[1], 
        SubSite   = sub, 
        Year      = yr, 
        DOY       = syn_doys,
        N         = scaled_template[syn_doys] * syn_eff, # Continuous Tweedie count
        Effort    = syn_eff, 
        LogEffort = log(syn_eff), 
        Weight    = 0.5,                                 # Regularization prior weight
        Latitude  = df_sy$Latitude[1]
      )
    }
  }
  
  # Join dynamic daily thermal drivers to synthetic records by (Site, Year, DOY)
  df_syn <- bind_rows(synthetic_list) %>%
    left_join(therm_daily, by = c("Site", "Year", "DOY")) %>%
    filter(!is.na(T_roll30) & !is.na(GDD5_roll60) & !is.na(GDD5_YTD))
  
  cols_keep <- c("Site", "SubSite", "Year", "DOY", "N", "Effort", "LogEffort", 
                 "Weight", "Latitude", "T_roll30", "GDD5_roll60", "GDD5_YTD")
  
  df_final <- bind_rows(
    df_filtered %>% mutate(Weight = 1.0) %>% dplyr::select(all_of(cols_keep)),
    df_syn %>% dplyr::select(all_of(cols_keep))
  ) %>%
    mutate(SubSite_fac = as.factor(SubSite), 
           Year_fac    = as.factor(Year))
  
  cat(sprintf("3. Augmented dataset ready: %d empirical (w = 1.0) + %d synthetic (w = 0.5) records.\n",
              sum(df_final$Weight == 1.0), sum(df_final$Weight == 0.5)))
  
  # ----------------------------------------------------------------------------
  # STAGE 4: Candidate Driver Selection & Functional ANOVA GAMM (Sections 2.3 - 2.5)
  # ----------------------------------------------------------------------------
  candidate_drivers <- c("T_roll30", "GDD5_roll60", "GDD5_YTD")
  model_candidates  <- list()
  aic_records       <- list()
  
  for (drv in candidate_drivers) {
    if (model_type == "dual") {
      # Full dual-tensor functional ANOVA decomposition (Section 2.4)
      f <- as.formula(paste0(
        "N ~ s(DOY, bs = 'cc', k = 15) + s(Latitude, k = 3) + s(", drv, ", k = 4) + ",
        "ti(DOY, Latitude, bs = c('cc', 'cr'), k = c(12, 3)) + ",
        "ti(DOY, ", drv, ", bs = c('cc', 'cr'), k = c(12, 4)) + ",
        "s(SubSite_fac, bs = 're') + s(Year_fac, bs = 're') + offset(LogEffort)"
      ))
    } else {
      # Modular single-tensor thermal decomposition for spatially restricted taxa (Section 2.5)
      f <- as.formula(paste0(
        "N ~ s(DOY, bs = 'cc', k = 15) + s(", drv, ", k = 4) + ",
        "ti(DOY, ", drv, ", bs = c('cc', 'cr'), k = c(12, 4)) + ",
        "s(SubSite_fac, bs = 're') + s(Year_fac, bs = 're') + offset(LogEffort)"
      ))
    }
    
    m <- bam(f, data = df_final, weights = Weight, family = tw(),
             method = "fREML", discrete = TRUE, gamma = 1.4,
             knots = list(DOY = c(0.5, 365.5)))
    
    sm <- summary(m)
    model_candidates[[drv]] <- m
    aic_records[[drv]] <- data.frame(
      Species           = sp_name,
      ModelStructure    = paste0( tools::toTitleCase(model_type), "-tensor (n = ", sum(sy_grid$FilterPass), ")"),
      Driver            = drv,
      EDF               = round(sum(m$edf), 2),
      LogLik            = round(as.numeric(logLik(m)), 2),
      AIC               = round(AIC(m), 2),
      DevianceExplained = round(sm$dev.expl * 100, 2)
    )
  }
  
  aic_table <- bind_rows(aic_records) %>%
    mutate(DeltaAIC = round(AIC - min(AIC), 2)) %>%
    arrange(DeltaAIC)
  
  champ_driver <- aic_table$Driver[1]
  mod_final    <- model_candidates[[champ_driver]]
  
  # Evaluate estimate concurvity for shape-modifying tensor interactions (Section 3.2)
  conc <- tryCatch(concurvity(mod_final, full = TRUE)["estimate", ], error = function(e) NULL)
  aic_table$Concurvity_ti_Lat <- NA_real_
  aic_table$Concurvity_ti_Drv <- NA_real_
  if (!is.null(conc)) {
    lat_idx <- grep("ti\\(.*Latitude", names(conc))
    drv_idx <- grep(paste0("ti\\(.*", champ_driver), names(conc))
    if (length(lat_idx) > 0) aic_table$Concurvity_ti_Lat[1] <- round(conc[lat_idx[1]], 3)
    if (length(drv_idx) > 0) aic_table$Concurvity_ti_Drv[1] <- round(conc[drv_idx[1]], 3)
  }
  
  cat("\n4. Candidate Dynamic Driver Comparison (Table 1):\n")
  print(aic_table)
  
  saveRDS(mod_final, paste0("Champion_GAMM_", sp_short, ".rds"))
  
  # ----------------------------------------------------------------------------
  # STAGE 5: Posterior Simulation on the Unit Simplex & Visualization (Section 2.4)
  # ----------------------------------------------------------------------------
  n_sims <- 1000
  set.seed(11)
  beta <- coef(mod_final)
  V    <- vcov(mod_final, unconditional = TRUE)
  br   <- mgcv::rmvn(n_sims, beta, V)
  
  # Helper to simulate 1,000 posterior density curves marginalized over random effects
  extract_posterior_density <- function(drv_val, lat_val) {
    nd <- data.frame(DOY = 1:365, LogEffort = 0, 
                     SubSite_fac = df_final$SubSite_fac[1], 
                     Year_fac    = df_final$Year_fac[1], 
                     Latitude    = lat_val)
    nd[[champ_driver]] <- drv_val
    
    X <- predict(mod_final, newdata = nd, type = "lpmatrix")
    re_cols <- grep("s\\(Year_fac\\)|s\\(SubSite_fac\\)", colnames(X))
    if (length(re_cols) > 0) X[, re_cols] <- 0
    
    sim_val  <- exp(X %*% t(br))
    # Normalize each posterior draw onto the unit simplex (sums to 1.0 across DOY 1:365)
    sim_dens <- apply(sim_val, 2, function(x) if (sum(x) > 0) x / sum(x) else rep(1/365, 365))
    
    data.frame(
      DOY = 1:365,
      fit = rowMeans(sim_dens),
      lwr = apply(sim_dens, 1, quantile, probs = 0.025),
      upr = apply(sim_dens, 1, quantile, probs = 0.975)
    )
  }
  
  # Predict thermal plasticity across empirical envelope (5th vs. 95th percentile)
  lat_mean <- mean(df_final$Latitude, na.rm = TRUE)
  drv_5    <- as.numeric(quantile(df_final[[champ_driver]], 0.05, na.rm = TRUE))
  drv_95   <- as.numeric(quantile(df_final[[champ_driver]], 0.95, na.rm = TRUE))
  
  d_cold <- extract_posterior_density(drv_5, lat_mean)  %>% mutate(Scenario = "Cold Regime")
  d_warm <- extract_posterior_density(drv_95, lat_mean) %>% mutate(Scenario = "Warm Regime")
  
  p_thermal <- ggplot(bind_rows(d_cold, d_warm), aes(x = DOY)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr, fill = Scenario), alpha = 0.2) +
    geom_line(aes(y = fit, color = Scenario), linewidth = 0.8) +
    scale_color_manual(values = c("Cold Regime" = "#0072B2", "Warm Regime" = "#D55E00")) +
    scale_fill_manual(values = c("Cold Regime" = "#0072B2", "Warm Regime" = "#D55E00")) +
    labs(title = if (model_type == "dual") "b)" else NULL, 
         x = "Day of Year", y = "Probability Density") +
    theme_classic(base_size = 10) +
    theme(legend.position = "bottom", legend.title = element_blank(),
          plot.title = element_text(size = 10, face = "bold"))
  
  if (model_type == "dual") {
    # Predict latitudinal plasticity under isothermal conditions (5th vs. 95th percentile)
    lat_5    <- as.numeric(quantile(df_final$Latitude, 0.05, na.rm = TRUE))
    lat_95   <- as.numeric(quantile(df_final$Latitude, 0.95, na.rm = TRUE))
    drv_mean <- mean(df_final[[champ_driver]], na.rm = TRUE)
    
    d_south <- extract_posterior_density(drv_mean, lat_5)  %>% mutate(Scenario = "South")
    d_north <- extract_posterior_density(drv_mean, lat_95) %>% mutate(Scenario = "North")
    
    p_lat <- ggplot(bind_rows(d_south, d_north), aes(x = DOY)) +
      geom_ribbon(aes(ymin = lwr, ymax = upr, fill = Scenario), alpha = 0.2) +
      geom_line(aes(y = fit, color = Scenario), linewidth = 0.8) +
      scale_color_manual(values = c("South" = "#D55E00", "North" = "#0072B2")) +
      scale_fill_manual(values = c("South" = "#D55E00", "North" = "#0072B2")) +
      labs(title = "a)", x = "Day of Year", y = "Probability Density") +
      theme_classic(base_size = 10) +
      theme(legend.position = "bottom", legend.title = element_blank(),
            plot.title = element_text(size = 10, face = "bold"))
    
    p_out <- p_lat + p_thermal + plot_layout(ncol = 2)
    ggsave(out_fig_name, p_out, width = 16, height = 8.5, units = "cm", dpi = 300)
  } else {
    ggsave(out_fig_name, p_thermal, width = 8, height = 7.5, units = "cm", dpi = 300)
  }
  
  cat(sprintf("5. Saved figure: %s\n", out_fig_name))
  return(aic_table)
}

# ==============================================================================
# EXECUTE PIPELINE FOR ALL THREE CASE-STUDY SPECIES (SECTION 3.3)
# ==============================================================================

# 3.3.1 Pterostichus madidus (Full dual-tensor pipeline with internal ECN template)
res_pm <- run_phenology_workflow(
  sp_name      = "Pterostichus madidus",
  sp_short     = "PM",
  raw_ecn_df   = ecn_species_data$Pterostichus_madidus,
  model_type   = "dual",
  out_fig_name = "Figure_2_Pterostichus_madidus_Plasticity.png"
)

# 3.3.2 Notiophilus biguttatus (Dual-tensor pipeline with external Wijster template)
res_nb <- run_phenology_workflow(
  sp_name      = "Notiophilus biguttatus",
  sp_short     = "NB",
  raw_ecn_df   = ecn_species_data$Notiophilus_biguttatus,
  external_df  = wijster_nb,
  model_type   = "dual",
  out_fig_name = "Figure_4_Notiophilus_biguttatus_Plasticity.png"
)

# 3.3.3 Calathus micropterus (Modular single-tensor pipeline for localized range)
res_cm <- run_phenology_workflow(
  sp_name      = "Calathus micropterus",
  sp_short     = "CM",
  raw_ecn_df   = ecn_species_data$Calathus_micropterus,
  model_type   = "single",
  out_fig_name = "Figure_5_Calathus_micropterus_Plasticity.png"
)

# Combine and export Table 1
table_1_all <- bind_rows(res_pm, res_nb, res_cm)
write.csv(table_1_all, "Table_1_Model_Selection_Summary.csv", row.names = FALSE)
cat("\n==================================================================\n")
cat("All models completed! Exported 'Table_1_Model_Selection_Summary.csv'.\n")
cat("==================================================================\n")