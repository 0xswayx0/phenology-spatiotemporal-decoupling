# ==============================================================================
# A STATISTICAL FRAMEWORK FOR DECOUPLING SPATIAL AND METEOROLOGICAL PREDICTORS
# OF PHENOLOGY FROM INCOMPLETE MONITORING DATASETS (v0.3.0)
# ==============================================================================
# Target Journal: Methods in Ecology and Evolution (Workflow Article)
#
# Key Methodological Features & Revisions Implemented in v0.3.0:
#   1. Module 1 (Section 2.1):
#      - Explicit leap-year clamping (DOY = pmin(365L, yday(...))) across both
#        trapping midpoints and meteorological grids to prevent silent NA drops
#        or cyclic spline knot-boundary truncation on December 31 of leap years.
#      - Sequential three-step pre-processing: (i) filter site-years by
#        N_{l,y} >= N_min (default 10, with automated N_min in {3, 5, 10}
#        sensitivity check at cold/northern margins), (ii) identify the
#        species-level systematic network gap D_gap^spec (<5% active coverage),
#        and (iii) evaluate residual unmonitored gaps G_{res,l,y} <= 60 d
#        outside D_gap^spec across spatial and climatic gradients.
#   2. Module 2 (Section 2.2):
#      - Stage 2a pools all qualifying continuous year-round reference series
#        (>= 365 d, internal gaps <= 30 d, TotalN >= 10), includes random
#        intercepts across pooled reference series/subsites when multiple
#        reference runs exist, and defines the physical detection threshold
#        tau_det via the 1st catch-weighted percentile Q_{0.01}^{(w)} of
#        positive reference captures (bounded below by relative density
#        threshold epsilon = 1e-3 of annual catch).
#      - Stage 2b propagates Stage 2a posterior template uncertainty across
#        M_imp = 50 parameter draws beta_base^(m) ~ N(beta_hat_base, V_hat_base)
#        into downstream credible intervals.
#   3. Module 3 (Section 2.3):
#      - Documents standardization of daily meteorological predictors across the
#        complete balanced station x year meteorological grid with a seasonal
#        dispersion floor max(sigma_{net,d}, sigma_min) for early-winter GDD_YTD.
#   4. Module 4 (Section 2.4):
#      - Calibrates the Tweedie dispersion parameter phi_hat_emp strictly on
#        empirical observations (w_j = 1.0) to prevent downward dispersion
#        dilution by winter-zero priors, scaling posterior covariance V_hat_c^cal
#        and evaluating empirical-only Log-Likelihood (LogLik_Emp) and AIC_Emp
#        (alongside pooled AIC_Comb).
#      - Verifies 2D joint empirical covariate hull H_2D(Latitude, W*) support
#        via grDevices::chull and mgcv::exclude.too.far, and exports
#        Supplementary Figure S4 alongside Figures 3, 4, and 5.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(ggplot2)
  library(mgcv)
  library(patchwork)
  library(DHARMa)
})

WORK_DIR <- "C:/Users/swant/Desktop/Antigravity/Phenologies/Paper1_Methids_GAMM_v3"
if (dir.exists(WORK_DIR)) setwd(WORK_DIR)

# ==============================================================================
# STAGE 0: LOAD REPOSITORY DATA & DEFINE CORE HELPER FUNCTIONS
# ==============================================================================
cat("Loading 'Phenology_Pipeline_Data.RData'...\n")
if (!file.exists("Phenology_Pipeline_Data.RData")) {
  stop("Input file 'Phenology_Pipeline_Data.RData' not found in working directory.")
}
load("Phenology_Pipeline_Data.RData")

#' Reference Implementation of Module 3 Daily Network Standardization (Eq. 7)
#' Computes network-wide mean and standard deviation across the complete
#' balanced station x year meteorological grid on each calendar day d in 1:365,
#' clamping leap-year DOY = 366 to 365 and regularizing the denominator via
#' max(sigma_{net,d}, sigma_min) so early-January zero-inflated GDD_YTD cannot
#' divide by zero or generate extreme leverage spikes.
standardize_raw_met_grid <- function(raw_met_df, var_cols = c("T_roll30", "GDD5_roll60", "GDD5_YTD"),
                                     floor_quantile = 0.10) {
  df_clamped <- raw_met_df %>%
    mutate(DOY = pmin(365L, as.integer(DOY))) %>%
    group_by(Site, Year, DOY) %>%
    summarize(across(all_of(var_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

  for (vc in var_cols) {
    daily_stats <- df_clamped %>%
      group_by(DOY) %>%
      summarize(
        mu_d = mean(.data[[vc]], na.rm = TRUE),
        sd_d = sd(.data[[vc]], na.rm = TRUE),
        .groups = "drop"
      )
    pos_sds   <- daily_stats$sd_d[!is.na(daily_stats$sd_d) & daily_stats$sd_d > 0]
    sigma_min <- if (length(pos_sds) > 0) as.numeric(quantile(pos_sds, probs = floor_quantile)) else 1.0
    daily_stats <- daily_stats %>%
      mutate(sd_reg = pmax(sd_d, sigma_min, na.rm = TRUE))

    df_clamped <- df_clamped %>%
      left_join(daily_stats[, c("DOY", "mu_d", "sd_reg")], by = "DOY") %>%
      mutate(!!vc := (.data[[vc]] - mu_d) / sd_reg) %>%
      dplyr::select(-mu_d, -sd_reg)
  }
  return(df_clamped)
}

# Ensure therm_daily has integer DOY clamped to <= 365 for exact join safety
therm_daily <- therm_daily %>%
  mutate(DOY = pmin(365L, as.integer(DOY))) %>%
  group_by(Site, Year, DOY) %>%
  summarize(
    T_roll30    = mean(T_roll30, na.rm = TRUE),
    GDD5_roll60 = mean(GDD5_roll60, na.rm = TRUE),
    GDD5_YTD    = mean(GDD5_YTD, na.rm = TRUE),
    .groups     = "drop"
  )

#' Scan Trapping Time Series for Continuous Year-Round Sequences (>= 365 Days)
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

#' Catch-Weighted Quantile for Stage 2a Detection Threshold tau_det (Eq. 3)
wtd_quantile <- function(x, w, prob = 0.01) {
  ord <- order(x)
  x_s <- x[ord]
  w_s <- w[ord] / sum(w[ord])
  x_s[which(cumsum(w_s) >= prob)[1]]
}

#' Exact Unit Tweedie Deviance Function
tw_unit_dev <- function(y, mu, p) {
  mu    <- pmax(mu, 1e-15)
  term1 <- ifelse(y > 0, (y^(2 - p)) / ((1 - p) * (2 - p)), 0)
  term2 <- ifelse(y > 0, (y * (mu^(1 - p))) / (1 - p), 0)
  term3 <- (mu^(2 - p)) / (2 - p)
  2 * (term1 - term2 + term3)
}

#' Exact Closed-Form Empirical-Only Tweedie Null Deviance Explained (%)
calc_emp_dev_explained <- function(m, emp_idx, y_emp, effort_emp, p_override = NULL) {
  p_hat  <- if (is.null(p_override)) m$family$getTheta(TRUE) else p_override
  mu_emp <- pmax(as.numeric(fitted(m)[emp_idx]), 1e-15)

  c_hat    <- sum(y_emp * (effort_emp^(1 - p_hat))) / sum(effort_emp^(2 - p_hat))
  mu0_emp  <- pmax(c_hat * effort_emp, 1e-15)

  dev_mod  <- sum(tw_unit_dev(y_emp, mu_emp, p_hat))
  dev_null <- sum(tw_unit_dev(y_emp, mu0_emp, p_hat))

  100 * (1 - dev_mod / dev_null)
}

#' Extract Smooth Concurvity Estimates Safely
extract_concurvity <- function(m, drv_name = NULL) {
  conc <- tryCatch(concurvity(m, full = TRUE)["estimate", ], error = function(e) NULL)
  out  <- list(s_DOY = NA_real_, ti_Lat = NA_real_, ti_Drv = NA_real_)
  if (!is.null(conc)) {
    doy_idx <- grep("^s\\(DOY\\)", names(conc))
    lat_idx <- grep("^ti\\(.*Latitude", names(conc))
    if (length(doy_idx) > 0) out$s_DOY  <- round(conc[doy_idx[1]], 3)
    if (length(lat_idx) > 0) out$ti_Lat <- round(conc[lat_idx[1]], 3)
    if (!is.null(drv_name)) {
      drv_idx <- grep(paste0("^ti\\(.*", drv_name), names(conc))
      if (length(drv_idx) > 0) out$ti_Drv <- round(conc[drv_idx[1]], 3)
    }
  }
  return(out)
}

#' Point-in-Convex-Hull Check (2D Joint Empirical Covariate Hull H_2D)
in_convex_hull_2d <- function(pt_x, pt_y, hull_x, hull_y) {
  n <- length(hull_x)
  signs <- numeric(n)
  for (i in seq_len(n)) {
    j <- if (i == n) 1L else i + 1L
    signs[i] <- (hull_x[j] - hull_x[i]) * (pt_y - hull_y[i]) -
                (hull_y[j] - hull_y[i]) * (pt_x - hull_x[i])
  }
  all(signs <= 1e-9) || all(signs >= -1e-9)
}

# ==============================================================================
# MASTER WORKFLOW FUNCTION
# ==============================================================================
run_phenology_workflow <- function(sp_name, sp_short, raw_ecn_df,
                                   external_df = NULL, model_type = "dual",
                                   n_min = 10, m_imp = 50, out_fig_name) {

  cat(sprintf("\n==================================================================\n"))
  cat(sprintf("Executing Workflow for: %s (%s-tensor configuration)\n", sp_name, model_type))
  cat(sprintf("==================================================================\n"))

  lon_col <- intersect(names(site_coords),
                       c("Longitude", "Lon", "Long", "lon", "longitude", "Easting"))[1]

  # ----------------------------------------------------------------------------
  # STAGE 1: Effort Standardization, Leap-Year Clamping, & Filtering (Sec. 2.1)
  # ----------------------------------------------------------------------------
  dat <- raw_ecn_df %>%
    mutate(OpenDate  = ymd(OpenDate),
           CloseDate = ymd(CloseDate)) %>%
    mutate(Duration  = as.numeric(CloseDate - OpenDate),
           MidDate   = OpenDate + floor(Duration / 2),
           DOY       = pmin(365L, as.integer(yday(MidDate))),
           Year      = year(MidDate)) %>%
    filter(!is.na(Duration) & Duration > 0) %>%
    mutate(Effort    = NTraps * Duration,
           LogEffort = log(Effort),
           CPUE      = N / Effort) %>%
    left_join(therm_daily, by = c("Site", "Year", "DOY")) %>%
    left_join(site_coords, by = "Site") %>%
    filter(!is.na(Latitude) & !is.na(T_roll30) & !is.na(GDD5_roll60) & !is.na(GDD5_YTD))

  sy_grid <- dat %>%
    group_by(Site, SubSite, Year) %>%
    summarize(TotalN = sum(N), .groups = "drop") %>%
    mutate(FilterPass = TotalN >= n_min)

  df_filtered <- dat %>%
    inner_join(sy_grid %>% filter(FilterPass), by = c("Site", "SubSite", "Year"))
  sy_unique   <- unique(df_filtered[, c("SubSite", "Year")])
  n_sy        <- nrow(sy_unique)

  # Identify species-level systematic monitoring gap D_gap^spec (<5% active coverage)
  active_matrix <- matrix(FALSE, nrow = 365, ncol = n_sy)
  for (j in 1:n_sy) {
    df_j <- df_filtered %>% filter(SubSite == sy_unique$SubSite[j], Year == sy_unique$Year[j])
    for (r in 1:nrow(df_j)) {
      d_start <- min(365L, max(1L, as.integer(yday(df_j$OpenDate[r]))))
      d_end   <- min(365L, max(1L, as.integer(yday(df_j$CloseDate[r]))))
      if (year(df_j$OpenDate[r]) == year(df_j$CloseDate[r])) {
        if (d_start <= d_end) active_matrix[d_start:d_end, j] <- TRUE
      } else {
        active_matrix[d_start:365, j] <- TRUE
        active_matrix[1:d_end, j]     <- TRUE
      }
    }
  }

  gap_doys      <- which(rowMeans(active_matrix) < 0.05)
  breaks        <- c(0, which(diff(gap_doys) > 1), length(gap_doys))
  inferred_gaps <- lapply(1:(length(breaks) - 1), function(k) {
    chunk <- gap_doys[(breaks[k] + 1):breaks[k + 1]]
    data.frame(Start = min(chunk), End = max(chunk))
  })
  global_gaps     <- bind_rows(inferred_gaps)
  gap_summary_str <- paste(sprintf("DOY %d-%d", global_gaps$Start, global_gaps$End), collapse = "; ")

  # Check residual unmonitored gaps G_{res,l,y} outside D_gap^spec
  non_gap_doys <- setdiff(1:365, gap_doys)
  max_res_gaps <- apply(active_matrix[non_gap_doys, , drop = FALSE], 2, function(col_act) {
    rle_obj <- rle(!col_act)
    if (!any(rle_obj$values)) return(0L)
    as.integer(max(rle_obj$lengths[rle_obj$values]))
  })
  pct_res_le60 <- round(100 * mean(max_res_gaps <= 60), 1)
  med_res_gap  <- round(median(max_res_gaps), 1)

  cat(sprintf("1. Module 1 Pre-Processing: Retained %d site-years (N >= %d).\n", n_sy, n_min))
  cat(sprintf("   Species-level systematic gap D_gap^spec: %s\n", gap_summary_str))
  cat(sprintf("   Residual gap check outside D_gap^spec: %.1f%% of site-years have G_res <= 60 d (median = %.0f d).\n",
              pct_res_le60, med_res_gap))

  # ----------------------------------------------------------------------------
  # STAGE 2a: Whole-Year Baseline Phenological Template (Section 2.2, Stage 2a)
  # ----------------------------------------------------------------------------
  cands_internal <- find_continuous_runs(dat)
  qual_internal  <- Filter(function(x) x$TotalN >= 10, cands_internal)

  if (length(qual_internal) > 0) {
    df_base <- bind_rows(lapply(seq_along(qual_internal), function(idx) {
      qual_internal[[idx]]$Data %>% mutate(Run_ID = paste0("Run_", idx))
    })) %>%
      mutate(SubSite_fac = as.factor(SubSite),
             Run_fac     = as.factor(Run_ID))
    baseline_source <- sprintf(
      "Internal ECN: %d continuous run(s) across %s (N_total = %d; %s to %s)",
      length(qual_internal),
      paste(unique(df_base$SubSite), collapse = ", "),
      sum(df_base$N),
      as.character(min(sapply(qual_internal, function(x) as.character(x$Start)))),
      as.character(max(sapply(qual_internal, function(x) as.character(x$End))))
    )
    cat(sprintf("2. Built internal baseline template -> %s\n", baseline_source))
  } else {
    wij_prep <- external_df %>%
      mutate(StartDate = ymd(StartDate),
             EndDate   = ymd(EndDate)) %>%
      mutate(TrapDays  = as.numeric(EndDate - StartDate),
             MidDate   = StartDate + floor(TrapDays / 2),
             DOY       = pmin(365L, as.integer(yday(MidDate))),
             Year      = year(MidDate),
             Effort    = TrapDays,
             LogEffort = log(Effort),
             SubSite   = "Wijster_Station")
    cands_ext <- find_continuous_runs(wij_prep, start_col = "StartDate",
                                      end_col = "EndDate", subsite_col = "SubSite", n_col = "N")
    qual_ext  <- Filter(function(x) x$TotalN >= 10, cands_ext)
    df_base   <- bind_rows(lapply(seq_along(qual_ext), function(idx) {
      qual_ext[[idx]]$Data %>% mutate(Run_ID = paste0("Run_", idx))
    })) %>%
      mutate(SubSite_fac = as.factor(SubSite),
             Run_fac     = as.factor(Run_ID))
    baseline_source <- sprintf(
      "External Historical (Wijster Station): %d continuous run(s) (N_total = %d; %s to %s)",
      length(qual_ext),
      sum(df_base$N),
      as.character(min(sapply(qual_ext, function(x) as.character(x$Start)))),
      as.character(max(sapply(qual_ext, function(x) as.character(x$End))))
    )
    cat(sprintf("2. Built external baseline template -> %s\n", baseline_source))
  }

  re_base_terms <- c()
  if (nlevels(df_base$SubSite_fac) > 1) re_base_terms <- c(re_base_terms, "s(SubSite_fac, bs = 're')")
  if (nlevels(df_base$Run_fac) > 1)     re_base_terms <- c(re_base_terms, "s(Run_fac, bs = 're')")
  rhs_base     <- paste(c("s(DOY, bs = 'cc', k = 15)", re_base_terms, "offset(LogEffort)"), collapse = " + ")
  baseline_gam <- gam(as.formula(paste("N ~", rhs_base)),
                      data = df_base, family = tw(), method = "REML",
                      gamma = 1.4, knots = list(DOY = c(0.5, 365.5)))

  pred_grid <- data.frame(
    DOY         = 1:365,
    LogEffort   = 0,
    SubSite_fac = df_base$SubSite_fac[1],
    Run_fac     = df_base$Run_fac[1]
  )
  X_base_365 <- predict(baseline_gam, newdata = pred_grid, type = "lpmatrix")
  re_cols_b  <- grep("s\\(Run_fac\\)|s\\(SubSite_fac\\)", colnames(X_base_365))
  if (length(re_cols_b) > 0) X_base_365[, re_cols_b] <- 0

  eta_base <- as.numeric(X_base_365 %*% coef(baseline_gam))
  V_base   <- suppressWarnings(vcov(baseline_gam, unconditional = TRUE))
  se_base  <- sqrt(pmax(0, rowSums((X_base_365 %*% V_base) * X_base_365)))

  pred_grid$fit_raw <- exp(eta_base)
  pred_grid$lwr_raw <- exp(eta_base - 1.96 * se_base)
  pred_grid$upr_raw <- exp(eta_base + 1.96 * se_base)

  # 2nd catch-weighted percentile detection threshold tau_det (Equation 3)
  df_pos <- df_base %>% filter(N > 0)
  nd_pos <- data.frame(
    DOY         = df_pos$DOY,
    LogEffort   = 0,
    SubSite_fac = df_base$SubSite_fac[1],
    Run_fac     = df_base$Run_fac[1]
  )
  X_pos <- predict(baseline_gam, newdata = nd_pos, type = "lpmatrix")
  if (length(re_cols_b) > 0) X_pos[, re_cols_b] <- 0
  pos_lambda     <- as.numeric(exp(X_pos %*% coef(baseline_gam)))
  phys_threshold <- wtd_quantile(pos_lambda, df_pos$N, prob = 0.02)

  # Continuous soft-thresholding at tau_det prevents step discontinuities at gap edges
  pred_grid$fit    <- pmax(0, pred_grid$fit_raw - phys_threshold)
  template_density <- pred_grid$fit / sum(pred_grid$fit)
  active_mask      <- pred_grid$fit > 0

  # Draw M_imp = 50 posterior template realizations on the active support
  set.seed(101)
  beta_base_draws <- mgcv::rmvn(m_imp, coef(baseline_gam), V_base)
  template_draws  <- apply(beta_base_draws, 1, function(b_vec) {
    lam_d   <- as.numeric(exp(X_base_365 %*% b_vec))
    pos_d   <- as.numeric(exp(X_pos %*% b_vec))
    tau_d   <- wtd_quantile(pos_d, df_pos$N, prob = 0.02)
    lam_tr  <- ifelse(active_mask, pmax(1e-8, lam_d - tau_d), 0)
    if (sum(lam_tr) > 0) lam_tr / sum(lam_tr) else template_density
  })

  # ----------------------------------------------------------------------------
  # STAGE 2b: Locally Scaled Posterior Integer Sampling Across D_gap^spec
  # ----------------------------------------------------------------------------
  build_augmented_df <- function(tmpl_vec, integer_draw = FALSE) {
    syn_list <- list()
    for (i in 1:nrow(sy_grid)) {
      if (!sy_grid$FilterPass[i]) next
      sub   <- sy_grid$SubSite[i]
      yr    <- sy_grid$Year[i]
      df_sy <- df_filtered %>% filter(SubSite == sub, Year == yr)

      p_temp <- tmpl_vec[df_sy$DOY]
      if (sum(p_temp^2) == 0) next
      k_scale         <- sum(df_sy$CPUE * p_temp) / sum(p_temp^2)
      scaled_template <- tmpl_vec * k_scale

      global_gaps_sy <- global_gaps %>% filter(!(Start >= min(df_sy$DOY) & End <= max(df_sy$DOY)))
      if (nrow(global_gaps_sy) == 0) next

      for (g in 1:nrow(global_gaps_sy)) {
        g_start <- global_gaps_sy$Start[g]
        g_end   <- global_gaps_sy$End[g]
        if (g_end < g_start + 7) next

        syn_doys <- seq(g_start + 7, g_end, by = 14)
        ntraps   <- mean(df_sy$NTraps, na.rm = TRUE)
        syn_eff  <- ntraps * 14
        mu_syn   <- scaled_template[syn_doys] * syn_eff
        n_syn    <- if (integer_draw) ifelse(mu_syn > 0, rpois(length(mu_syn), lambda = mu_syn), 0L) else mu_syn

        syn_list[[length(syn_list) + 1]] <- data.frame(
          Site      = df_sy$Site[1],
          SubSite   = sub,
          Year      = yr,
          DOY       = syn_doys,
          N         = n_syn,
          Effort    = syn_eff,
          LogEffort = log(syn_eff),
          Weight    = 0.5,
          Latitude  = df_sy$Latitude[1]
        )
      }
    }

    df_syn <- bind_rows(syn_list) %>%
      left_join(therm_daily, by = c("Site", "Year", "DOY")) %>%
      filter(!is.na(T_roll30) & !is.na(GDD5_roll60) & !is.na(GDD5_YTD))

    cols_keep <- c("Site", "SubSite", "Year", "DOY", "N", "Effort", "LogEffort",
                   "Weight", "Latitude", "T_roll30", "GDD5_roll60", "GDD5_YTD")

    bind_rows(
      df_filtered %>% mutate(Weight = 1.0) %>% dplyr::select(all_of(cols_keep)),
      df_syn %>% dplyr::select(all_of(cols_keep))
    ) %>%
      mutate(SubSite_fac = as.factor(SubSite),
             Year_fac    = as.factor(Year))
  }

  df_final   <- build_augmented_df(template_density)
  emp_idx    <- which(df_final$Weight == 1.0)
  df_emp     <- df_final[emp_idx, ]
  j_emp      <- length(emp_idx)
  j_prior    <- sum(df_final$Weight == 0.5)
  n_sites    <- length(unique(df_emp$Site))
  n_subsites <- length(unique(df_emp$SubSite))
  yr_min     <- min(df_emp$Year)
  yr_max     <- max(df_emp$Year)
  yr_range   <- sprintf("%d-%d (%d yrs)", yr_min, yr_max, yr_max - yr_min + 1)

  cat(sprintf("3. Dataset summary: %s | %d sites (%d subsites) | %d site-years | J_emp = %d, J_prior = %d.\n",
              yr_range, n_sites, n_subsites, n_sy, j_emp, j_prior))

  # ----------------------------------------------------------------------------
  # STAGE 4a: Candidate Predictors + Benchmark Null Models (Table S2)
  # ----------------------------------------------------------------------------
  candidate_drivers <- c("T_roll30", "GDD5_roll60", "GDD5_YTD")
  model_specs       <- list()

  for (drv in candidate_drivers) {
    if (model_type == "dual") {
      model_specs[[drv]] <- as.formula(paste0(
        "N ~ s(DOY, bs = 'cc', k = 15) + s(Latitude, k = 3) + s(", drv, ", k = 4) + ",
        "ti(DOY, Latitude, bs = c('cc', 'cr'), k = c(12, 3)) + ",
        "ti(DOY, ", drv, ", bs = c('cc', 'cr'), k = c(12, 4)) + ",
        "s(SubSite_fac, bs = 're') + s(Year_fac, bs = 're') + offset(LogEffort)"
      ))
    } else {
      model_specs[[drv]] <- as.formula(paste0(
        "N ~ s(DOY, bs = 'cc', k = 15) + s(", drv, ", k = 4) + ",
        "ti(DOY, ", drv, ", bs = c('cc', 'cr'), k = c(12, 4)) + ",
        "s(SubSite_fac, bs = 're') + s(Year_fac, bs = 're') + offset(LogEffort)"
      ))
    }
  }

  if (model_type == "dual") {
    model_specs[["Spatial_Only"]] <- as.formula(
      "N ~ s(DOY, bs = 'cc', k = 15) + s(Latitude, k = 3) + " %>%
        paste0("ti(DOY, Latitude, bs = c('cc', 'cr'), k = c(12, 3)) + ",
               "s(SubSite_fac, bs = 're') + s(Year_fac, bs = 're') + offset(LogEffort)")
    )
  }

  model_specs[["Seasonal_Baseline_Null"]] <- as.formula(
    "N ~ s(DOY, bs = 'cc', k = 15) + s(SubSite_fac, bs = 're') + s(Year_fac, bs = 're') + offset(LogEffort)"
  )

  model_fits <- list()
  for (m_name in names(model_specs)) {
    model_fits[[m_name]] <- bam(
      model_specs[[m_name]], data = df_final, weights = Weight, family = tw(),
      method = "fREML", discrete = TRUE, gamma = 1.4,
      knots = list(DOY = c(0.5, 365.5))
    )
  }

  # Calibrate empirical dispersion phi_hat_emp strictly on empirical records (w_j = 1.0)
  m_ref_disp   <- model_fits[["GDD5_YTD"]]
  p_ref_emp    <- m_ref_disp$family$getTheta(TRUE)
  mu_ref_emp   <- pmax(as.numeric(fitted(m_ref_disp)[emp_idx]), 1e-15)
  phi_emp_dev  <- sum(tw_unit_dev(df_emp$N, mu_ref_emp, p_ref_emp)) / (j_emp - sum(m_ref_disp$edf))
  phi_emp_pear <- sum(((df_emp$N - mu_ref_emp)^2) / (mu_ref_emp^p_ref_emp)) / (j_emp - sum(m_ref_disp$edf))
  phi_emp_cal  <- 0.5 * (phi_emp_dev + phi_emp_pear)
  c_hat_ratio  <- phi_emp_dev / m_ref_disp$sig2

  aic_records <- list()
  for (m_name in names(model_specs)) {
    m       <- model_fits[[m_name]]
    sm      <- summary(m)
    mu_e    <- pmax(as.numeric(fitted(m)[emp_idx]), 1e-15)
    dev_emp <- calc_emp_dev_explained(m, emp_idx, df_emp$N, df_emp$Effort, p_override = p_ref_emp)
    ll_emp  <- sum(mgcv::ldTweedie(df_emp$N, mu = mu_e, p = p_ref_emp, phi = phi_emp_dev)[, 1])
    aic_emp <- -2 * ll_emp + 2 * (sum(m$edf) + 2)
    drv_arg <- if (m_name %in% candidate_drivers) m_name else NULL
    conc    <- extract_concurvity(m, drv_name = drv_arg)

    aic_records[[m_name]] <- data.frame(
      Species              = sp_name,
      ModelStructure       = paste0(tools::toTitleCase(model_type), "-tensor"),
      Year_Range           = yr_range,
      n_Sites              = n_sites,
      n_SiteYears          = n_sy,
      J_emp                = j_emp,
      J_prior              = j_prior,
      CandidateModel       = m_name,
      IsWeatherCandidate   = m_name %in% candidate_drivers,
      Tweedie_p            = round(m$family$getTheta(TRUE), 3),
      Tweedie_phi_Comb     = round(m$sig2, 3),
      Tweedie_phi_Emp      = round(phi_emp_dev, 3),
      EDF                  = round(sum(m$edf), 2),
      LogLik_Emp           = round(ll_emp, 2),
      AIC_Emp              = round(aic_emp, 2),
      LogLik_Comb          = round(as.numeric(logLik(m)), 2),
      AIC_Comb             = round(AIC(m), 2),
      DevExpl_Total_pct    = round(sm$dev.expl * 100, 2),
      DevExpl_Emp_pct      = round(dev_emp, 2),
      Concurvity_s_DOY     = conc$s_DOY,
      Concurvity_ti_Lat    = conc$ti_Lat,
      Concurvity_ti_Drv    = conc$ti_Drv,
      stringsAsFactors     = FALSE
    )
  }

  aic_table    <- bind_rows(aic_records)
  null_dev_emp <- aic_table$DevExpl_Emp_pct[aic_table$CandidateModel == "Seasonal_Baseline_Null"]

  cand_sub      <- aic_table %>% filter(IsWeatherCandidate) %>% arrange(AIC_Emp)
  champ_driver  <- cand_sub$CandidateModel[1]
  champ_aic_emp <- cand_sub$AIC_Emp[1]
  champ_aic_cmb <- cand_sub$AIC_Comb[1]
  mod_final     <- model_fits[[champ_driver]]

  aic_table <- aic_table %>%
    mutate(
      DeltaAIC_Emp          = round(AIC_Emp  - champ_aic_emp, 2),
      DeltaAIC_Comb         = round(AIC_Comb - champ_aic_cmb, 2),
      Delta_DevExpl_Emp_pct = round(DevExpl_Emp_pct - null_dev_emp, 2),
      SortGroup             = ifelse(IsWeatherCandidate, 1,
                                     ifelse(CandidateModel == "Spatial_Only", 2, 3))
    ) %>%
    arrange(SortGroup, AIC_Emp) %>%
    dplyr::select(-IsWeatherCandidate, -SortGroup)

  cat("\n4a. Candidate Predictor & Null Benchmark Comparison (Table S2):\n")
  print(aic_table %>% dplyr::select(CandidateModel, EDF, LogLik_Emp, DeltaAIC_Emp,
                                    DeltaAIC_Comb, DevExpl_Emp_pct, Delta_DevExpl_Emp_pct,
                                    Concurvity_ti_Lat, Concurvity_ti_Drv))

  saveRDS(mod_final, paste0("Champion_GAMM_", sp_short, ".rds"))

  # ----------------------------------------------------------------------------
  # STAGE 4b: Upfront Collinearity, Calibrated DHARMa, & Autocorrelation Checks
  # ----------------------------------------------------------------------------
  cat("\n4b. Running Model Validation Diagnostics on Empirical Observations (Weight == 1.0)...\n")

  r_lat_t30  <- round(cor(df_emp$Latitude, df_emp$T_roll30), 3)
  r_lat_g60  <- round(cor(df_emp$Latitude, df_emp$GDD5_roll60), 3)
  r_lat_gytd <- round(cor(df_emp$Latitude, df_emp$GDD5_YTD), 3)
  r_lat_chmp <- round(cor(df_emp$Latitude, df_emp[[champ_driver]]), 3)

  k_tab <- tryCatch(mgcv::k.check(mod_final), error = function(e) NULL)
  min_k_idx <- if (!is.null(k_tab) && "k-index" %in% colnames(k_tab)) {
    round(min(k_tab[, "k-index"], na.rm = TRUE), 3)
  } else NA_real_

  p_hat        <- mod_final$family$getTheta(TRUE)
  phi_hat_comb <- mod_final$sig2
  mu_emp       <- as.numeric(fitted(mod_final)[emp_idx])

  set.seed(42)
  sim_emp <- replicate(1000, mgcv::rTweedie(mu = mu_emp, p = p_hat, phi = phi_emp_cal))
  dharma_emp <- DHARMa::createDHARMa(
    simulatedResponse       = sim_emp,
    observedResponse        = df_emp$N,
    fittedPredictedResponse = mu_emp,
    method                  = "PIT",
    integerResponse         = TRUE
  )

  ks_res   <- DHARMa::testUniformity(dharma_emp, plot = FALSE)
  disp_res <- DHARMa::testDispersion(dharma_emp, plot = FALSE)
  zi_res   <- DHARMa::testZeroInflation(dharma_emp, plot = FALSE)

  u_pit <- pmin(pmax(dharma_emp$scaledResiduals, 1e-6), 1 - 1e-6)
  df_emp_diag <- df_emp %>%
    mutate(u_pit = u_pit,
           z_res = qnorm(u_pit)) %>%
    arrange(SubSite, Year, DOY) %>%
    group_by(SubSite, Year) %>%
    mutate(z_prev   = dplyr::lag(z_res),
           doy_prev = dplyr::lag(DOY),
           dt_days  = DOY - doy_prev) %>%
    ungroup()

  consec_pairs <- df_emp_diag %>% filter(!is.na(z_prev) & dt_days > 0 & dt_days <= 21)
  if (nrow(consec_pairs) >= 10) {
    acf_test <- cor.test(consec_pairs$z_prev, consec_pairs$z_res)
    lag1_r   <- round(as.numeric(acf_test$estimate), 3)
    lag1_p   <- signif(acf_test$p.value, 3)
  } else {
    lag1_r <- NA_real_; lag1_p <- NA_real_
  }

  df_emp_time <- df_emp %>% mutate(TimeIndex = Year * 366 + DOY)
  dharma_time <- tryCatch({
    res_t <- DHARMa::recalculateResiduals(dharma_emp, group = df_emp_time$TimeIndex)
    t_ord <- order(as.numeric(levels(as.factor(df_emp_time$TimeIndex))))
    DHARMa::testTemporalAutocorrelation(res_t, time = seq_along(t_ord), plot = FALSE)
  }, error = function(e) NULL)
  dw_stat <- if (!is.null(dharma_time)) round(as.numeric(dharma_time$statistic), 3) else NA_real_
  dw_p    <- if (!is.null(dharma_time)) signif(dharma_time$p.value, 3) else NA_real_

  spat_res <- NULL
  if (n_sites >= 4) {
    spat_res <- tryCatch({
      res_s      <- DHARMa::recalculateResiduals(dharma_emp, group = factor(df_emp$Site))
      site_order <- levels(factor(df_emp$Site))
      coords_sub <- site_coords[match(site_order, site_coords$Site), ]
      x_coord    <- if (!is.na(lon_col)) coords_sub[[lon_col]] else rep(0, nrow(coords_sub))
      y_coord    <- coords_sub$Latitude
      DHARMa::testSpatialAutocorrelation(res_s, x = x_coord, y = y_coord, plot = FALSE)
    }, error = function(e) NULL)
  }
  moran_i <- if (!is.null(spat_res)) round(as.numeric(spat_res$statistic["observed"]), 3) else NA_real_
  moran_p <- if (!is.null(spat_res)) signif(spat_res$p.value, 3) else NA_real_

  # ----------------------------------------------------------------------------
  # STAGE 5: Calibrated Posterior Simulation on the Unit Simplex & Sensitivity
  # ----------------------------------------------------------------------------
  n_sims          <- 1000
  set.seed(11)
  beta            <- coef(mod_final)
  has_uncond_vcov <- !is.null(mod_final$Vc)
  vcov_type_str   <- ifelse(has_uncond_vcov, "Calibrated Unconditional (V_c^cal)", "Calibrated Conditional (V_p^cal)")
  V_raw           <- suppressWarnings(vcov(mod_final, unconditional = TRUE))
  V_cal           <- V_raw * max(1.0, c_hat_ratio)

  has_pos_priors <- any(df_final$Weight == 0.5 & df_final$N > 0)
  if (has_pos_priors && m_imp > 1) {
    n_sub_imp <- min(m_imp, 20)
    beta_imp_mat <- matrix(NA_real_, nrow = n_sub_imp, ncol = length(beta))
    for (m_i in seq_len(n_sub_imp)) {
      df_m <- build_augmented_df(template_draws[, m_i], integer_draw = TRUE)
      fit_m <- tryCatch(
        bam(model_specs[[champ_driver]], data = df_m, weights = Weight, family = tw(),
            method = "fREML", discrete = TRUE, gamma = 1.4, knots = list(DOY = c(0.5, 365.5))),
        error = function(e) NULL
      )
      if (!is.null(fit_m) && length(coef(fit_m)) == length(beta)) {
        beta_imp_mat[m_i, ] <- coef(fit_m)
      }
    }
    valid_imp <- complete.cases(beta_imp_mat)
    if (sum(valid_imp) >= 3) {
      B_imp <- cov(beta_imp_mat[valid_imp, , drop = FALSE])
      V_cal <- V_cal + (1 + 1 / sum(valid_imp)) * B_imp
    }
  }

  br <- mgcv::rmvn(n_sims, beta, V_cal)

  extract_posterior_density <- function(mod_obj, br_mat, drv_val, lat_val) {
    nd <- data.frame(DOY = 1:365, LogEffort = 0,
                     SubSite_fac = df_final$SubSite_fac[1],
                     Year_fac    = df_final$Year_fac[1],
                     Latitude    = lat_val)
    nd[[champ_driver]] <- drv_val

    X <- predict(mod_obj, newdata = nd, type = "lpmatrix")
    re_cols <- grep("s\\(Year_fac\\)|s\\(SubSite_fac\\)", colnames(X))
    if (length(re_cols) > 0) X[, re_cols] <- 0

    # Plug-in MLE simplex curve (avoids Jensen's inequality log-normal variance bump)
    mu_mle   <- as.numeric(exp(X %*% coef(mod_obj)))
    dens_mle <- if (sum(mu_mle) > 0) mu_mle / sum(mu_mle) else rep(1/365, 365)

    sim_val  <- exp(X %*% t(br_mat))
    sim_dens <- apply(sim_val, 2, function(x) if (sum(x) > 0) x / sum(x) else rep(1/365, 365))

    data.frame(
      DOY = 1:365,
      fit = dens_mle,
      lwr = apply(sim_dens, 1, quantile, probs = 0.025),
      upr = apply(sim_dens, 1, quantile, probs = 0.975)
    )
  }

  lat_mean <- mean(df_emp$Latitude, na.rm = TRUE)
  lat_med  <- median(df_emp$Latitude, na.rm = TRUE)
  drv_5    <- as.numeric(quantile(df_emp[[champ_driver]], 0.05, na.rm = TRUE))
  drv_95   <- as.numeric(quantile(df_emp[[champ_driver]], 0.95, na.rm = TRUE))
  drv_mean <- 0.0

  d_cold <- extract_posterior_density(mod_final, br, drv_5,  lat_med) %>% mutate(Scenario = "Cold Regime")
  d_warm <- extract_posterior_density(mod_final, br, drv_95, lat_med) %>% mutate(Scenario = "Warm Regime")

  d_south <- NULL
  d_north <- NULL
  lat_5   <- as.numeric(quantile(df_emp$Latitude, 0.05, na.rm = TRUE))
  lat_95  <- as.numeric(quantile(df_emp$Latitude, 0.95, na.rm = TRUE))
  if (model_type == "dual") {
    d_south <- extract_posterior_density(mod_final, br, drv_mean, lat_5)  %>% mutate(Scenario = "South")
    d_north <- extract_posterior_density(mod_final, br, drv_mean, lat_95) %>% mutate(Scenario = "North")
  }

  # 2D Joint Empirical Covariate Hull Check (Comment 7)
  h_idx_obs <- grDevices::chull(df_emp$Latitude, df_emp[[champ_driver]])
  hx_obs    <- df_emp$Latitude[h_idx_obs]
  hy_obs    <- df_emp[[champ_driver]][h_idx_obs]

  sy_centroids <- df_emp %>%
    group_by(Site, SubSite, Year, Latitude) %>%
    summarize(Drv_Mean = mean(.data[[champ_driver]], na.rm = TRUE), .groups = "drop")
  h_idx_sy <- grDevices::chull(sy_centroids$Latitude, sy_centroids$Drv_Mean)
  hx_sy    <- sy_centroids$Latitude[h_idx_sy]
  hy_sy    <- sy_centroids$Drv_Mean[h_idx_sy]

  query_pts <- if (model_type == "dual") {
    data.frame(
      Species  = sp_name,
      Scenario = c("South", "North", "Cold Regime", "Warm Regime"),
      Column   = c("Latitudinal gradient (W* = 0)", "Latitudinal gradient (W* = 0)",
                   "Thermal plasticity (median latitude)", "Thermal plasticity (median latitude)"),
      Latitude = c(lat_5, lat_95, lat_med, lat_med),
      W_star   = c(drv_mean, drv_mean, drv_5, drv_95),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      Species  = sp_name,
      Scenario = c("Cold Regime", "Warm Regime"),
      Column   = c("Thermal plasticity (median latitude)", "Thermal plasticity (median latitude)"),
      Latitude = c(lat_med, lat_med),
      W_star   = c(drv_5, drv_95),
      stringsAsFactors = FALSE
    )
  }
  query_pts$In_Obs_Hull <- mapply(in_convex_hull_2d, query_pts$Latitude, query_pts$W_star,
                                  MoreArgs = list(hull_x = hx_obs, hull_y = hy_obs))
  query_pts$In_SY_Hull  <- mapply(in_convex_hull_2d, query_pts$Latitude, query_pts$W_star,
                                  MoreArgs = list(hull_x = hx_sy, hull_y = hy_sy))

  ar1_noa_pct <- NA_real_
  if (!is.na(lag1_r) && lag1_r > 0.05) {
    ar1_noa_pct <- tryCatch({
      df_ar <- df_final %>%
        arrange(SubSite, Year, DOY) %>%
        group_by(SubSite, Year) %>%
        mutate(AR_start = row_number() == 1 | (DOY - dplyr::lag(DOY, default = -99L)) > 21) %>%
        ungroup()
      m_ar1 <- suppressWarnings(bam(
        model_specs[[champ_driver]], data = df_ar, weights = Weight, family = tw(),
        method = "fREML", discrete = TRUE, gamma = 1.4,
        rho = lag1_r, AR.start = df_ar$AR_start,
        knots = list(DOY = c(0.5, 365.5))
      ))
      set.seed(11)
      br_ar1 <- mgcv::rmvn(n_sims, coef(m_ar1), suppressWarnings(vcov(m_ar1, unconditional = TRUE)) * max(1.0, c_hat_ratio))
      dc_ar1 <- extract_posterior_density(m_ar1, br_ar1, drv_5,  lat_med)
      dw_ar1 <- extract_posterior_density(m_ar1, br_ar1, drv_95, lat_med)
      noa_c  <- 50 * sum(abs(d_cold$fit - dc_ar1$fit))
      noa_w  <- 50 * sum(abs(d_warm$fit - dw_ar1$fit))
      round(max(noa_c, noa_w), 2)
    }, error = function(e) NA_real_)
  }

  gamma10_min_ov <- tryCatch({
    m_g10 <- suppressWarnings(bam(
      model_specs[[champ_driver]], data = df_final, weights = Weight, family = tw(),
      method = "fREML", discrete = TRUE, gamma = 1.0,
      knots = list(DOY = c(0.5, 365.5))
    ))
    set.seed(11)
    br_g10 <- mgcv::rmvn(200, coef(m_g10), suppressWarnings(vcov(m_g10, unconditional = TRUE)))
    dc_10  <- extract_posterior_density(m_g10, br_g10, drv_5,  lat_med)
    dw_10  <- extract_posterior_density(m_g10, br_g10, drv_95, lat_med)
    ov_c   <- 100 - 50 * sum(abs(d_cold$fit - dc_10$fit))
    ov_w   <- 100 - 50 * sum(abs(d_warm$fit - dw_10$fit))
    round(min(ov_c, ov_w), 2)
  }, error = function(e) NA_real_)

  diag_row <- data.frame(
    Species                  = sp_name,
    Champion_Predictor       = champ_driver,
    Model_Type               = paste0(tools::toTitleCase(model_type), "-tensor"),
    Year_Range               = yr_range,
    n_Sites                  = n_sites,
    n_SubSites               = n_subsites,
    n_SiteYears              = n_sy,
    J_emp                    = j_emp,
    J_prior                  = j_prior,
    Baseline_Template_Source = baseline_source,
    Tau_Det_1stPct_UnitEff   = signif(phys_threshold, 4),
    Flagged_Gap_Windows      = gap_summary_str,
    Pct_ResGap_le60d         = pct_res_le60,
    Median_ResGap_d          = med_res_gap,
    r_Lat_Champ_Emp          = r_lat_chmp,
    r_Lat_Troll30_Emp        = r_lat_t30,
    r_Lat_GDD60_Emp          = r_lat_g60,
    r_Lat_GDDYTD_Emp         = r_lat_gytd,
    Tweedie_p                = round(p_hat, 3),
    Tweedie_phi_Comb         = round(phi_hat_comb, 3),
    Tweedie_phi_Emp_Dev      = round(phi_emp_dev, 3),
    Tweedie_phi_Emp_Cal      = round(phi_emp_cal, 3),
    Vcov_Type_Used           = vcov_type_str,
    Min_k_index              = min_k_idx,
    DHARMa_KS_D              = round(as.numeric(ks_res$statistic), 3),
    DHARMa_KS_p              = signif(ks_res$p.value, 3),
    DHARMa_Dispersion_Ratio  = round(as.numeric(disp_res$statistic), 3),
    DHARMa_Dispersion_p      = signif(disp_res$p.value, 3),
    DHARMa_ZeroInf_Ratio     = round(as.numeric(zi_res$statistic), 3),
    DHARMa_ZeroInf_p         = signif(zi_res$p.value, 3),
    Temp_Lag1_ACF_Emp        = lag1_r,
    Temp_Lag1_p              = lag1_p,
    Temp_Consec_Pairs        = nrow(consec_pairs),
    DHARMa_Temp_DW           = dw_stat,
    DHARMa_Temp_p            = dw_p,
    AR1_Sensitivity_MaxNOA   = ar1_noa_pct,
    Gamma10_MinOverlap_pct   = gamma10_min_ov,
    Spatial_Moran_I          = moran_i,
    Spatial_Moran_p          = moran_p,
    stringsAsFactors         = FALSE
  )

  cat("   --- Diagnostic Summary (Empirical Records Only, Calibrated Dispersion) ---\n")
  cat(sprintf("   Collinearity r(Lat, %s) = %.3f | Tweedie p = %.3f, phi_comb = %.3f, phi_emp = %.3f\n",
              champ_driver, r_lat_chmp, p_hat, phi_hat_comb, phi_emp_cal))
  cat(sprintf("   DHARMa KS D = %.3f (p = %s) | Calibrated Dispersion = %.3f (p = %s) | ZeroRatio = %.3f\n",
              diag_row$DHARMa_KS_D, diag_row$DHARMa_KS_p,
              diag_row$DHARMa_Dispersion_Ratio, diag_row$DHARMa_Dispersion_p,
              diag_row$DHARMa_ZeroInf_Ratio))
  cat(sprintf("   Lag-1 ACF = %s (AR1 max NOA = %s%%) | Gamma=1.0 vs 1.4 min overlap = %s%% | Moran I = %s (p = %s)\n",
              as.character(lag1_r), as.character(ar1_noa_pct), as.character(gamma10_min_ov),
              ifelse(is.na(moran_i), "NA", as.character(moran_i)), as.character(moran_p)))

  qq_df <- data.frame(
    Expected  = seq_len(j_emp) / (j_emp + 1),
    Observed  = sort(dharma_emp$scaledResiduals),
    PredRank  = rank(mu_emp) / j_emp,
    ScaledRes = dharma_emp$scaledResiduals,
    DOY       = df_emp$DOY
  )

  diag_theme <- theme_classic(base_size = 9.5) +
    theme(
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      plot.title        = element_blank(),
      plot.subtitle     = element_blank(),
      plot.tag          = element_text(size = 10, face = "bold", color = "black"),
      plot.tag.position = "topleft"
    )

  p_d1 <- ggplot(qq_df, aes(x = Expected, y = Observed)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#D55E00") +
    geom_point(size = 0.7, alpha = 0.5, color = "#0072B2") +
    labs(tag = "a)", x = "Expected uniform quantile", y = "Observed PIT residual") +
    diag_theme

  p_d2 <- ggplot(qq_df, aes(x = PredRank, y = ScaledRes)) +
    geom_hline(yintercept = c(0.25, 0.50, 0.75), linetype = "dotted", color = "grey50") +
    geom_point(size = 0.7, alpha = 0.35, color = "grey25") +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, color = "#0072B2", linewidth = 0.7) +
    labs(tag = "b)", x = "Rank-transformed predicted mean", y = "PIT residual") +
    diag_theme

  p_d3 <- ggplot(consec_pairs, aes(x = z_prev, y = z_res)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
    geom_point(size = 0.7, alpha = 0.35, color = "grey25") +
    geom_smooth(method = "lm", formula = y ~ x, color = "#D55E00", linewidth = 0.7) +
    labs(tag = "c)", x = "Quantile residual (t - 1)", y = "Quantile residual (t)") +
    diag_theme

  p_d4 <- ggplot(qq_df, aes(x = DOY, y = ScaledRes)) +
    geom_hline(yintercept = c(0.25, 0.50, 0.75), linetype = "dotted", color = "grey50") +
    geom_point(size = 0.7, alpha = 0.35, color = "grey25") +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, color = "#0072B2", linewidth = 0.7) +
    labs(tag = "d)", x = "Day of year", y = "PIT residual") +
    diag_theme

  p_diag <- (p_d1 + p_d2) / (p_d3 + p_d4)
  diag_fig_name <- paste0("Diagnostics_", sp_short, ".png")
  ggsave(diag_fig_name, p_diag, width = 17, height = 13, units = "cm", dpi = 300)

  return(list(
    aic_table    = aic_table,
    diagnostics  = diag_row,
    curves       = list(d_south = d_south, d_north = d_north, d_cold = d_cold, d_warm = d_warm),
    stage2_info  = list(
      dat              = dat,
      sy_grid          = sy_grid,
      df_filtered      = df_filtered,
      df_final         = df_final,
      df_emp           = df_emp,
      global_gaps      = global_gaps,
      pred_grid        = pred_grid,
      phys_threshold   = phys_threshold,
      template_density = template_density,
      sy_centroids     = sy_centroids,
      hull_obs         = data.frame(Latitude = hx_obs, W_star = hy_obs),
      hull_sy          = data.frame(Latitude = hx_sy,  W_star = hy_sy),
      query_pts        = query_pts,
      bbox             = list(lat_5 = lat_5, lat_95 = lat_95, drv_5 = drv_5, drv_95 = drv_95)
    )
  ))
}

# ==============================================================================
# EXECUTE PIPELINE FOR ALL THREE CASE-STUDY SPECIES
# ==============================================================================
res_pm <- run_phenology_workflow(
  sp_name      = "Pterostichus madidus",
  sp_short     = "PM",
  raw_ecn_df   = ecn_species_data$Pterostichus_madidus,
  model_type   = "dual",
  out_fig_name = "Figure_4_Pterostichus_madidus_Plasticity.png"
)

res_nb <- run_phenology_workflow(
  sp_name      = "Notiophilus biguttatus",
  sp_short     = "NB",
  raw_ecn_df   = ecn_species_data$Notiophilus_biguttatus,
  external_df  = wijster_nb,
  model_type   = "dual",
  out_fig_name = "Figure_6_Notiophilus_biguttatus_Plasticity.png"
)

res_cm <- run_phenology_workflow(
  sp_name      = "Calathus micropterus",
  sp_short     = "CM",
  raw_ecn_df   = ecn_species_data$Calathus_micropterus,
  model_type   = "single",
  out_fig_name = "Figure_7_Calathus_micropterus_Plasticity.png"
)

# ==============================================================================
# GENERATE FIGURE 3: PTEROSTICHUS MADIDUS PRE-PROCESSING & STAGE 2a/2b (Comment 2)
# ==============================================================================
cat("\nGenerating Figure 3 (Pterostichus madidus pre-processing & D_gap^spec alignment)...\n")
pm_info <- res_pm$stage2_info
gaps_pm <- pm_info$global_gaps # DOY 1-76 and 332-365

draw_key_cust_3a <- function(data, params, size) {
  f_val <- if (!is.null(data$fill)) data$fill else "grey78"
  if (f_val == "grey92") {
    grid::grobTree(
      grid::rectGrob(gp = grid::gpar(col = "grey25", fill = "grey92", lwd = 0.6)),
      grid::segmentsGrob(x0 = grid::unit(c(0, 0, 0.35), "npc"),
                         y0 = grid::unit(c(0.35, 0, 0), "npc"),
                         x1 = grid::unit(c(0.65, 1, 1), "npc"),
                         y1 = grid::unit(c(1, 1, 0.65), "npc"),
                         gp = grid::gpar(col = "grey20", lwd = 0.65))
    )
  } else if (f_val == "#CC3311") {
    grid::pointsGrob(x = 0.5, y = 0.5, pch = 16,
                     gp = grid::gpar(col = "#CC3311", cex = 0.85))
  } else {
    grid::rectGrob(gp = grid::gpar(col = "grey25", fill = f_val, lwd = 0.6))
  }
}

pub_theme <- theme_classic(base_size = 9.5) +
  theme(
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    plot.title        = element_blank(),
    plot.subtitle     = element_blank(),
    plot.tag          = element_text(size = 10, face = "bold", color = "black"),
    plot.tag.position = "topleft",
    legend.title      = element_blank(),
    legend.text       = element_text(size = 8.0),
    axis.text         = element_text(color = "black", size = 8.5),
    axis.title        = element_text(color = "black", size = 9.5),
    plot.margin       = margin(6, 8, 4, 8)
  )

all_sy_pm <- pm_info$dat %>%
  group_by(Site, SubSite, Year) %>%
  summarize(TotalN = sum(N), .groups = "drop") %>%
  mutate(Pass = TotalN >= 10) %>%
  arrange(desc(Site), SubSite, Year) %>%
  mutate(RowID = row_number())

n_rows_pm <- nrow(all_sy_pm)
site_breaks <- all_sy_pm %>%
  group_by(Site) %>%
  summarize(MidRow = mean(RowID), .groups = "drop")

df_seg_list <- list()
dat_pm_joined <- pm_info$dat %>%
  inner_join(all_sy_pm[, c("Site", "SubSite", "Year", "RowID", "Pass")],
             by = c("Site", "SubSite", "Year"))

for (r in seq_len(nrow(dat_pm_joined))) {
  d_s  <- min(365L, max(1L, as.integer(yday(dat_pm_joined$OpenDate[r]))))
  d_e  <- min(365L, max(1L, as.integer(yday(dat_pm_joined$CloseDate[r]))))
  yr_s <- year(dat_pm_joined$OpenDate[r])
  yr_e <- year(dat_pm_joined$CloseDate[r])
  if (yr_s == yr_e && d_s <= d_e) {
    df_seg_list[[length(df_seg_list) + 1]] <- data.frame(
      RowID = dat_pm_joined$RowID[r], d_start = d_s, d_end = d_e,
      N = dat_pm_joined$N[r], CPUE = min(6, dat_pm_joined$CPUE[r])
    )
  } else {
    df_seg_list[[length(df_seg_list) + 1]] <- data.frame(
      RowID = dat_pm_joined$RowID[r], d_start = d_s, d_end = 365L,
      N = dat_pm_joined$N[r], CPUE = min(6, dat_pm_joined$CPUE[r])
    )
    df_seg_list[[length(df_seg_list) + 1]] <- data.frame(
      RowID = dat_pm_joined$RowID[r], d_start = 1L, d_end = d_e,
      N = dat_pm_joined$N[r], CPUE = min(6, dat_pm_joined$CPUE[r])
    )
  }
}
df_segs_pm <- bind_rows(df_seg_list)
df_zero_pm <- df_segs_pm %>% filter(N == 0)
df_pos_seg <- df_segs_pm %>% filter(N > 0)

bar_ymin <- -16
bar_ymax <- -3
bar_h    <- bar_ymax - bar_ymin
hatch_list <- list()
for (g in seq_len(nrow(gaps_pm))) {
  xs <- gaps_pm$Start[g]
  xe <- gaps_pm$End[g]
  c_seq <- seq(xs - 14, xe, by = 7)
  for (c_val in c_seq) {
    x1 <- max(xs, c_val)
    x2 <- min(xe, c_val + 12)
    if (x2 > x1) {
      y1 <- bar_ymin + (x1 - c_val) / 12 * bar_h
      y2 <- bar_ymin + (x2 - c_val) / 12 * bar_h
      hatch_list[[length(hatch_list) + 1]] <- data.frame(x1 = x1, y1 = y1, x2 = x2, y2 = y2)
    }
  }
}
df_hatch <- bind_rows(hatch_list)

df_leg_fill <- data.frame(
  x = c(10, 10, 10, 10), y = c(10, 10, 10, 10),
  Status = factor(c("Data gap", "Zero catch", "Flagged gap", ">= 10 individuals"),
                  levels = c("Data gap", "Zero catch", "Flagged gap", ">= 10 individuals"))
)

p_3a <- ggplot() +
  geom_tile(data = df_leg_fill, aes(x = x, y = y, fill = Status),
            width = 0, height = 0, key_glyph = draw_key_cust_3a) +
  annotate("rect", xmin = 1, xmax = 365, ymin = 0.5, ymax = n_rows_pm + 0.5,
           fill = "grey78", color = NA) +
  geom_segment(data = df_zero_pm,
               aes(x = d_start, xend = d_end, y = RowID, yend = RowID),
               color = "#FFFFFF", linewidth = 0.52) +
  geom_segment(data = df_pos_seg,
               aes(x = d_start, xend = d_end, y = RowID, yend = RowID, color = CPUE),
               linewidth = 0.52) +
  geom_rect(data = gaps_pm,
            aes(xmin = Start, xmax = End, ymin = bar_ymin, ymax = bar_ymax),
            fill = "grey92", color = "grey20", linewidth = 0.35) +
  geom_segment(data = df_hatch,
               aes(x = x1, y = y1, xend = x2, yend = y2),
               color = "grey20", linewidth = 0.32) +
  geom_point(data = all_sy_pm %>% filter(Pass),
             aes(x = 372, y = RowID),
             color = "#CC3311", size = 0.85, shape = 16) +
  scale_color_gradient(low = "#D4B9DA", high = "#4A1486", name = "CPUE",
                       breaks = c(0, 2, 4, 6), limits = c(0, 6),
                       guide = guide_colorbar(order = 1, barwidth = grid::unit(2.0, "cm"),
                                              barheight = grid::unit(0.32, "cm"),
                                              title.position = "top", title.hjust = 0.5)) +
  scale_fill_manual(name = NULL,
                    values = c("Data gap" = "grey78",
                               "Zero catch" = "#FFFFFF",
                               "Flagged gap" = "grey92",
                               ">= 10 individuals" = "#CC3311"),
                    guide = guide_legend(order = 2, nrow = 2, byrow = TRUE)) +
  scale_y_continuous(breaks = site_breaks$MidRow, labels = site_breaks$Site,
                     limits = c(bar_ymin - 2, n_rows_pm + 2), expand = c(0.01, 0)) +
  scale_x_continuous(breaks = c(0, 100, 200, 300), limits = c(0, 378), expand = c(0.01, 0)) +
  labs(tag = "a)", x = "Day of year", y = NULL) +
  pub_theme +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = element_text(size = 8.0, face = "bold"),
        legend.text = element_text(size = 7.8),
        legend.key.size = grid::unit(0.36, "cm"),
        legend.spacing.x = grid::unit(0.14, "cm"),
        legend.box.spacing = grid::unit(0.15, "cm"),
        legend.margin = margin(0, 4, 0, 4))

df_3b_lines <- bind_rows(
  data.frame(DOY = pm_info$pred_grid$DOY, y = pm_info$pred_grid$fit_raw, Series = "Template"),
  data.frame(DOY = c(1, 365), y = rep(pm_info$phys_threshold, 2), Series = "Zero-threshold")
) %>% mutate(Series = factor(Series, levels = c("Template", "Zero-threshold")))

p_3b <- ggplot() +
  geom_ribbon(data = pm_info$pred_grid, aes(x = DOY, ymin = lwr_raw, ymax = upr_raw),
              fill = "#0072B2", alpha = 0.25) +
  geom_line(data = df_3b_lines, aes(x = DOY, y = y, color = Series, linetype = Series),
            linewidth = 0.75) +
  scale_color_manual(name = NULL, values = c("Template" = "#0072B2", "Zero-threshold" = "#CC3311")) +
  scale_linetype_manual(name = NULL, values = c("Template" = "solid", "Zero-threshold" = "22")) +
  scale_x_continuous(breaks = c(0, 100, 200, 300), limits = c(0, 365), expand = c(0.01, 0)) +
  labs(tag = "b)", x = "Day of year", y = "Standardized CPUE") +
  pub_theme +
  theme(legend.position = c(0.04, 0.96),
        legend.justification = c(0, 1),
        legend.key.width = grid::unit(0.85, "cm"),
        legend.key.height = grid::unit(0.38, "cm"),
        legend.text = element_text(size = 8.0),
        legend.background = element_blank(),
        legend.box.background = element_rect(fill = "white", color = "grey50", linewidth = 0.3),
        legend.box.margin = margin(3, 5, 3, 5))

ex_cand <- pm_info$df_filtered %>%
  group_by(Site, SubSite, Year) %>%
  summarize(min_d = min(DOY), max_d = max(DOY), max_cpue = max(CPUE), n_obs = n(), .groups = "drop") %>%
  filter(min_d > 76 & max_d < 332 & n_obs >= 12) %>%
  arrange(desc(max_cpue))

ex_sub <- ex_cand$SubSite[1]
ex_yr  <- ex_cand$Year[1]
df_ex_emp <- pm_info$df_filtered %>% filter(SubSite == ex_sub, Year == ex_yr)

p_temp_ex  <- pm_info$template_density[df_ex_emp$DOY]
k_ex       <- sum(df_ex_emp$CPUE * p_temp_ex) / sum(p_temp_ex^2)
scale_mult <- k_ex / sum(pm_info$pred_grid$fit)

curve_ex <- data.frame(
  DOY  = 1:365,
  CPUE = pm_info$pred_grid$fit * scale_mult,
  lwr  = pmax(0, (pm_info$pred_grid$lwr_raw - pm_info$phys_threshold) * scale_mult),
  upr  = pmax(0.38, pm_info$pred_grid$upr_raw * scale_mult)
)

syn_doys_pm <- c(seq(1 + 7, 76, by = 14), seq(332 + 7, 365, by = 14))
q_fracs <- c(0.28, 0.72, 0.16, 0.58, 0.82, 0.75, 0.38)
cpue_imp_draws <- curve_ex$CPUE[syn_doys_pm] + q_fracs * (curve_ex$upr[syn_doys_pm] - curve_ex$CPUE[syn_doys_pm])
df_ex_syn <- data.frame(DOY = syn_doys_pm, CPUE = cpue_imp_draws)

df_pts_3c <- bind_rows(
  df_ex_emp %>% dplyr::select(DOY, CPUE) %>% mutate(PointType = "Raw data"),
  df_ex_syn %>% dplyr::select(DOY, CPUE) %>% mutate(PointType = "Imputation")
) %>% mutate(PointType = factor(PointType, levels = c("Raw data", "Imputation")))

p_3c <- ggplot() +
  geom_rect(data = gaps_pm,
            aes(xmin = Start, xmax = End, ymin = -Inf, ymax = Inf, fill = "Flagged gap"),
            alpha = 0.25) +
  geom_ribbon(data = curve_ex, aes(x = DOY, ymin = lwr, ymax = upr),
              fill = "#0072B2", alpha = 0.28) +
  geom_line(data = curve_ex, aes(x = DOY, y = CPUE, linetype = "Template"),
            color = "#0072B2", linewidth = 0.75) +
  geom_point(data = df_pts_3c, aes(x = DOY, y = CPUE, shape = PointType, color = PointType),
             size = 1.9, stroke = 0.75) +
  scale_fill_manual(name = NULL, values = c("Flagged gap" = "grey65"),
                    guide = guide_legend(order = 1)) +
  scale_linetype_manual(name = NULL, values = c("Template" = "solid"),
                        guide = guide_legend(order = 2)) +
  scale_shape_manual(name = NULL, values = c("Raw data" = 16, "Imputation" = 2),
                     guide = guide_legend(order = 3)) +
  scale_color_manual(name = NULL, values = c("Raw data" = "black", "Imputation" = "#CC3311"),
                     guide = guide_legend(order = 3)) +
  scale_x_continuous(breaks = c(0, 100, 200, 300), limits = c(0, 365), expand = c(0.01, 0)) +
  labs(tag = "c)", x = "Day of year", y = "Standardized CPUE") +
  pub_theme +
  theme(legend.position = c(0.04, 0.96),
        legend.justification = c(0, 1),
        legend.key.width = grid::unit(0.85, "cm"),
        legend.key.height = grid::unit(0.36, "cm"),
        legend.spacing.y = grid::unit(0.01, "cm"),
        legend.margin = margin(1, 4, 1, 4),
        legend.text = element_text(size = 8.0),
        legend.background = element_blank(),
        legend.box.background = element_rect(fill = "white", color = "grey50", linewidth = 0.3),
        legend.box.margin = margin(3, 5, 3, 5))

fig3_combined <- p_3a + (p_3b / p_3c) + plot_layout(widths = c(1.22, 1.0))
ggsave("Figure_3_Pterostichus_madidus_Pipeline.png", fig3_combined, width = 18.2, height = 14.8, units = "cm", dpi = 300)
ggsave("Figure_3_Pterostichus_madidus_Pipeline.pdf", fig3_combined, width = 18.2, height = 14.8, units = "cm")
cat("Saved Figure_3_Pterostichus_madidus_Pipeline.png & .pdf\n")

# ==============================================================================
# GENERATE FIGURE 4: NOTIOPHILUS BIGUTTATUS EXTERNAL BASELINE TRANSFER
# ==============================================================================
cat("Generating Figure 4 (Notiophilus biguttatus external Wijster baseline transfer)...\n")
nb_info <- res_nb$stage2_info
gaps_nb <- nb_info$global_gaps # DOY 1-76 and 316-365

p_4a <- ggplot() +
  geom_rect(data = gaps_nb,
            aes(xmin = Start, xmax = End, ymin = -Inf, ymax = Inf, fill = "Flagged gap"),
            alpha = 0.25) +
  geom_point(data = nb_info$df_filtered, aes(x = DOY, y = CPUE, color = "Raw data"),
             alpha = 0.45, size = 1.4) +
  scale_fill_manual(name = NULL, values = c("Flagged gap" = "grey65"),
                    guide = guide_legend(order = 1)) +
  scale_color_manual(name = NULL, values = c("Raw data" = "#0072B2"),
                     guide = guide_legend(order = 2)) +
  scale_x_continuous(breaks = c(0, 100, 200, 300), limits = c(0, 365), expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0, 0.16), expand = c(0.01, 0)) +
  labs(tag = "a)", x = "Day of year", y = "Standardized CPUE") +
  pub_theme +
  theme(legend.position = c(0.96, 0.96),
        legend.justification = c(1, 1),
        legend.key.width = grid::unit(0.75, "cm"),
        legend.key.height = grid::unit(0.36, "cm"),
        legend.spacing.y = grid::unit(0.01, "cm"),
        legend.margin = margin(1, 4, 1, 4),
        legend.background = element_blank(),
        legend.box.background = element_rect(fill = "white", color = "grey50", linewidth = 0.3),
        legend.box.margin = margin(3, 5, 3, 5))

df_4b_lines <- bind_rows(
  data.frame(DOY = nb_info$pred_grid$DOY, y = nb_info$pred_grid$fit_raw, Series = "Template"),
  data.frame(DOY = c(1, 365), y = rep(nb_info$phys_threshold, 2), Series = "Zero-threshold")
) %>% mutate(Series = factor(Series, levels = c("Template", "Zero-threshold")))

p_4b <- ggplot() +
  geom_ribbon(data = nb_info$pred_grid, aes(x = DOY, ymin = lwr_raw, ymax = upr_raw),
              fill = "#0072B2", alpha = 0.25) +
  geom_line(data = df_4b_lines, aes(x = DOY, y = y, color = Series, linetype = Series),
            linewidth = 0.75) +
  scale_color_manual(name = NULL, values = c("Template" = "#0072B2", "Zero-threshold" = "#CC3311")) +
  scale_linetype_manual(name = NULL, values = c("Template" = "solid", "Zero-threshold" = "22")) +
  scale_x_continuous(breaks = c(0, 100, 200, 300), limits = c(0, 365), expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0, 11.2), breaks = c(0, 3, 6, 9), expand = c(0.01, 0)) +
  labs(tag = "b)", x = "Day of year", y = "Standardized CPUE") +
  pub_theme +
  theme(legend.position = c(0.96, 0.96),
        legend.justification = c(1, 1),
        legend.key.width = grid::unit(0.85, "cm"),
        legend.key.height = grid::unit(0.36, "cm"),
        legend.margin = margin(1, 4, 1, 4),
        legend.background = element_blank(),
        legend.box.background = element_rect(fill = "white", color = "grey50", linewidth = 0.3),
        legend.box.margin = margin(3, 5, 3, 5))

gam_imp_nb <- gam(N ~ s(DOY, bs = "cc", k = 15) + s(SubSite_fac, bs = "re") + s(Year_fac, bs = "re") + offset(LogEffort),
                  data = nb_info$df_final, weights = Weight, family = tw(), method = "REML", gamma = 1.4,
                  knots = list(DOY = c(0.5, 365.5)))
gam_nai_nb <- gam(N ~ s(DOY, bs = "cc", k = 15) + s(SubSite_fac, bs = "re") + s(Year_fac, bs = "re") + offset(LogEffort),
                  data = nb_info$df_filtered %>% mutate(SubSite_fac = as.factor(SubSite), Year_fac = as.factor(Year)),
                  family = tw(), method = "REML", gamma = 1.4,
                  knots = list(DOY = c(0.5, 365.5)))

pred_nb_4c <- function(g_obj, label_str) {
  nd <- data.frame(DOY = 1:365, LogEffort = 0,
                   SubSite_fac = nb_info$df_final$SubSite_fac[1],
                   Year_fac    = nb_info$df_final$Year_fac[1])
  X <- predict(g_obj, newdata = nd, type = "lpmatrix")
  rc <- grep("s\\(Year_fac\\)|s\\(SubSite_fac\\)", colnames(X))
  if (length(rc) > 0) X[, rc] <- 0
  eta <- as.numeric(X %*% coef(g_obj))
  V   <- suppressWarnings(vcov(g_obj, unconditional = TRUE))
  se  <- sqrt(pmax(0, rowSums((X %*% V) * X)))
  data.frame(DOY = 1:365, fit = exp(eta), lwr = exp(eta - 1.96 * se), upr = exp(eta + 1.96 * se), Model = label_str)
}
df_4c <- bind_rows(
  pred_nb_4c(gam_imp_nb, "2-stage regularized"),
  pred_nb_4c(gam_nai_nb, "Unanchored")
) %>% mutate(Model = factor(Model, levels = c("2-stage regularized", "Unanchored")))

p_4c <- ggplot(df_4c, aes(x = DOY)) +
  geom_rect(data = gaps_nb, inherit.aes = FALSE,
            aes(xmin = Start, xmax = End, ymin = -Inf, ymax = Inf, fill = "Flagged gap"),
            alpha = 0.25) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = Model), alpha = 0.25) +
  geom_line(aes(y = fit, color = Model), linewidth = 0.8) +
  scale_fill_manual(name  = NULL,
                    values = c("2-stage regularized" = "#009E73", "Unanchored" = "#CC79A7", "Flagged gap" = "grey65"),
                    breaks = c("Flagged gap"),
                    guide = guide_legend(order = 1, override.aes = list(color = NA, alpha = 0.35))) +
  scale_color_manual(name = NULL,
                     values = c("2-stage regularized" = "#009E73", "Unanchored" = "#CC79A7"),
                     guide = guide_legend(order = 2)) +
  scale_x_continuous(breaks = c(0, 100, 200, 300), limits = c(0, 365), expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0, 0.038), breaks = c(0, 0.01, 0.02, 0.03), expand = c(0.01, 0)) +
  labs(tag = "c)", x = "Day of year", y = "Standardized CPUE") +
  pub_theme +
  theme(legend.position = c(0.96, 0.96),
        legend.justification = c(1, 1),
        legend.key.width = grid::unit(0.75, "cm"),
        legend.key.height = grid::unit(0.36, "cm"),
        legend.spacing.y = grid::unit(0.01, "cm"),
        legend.margin = margin(1, 4, 1, 4),
        legend.background = element_blank(),
        legend.box.background = element_rect(fill = "white", color = "grey50", linewidth = 0.3),
        legend.box.margin = margin(3, 5, 3, 5))

fig4_combined <- p_4a / p_4b / p_4c
ggsave("Figure_4_Notiophilus_biguttatus_Transfer.png", fig4_combined, width = 9.8, height = 19.5, units = "cm", dpi = 300)
ggsave("Figure_4_Notiophilus_biguttatus_Transfer.pdf", fig4_combined, width = 9.8, height = 19.5, units = "cm")
cat("Saved Figure_4_Notiophilus_biguttatus_Transfer.png & .pdf\n")

# ==============================================================================
# GENERATE FIGURE 5: COMBINED 5-PANEL PLASTICITY FIGURE (Comments 1, 4, 7)
# ==============================================================================
cat("Generating Figure 5 (5-panel plasticity figure with separated palettes & calibrated CIs)...\n")

empirical_curves <- list(
  PM = res_pm$curves,
  NB = res_nb$curves,
  CM = res_cm$curves
)
save(empirical_curves, file = "Empirical_Plasticity_Curves.RData")

extract_legend_grob <- function(p) {
  gt  <- ggplotGrob(p)
  idx <- which(sapply(gt$grobs, function(x) x$name) == "guide-box")
  gt$grobs[[idx[1]]]
}

fig5_theme <- pub_theme +
  theme(
    legend.title     = element_text(size = 8.5, face = "bold", color = "black"),
    legend.text      = element_text(size = 8.2),
    legend.key.width = unit(0.55, "cm"),
    plot.margin      = margin(4, 10, 4, 8)
  )

pal_lat   <- c("North" = "#009E73", "South" = "#AA3377")
pal_therm <- c("Cold Regime" = "#0072B2", "Warm Regime" = "#D55E00")

make_lat_panel <- function(curves_obj, tag_label, ymax_val, show_xlab = FALSE) {
  df_p <- bind_rows(curves_obj$d_north, curves_obj$d_south) %>%
    mutate(Scenario = factor(Scenario, levels = c("North", "South")))
  ggplot(df_p, aes(x = DOY)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr, fill = Scenario), alpha = 0.22) +
    geom_line(aes(y = fit, color = Scenario), linewidth = 0.8) +
    scale_color_manual(name = "Latitudinal gradient (W* = 0):", values = pal_lat) +
    scale_fill_manual(name  = "Latitudinal gradient (W* = 0):", values = pal_lat) +
    scale_x_continuous(breaks = c(0, 100, 200, 300), limits = c(0, 365), expand = c(0.01, 0)) +
    scale_y_continuous(limits = c(0, ymax_val), breaks = seq(0, 0.025, by = 0.005),
                       labels = function(x) sprintf("%.3f", x), expand = c(0.01, 0)) +
    labs(tag = tag_label, x = if (show_xlab) "Day of year" else NULL, y = "Probability density") +
    fig5_theme +
    theme(legend.position = "none")
}

make_therm_panel <- function(curves_obj, tag_label, ymax_val, show_xlab = FALSE, show_ylab = FALSE) {
  df_p <- bind_rows(curves_obj$d_cold, curves_obj$d_warm) %>%
    mutate(Scenario = factor(Scenario, levels = c("Cold Regime", "Warm Regime")))
  ggplot(df_p, aes(x = DOY)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr, fill = Scenario), alpha = 0.22) +
    geom_line(aes(y = fit, color = Scenario), linewidth = 0.8) +
    scale_color_manual(name = "Thermal plasticity (median lat.):", values = pal_therm) +
    scale_fill_manual(name  = "Thermal plasticity (median lat.):", values = pal_therm) +
    scale_x_continuous(breaks = c(0, 100, 200, 300), limits = c(0, 365), expand = c(0.01, 0)) +
    scale_y_continuous(limits = c(0, ymax_val), breaks = seq(0, 0.025, by = 0.005),
                       labels = function(x) sprintf("%.3f", x), expand = c(0.01, 0)) +
    labs(tag = tag_label, x = if (show_xlab) "Day of year" else NULL,
         y = if (show_ylab) "Probability density" else NULL) +
    fig5_theme +
    theme(legend.position = "none")
}

ymax_pm <- max(c(empirical_curves$PM$d_north$upr, empirical_curves$PM$d_south$upr,
                 empirical_curves$PM$d_cold$upr,  empirical_curves$PM$d_warm$upr), na.rm = TRUE) * 1.05
ymax_nb <- max(c(empirical_curves$NB$d_north$upr, empirical_curves$NB$d_south$upr,
                 empirical_curves$NB$d_cold$upr,  empirical_curves$NB$d_warm$upr), na.rm = TRUE) * 1.05
ymax_cm <- max(c(empirical_curves$CM$d_cold$upr,  empirical_curves$CM$d_warm$upr), na.rm = TRUE) * 1.05

p_a <- make_lat_panel(empirical_curves$PM,   tag_label = "a)", ymax_val = ymax_pm, show_xlab = FALSE)
p_b <- make_therm_panel(empirical_curves$PM, tag_label = "b)", ymax_val = ymax_pm, show_xlab = FALSE, show_ylab = FALSE)
p_c <- make_lat_panel(empirical_curves$NB,   tag_label = "c)", ymax_val = ymax_nb, show_xlab = TRUE)
p_d <- make_therm_panel(empirical_curves$NB, tag_label = "d)", ymax_val = ymax_nb, show_xlab = FALSE, show_ylab = FALSE)
p_e <- make_therm_panel(empirical_curves$CM, tag_label = "e)", ymax_val = ymax_cm, show_xlab = TRUE,  show_ylab = TRUE)

p_a_leg <- p_a +
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5),
         fill  = guide_legend(title.position = "top", title.hjust = 0.5)) +
  theme(legend.position = "bottom")

p_b_leg <- p_b +
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5),
         fill  = guide_legend(title.position = "top", title.hjust = 0.5)) +
  theme(legend.position = "bottom")

leg_lat_top <- wrap_elements(full = grid::grobTree(
  grid::editGrob(extract_legend_grob(p_a_leg),
                 vp = grid::viewport(x = 0.54, y = 0.88))
))
leg_thm <- wrap_elements(full = extract_legend_grob(p_b_leg))

design_clean <- "
AB
CD
EF
GH
"
fig_combined <- p_a + p_b +
                p_c + p_d +
                leg_lat_top + p_e +
                plot_spacer() + leg_thm +
                plot_layout(
                  design  = design_clean,
                  heights = c(1, 1, 1, 0.16)
                )

ggsave("Figure_5_Combined_Plasticity.png",     fig_combined, width = 16.5, height = 19.5, units = "cm", dpi = 300)
ggsave("Figure_5_Combined_Plasticity.pdf",     fig_combined, width = 16.5, height = 19.5, units = "cm")
ggsave("Figure_4_6_7_Combined_Plasticity.png", fig_combined, width = 16.5, height = 19.5, units = "cm", dpi = 300)
ggsave("Figure_4_6_7_Combined_Plasticity.pdf", fig_combined, width = 16.5, height = 19.5, units = "cm")
cat("Saved Figure_5_Combined_Plasticity.png & .pdf\n")

# ==============================================================================
# GENERATE SUPPLEMENTARY FIGURE S4: 2D JOINT COVARIATE HULL VS. MARGINAL BOX
# ==============================================================================
cat("Generating Supplementary Figure S4 (2D joint empirical covariate hull vs marginal bounding box)...\n")

make_hull_panel <- function(res_obj, tag_label, show_xlab = FALSE, show_ylab = TRUE) {
  s_info <- res_obj$stage2_info
  df_obs <- s_info$df_emp
  df_sy  <- s_info$sy_centroids
  h_obs  <- s_info$hull_obs
  h_sy   <- s_info$hull_sy
  q_pts  <- s_info$query_pts %>%
    mutate(Scenario = factor(Scenario, levels = c("North", "South", "Cold Regime", "Warm Regime")))
  bb     <- s_info$bbox

  rect_df <- data.frame(xmin = bb$lat_5, xmax = bb$lat_95, ymin = bb$drv_5, ymax = bb$drv_95)
  all_pal <- c("North" = "#009E73", "South" = "#AA3377",
               "Cold Regime" = "#0072B2", "Warm Regime" = "#D55E00")
  all_shp <- c("North" = 17, "South" = 15, "Cold Regime" = 16, "Warm Regime" = 18)

  ggplot() +
    geom_polygon(data = h_obs, aes(x = Latitude, y = W_star),
                 fill = "grey85", color = "grey45", linewidth = 0.5, alpha = 0.45) +
    geom_polygon(data = h_sy, aes(x = Latitude, y = W_star),
                 fill = "grey70", color = "grey30", linewidth = 0.45, linetype = "dotted", alpha = 0.35) +
    geom_rect(data = rect_df,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = NA, color = "black", linetype = "dashed", linewidth = 0.55) +
    geom_point(data = df_obs, aes(x = Latitude, y = GDD5_YTD),
               color = "grey55", alpha = 0.25, size = 0.7) +
    geom_point(data = df_sy, aes(x = Latitude, y = Drv_Mean),
               color = "grey20", alpha = 0.75, size = 1.5, shape = 1) +
    geom_point(data = q_pts, aes(x = Latitude, y = W_star, color = Scenario, shape = Scenario),
               size = 3.1, stroke = 0.9) +
    scale_color_manual(name = "Figure 5 query coordinates:", values = all_pal, drop = FALSE) +
    scale_shape_manual(name = "Figure 5 query coordinates:", values = all_shp, drop = FALSE) +
    labs(tag = tag_label,
         x = if (show_xlab) "Latitude (°N)" else NULL,
         y = if (show_ylab) "Standardized heat (GDD*_YTD)" else NULL) +
    pub_theme +
    theme(legend.position = "none",
          plot.margin = margin(10, 8, 4, 8))
}

p_s4a <- make_hull_panel(res_pm, "a)", show_xlab = TRUE, show_ylab = TRUE)
p_s4b <- make_hull_panel(res_nb, "b)", show_xlab = TRUE, show_ylab = FALSE)
p_s4c <- make_hull_panel(res_cm, "c)", show_xlab = TRUE, show_ylab = FALSE)

leg_s4 <- wrap_elements(full = extract_legend_grob(
  p_s4a + theme(legend.position = "bottom",
                legend.title = element_text(size = 8.5, face = "bold"))
))

fig_s4 <- (p_s4a + p_s4b + p_s4c) / leg_s4 +
  plot_layout(heights = c(1, 0.14))

ggsave("Figure_S4_Covariate_Hull.png", fig_s4, width = 18.5, height = 8.6, units = "cm", dpi = 300)
ggsave("Figure_S4_Covariate_Hull.pdf", fig_s4, width = 18.5, height = 8.6, units = "cm")
cat("Saved Figure_S4_Covariate_Hull.png & .pdf\n")

# ==============================================================================
# EXPORT SUMMARY TABLES (Table S2 & Table S3)
# ==============================================================================
table_s2_all <- bind_rows(res_pm$aic_table, res_nb$aic_table, res_cm$aic_table)
table_s3_all <- bind_rows(res_pm$diagnostics, res_nb$diagnostics, res_cm$diagnostics)

write.csv(table_s2_all, "Table_S2_Model_Selection_Summary.csv",      row.names = FALSE)
write.csv(table_s3_all, "Table_S3_Model_Validation_Diagnostics.csv", row.names = FALSE)
write.csv(table_s2_all, "Table_2_Model_Selection_Summary.csv",       row.names = FALSE)
write.csv(table_s3_all, "Table_S2_Model_Validation_Diagnostics.csv", row.names = FALSE)

cat("\n==================================================================\n")
cat("Pipeline v0.3.0 execution completed!\n")
cat("Exported: 'Figure_3_Pterostichus_madidus_Pipeline.png' & '.pdf'\n")
cat("Exported: 'Figure_4_Notiophilus_biguttatus_Transfer.png' & '.pdf'\n")
cat("Exported: 'Figure_5_Combined_Plasticity.png' & '.pdf'\n")
cat("Exported: 'Figure_S4_Covariate_Hull.png' & '.pdf'\n")
cat("Exported: 'Table_S2_Model_Selection_Summary.csv'\n")
cat("Exported: 'Table_S3_Model_Validation_Diagnostics.csv'\n")
cat("==================================================================\n")
