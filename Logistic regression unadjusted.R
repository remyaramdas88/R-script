
# Packages

library(readr)
library(dplyr)
library(lubridate)
library(janitor)

# Helpers

virus_01 <- function(x) {
  x <- suppressWarnings(as.integer(x))
  x[x == 999] <- NA_integer_
  ifelse(x %in% c(0L, 1L), x, NA_integer_)
}

parse_date_any <- function(x) {
  if (inherits(x, "Date")) return(as.Date(x))
  y <- suppressWarnings(dmy(x))
  if (all(is.na(y))) y <- suppressWarnings(ymd(x))
  as.Date(y)
}

ppm_to_ppb <- function(x) {
  # If value < 1, treat as ppm and convert; otherwise assume already ppb
  ifelse(is.na(x), NA_real_, ifelse(x < 1, x * 1000, x))
}

fit_and_show <- function(formula, data, label) {
  cat("\n====================\n", label, "\n====================\n", sep = "")
  m <- glm(formula, data = data, family = binomial, na.action = na.exclude)
  
  # GLM coefficient table
  print(summary(m)$coefficients)
  
  # OR + 95% Wald CI (no profiling)
  cf <- coef(m)
  se <- sqrt(diag(vcov(m)))
  OR  <- exp(cf)
  LCL <- exp(cf - 1.96 * se)
  UCL <- exp(cf + 1.96 * se)
  
  out <- data.frame(OR = OR, `2.5 %` = LCL, `97.5 %` = UCL, check.names = FALSE)
  print(out)
  
  invisible(m)
}


climate_path   <- "R studio/Patient climate data 4dayAvg new.csv"
pollution_path <- "R studio/Patient pollution data New.csv"

# 1) READ + CLEAN

climate <- read.csv(climate_path, check.names = FALSE) %>%
  clean_names() %>%
  mutate(
    # standardize outcomes to lowercase
    rsv      = virus_01(rsv),
    hrv      = virus_01(hrv),
    adm_date = parse_date_any(admission_date),
    
    tmax_daily = suppressWarnings(as.numeric(max_temp)),
    tmin_daily = suppressWarnings(as.numeric(min_temp)),
    rain_daily = suppressWarnings(as.numeric(rainfall)),
    
    tmax_4day = suppressWarnings(as.numeric(mean_max_temp_4_day_average)),
    tmin_4day = suppressWarnings(as.numeric(mean_min_temp_4_day_average)),
    rain_4day = suppressWarnings(as.numeric(mean_rainfall_4_day_average))
  ) %>%
  mutate(across(
    c(tmax_daily, tmin_daily, rain_daily, tmax_4day, tmin_4day, rain_4day),
    ~ ifelse(. == 999, NA, .)
  )) %>%
  filter(
    adm_date >= as.Date("2002-01-01"),
    adm_date <= as.Date("2022-12-31")
  )

# --- POLLUTION ---

pollution <- read.csv(pollution_path, check.names = FALSE) %>%
  clean_names() %>%
  mutate(
    # standardize outcomes to lowercase to match climate
    rsv      = virus_01(rsv),
    hrv      = virus_01(hrv),
    adm_date = parse_date_any(admission_date),
    
    # daily raw
    pm25_daily = suppressWarnings(as.numeric(pm2_5)),
    no_daily   = suppressWarnings(as.numeric(no)),   # ppm in your file
    no2_daily  = suppressWarnings(as.numeric(no2)),  # ppm
    
    # 4-day raw
    pm25_4 = suppressWarnings(as.numeric(pm2_5_4_day_average)),
    no_4   = suppressWarnings(as.numeric(no_4_day_average)),   # ppm
    no2_4  = suppressWarnings(as.numeric(no2_4_day_average))   # ppm
  ) %>%
  mutate(
    # ppm -> ppb
    no_daily_ppb  = ppm_to_ppb(no_daily),
    no2_daily_ppb = ppm_to_ppb(no2_daily),
    no_4_ppb      = ppm_to_ppb(no_4),
    no2_4_ppb     = ppm_to_ppb(no2_4),
    
    # logs
    log_pm25_daily    = log10(pm25_daily + 1),
    log_pm25_4        = log10(pm25_4 + 1),
    log_no_daily_ppb  = log10(no_daily_ppb + 1),
    log_no2_daily_ppb = log10(no2_daily_ppb + 1),
    log_no_4_ppb      = log10(no_4_ppb + 1),
    log_no2_4_ppb     = log10(no2_4_ppb + 1)
  ) %>%
  filter(
    adm_date >= as.Date("2002-01-01"),
    adm_date <= as.Date("2022-12-31")
  )

# (Optional) quick guardrails
stopifnot(
  all(c("rsv", "hrv") %in% names(climate)),
  all(c("rsv", "hrv") %in% names(pollution))
)

# 2) MODELS (UNADJUSTED)

run_for_outcome <- function(outcome = c("rsv", "hrv")) {
  outcome_var <- match.arg(tolower(outcome), c("rsv", "hrv"))
  
  # ----- CLIMATE: daily -----
  fit_and_show(as.formula(paste(outcome_var, "~ rain_daily")), climate,
               paste("Climate daily: Rainfall (", outcome_var, ")", sep = ""))
  fit_and_show(as.formula(paste(outcome_var, "~ tmax_daily")), climate,
               paste("Climate daily: Max Temp (", outcome_var, ")", sep = ""))
  fit_and_show(as.formula(paste(outcome_var, "~ tmin_daily")), climate,
               paste("Climate daily: Min Temp (", outcome_var, ")", sep = ""))
  
  # ----- CLIMATE: 4-day -----
  fit_and_show(as.formula(paste(outcome_var, "~ rain_4day")), climate,
               paste("Climate 4-day: Rainfall (", outcome_var, ")", sep = ""))
  fit_and_show(as.formula(paste(outcome_var, "~ tmax_4day")), climate,
               paste("Climate 4-day: Max Temp (", outcome_var, ")", sep = ""))
  fit_and_show(as.formula(paste(outcome_var, "~ tmin_4day")), climate,
               paste("Climate 4-day: Min Temp (", outcome_var, ")", sep = ""))
  
  # ----- POLLUTION: daily (raw + log) -----
  fit_and_show(as.formula(paste(outcome_var, "~ pm25_daily")), pollution,
               paste("Pollution daily: PM2.5 (", outcome_var, ")", sep = ""))
  fit_and_show(as.formula(paste(outcome_var, "~ no_daily_ppb")), pollution,
               paste("Pollution daily: NO (ppb) (", outcome_var, ")", sep = ""))
  fit_and_show(as.formula(paste(outcome_var, "~ no2_daily_ppb")), pollution,
               paste("Pollution daily: NO2 (ppb) (", outcome_var, ")", sep = ""))
  fit_and_show(as.formula(paste(outcome_var, "~ log_pm25_daily")), pollution,
               paste("Pollution daily (log10+1): PM2.5 (", outcome_var, ")", sep = ""))
  fit_and_show(as.formula(paste(outcome_var, "~ log_no_daily_ppb")), pollution,
               paste("Pollution daily (log10+1): NO (ppb) (", outcome_var, ")", sep = ""))
  fit_and_show(as.formula(paste(outcome_var, "~ log_no2_daily_ppb")), pollution,
               paste("Pollution daily (log10+1): NO2 (ppb) (", outcome_var, ")", sep = ""))
  
  # ----- POLLUTION: 4-day (raw + log) -----
  fit_and_show(as.formula(paste(outcome_var, "~ pm25_4")), pollution,
               paste("Pollution 4-day: PM2.5 (", outcome_var, ")", sep = ""))
  fit_and_show(as.formula(paste(outcome_var, "~ no_4_ppb")), pollution,
               paste("Pollution 4-day: NO (ppb) (", outcome_var, ")", sep = ""))
  fit_and_show(as.formula(paste(outcome_var, "~ no2_4_ppb")), pollution,
               paste("Pollution 4-day: NO2 (ppb) (", outcome_var, ")", sep = ""))
  fit_and_show(as.formula(paste(outcome_var, "~ log_pm25_4")), pollution,
               paste("Pollution 4-day (log10+1): PM2.5 (", outcome_var, ")", sep = ""))
  fit_and_show(as.formula(paste(outcome_var, "~ log_no_4_ppb")), pollution,
               paste("Pollution 4-day (log10+1): NO (ppb) (", outcome_var, ")", sep = ""))
  fit_and_show(as.formula(paste(outcome_var, "~ log_no2_4_ppb")), pollution,
               paste("Pollution 4-day (log10+1): NO2 (ppb) (", outcome_var, ")", sep = ""))
}

# Run
run_for_outcome("rsv")
run_for_outcome("hrv")
