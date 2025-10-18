# T-tests + Mann–Whitney with 95% CIs (2002–2022)

library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(lubridate)


num_safely <- function(x) {
  if (is.numeric(x)) return(x)
  suppressWarnings(as.numeric(x))
}

parse_date_any <- function(df, candidates) {
  nm <- candidates[candidates %in% names(df)][1]
  if (is.na(nm)) return(rep(as.Date(NA), nrow(df)))
  x <- df[[nm]]
  if (inherits(x, "Date")) return(as.Date(x))
  out <- suppressWarnings(dmy(x))
  if (all(is.na(out))) out <- suppressWarnings(ymd(x))
  if (all(is.na(out))) out <- suppressWarnings(parse_date_time(x, orders = c("dmy","ymd","mdy")))
  as.Date(out)
}

as_bin01 <- function(x) ifelse(x %in% c(0,1), as.integer(x), NA_integer_)

fmt_ci <- function(lo, hi, digits = 4) {
  if (is.na(lo) | is.na(hi)) return("[NA, NA]")
  paste0("[", signif(lo, digits), ", ", signif(hi, digits), "]")
}


mw_and_t <- function(data, group_col, var, log_transform = FALSE) {
  g <- data[[group_col]]
  x <- num_safely(data[[var]])
  if (log_transform) x <- log10(x + 1)
  
  keep <- !is.na(g) & !is.na(x)
  g <- g[keep]; x <- x[keep]
  
  # if only one group, return skeleton
  if (length(unique(g)) < 2L) {
    return(tibble(
      variable = var,
      n_0 = sum(g == 0), mean_0 = NA_real_, sd_0 = NA_real_,
      n_1 = sum(g == 1), mean_1 = NA_real_, sd_1 = NA_real_,
      t_stat = NA_real_, t_p = NA_real_,
      mw_W = NA_real_, mw_p = NA_real_,
      ci_low = NA_real_, ci_high = NA_real_
    ))
  }
  
  x0 <- x[g == 0]; x1 <- x[g == 1]
  mean_0 <- mean(x0); sd_0 <- sd(x0)
  mean_1 <- mean(x1); sd_1 <- sd(x1)
  
  
  t_res <- tryCatch(t.test(x ~ g, var.equal = TRUE, conf.level = 0.95),
                    error = function(e) NULL)
  
  mw_res <- tryCatch(wilcox.test(x ~ g, exact = FALSE),
                     error = function(e) NULL)
  
  tibble(
    variable = var,
    n_0 = sum(g == 0), mean_0 = mean_0, sd_0 = sd_0,
    n_1 = sum(g == 1), mean_1 = mean_1, sd_1 = sd_1,
    t_stat = ifelse(is.null(t_res), NA_real_, unname(t_res$statistic)),
    t_p    = ifelse(is.null(t_res), NA_real_, t_res$p.value),
    mw_W   = ifelse(is.null(mw_res), NA_real_, unname(mw_res$statistic)),
    mw_p   = ifelse(is.null(mw_res), NA_real_, mw_res$p.value),
    ci_low = ifelse(is.null(t_res), NA_real_, t_res$conf.int[1]),
    ci_high= ifelse(is.null(t_res), NA_real_, t_res$conf.int[2])
  )
}


run_block <- function(df, group_col, vars, log_transform = FALSE) {
  vars %>% map_dfr(~ mw_and_t(df, group_col, ., log_transform))
}

climate_path <- "R studio/Patient climate data 4dayAvg new.csv"
pollut_path  <- "R studio/Patient pollution data New.csv"

climate_raw   <- read_csv(climate_path, show_col_types = FALSE)
pollutant_raw <- read_csv(pollut_path,  show_col_types = FALSE)

# (a) Climate
climate <- climate_raw %>%
  mutate(across(where(is.numeric), ~ ifelse(. == 999, NA, .))) %>%
  mutate(
    RSV = as_bin01(RSV),
    HRV = as_bin01(HRV),
    AdmissionDate_any = parse_date_any(., c("AdmissionDate","Admission date","admission_date"))
  ) %>%
  select(RSV, HRV,
         `Max Temp`, `Min Temp`, `Rainfall`,
         `Mean Max Temp 4 day average`, `Mean Min Temp 4 day average`, `Mean Rainfall 4 day average`,
         AdmissionDate_any) %>%
  filter(AdmissionDate_any >= as.Date("2002-01-01"),
         AdmissionDate_any <= as.Date("2022-12-31")) %>%
  rename(
    tmax_daily = `Max Temp`,
    tmin_daily = `Min Temp`,
    rain_daily = `Rainfall`,
    tmax_4day  = `Mean Max Temp 4 day average`,
    tmin_4day  = `Mean Min Temp 4 day average`,
    rain_4day  = `Mean Rainfall 4 day average`,
    adm_date   = AdmissionDate_any
  )

# (b) Pollutants
pollutant <- pollutant_raw %>%
  mutate(
    RSV = as_bin01(rsv),
    HRV = as_bin01(hrv),
    admission_any = parse_date_any(., c("admission_date","AdmissionDate","Admission date"))
  ) %>%
  mutate(
    PM25_daily = num_safely(`PM2.5`),
    PM25_4     = num_safely(`PM2.5 4 day average`),
    NO_daily   = num_safely(NO),
    NO2_daily  = num_safely(NO2),
    NO_4       = num_safely(`NO 4 day average`),
    NO2_4      = num_safely(`NO2 4 day average`)
  ) %>%
  mutate(
    NO_daily_ppb  = NO_daily  * 1000,
    NO2_daily_ppb = NO2_daily * 1000,
    NO_4_ppb      = NO_4      * 1000,
    NO2_4_ppb     = NO2_4     * 1000
  ) %>%
  filter(admission_any >= as.Date("2002-01-01"),
         admission_any <= as.Date("2022-12-31")) %>%
  select(RSV, HRV,
         PM25_daily, NO_daily_ppb, NO2_daily_ppb,
         PM25_4,     NO_4_ppb,    NO2_4_ppb,
         admission_any) %>%
  rename(adm_date = admission_any)

clim_daily_vars <- c("tmax_daily", "tmin_daily", "rain_daily")
clim_4day_vars  <- c("tmax_4day",  "tmin_4day",  "rain_4day")

poll_daily_vars <- c("PM25_daily", "NO_daily_ppb", "NO2_daily_ppb")
poll_4day_vars  <- c("PM25_4",     "NO_4_ppb",     "NO2_4_ppb")

## RSV
tbl_clim_RSV_daily      <- run_block(climate,   "RSV", clim_daily_vars)
tbl_clim_RSV_4day       <- run_block(climate,   "RSV", clim_4day_vars)

tbl_poll_RSV_daily      <- run_block(pollutant, "RSV", poll_daily_vars)
tbl_poll_RSV_daily_log  <- run_block(pollutant, "RSV", poll_daily_vars,  log_transform = TRUE)

tbl_poll_RSV_4day       <- run_block(pollutant, "RSV", poll_4day_vars)
tbl_poll_RSV_4day_log   <- run_block(pollutant, "RSV", poll_4day_vars,   log_transform = TRUE)

## HRV
tbl_clim_HRV_daily      <- run_block(climate,   "HRV", clim_daily_vars)
tbl_clim_HRV_4day       <- run_block(climate,   "HRV", clim_4day_vars)

tbl_poll_HRV_daily      <- run_block(pollutant, "HRV", poll_daily_vars)
tbl_poll_HRV_daily_log  <- run_block(pollutant, "HRV", poll_daily_vars,  log_transform = TRUE)

tbl_poll_HRV_4day       <- run_block(pollutant, "HRV", poll_4day_vars)
tbl_poll_HRV_4day_log   <- run_block(pollutant, "HRV", poll_4day_vars,   log_transform = TRUE)

add_block <- function(df, virus, window, domain, logflag = FALSE) {
  df %>%
    mutate(
      Virus = virus,
      Window = window,
      Domain = domain,
      Transform = ifelse(logflag, "log10+1", "raw")
    ) %>%
    relocate(Virus, Window, Domain, Transform, variable)
}

out <- bind_rows(
  add_block(tbl_clim_RSV_daily,     "RSV", "Daily", "Climate"),
  add_block(tbl_clim_RSV_4day,      "RSV", "4-day", "Climate"),
  add_block(tbl_poll_RSV_daily,     "RSV", "Daily", "Pollution"),
  add_block(tbl_poll_RSV_daily_log, "RSV", "Daily", "Pollution", TRUE),
  add_block(tbl_poll_RSV_4day,      "RSV", "4-day", "Pollution"),
  add_block(tbl_poll_RSV_4day_log,  "RSV", "4-day", "Pollution", TRUE),
  
  add_block(tbl_clim_HRV_daily,     "HRV", "Daily", "Climate"),
  add_block(tbl_clim_HRV_4day,      "HRV", "4-day", "Climate"),
  add_block(tbl_poll_HRV_daily,     "HRV", "Daily", "Pollution"),
  add_block(tbl_poll_HRV_daily_log, "HRV", "Daily", "Pollution", TRUE),
  add_block(tbl_poll_HRV_4day,      "HRV", "4-day", "Pollution"),
  add_block(tbl_poll_HRV_4day_log,  "HRV", "4-day", "Pollution", TRUE)
)

# View (now includes 95% CIs for both tests)
print(out, n = 60, width = Inf)