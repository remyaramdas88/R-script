# Simple adjusted logistic regression (2002–2022)
#adjusted age and gender
# ---------------------------
library(readr)
library(dplyr)
library(lubridate)

# --- helpers ---
virus_01 <- function(x) {
  x <- suppressWarnings(as.integer(x))
  x[x == 999] <- NA_integer_
  ifelse(x %in% c(0L,1L), x, NA_integer_)
}
parse_date_any <- function(x) {
  if (inherits(x, "Date")) return(as.Date(x))
  out <- suppressWarnings(dmy(x))
  if (all(is.na(out))) out <- suppressWarnings(ymd(x))
  if (all(is.na(out))) out <- suppressWarnings(parse_date_time(x, orders = c("dmy","ymd","mdy")))
  as.Date(out)
}
gfac <- function(x) factor(ifelse(toupper(x)=="M","M","F"), levels = c("F","M"))

# CLIMATE 
climate <- read.csv("R studio/Patient climate data 4dayAvg new.csv",
                    header = TRUE, check.names = FALSE)
climate <- climate %>%
  mutate(
    RSV = virus_01(RSV),
    HRV = virus_01(HRV),
    age = suppressWarnings(as.numeric(Age)),
    gender = gfac(gender),
    adm_date = parse_date_any(`Admission date`),
    `Max Temp` = suppressWarnings(as.numeric(`Max Temp`)),
    `Min Temp` = suppressWarnings(as.numeric(`Min Temp`)),
    Rainfall   = suppressWarnings(as.numeric(Rainfall)),
    `Mean Max Temp 4 day average` = suppressWarnings(as.numeric(`Mean Max Temp 4 day average`)),
    `Mean Min Temp 4 day average` = suppressWarnings(as.numeric(`Mean Min Temp 4 day average`)),
    `Mean Rainfall 4 day average` = suppressWarnings(as.numeric(`Mean Rainfall 4 day average`))
  ) %>%
  mutate(across(c(`Max Temp`,`Min Temp`,Rainfall,
                  `Mean Max Temp 4 day average`,
                  `Mean Min Temp 4 day average`,
                  `Mean Rainfall 4 day average`),
                ~ ifelse(. == 999, NA, .))) %>%
  filter(adm_date >= as.Date("2002-01-01"),
         adm_date <= as.Date("2022-12-31")) %>%
  rename(
    tmax_daily = `Max Temp`,
    tmin_daily = `Min Temp`,
    rain_daily = Rainfall,
    tmax_4day  = `Mean Max Temp 4 day average`,
    tmin_4day  = `Mean Min Temp 4 day average`,
    rain_4day  = `Mean Rainfall 4 day average`
  )

# POLLUTION
pollution <- read.csv("R studio/Patient pollution data New.csv",
                      header = TRUE, check.names = FALSE)
pollution <- pollution %>%
  mutate(
    RSV = virus_01(coalesce(rsv)),
    HRV = virus_01(coalesce(hrv)),
    age = suppressWarnings(as.numeric(coalesce(age))),
    gender = gfac(gender),
    adm_date = parse_date_any(admission_date),
    PM25_daily = suppressWarnings(as.numeric(`PM2.5`)),
    PM25_4     = suppressWarnings(as.numeric(`PM2.5 4 day average`)),
    NO_daily   = suppressWarnings(as.numeric(NO)),
    NO2_daily  = suppressWarnings(as.numeric(NO2)),
    NO_4       = suppressWarnings(as.numeric(`NO 4 day average`)),
    NO2_4      = suppressWarnings(as.numeric(`NO2 4 day average`))
  ) %>%
  # convert ppm -> ppb
  mutate(
    NO_daily_ppb  = NO_daily  * 1000,
    NO2_daily_ppb = NO2_daily * 1000,
    NO_4_ppb      = NO_4      * 1000,
    NO2_4_ppb     = NO2_4     * 1000
  ) %>%
  # log10(+1)
  mutate(
    log_PM25_daily    = log10(PM25_daily + 1),
    log_PM25_4        = log10(PM25_4 + 1),
    log_NO_daily_ppb  = log10(NO_daily_ppb + 1),
    log_NO2_daily_ppb = log10(NO2_daily_ppb + 1),
    log_NO_4_ppb      = log10(NO_4_ppb + 1),
    log_NO2_4_ppb     = log10(NO2_4_ppb + 1)
  ) %>%
  filter(adm_date >= as.Date("2002-01-01"),
         adm_date <= as.Date("2022-12-31"))


# MODELS (adjusted for age + gender)

fit_and_show <- function(formula, data, label) {
  cat("\n====================\n", label, "\n====================\n", sep = "")
  m <- glm(formula, data = data, family = binomial, na.action = na.exclude)
  print(summary(m))
  # OR and 95% CI (exponentiated)
  print(exp(cbind(OR = coef(m), confint(m))))
  invisible(m)
}

run_for_outcome <- function(outcome = "RSV") {
  
  # CLIMATE: daily (rain, tmax, tmin)
  fit_and_show(as.formula(paste(outcome, "~ rain_daily + age + gender")), climate,
               paste("Climate daily: Rainfall + age + gender (Outcome:", outcome, ")"))
  fit_and_show(as.formula(paste(outcome, "~ tmax_daily + age + gender")), climate,
               paste("Climate daily: Max Temp + age + gender (Outcome:", outcome, ")"))
  fit_and_show(as.formula(paste(outcome, "~ tmin_daily + age + gender")), climate,
               paste("Climate daily: Min Temp + age + gender (Outcome:", outcome, ")"))
  
  # CLIMATE: 4-day means 
  fit_and_show(as.formula(paste(outcome, "~ rain_4day + age + gender")), climate,
               paste("Climate 4-day: Rainfall + age + gender (Outcome:", outcome, ")"))
  fit_and_show(as.formula(paste(outcome, "~ tmax_4day + age + gender")), climate,
               paste("Climate 4-day: Max Temp + age + gender (Outcome:", outcome, ")"))
  fit_and_show(as.formula(paste(outcome, "~ tmin_4day + age + gender")), climate,
               paste("Climate 4-day: Min Temp + age + gender (Outcome:", outcome, ")"))
  
  # POLLUTANTS: daily raw (PM2.5 ug/m3; NO/NO2 in ppb) 
  fit_and_show(as.formula(paste(outcome, "~ PM25_daily + age + gender")), pollution,
               paste("Pollution daily: PM2.5 + age + gender (Outcome:", outcome, ")"))
  fit_and_show(as.formula(paste(outcome, "~ NO_daily_ppb + age + gender")), pollution,
               paste("Pollution daily: NO(ppb) + age + gender (Outcome:", outcome, ")"))
  fit_and_show(as.formula(paste(outcome, "~ NO2_daily_ppb + age + gender")), pollution,
               paste("Pollution daily: NO2(ppb) + age + gender (Outcome:", outcome, ")"))
  
  # POLLUTANTS: 4-day raw 
  fit_and_show(as.formula(paste(outcome, "~ PM25_4 + age + gender")), pollution,
               paste("Pollution 4-day: PM2.5 + age + gender (Outcome:", outcome, ")"))
  fit_and_show(as.formula(paste(outcome, "~ NO_4_ppb + age + gender")), pollution,
               paste("Pollution 4-day: NO(ppb) + age + gender (Outcome:", outcome, ")"))
  fit_and_show(as.formula(paste(outcome, "~ NO2_4_ppb + age + gender")), pollution,
               paste("Pollution 4-day: NO2(ppb) + age + gender (Outcome:", outcome, ")"))
  
  # POLLUTANTS: daily log10(+1)
  fit_and_show(as.formula(paste(outcome, "~ log_PM25_daily + age + gender")), pollution,
               paste("Pollution daily (log10+1): PM2.5 + age + gender (Outcome:", outcome, ")"))
  fit_and_show(as.formula(paste(outcome, "~ log_NO_daily_ppb + age + gender")), pollution,
               paste("Pollution daily (log10+1): NO(ppb) + age + gender (Outcome:", outcome, ")"))
  fit_and_show(as.formula(paste(outcome, "~ log_NO2_daily_ppb + age + gender")), pollution,
               paste("Pollution daily (log10+1): NO2(ppb) + age + gender (Outcome:", outcome, ")"))
  
  # POLLUTANTS: 4-day log10(+1)
  fit_and_show(as.formula(paste(outcome, "~ log_PM25_4 + age + gender")), pollution,
               paste("Pollution 4-day (log10+1): PM2.5 + age + gender (Outcome:", outcome, ")"))
  fit_and_show(as.formula(paste(outcome, "~ log_NO_4_ppb + age + gender")), pollution,
               paste("Pollution 4-day (log10+1): NO(ppb) + age + gender (Outcome:", outcome, ")"))
  fit_and_show(as.formula(paste(outcome, "~ log_NO2_4_ppb + age + gender")), pollution,
               paste("Pollution 4-day (log10+1): NO2(ppb) + age + gender (Outcome:", outcome, ")"))
}

run_for_outcome("RSV")  
run_for_outcome("HRV")
