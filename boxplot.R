
library(readr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(patchwork)

climate_path <- "R studio/Patient climate data 4dayAvg new.csv"
pollut_path  <- "R studio/Patient pollution data New.csv"

num_safely <- function(x) { if (is.numeric(x)) x else suppressWarnings(as.numeric(x)) }
as_bin01   <- function(x) ifelse(x %in% c(0,1), as.integer(x), NA_integer_)
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

# Climate
climate_raw <- read_csv(climate_path, show_col_types = FALSE)
climate <- climate_raw %>%
  mutate(across(where(is.numeric), ~ ifelse(. == 999, NA, .))) %>%
  mutate(
    RSV = as_bin01(RSV),
    HRV = as_bin01(HRV),
    adm_any = parse_date_any(., c("AdmissionDate","Admission date","admission_date"))
  ) %>%
  select(
    RSV, HRV,
    `Max Temp`, `Min Temp`, `Rainfall`,
    `Mean Max Temp 4 day average`, `Mean Min Temp 4 day average`, `Mean Rainfall 4 day average`,
    adm_any
  ) %>%
  rename(
    tmax_daily = `Max Temp`,
    tmin_daily = `Min Temp`,
    rain_daily = `Rainfall`,
    tmax_4day  = `Mean Max Temp 4 day average`,
    tmin_4day  = `Mean Min Temp 4 day average`,
    rain_4day  = `Mean Rainfall 4 day average`,
    adm_date   = adm_any
  ) %>%
  filter(adm_date >= as.Date("2002-01-01"), adm_date <= as.Date("2022-12-31"))

# Pollution
pollutant_raw <- read_csv(pollut_path, show_col_types = FALSE)
pollutant <- pollutant_raw %>%
  mutate(
    RSV = as_bin01(rsv),
    HRV = as_bin01(hrv),
    adm_any = parse_date_any(., c("admission_date","AdmissionDate","Admission date")),
    PM25_daily = num_safely(`PM2.5`),
    PM25_4     = num_safely(`PM2.5 4 day average`),
    NO_daily   = num_safely(NO),
    NO2_daily  = num_safely(NO2),
    NO_4       = num_safely(`NO 4 day average`),
    NO2_4      = num_safely(`NO2 4 day average`),
    # convert NO/NO2 (ppm) -> ppb if needed
    NO_daily_ppb  = NO_daily  * 1000,
    NO2_daily_ppb = NO2_daily * 1000,
    NO_4_ppb      = NO_4      * 1000,
    NO2_4_ppb     = NO2_4     * 1000
  ) %>%
  select(
    RSV, HRV,
    PM25_daily, NO_daily_ppb, NO2_daily_ppb,
    PM25_4,     NO_4_ppb,    NO2_4_ppb,
    adm_any
  ) %>%
  rename(adm_date = adm_any) %>%
  filter(adm_date >= as.Date("2002-01-01"), adm_date <= as.Date("2022-12-31"))

clim_daily_vars <- c("tmax_daily","tmin_daily","rain_daily")
clim_4day_vars  <- c("tmax_4day","tmin_4day","rain_4day")
poll_daily_vars <- c("PM25_daily","NO_daily_ppb","NO2_daily_ppb")
poll_4day_vars  <- c("PM25_4","NO_4_ppb","NO2_4_ppb")

lab_var <- c(
  tmax_daily   = "Maximum temperature (°C) / daily",
  tmin_daily   = "Minimum temperature (°C) / daily",
  rain_daily   = "Rainfall (mm) / daily",
  tmax_4day    = "Maximum temperature (°C) / 4-day mean",
  tmin_4day    = "Minimum temperature (°C) / 4-day mean",
  rain_4day    = "Rainfall (mm) / 4-day mean",
  PM25_daily   = "PM2.5 (µg/m³) / daily",
  NO_daily_ppb = "NO (ppb) / daily",
  NO2_daily_ppb= "NO2 (ppb) / daily",
  PM25_4       = "PM2.5 (µg/m³) / 4-day mean",
  NO_4_ppb     = "NO (ppb) / 4-day mean",
  NO2_4_ppb    = "NO2 (ppb) / 4-day mean"
)

p_to_stars <- function(p){
  if (is.na(p)) ""
  else if (p < 1e-3) "***"
  else if (p < 1e-2) "**"
  else if (p < 0.05) "*"
}
is_4day_var <- function(v) grepl("(_4day$|_4$|4[_-]?day|_4_ppb$|_4$)", v, perl = TRUE)

make_box_panel <- function(df, outcome_var, value_var, panel_title) {
  d <- df %>%
    select(all_of(c(outcome_var, value_var))) %>%
    filter(!is.na(.data[[outcome_var]]), !is.na(.data[[value_var]])) %>%
    mutate(group = factor(ifelse(.data[[outcome_var]] == 1, "Positive", "Negative"),
                          levels = c("Negative","Positive")))
  if (nrow(d) < 3L) return(NULL)
  
  # Mann–Whitney p-value
  fml <- reformulate("group", response = value_var)
  mwp <- tryCatch(wilcox.test(fml, data = d, exact = FALSE)$p.value,
                  error = function(e) NA_real_)
  
  y_max <- max(d[[value_var]], na.rm = TRUE)
  y_min <- min(d[[value_var]], na.rm = TRUE)
  y_pad <- (y_max - y_min) * 0.08
  y_bar <- y_max + y_pad
  
  ggplot(d, aes(x = group, y = .data[[value_var]], fill = group)) +
    geom_boxplot(outlier.alpha = 0.35, width = 0.65, linewidth = 0.3) +
    scale_fill_manual(values = c(Negative = "blue", Positive = "orange")) +
    geom_segment(aes(x = 1, xend = 2, y = y_bar, yend = y_bar), linewidth = 0.3) +
    geom_segment(aes(x = 1, xend = 1, y = y_bar, yend = y_bar - y_pad*0.35), linewidth = 0.3) +
    geom_segment(aes(x = 2, xend = 2, y = y_bar, yend = y_bar - y_pad*0.35), linewidth = 0.3) +
    annotate("text", x = 1.5, y = y_bar + y_pad*0.55, label = p_to_stars(mwp),
             family = "Times New Roman", fontface = "bold", size = 3.6) +
    labs(
      title = panel_title,
      x = "Detection group",
      y = if (is_4day_var(value_var)) "4-day mean value" else "Daily value"
    ) +
    coord_cartesian(ylim = c(y_min, y_bar + y_pad*0.6)) +
    theme_bw(base_family = "Times New Roman") +
    theme(
      plot.title = element_text(size = 11.5, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 11.5),
      axis.title.y = element_text(size = 11.5),
      axis.text    = element_text(size = 11.5),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
}


make_6_box <- function(df, outcome_var, daily_vars, day4_vars, big_title) {
  p1 <- make_box_panel(df, outcome_var, daily_vars[1], lab_var[[ daily_vars[1] ]])
  p2 <- make_box_panel(df, outcome_var, daily_vars[2], lab_var[[ daily_vars[2] ]])
  p3 <- make_box_panel(df, outcome_var, daily_vars[3], lab_var[[ daily_vars[3] ]])
  p4 <- make_box_panel(df, outcome_var, day4_vars[1],  lab_var[[ day4_vars[1] ]])
  p5 <- make_box_panel(df, outcome_var, day4_vars[2],  lab_var[[ day4_vars[2] ]])
  p6 <- make_box_panel(df, outcome_var, day4_vars[3],  lab_var[[ day4_vars[3] ]])
  
  ((p1 | p2 | p3) / (p4 | p5 | p6)) +
    plot_annotation(
      title = big_title,
      tag_levels = "a", tag_prefix = "(", tag_suffix = ")",
      theme = theme(
        plot.title = element_text(family = "Times New Roman", size = 12.5, face = "bold", hjust = 0.5),
        plot.tag   = element_text(family = "Times New Roman", size = 11)
      )
    )
}

fig_rsv_clim <- make_6_box(
  climate, "RSV", clim_daily_vars, clim_4day_vars,
  "Distribution of climate variables by RSV detection (daily and 4-day)"
)
fig_rsv_poll <- make_6_box(
  pollutant, "RSV", poll_daily_vars, poll_4day_vars,
  "Distribution of pollutant exposures by RSV detection (daily and 4-day)"
)
fig_rv_clim <- make_6_box(
  climate, "HRV", clim_daily_vars, clim_4day_vars,
  "Distribution of climate variables by RV detection (daily and 4-day)"
)
fig_rv_poll <- make_6_box(
  pollutant, "HRV", poll_daily_vars, poll_4day_vars,
  "Distribution of pollutant exposures by RV detection (daily and 4-day)"
)

print(fig_rsv_clim); print(fig_rsv_poll); print(fig_rv_clim); print(fig_rv_poll)
