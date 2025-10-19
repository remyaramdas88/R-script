# 6-panel figures (daily + 4-day) for RSV/RV

library(patchwork)
climate_path   <- "R studio/Patient climate data 4dayAvg new.csv"
pollution_path <- "R studio/Patient pollution data New.csv"

out_dir <- "plots_all"
if (!dir.exists(out_dir)) dir.create(out_dir)

label_daily  <- c(
  tmax_daily = "Maximum temperature (°C)/daily ",
  tmin_daily = "Minimum temperature (°C)/daily",
  rain_daily = "Rainfall (mm)/daily",
  pm25_daily = "PM2.5 (µg/m³)/daily",
  no_daily_ppb = "NO (ppb)/daily ",
  no2_daily_ppb = "NO2 (ppb)/daily"
)
label_4day <- c(
  tmax_4day = "Maximum temperature (°C)/4-day mean ",
  tmin_4day = "Minimum temperature (°C)/4-day mean",
  rain_4day = "Rainfall (mm)/4-day mean ",
  pm25_4 = "PM2.5 (µg/m³)/4-day mean ",
  no_4_ppb = "NO (ppb)/4-day mean",
  no2_4_ppb = "NO2 (ppb)/4-day mean"
)

make_pp_plot_col <- function(df_pred, xvar, outcome_name, xlab, title) {
  if (toupper(outcome_name) == "RSV") {
    line_col <- "darkorange"
    fill_col <- "blue"
  } else {
    line_col <- "blue"
    fill_col <- "orange"
  }
  
  ggplot(df_pred, aes(x = .data[[xvar]], y = prob)) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), fill = fill_col, alpha = 0.25) +
    geom_line(size = 0.9, color = line_col) +
    scale_y_continuous(limits = c(0, 1),
                       labels = scales::percent_format(accuracy = 1)) +
    labs(
      title = title,
      x = xlab,
      y = "Predicted probability"
    ) +
    theme_bw(base_family = "Arial") +
    theme(
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 9),
      axis.text  = element_text(size = 9),
      panel.grid.minor = element_blank()
    )
}

pp_panel <- function(df, outcome_var, outcome_name, xvar, subtitle) {
  pd <- fit_predict_curve(df, outcome_var, xvar)
  if (is.null(pd)) return(NULL)
  xlab <- ifelse(xvar %in% names(pp_labels), pp_labels[[xvar]], xvar)
  make_pp_plot_col(pd, xvar = xvar, outcome_name = outcome_name,
                   xlab = xlab, title = subtitle)
}


make_6panel <- function(outcome_var, outcome_name, kind = c("climate","pollution")) {
  kind <- match.arg(kind)
  
  if (kind == "climate") {
    daily_vars <- c("tmax_daily","tmin_daily","rain_daily")
    day4_vars  <- c("tmax_4day","tmin_4day","rain_4day")
    df_use <- climate
    big_title <- paste0("Predicted probability of ", outcome_name,
                        " detection by climatic variables (daily and 4-day)")
  } else {
    daily_vars <- c("pm25_daily","no_daily_ppb","no2_daily_ppb")
    day4_vars  <- c("pm25_4","no_4_ppb","no2_4_ppb")
    df_use <- pollution
    big_title <- paste0("Predicted probability of ", outcome_name,
                        " detection by pollutant exposures (daily and 4-day)")
  }
  
  
  p1 <- pp_panel(df_use, outcome_var, outcome_name, daily_vars[1], label_daily[[daily_vars[1]]])
  p2 <- pp_panel(df_use, outcome_var, outcome_name, daily_vars[2], label_daily[[daily_vars[2]]])
  p3 <- pp_panel(df_use, outcome_var, outcome_name, daily_vars[3], label_daily[[daily_vars[3]]])
  p4 <- pp_panel(df_use, outcome_var, outcome_name, day4_vars[1],  label_4day[[day4_vars[1]]])
  p5 <- pp_panel(df_use, outcome_var, outcome_name, day4_vars[2],  label_4day[[day4_vars[2]]])
  p6 <- pp_panel(df_use, outcome_var, outcome_name, day4_vars[3],  label_4day[[day4_vars[3]]])
  
 
  top_row <- wrap_plots(Filter(Negate(is.null), list(p1,p2,p3)), nrow = 1)
  bot_row <- wrap_plots(Filter(Negate(is.null), list(p4,p5,p6)), nrow = 1)
  
  fig <- (top_row / bot_row) +
    plot_annotation(
      title = big_title,
      tag_levels = 'a',  # <-- adds a, b, c, d... automatically
      tag_prefix = "(",  # makes it (a), (b), etc.
      tag_suffix = ")",
      theme = theme(
        plot.title = element_text(family = "Arial", size = 11, face = "bold", hjust = 0.5),
        plot.tag = element_text(family = "Arial", size = 11, face = "bold")
      )
    )
  fig
}

# Generate all four figures (RSV/RV × climate/pollution) 
fig_rsv_clim <- make_6panel(outcome_var = "rsv", outcome_name = "RSV", kind = "climate")
fig_rsv_poll <- make_6panel(outcome_var = "rsv", outcome_name = "RSV", kind = "pollution")
fig_rv_clim  <- make_6panel(outcome_var = "hrv", outcome_name = "RV",  kind = "climate")
fig_rv_poll  <- make_6panel(outcome_var = "hrv", outcome_name = "RV",  kind = "pollution")

# Show in RStudio
print(fig_rsv_clim)
print(fig_rsv_poll)
print(fig_rv_clim)
print(fig_rv_poll)

ggsave(file.path(out_dir, "RSV_climate_daily_4day_6panel.png"),
       fig_rsv_clim, width = 14, height = 9, dpi = 600, device = grDevices::cairo_png)

ggsave(file.path(out_dir, "RSV_pollutants_daily_4day_6panel.png"),
       fig_rsv_poll, width = 14, height = 9, dpi = 600, device = grDevices::cairo_png)

ggsave(file.path(out_dir, "RV_climate_daily_4day_6panel.png"),
       fig_rv_clim, width = 14, height = 9, dpi = 600, device = grDevices::cairo_png)

ggsave(file.path(out_dir, "RV_pollutants_daily_4day_6panel.png"),
       fig_rv_poll, width = 14, height = 9, dpi = 600, device = grDevices::cairo_png)

