#RV subtypes

# Load libraries
library(dplyr)
library(ggplot2)

# 1. Filter for RV-A and RV-C only
rv_ac <- rv_pos %>%
  filter(subtype %in% c("RV-A", "RV-C")) %>%
  mutate(
    subtype_bin = if_else(subtype == "RV-C", 1, 0)  # 1 = RV-C, 0 = RV-A
  )

# 2. Choose the climate variable you want to test
# Options: "mean_max_temp_4_day_average", "mean_min_temp_4_day_average", "mean_rainfall_4_day_average"
temp_var <- "mean_rainfall_4_day_average"

# 3. Run logistic regression (RV-C vs RV-A)
model <- glm(
  formula = as.formula(paste("subtype_bin ~", temp_var, "+ age + gender")),
  data = rv_ac,
  family = binomial
)

# 4. View model summary
summary(model)


# Age summaries (2002–2022)
#--------------------------------

# Packages
library(dplyr)
library(lubridate)
library(writexl)

# -------- Helper: parse admission date safely ----------
parse_adm_date <- function(x) {
  # If numeric (Excel serial), convert from 1899-12-30
  x_num <- suppressWarnings(as.numeric(x))
  as_date <- ifelse(!is.na(x_num),
                    as.Date(x_num, origin = "1899-12-30"),
                    NA)
  as_date <- as.Date(as_date, origin = "1970-01-01")  # keep Date class if already Date
  
  # If not numeric, try common formats (your file is dd-mm-YYYY)
  if (all(is.na(as_date))) {
    as_date <- suppressWarnings(parse_date_time(x, orders = c("dmy","ymd","mdy")))
    as_date <- as.Date(as_date)
  }
  as_date
}

# -------- Read data  --------
df <- read.csv("R studio/Patient climate data 4dayAvg new.csv",
               header = TRUE, check.names = FALSE)

# Prefer "Admission date" if present; else use "AdmissionDate"
df <- df %>%
  mutate(
    Admission_date_raw = dplyr::coalesce(
      `Admission date` %||% NULL,
      AdmissionDate %||% NULL
    ),
    Admission_date = parse_adm_date(Admission_date_raw)
  )

# -------- Filter to 2002–2022 -----------------------
df_filtered <- df %>%
  filter(!is.na(Admission_date),
         year(Admission_date) >= 2002,
         year(Admission_date) <= 2022)

# -------- Clean virus flags (0/1; 999 -> NA) -------
clean_virus <- function(v) {
  v <- suppressWarnings(as.integer(v))
  v[v == 999] <- NA_integer_
  ifelse(v %in% c(0L, 1L), v, NA_integer_)
}

df_filtered <- df_filtered %>%
  mutate(
    RSV = clean_virus(RSV),
    HRV = clean_virus(HRV)
  )

# -------- Gender summary ----------------------------
# Map M/F to labelled categories used in your table
gender_summary <- df_filtered %>%
  mutate(Category = ifelse(toupper(gender) == "F", "Female (F)",
                           ifelse(toupper(gender) == "M", "Male (M)", NA_character_))) %>%
  filter(!is.na(Category)) %>%
  summarise(
    Mean_Age   = mean(Age, na.rm = TRUE),
    Median_Age = median(Age, na.rm = TRUE),
    SD_Age     = sd(Age, na.rm = TRUE),
    Total_Case = n(),
    .by = Category
  )

# -------- RSV summary -------------------------------
rsv_summary <- df_filtered %>%
  filter(!is.na(RSV)) %>%
  mutate(Category = ifelse(RSV == 1, "RSV Positive (1)", "RSV Negative (0)")) %>%
  summarise(
    Mean_Age   = mean(Age, na.rm = TRUE),
    Median_Age = median(Age, na.rm = TRUE),
    SD_Age     = sd(Age, na.rm = TRUE),
    Total_Case = n(),
    .by = Category
  )

# -------- RV summary -------------------------------
hrv_summary <- df_filtered %>%
  filter(!is.na(HRV)) %>%
  mutate(Category = ifelse(HRV == 1, "HRV Positive (1)", "HRV Negative (0)")) %>%
  summarise(
    Mean_Age   = mean(Age, na.rm = TRUE),
    Median_Age = median(Age, na.rm = TRUE),
    SD_Age     = sd(Age, na.rm = TRUE),
    Total_Case = n(),
    .by = Category
  )

# -------- Overall summary ---------------------------
overall_summary <- df_filtered %>%
  summarise(
    Category   = "Overall",
    Mean_Age   = mean(Age, na.rm = TRUE),
    Median_Age = median(Age, na.rm = TRUE),
    SD_Age     = sd(Age, na.rm = TRUE),
    Total_Case = n()
  )

# -------- 9) Combine all tables ------------------------
final_summary <- bind_rows(
  gender_summary %>% mutate(Category = as.character(Category)),
  rsv_summary    %>% mutate(Category = as.character(Category)),
  hrv_summary    %>% mutate(Category = as.character(Category)),
  overall_summary%>% mutate(Category = as.character(Category))
) %>%
  select(Category, Mean_Age, Median_Age, SD_Age, Total_Case)

# -------- 10) View & Save ------------------------------
print(final_summary)


