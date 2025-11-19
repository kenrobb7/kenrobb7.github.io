# charls_multimorbidity_analysis.R ----------------------------------------------------
#
# This script implements the analytic plan that links lifetime adversity to
# incident cardiometabolic multimorbidity (CMM) using CHARLS 2014-2020 data.
# It is written in a modular way so each step (data import, cohort construction,
# adversity scoring, outcome derivation, modeling, and reporting) can be executed
# independently or end-to-end. Replace the placeholder file names with the actual
# CHARLS datasets available in your environment.
# -------------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(broom)
  library(gt)
  library(glue)
})

# ---------------------------------------------------------------------------- #
# 0. User inputs
# ---------------------------------------------------------------------------- #
data_dir <- "data"  # <-- change to the folder where the cleaned CHARLS files live
output_dir <- "outputs"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

life_history_path <- file.path(data_dir, "charls2014_life_history.csv")
w2015_path <- file.path(data_dir, "charls2015_wave3.csv")
w2020_path <- file.path(data_dir, "charls2020_wave5.csv")

# ---------------------------------------------------------------------------- #
# 1. Helper functions
# ---------------------------------------------------------------------------- #

affirm_binary <- function(x) {
  # CHARLS stores many items as yes/no, 1/2, or other encodings.
  # This helper coerces them to 0/1 while keeping NAs intact.
  case_when(
    is.na(x) ~ NA_real_,
    x %in% c(1, "yes", "Yes", "Y") ~ 1,
    x %in% c(0, 2, "no", "No", "N") ~ 0,
    TRUE ~ NA_real_
  )
}

construct_adversity_scores <- function(df, prefix, threat_items, deprivation_items) {
  # prefix: a string ("child" or "adult") used for naming the resulting columns
  df %>%
    mutate(
      across(all_of(threat_items), affirm_binary),
      across(all_of(deprivation_items), affirm_binary)
    ) %>%
    mutate(
      !!glue::glue("{prefix}_threat") := rowSums(across(all_of(threat_items)), na.rm = TRUE),
      !!glue::glue("{prefix}_deprivation") := rowSums(across(all_of(deprivation_items)), na.rm = TRUE)
    ) %>%
    mutate(
      across(
        .cols = c(!!glue::glue("{prefix}_threat"), !!glue::glue("{prefix}_deprivation")),
        .fns = ~ if_else(is.na(.x), NA_real_, .x)
      ),
      !!glue::glue("{prefix}_total") := !!sym(glue::glue("{prefix}_threat")) +
        !!sym(glue::glue("{prefix}_deprivation"))
    )
}

categorise_score <- function(score, breaks = c(-Inf, 0, 2, 4, Inf), labels = c("0", "1-2", "3-4", "5+")) {
  cut(score, breaks = breaks, labels = labels, right = TRUE)
}

summarise_for_table <- function(df, group_var) {
  group_var <- enquo(group_var)
  df %>%
    group_by(!!group_var) %>%
    summarise(
      n = n(),
      age_mean = mean(age, na.rm = TRUE),
      female_pct = mean(sex == "Female", na.rm = TRUE) * 100,
      bmi_mean = mean(bmi, na.rm = TRUE),
      current_smoker_pct = mean(smoke_status == "Current", na.rm = TRUE) * 100
    ) %>%
    ungroup()
}

fit_logistic_models <- function(df, exposure_var) {
  exposure_var <- enquo(exposure_var)

  model1 <- glm(
    incident_cmm ~ !!exposure_var + age + sex,
    data = df, family = binomial(link = "logit")
  )

  model2 <- glm(
    incident_cmm ~ !!exposure_var + age + sex + residence + education + marital_status + childhood_econ + hh_income_quintile,
    data = df, family = binomial(link = "logit")
  )

  model3 <- glm(
    incident_cmm ~ !!exposure_var + age + sex + residence + education + marital_status + childhood_econ + hh_income_quintile +
      smoke_status + drink_status + physical_activity + bmi + baseline_single_cvd,
    data = df, family = binomial(link = "logit")
  )

  bind_rows(
    tidy(model1, exponentiate = TRUE, conf.int = TRUE) %>% mutate(model = "Model 1"),
    tidy(model2, exponentiate = TRUE, conf.int = TRUE) %>% mutate(model = "Model 2"),
    tidy(model3, exponentiate = TRUE, conf.int = TRUE) %>% mutate(model = "Model 3")
  ) %>%
    filter(term == as_name(exposure_var))
}

# ---------------------------------------------------------------------------- #
# 2. Data import
# ---------------------------------------------------------------------------- #
message("Reading raw data ...")
life_history <- read_csv(life_history_path, show_col_types = FALSE)
wave2015 <- read_csv(w2015_path, show_col_types = FALSE)
wave2020 <- read_csv(w2020_path, show_col_types = FALSE)

# ---------------------------------------------------------------------------- #
# 3. Adversity measures
# ---------------------------------------------------------------------------- #
child_threat <- c("ch_phys_abuse", "ch_emotional_abuse", "ch_witness_violence", "ch_parent_loss")
child_deprivation <- c("ch_food_insecure", "ch_poor_housing", "ch_school_dropout")
adult_threat <- c("ad_spouse_loss", "ad_child_loss", "ad_violence")
adult_deprivation <- c("ad_poverty", "ad_longterm_unemployment", "ad_functional_loss")

life_history_scores <- life_history %>%
  construct_adversity_scores("child", child_threat, child_deprivation) %>%
  construct_adversity_scores("adult", adult_threat, adult_deprivation) %>%
  mutate(
    lifetime_total = child_total + adult_total,
    lifetime_threat = child_threat + adult_threat,
    lifetime_deprivation = child_deprivation + adult_deprivation,
    lifetime_total_cat = categorise_score(lifetime_total),
    lifetime_threat_cat = categorise_score(lifetime_threat, breaks = c(-Inf, 0, 1, 2, Inf), labels = c("0", "1", "2", "3+")),
    lifetime_deprivation_cat = categorise_score(lifetime_deprivation, breaks = c(-Inf, 0, 1, 2, Inf), labels = c("0", "1", "2", "3+"))
  ) %>%
  select(id, starts_with("child_"), starts_with("adult_"), starts_with("lifetime_"))

# ---------------------------------------------------------------------------- #
# 4. Cardiometabolic disease counts (2015 & 2020)
# ---------------------------------------------------------------------------- #
cm_vars <- c("hypertension", "diabetes", "dyslipidemia", "heart_disease", "stroke")

derive_cm_counts <- function(df, suffix) {
  df %>%
    mutate(across(all_of(cm_vars), affirm_binary)) %>%
    mutate(!!glue::glue("cm_count_{suffix}") := rowSums(across(all_of(cm_vars)), na.rm = TRUE)) %>%
    rename_with(~ glue::glue("{.x}_{suffix}"), all_of(cm_vars))
}

wave2015_cm <- derive_cm_counts(wave2015, "2015")
wave2020_cm <- derive_cm_counts(wave2020, "2020")

# ---------------------------------------------------------------------------- #
# 5. Cohort construction
# ---------------------------------------------------------------------------- #
message("Constructing analytic cohort ...")

cohort <- wave2015_cm %>%
  inner_join(life_history_scores, by = "id") %>%
  inner_join(wave2020_cm %>% select(id, cm_count_2020, ends_with("_2020")), by = "id") %>%
  filter(age >= 45) %>%
  mutate(
    baseline_cmm = cm_count_2015 >= 2,
    baseline_single_cvd = cm_count_2015 == 1,
    incident_cmm = case_when(
      baseline_cmm ~ NA,  # drop those with existing multimorbidity
      cm_count_2020 >= 2 ~ 1,
      cm_count_2020 < 2 ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!baseline_cmm, !is.na(incident_cmm))

# apply covariate encodings ---------------------------------------------------
cohort <- cohort %>%
  mutate(
    sex = factor(if_else(sex %in% c(1, "Male"), "Male", "Female")),
    residence = factor(residence, levels = c("Rural", "Urban")),
    education = factor(education, levels = c("No schooling", "Primary", "Middle", "High+"), ordered = TRUE),
    marital_status = factor(marital_status, levels = c("Married", "Other")),
    smoke_status = factor(smoke_status, levels = c("Never", "Former", "Current"), ordered = TRUE),
    drink_status = factor(drink_status, levels = c("Never", "Occasional", "Frequent"), ordered = TRUE),
    physical_activity = factor(physical_activity, levels = c("Low", "Moderate", "High"), ordered = TRUE),
    childhood_econ = factor(childhood_econ, levels = c("Good", "Average", "Poor"), ordered = TRUE),
    hh_income_quintile = factor(hh_income_quintile, levels = 1:5, ordered = TRUE)
  )

# ---------------------------------------------------------------------------- #
# 6. Descriptive table
# ---------------------------------------------------------------------------- #
message("Producing descriptive statistics ...")

descriptive_table <- summarise_for_table(cohort, lifetime_total_cat) %>%
  gt() %>%
  fmt_number(columns = c(age_mean, bmi_mean), decimals = 1) %>%
  fmt_number(columns = c(female_pct, current_smoker_pct), decimals = 1, pattern = "{x}%")

gtsave(descriptive_table, file = file.path(output_dir, "table1_lifetime_adversity.png"))

# ---------------------------------------------------------------------------- #
# 7. Regression models
# ---------------------------------------------------------------------------- #
message("Running logistic regression models ...")

model_total <- fit_logistic_models(cohort, lifetime_total_cat)
model_threat <- fit_logistic_models(cohort, lifetime_threat_cat)
model_deprivation <- fit_logistic_models(cohort, lifetime_deprivation_cat)

model_results <- bind_rows(
  model_total %>% mutate(exposure = "Lifetime total adversity"),
  model_threat %>% mutate(exposure = "Lifetime threat adversity"),
  model_deprivation %>% mutate(exposure = "Lifetime deprivation adversity")
)

write_csv(model_results, file.path(output_dir, "logistic_models_results.csv"))

# Optional: formatted table ---------------------------------------------------
model_table <- model_results %>%
  mutate(
    estimate_fmt = sprintf("%.2f (%.2f, %.2f)", estimate, conf.low, conf.high)
  ) %>%
  select(exposure, model, term, estimate_fmt) %>%
  pivot_wider(names_from = model, values_from = estimate_fmt)

model_table_gt <- model_table %>%
  gt(rowname_col = "term", groupname_col = "exposure") %>%
  tab_header(
    title = "Lifetime adversity and incident cardiometabolic multimorbidity",
    subtitle = "Odds ratios (95% CI) from logistic regression models"
  )

gtsave(model_table_gt, file = file.path(output_dir, "table2_logistic_models.png"))

message("All done! Outputs saved to: ", normalizePath(output_dir))
