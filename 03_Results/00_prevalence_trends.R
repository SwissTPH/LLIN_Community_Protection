#################################
# Figure 2: Prevalence trends
#
# Created: September 2026
#
# Purpose:
#   Generate the prevalence trend figure showing malaria
#   parasite prevalence among LLIN users and non-users.
#
# Data:
#   OpenMalaria simulation outputs stored in:
#   02_Simulation_Data/OpenMalaria_simulation_data.sqlite
#
#################################

# Clear workspace
rm(list = ls())


# ------------------------------------------------------------
# Load required packages
# ------------------------------------------------------------

library(tidyverse)
library(RSQLite)
library(DBI)
library(here)
library(grid)


# ------------------------------------------------------------
# Load OpenMalaria simulation results
# ------------------------------------------------------------

# Define the path to the SQLite database relative to the
# repository root. This makes the code portable across users.
db_file <- here(
  "02_Simulation_Data",
  "OpenMalaria_simulation_data.sqlite"
)

# Check that the database exists
if (!file.exists(db_file)) {
  stop(
    "SQLite database not found: ",
    db_file
  )
}

# Open database connection
conn <- DBI::dbConnect(
  RSQLite::SQLite(),
  db_file
)

# List database tables
DBI::dbListTables(conn)

# Extract the OpenMalaria results table
simulation_results <- DBI::dbReadTable(
  conn,
  "om_results"
)

# Always disconnect from the database
DBI::dbDisconnect(conn)


# ------------------------------------------------------------
# Load scenario information
# ------------------------------------------------------------

# Scenario information is required to link simulation results
# to the corresponding intervention and transmission settings.
scenario_file <- here(
  "02_Simulation_Data",
  "scenarios.rds"
)

if (!file.exists(scenario_file)) {
  stop(
    "Scenario file not found: ",
    scenario_file
  )
}

scenario_data <- readRDS(
  scenario_file
)


# ------------------------------------------------------------
# Merge simulation results with scenario information
# ------------------------------------------------------------

scenario_data <- scenario_data %>%
  mutate(
    scenario_id = ID
  )

openmalaria_data <- simulation_results %>%
  left_join(
    scenario_data,
    by = "scenario_id"
  ) %>%
  mutate(
    year = lubridate::year(
      as.Date(date)
    )
  )


# ------------------------------------------------------------
# Basic data checks
# ------------------------------------------------------------

message(
  "Number of simulation records: ",
  nrow(openmalaria_data)
)

message(
  "Number of scenarios: ",
  n_distinct(openmalaria_data$scenario_id)
)

message(
  "Years covered: ",
  min(openmalaria_data$year, na.rm = TRUE),
  "–",
  max(openmalaria_data$year, na.rm = TRUE)
)


# ------------------------------------------------------------
# 1. Define population groups and EIR categories
# ------------------------------------------------------------

prevalence_data <- openmalaria_data %>%
  select(
    year,
    age_group,
    seed,
    EIR,
    futNetcovstart2023,
    prevalenceRate
  ) %>%
  mutate(
    population = case_when(
      age_group == "LLINusers_0-100" ~ "LLIN users",
      age_group == "0-100"          ~ "LLIN non-users",
      TRUE                          ~ NA_character_
    ),
    EIR_cat = case_when(
      EIR <= 2             ~ "Low",
      EIR > 2 & EIR <= 8   ~ "Moderate",
      EIR > 8              ~ "High",
      TRUE                 ~ NA_character_
    ),
    EIR_cat = factor(
      EIR_cat,
      levels = c("Low", "Moderate", "High")
    )
  ) %>%
  filter(
    year > 2021,
    year < 2032,
    !is.na(population),
    !is.na(EIR_cat)
  ) %>%
  group_by(
    year,
    seed,
    population,
    EIR_cat,
    futNetcovstart2023,
    EIR
  ) %>%
  summarise(
    prevalenceRate = mean(
      prevalenceRate,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# ------------------------------------------------------------
# 2. Calculate baseline prevalence for LLIN non-users
# ------------------------------------------------------------

baseline_prevalence <- prevalence_data %>%
  filter(
    futNetcovstart2023 == 0,
    population == "LLIN non-users"
  ) %>%
  group_by(
    year,
    seed,
    EIR_cat,
    EIR
  ) %>%
  summarise(
    baseline_prev = mean(
      prevalenceRate,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# ------------------------------------------------------------
# 3. Calculate prevalence reduction for users and non-users
# ------------------------------------------------------------

prevalence_by_population <- prevalence_data %>%
  left_join(
    baseline_prevalence,
    by = c(
      "year",
      "seed",
      "EIR_cat",
      "EIR"
    )
  ) %>%
  mutate(
    prev_reduction =
      100 * (baseline_prev - prevalenceRate) /
      baseline_prev
  )


# ------------------------------------------------------------
# 4. Calculate overall population prevalence
# ------------------------------------------------------------

overall_prevalence_by_seed <- openmalaria_data %>%
  select(
    year,
    EIR,
    age_group,
    seed,
    futNetcovstart2023,
    nHost,
    nPatent
  ) %>%
  mutate(
    population = case_when(
      age_group == "0-100"          ~ "LLIN non-users",
      age_group == "LLINusers_0-100" ~ "LLIN users",
      TRUE                           ~ NA_character_
    ),
    EIR_cat = case_when(
      EIR <= 2             ~ "Low",
      EIR > 2 & EIR <= 8   ~ "Moderate",
      EIR > 8              ~ "High",
      TRUE                 ~ NA_character_
    ),
    EIR_cat = factor(
      EIR_cat,
      levels = c("Low", "Moderate", "High")
    )
  ) %>%
  filter(
    year > 2021,
    year < 2032,
    !is.na(population),
    !is.na(EIR_cat)
  ) %>%
  group_by(
    year,
    seed,
    EIR_cat,
    futNetcovstart2023,
    EIR
  ) %>%
  summarise(
    total_pop = sum(nHost, na.rm = TRUE),
    total_patent = sum(nPatent, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    prevalenceRate = total_patent / total_pop
  )


# ------------------------------------------------------------
# 5. Calculate baseline overall population prevalence
# ------------------------------------------------------------

overall_baseline_prevalence <- overall_prevalence_by_seed %>%
  filter(
    futNetcovstart2023 == 0
  ) %>%
  select(
    year,
    seed,
    EIR,
    EIR_cat,
    baseline_prev = prevalenceRate
  )


# ------------------------------------------------------------
# 6. Calculate overall population prevalence reduction
# ------------------------------------------------------------

overall_prevalence_by_seed <- overall_prevalence_by_seed %>%
  left_join(
    overall_baseline_prevalence,
    by = c(
      "year",
      "seed",
      "EIR",
      "EIR_cat"
    )
  ) %>%
  mutate(
    prev_reduction =
      100 * (baseline_prev - prevalenceRate) /
      baseline_prev,
    population = "Overall population"
  )


# ------------------------------------------------------------
# 7. Combine population-specific and overall estimates
# ------------------------------------------------------------

figure2_seed_data <- bind_rows(
  prevalence_by_population %>%
    select(
      year,
      seed,
      EIR,
      EIR_cat,
      futNetcovstart2023,
      population,
      prevalenceRate,
      prev_reduction
    ),
  
  overall_prevalence_by_seed %>%
    select(
      year,
      seed,
      EIR,
      EIR_cat,
      futNetcovstart2023,
      population,
      prevalenceRate,
      prev_reduction
    )
) %>%
  mutate(
    population = factor(
      population,
      levels = c(
        "LLIN users",
        "LLIN non-users",
        "Overall population"
      )
    ),
    futNetcovstart2023 = round(
      futNetcovstart2023,
      3
    )
  )


# ------------------------------------------------------------
# 8. Summarise prevalence across simulation seeds
# ------------------------------------------------------------

figure2_summary <- figure2_seed_data %>%
  group_by(
    year,
    EIR_cat,
    futNetcovstart2023,
    population
  ) %>%
  summarise(
    med = median(
      prevalenceRate,
      na.rm = TRUE
    ),
    q25 = quantile(
      prevalenceRate,
      0.25,
      na.rm = TRUE
    ),
    q75 = quantile(
      prevalenceRate,
      0.75,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# ------------------------------------------------------------
# 9. Add 2022 prevalence as the counterfactual reference
# ------------------------------------------------------------

counterfactual_2022 <- figure2_summary %>%
  filter(
    futNetcovstart2023 == 0,
    population == "Overall population",
    year == 2022
  ) %>%
  select(
    EIR_cat,
    med
  ) %>%
  mutate(
    population = "Counterfactual (2022 PfPR)"
  )


counterfactual_reference <- expand.grid(
  year = unique(figure2_summary$year),
  futNetcovstart2023 = unique(
    figure2_summary$futNetcovstart2023
  ),
  EIR_cat = unique(
    figure2_summary$EIR_cat
  )
) %>%
  left_join(
    counterfactual_2022,
    by = "EIR_cat"
  )


# ------------------------------------------------------------
# 10. Prepare data for Figure 2
# ------------------------------------------------------------

figure2_plot_data <- bind_rows(
  figure2_summary,
  counterfactual_reference
) %>%
  filter(
    futNetcovstart2023 == 0.5,
    EIR_cat == "High",
    population != "Overall population"
  )