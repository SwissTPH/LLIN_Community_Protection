#################################
# Figure 00: Prevalence trends
#
# Created: September 2026
#
# Purpose:
#   Generate the prevalence trend figure showing malaria
#   parasite prevalence among ITN users and non-users at 50% ITN usage
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
      age_group == "LLINusers_0-100" ~ "ITN users",
      age_group == "0-100"          ~ "ITN non-users",
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
# 2. Calculate baseline prevalence for ITN non-users
# ------------------------------------------------------------

baseline_prevalence <- prevalence_data %>%
  filter(
    futNetcovstart2023 == 0,
    population == "ITN non-users"
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
      age_group == "0-100"          ~ "ITN non-users",
      age_group == "LLINusers_0-100" ~ "ITN users",
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

figure00_seed_data <- bind_rows(
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
        "ITN users",
        "ITN non-users",
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

figure00_summary <- figure00_seed_data %>%
  group_by(
    year,
    EIR_cat,
    futNetcovstart2023,
    population
  ) %>%
  summarise(
    mean_prevalence = mean(
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

counterfactual_2022 <- figure00_summary %>%
  filter(
    futNetcovstart2023 == 0,
    population == "Overall population",
    year == 2022
  ) %>%
  select(
    EIR_cat,
    mean_prevalence
  ) %>%
  mutate(
    population = "Counterfactual (2022 PfPR)"
  )


counterfactual_reference <- expand.grid(
  year = unique(figure00_summary$year),
  futNetcovstart2023 = unique(
    figure00_summary$futNetcovstart2023
  ),
  EIR_cat = unique(
    figure00_summary$EIR_cat
  )
) %>%
  left_join(
    counterfactual_2022,
    by = "EIR_cat"
  )


# ------------------------------------------------------------
# 10. Prepare data for Figure 2
# ------------------------------------------------------------

figure00_plot_data <- bind_rows(
  figure00_summary,
  counterfactual_reference
) %>%
  filter(
    futNetcovstart2023 == 0.5,
    #EIR_cat == "High",
    population != "Overall population"
  )

# ------------------------------------------------------------
# 10. Plot Figure00 as yearl prevalence trend at 50% ITN usage
# ------------------------------------------------------------

theme_pub <- function(base_size = 10) {
  
  theme_bw(
    base_size = base_size,
    base_family = "Arial"
  ) +
    theme(
      # Overall text
      text = element_text(
        family = "Arial",
        colour = "black"
      ),
      
      # Legend
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_blank(),
      legend.text = element_text(
        size = base_size - 1
      ),
      legend.key = element_blank(),
      legend.key.height = unit(0.5, "lines"),
      legend.spacing.x = unit(0.25, "cm"),
      
      # Facets
      strip.text = element_text(
        face = "bold",
        size = base_size,
        colour = "black"
      ),
      strip.background = element_rect(
        fill = "white",
        colour = NA
      ),
      
      # Grid
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(
        colour = "grey85",
        linewidth = 0.25
      ),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.8, "lines"),
      
      # Axes
      axis.title = element_text(
        face = "bold",
        size = base_size,
        colour = "black"
      ),
      axis.text = element_text(
        colour = "black",
        size = base_size - 1
      ),
      
      # Panel border
      panel.border = element_rect(
        colour = "black",
        linewidth = 0.5
      ),
      
      # Axis ticks
      axis.ticks = element_line(
        linewidth = 0.4,
        colour = "black"
      ),
      
      # Do not put the manuscript figure title inside the plot
      plot.title = element_blank(),
      
      # Useful for multi-panel labels such as A and B
      plot.tag = element_text(
        family = "Arial",
        face = "bold",
        size = base_size + 1,
        colour = "black"
      ),
      
      # Small outer border while avoiding excessive whitespace
      plot.margin = margin(5, 5, 5, 5),
      
      # White background
      plot.background = element_rect(
        fill = "white",
        colour = NA
      )
    )
}


cols_pop <- c(
  "ITN users"     = "#432CA1",
  "ITN non-users" = "#D55E00",
  "Counterfactual (2022 PfPR)" = "grey"
)

figure00_prevalence_trend <- ggplot(
  figure00_plot_data,
  aes(
    x = year,
    y = mean_prevalence,
    group = population,
    colour = population,
    fill = population
  )
) +
  
  # Interquartile range across simulation seeds
  geom_ribbon(
    aes(
      ymin = q25,
      ymax = q75
    ),
    alpha = 0.18,
    colour = NA,
    show.legend = FALSE
  ) +
  
  # Mean prevalence reduction
  geom_line(
    linewidth = 0.8
  ) +
  
  geom_point(
    size = 1.5
  ) +
  
  # Transmission intensity panels
  facet_wrap(
    ~ EIR_cat,
    nrow = 1
  ) +
  
  scale_colour_manual(
    values = cols_pop,
    breaks = c(
      "ITN non-users",
      "ITN users",
      "Counterfactual (2022 PfPR)"
    ),
    name = NULL
  ) +
  
  scale_fill_manual(
    values = cols_pop,
    breaks = c(
      "ITN non-users",
      "ITN users",
      "Counterfactual (2022 PfPR)"
    ),
    name = NULL
  ) +
  
  
  labs(
    x = "Year",
    y = "All-age prevalence"
  ) +
  
  theme_pub(
    base_size = 10
  ) +
  
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    axis.text.x = element_text(
      size = 8.5,
      hjust = 0.5
    )
  ) +
  
  guides(
    colour = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        linewidth = 0.8,
        size = 1.5
      )
    )
  )

figure00_prevalence_trend
