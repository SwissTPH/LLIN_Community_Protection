#################################
# Figure 2: Prevalence reduction
#
# Created: September 2026
#
# Purpose:
#   Generate the prevalence reduction figure showing malaria
#   parasite prevalence reduction among ITN users and non-users.
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
# 2. Calculate baseline prevalence for I non-users
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


# ============================================================
# Figure 2: Relative prevalence reduction
#
# Mean relative prevalence reduction across the intervention
# period (2023–2031) among ITN users and non-users.
# ============================================================

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
  "ITN non-users" = "#D55E00"
)


# ------------------------------------------------------------
# Summarize prevalence reduction
# ------------------------------------------------------------
#
# For each year, EIR category, ITN usage level and population
# group, calculate the mean and IQR across simulation seeds and 
# averaged across 2023–2031.
# ------------------------------------------------------------

figure2_annual_reduction <- figure2_seed_data %>%
  filter(
    year >= 2023,
    year <= 2031,
    population %in% c(
      "ITN non-users",
      "ITN users"
    )
  ) %>%
  mutate(
    EIR_cat = factor(
      EIR_cat,
      levels = c("Low", "Moderate", "High"),
      labels = c(
        "Low PfPR (<10%)",
        "Moderate PfPR (10–35%)",
        "High PfPR (>35%)"
      )
    )
  ) %>%
  group_by(
    EIR_cat,
    futNetcovstart2023,
    population
  ) %>%
  summarise(
    mean_reduction = mean(
      prev_reduction,
      na.rm = TRUE
    ),
    q25_reduction = quantile(
      prev_reduction,
      0.25,
      na.rm = TRUE
    ),
    q75_reduction = quantile(
      prev_reduction,
      0.75,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# ------------------------------------------------------------
# Plot Figure 2
# ------------------------------------------------------------

figure2_prevalence_reduction <- ggplot(
  figure2_annual_reduction,
  aes(
    x = futNetcovstart2023,
    y = mean_reduction,
    group = population,
    colour = population,
    fill = population
  )
) +
  
  # Interquartile range across simulation seeds
  geom_ribbon(
    aes(
      ymin = q25_reduction,
      ymax = q75_reduction
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
      "ITN users"
    ),
    name = NULL
  ) +
  
  scale_fill_manual(
    values = cols_pop,
    breaks = c(
      "ITN non-users",
      "ITN users"
    ),
    name = NULL
  ) +
  
  scale_x_continuous(
    breaks = seq(0, 0.8, by = 0.2),
    labels = scales::label_percent(
      accuracy = 1
    ),
    limits = c(0, 0.8),
    expand = expansion(
      mult = c(0.025, 0.04)
    )
  ) +
  
  scale_y_continuous(
    breaks = seq(0, 100, by = 25),
    limits = c(0, 102),
    expand = expansion(
      mult = c(0, 0.01)
    )
  ) +
  
  labs(
    x = "ITN coverage",
    y = "Relative prevalence reduction (%)"
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


# Display Figure 2
figure2_prevalence_reduction

#======================================
#save figure
#======================================

ggsave(
  filename = here::here("04_Figures", 
                        "Figure2_prevalence_reduction.tiff"),
  plot = figure2_prevalence_reduction,
  device = ragg::agg_tiff,
  width = 19.05,
  height = 10.5,
  units = "cm",
  res = 600,
  compression = "lzw",
  background = "white"
)
#====================================
#gap between users and non-users
#====================================

# gap btn user and non-user protection for fig2A
user_nonuser_gap <- figure2_annual_reduction|>
  pivot_wider(id_cols = c(EIR_cat, futNetcovstart2023),
              names_from = population, values_from = mean_reduction)|>
  mutate(gap = `ITN users` - `ITN non-users`)|>
  filter(futNetcovstart2023 %in% c(0.1, 0.8))


#################################
# Figure 4: Community protection and marginal gains
#
# Created: September 2026
#
# Purpose:
#   Quantify community-level protection experienced by ITN
#   non-users across transmission intensity and ITN usage.
#
#   Panel A:
#     Reduction in prevalence among ITN non-users.
#
#   Panel B:
#     Additional reduction in prevalence among ITN non-users
#     for each 10-percentage-point increase in ITN usage,
#     together with the usage level at which the marginal gain
#     is greatest.
#
# Starting data:
#   figure2_seed_data in figure2_prevalence_reduction.R script
#################################

library(tidyverse)
library(ggsci)
library(patchwork)
library(metR)
library(akima)


# ============================================================
# 1. PREPARE FIGURE 4 DATA
# ============================================================

# Figure 2 contains seed-level prevalence reduction estimates.
# For Figure 4, calculate the mean reduction over the
# 2023–2031 period and across simulation seeds for each
# EIR × ITN usage × population combination.

figure4_summary_data <- figure2_seed_data %>%
  filter(
    year >= 2023,
    year <= 2031
  ) %>%
  group_by(
    EIR,
    EIR_cat,
    population,
    futNetcovstart2023
  ) %>%
  summarise(
    mean_prev_reduction =
      mean(
        prev_reduction,
        na.rm = TRUE
      ),
    .groups = "drop"
  )


# ============================================================
# 2. NON-USER COMMUNITY PROTECTION DATA
# ============================================================

figure4_nonuser_data <- figure4_summary_data %>%
  filter(
    population == "ITN non-users"
  ) %>%
  rename(
    usage = futNetcovstart2023,
    reduction_nonuser = mean_prev_reduction
  ) %>%
  arrange(
    EIR,
    usage
  )


# ============================================================
# 3. MONOTONIC SMOOTHING
# ============================================================

# Prevalence reduction among non-users is expected to:
#
#   1. Increase with increasing ITN usage.
#   2. Decrease with increasing EIR at a fixed usage level.
#
# Two sequential isotonic regression passes are therefore
# applied to enforce these expected relationships.
#
# Pass 1:
#   Reduction is constrained to be non-decreasing with usage.
#
# Pass 2:
#   Reduction is constrained to be non-increasing with EIR.

figure4_nonuser_data <- figure4_nonuser_data %>%
  arrange(
    EIR,
    usage
  ) %>%
  group_by(
    EIR
  ) %>%
  mutate(
    reduction_nonuser =
      isoreg(
        usage,
        reduction_nonuser
      )$yf
  ) %>%
  ungroup() %>%
  arrange(
    usage,
    EIR
  ) %>%
  group_by(
    usage
  ) %>%
  mutate(
    reduction_nonuser =
      -isoreg(
        EIR,
        -reduction_nonuser
      )$yf
  ) %>%
  ungroup() %>%
  arrange(
    EIR,
    usage
  )


# ============================================================
# 4. CALCULATE MARGINAL COMMUNITY PROTECTION
# ============================================================

# Marginal gain is the additional prevalence reduction among
# non-users associated with the next 10-percentage-point
# increase in ITN usage.
#
# For example:
#   40% -> 50% usage
#
# marginal gain =
#   reduction at 50% - reduction at 40%
#
# There is no marginal gain calculated at 80% because a
# 90% usage level was not simulated.

figure4_nonuser_data <- figure4_nonuser_data %>%
  group_by(
    EIR
  ) %>%
  mutate(
    marginal_gain =
      lead(reduction_nonuser) -
      reduction_nonuser
  ) %>%
  ungroup()


# ============================================================
# 5. CREATE INTERPOLATION GRID
# ============================================================

figure4_usage_seq <- seq(
  min(figure4_nonuser_data$usage, na.rm = TRUE),
  max(figure4_nonuser_data$usage, na.rm = TRUE),
  length.out = 400
)

figure4_eir_seq <- 10^seq(
  log10(min(figure4_nonuser_data$EIR, na.rm = TRUE)),
  log10(max(figure4_nonuser_data$EIR, na.rm = TRUE)),
  length.out = 400
)


# Function to interpolate the observed data onto a fine grid

figure4_make_grid <- function(data, variable) {
  
  interpolation_data <- data %>%
    filter(
      !is.na(.data[[variable]])
    )
  
  interpolation_result <- akima::interp(
    x = interpolation_data$usage,
    y = interpolation_data$EIR,
    z = interpolation_data[[variable]],
    xo = figure4_usage_seq,
    yo = figure4_eir_seq,
    duplicate = "mean",
    linear = TRUE,
    extrap = FALSE
  )
  
  expand.grid(
    usage = interpolation_result$x,
    eir = interpolation_result$y
  ) %>%
    mutate(
      !!variable := as.vector(
        interpolation_result$z
      )
    )
}


# Interpolated grids for the two panels

figure4_reduction_grid <- figure4_make_grid(
  figure4_nonuser_data,
  "reduction_nonuser"
)

figure4_marginal_grid <- figure4_make_grid(
  figure4_nonuser_data,
  "marginal_gain"
)


# ============================================================
# 6. IDENTIFY USAGE LEVEL WITH MAXIMUM MARGINAL GAIN
# ============================================================

# Identify the usage level with the largest marginal gain
# separately for each EIR using the observed data.
#
# The loess curve is used only to provide a smooth visual
# representation of the peak usage across the EIR gradient.

figure4_peak_curve_raw <- figure4_nonuser_data %>%
  filter(
    !is.na(marginal_gain)
  ) %>%
  group_by(
    EIR
  ) %>%
  slice_max(
    marginal_gain,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  transmute(
    eir = EIR,
    usage
  )


figure4_peak_fit <- loess(
  usage ~ log10(eir),
  data = figure4_peak_curve_raw,
  span = 0.4
)


figure4_peak_curve <- data.frame(
  eir = figure4_eir_seq
) %>%
  mutate(
    usage =
      predict(
        figure4_peak_fit,
        newdata = data.frame(
          eir = figure4_eir_seq
        )
      ),
    
    # Prevent the loess curve from extending beyond
    # the observed usage range.
    usage =
      pmin(
        pmax(
          usage,
          min(
            figure4_peak_curve_raw$usage,
            na.rm = TRUE
          )
        ),
        max(
          figure4_peak_curve_raw$usage,
          na.rm = TRUE
        )
      )
  )


# ============================================================
# 7. TRANSMISSION THRESHOLDS
# ============================================================

# EIR thresholds corresponding to the PfPR categories:
#
#   Low:       PfPR <10%
#   Moderate:  PfPR 10–35%
#   High:      PfPR >35%

figure4_thresholds <- data.frame(
  eir = c(2, 8),
  label = c(
    "PfPR 10%",
    "PfPR 35%"
  )
)


# ============================================================
# 8. COMMON HEATMAP SETTINGS
# ============================================================

figure4_base_layers <- list(
  
  geom_hline(
    data = figure4_thresholds,
    aes(
      yintercept = eir
    ),
    colour = "white",
    linetype = "dashed",
    linewidth = 0.35
  ),
  
  scale_y_log10(
    breaks = c(
      1, 2, 8, 20, 50, 100
    ),
    labels = c(
      1, 2, 8, 20, 50, 100
    ),
    expand = c(0, 0)
  ),
  
  scale_x_continuous(
    breaks = seq(
      10,
      80,
      10
    ),
    labels = function(x) {
      paste0(
        x,
        "%"
      )
    },
    expand = c(0, 0)
  ),
  
  labs(
    x = "ITN usage"
  ),
  
  theme_pub(),
  
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.key.width = unit(
      1.4,
      "cm"
    ),
    legend.title = element_text(
      size = 9
    ),
    legend.title.position = "top"
  )
)


# ============================================================
# FIGURE 4A
# COMMUNITY PROTECTION LEVEL
# ============================================================

figure2B <- ggplot(
  figure4_reduction_grid,
  aes(
    usage * 100,
    eir
  )
) +
  
  geom_raster(
    aes(
      fill = reduction_nonuser
    ),
    interpolate = TRUE
  ) +
  
  geom_contour(
    aes(
      z = reduction_nonuser
    ),
    breaks = c(
      20,
      40,
      60,
      80
    ),
    colour = "black",
    linewidth = 0.3
  ) +
  
  geom_text_contour(
    aes(
      z = reduction_nonuser
    ),
    breaks = c(
      20,
      40,
      60,
      80
    ),
    size = 2.6,
    stroke = 0.2,
    colour = "black"
  ) +
  
  scale_fill_viridis_c(
    limits = c(
      0,
      100
    ),
    name = "Reduction in prevalence among non-users (%)"
  ) +
  
  figure4_base_layers +
  
  labs(
    x = "ITN coverage",
    y = "Annual EIR (ib/person/year)"
  )


# ============================================================
# FIGURE 2C
# MARGINAL COMMUNITY PROTECTION
# ============================================================
transition_labels <- figure4_marginal_grid %>%
  filter(usage <= 0.7) %>%
  distinct(usage) %>%
  arrange(usage) %>%
  mutate(
    transition = paste0(
      seq(0, 60, length.out = n()),
      "% → ",
      seq(10, 70, length.out = n()),
      "%"
    )
  )


# Actual usage values corresponding to the 10-percentage-point transitions
x_breaks <- figure4_marginal_grid %>%
  filter(usage <= 0.7) %>%
  pull(usage) %>%
  unique() %>%
  sort()

# Find the usage value closest to each desired transition endpoint
target_usage <- seq(0.1, 0.7, by = 0.1)

x_breaks_transition <- sapply(
  target_usage,
  function(x) x_breaks[which.min(abs(x_breaks - x))]
)

x_labels <- c(
  "0% → 10%",
  "10% → 20%",
  "20% → 30%",
  "30% → 40%",
  "40% → 50%",
  "50% → 60%",
  "60% → 70%"
)
# Plot
figure2C <- ggplot(
  figure4_marginal_grid %>%
    filter(usage <= 0.7),
  aes(
    x = usage * 100,
    y = eir
  )
) +
  
  geom_raster(
    aes(fill = marginal_gain),
    interpolate = TRUE
  ) +
  
  geom_contour(
    aes(z = marginal_gain),
    breaks = c(10, 15, 25),
    colour = "white",
    linewidth = 0.25
  ) +
  
  geom_text_contour(
    aes(z = marginal_gain),
    breaks = c(10, 15, 25),
    size = 2.6,
    stroke = 0.2,
    colour = "black"
  ) +
  
  scale_fill_viridis_c(
    option = "magma",
    name = "Additional non-user reduction per 10-percentage-point increase in coverage"
  ) +
  
  scale_x_continuous(
    breaks = x_breaks_transition * 100,
    labels = x_labels
  ) +
  
  figure4_base_layers +
  
  labs(
    x = "ITN coverage",
    y = "Annual EIR (ib/person/year)"
  ) +
  
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 1,
      vjust = 1
    )
  )

figure2C
# ============================================================
# COMBINE FIGURE 4
# ============================================================

# ============================================================
# COMBINE FIGURE 3
# ============================================================

figure2A <- figure2_prevalence_reduction +
  labs(tag = "A") +
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 14
    ),
    plot.tag.position = c(0, 1)
  )


figure2B <- figure2B +
  labs(tag = "B") +
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 14
    ),
    plot.tag.position = c(0, 1)
  )


figure2C <- figure2C +
  labs(tag = "C") +
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 14
    ),
    plot.tag.position = c(-0.01, 1)
  )


pbc <- (
  figure2B +
    plot_spacer() +
    figure2C
) +
  plot_layout(
    widths = c(1, 0.04, 1)
  ) &
  theme(
    legend.position = "bottom"
  )


figure2 <- figure2A / pbc +
  plot_layout(
    heights = c(1.3, 1)
  )


# Display Figure 3
figure2

# ============================================================
# SAVE FIGURE 4
# ============================================================

ggsave(
  filename = here::here(
    "04_Figures",
    "Fig2.tiff"
  ),
  plot = figure2,
  device = ragg::agg_tiff,
  width = 22,
  height = 19,
  units = "cm",
  res = 600,
  compression = "lzw",
  background = "white"
)

