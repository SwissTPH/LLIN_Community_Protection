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

figure4A <- ggplot(
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
    name = "Reduction in prevalence among\nnon-users (%)"
  ) +
  
  figure4_base_layers +
  
  labs(
    x = "ITN usage",
    y = "Annual EIR (ib/person/year)"
  )


# ============================================================
# FIGURE 4B
# MARGINAL COMMUNITY PROTECTION
# ============================================================

figure4B <- ggplot(
  figure4_marginal_grid %>%
    filter(
      usage <= 0.7
    ),
  aes(
    usage * 100,
    eir
  )
) +
  
  geom_raster(
    aes(
      fill = marginal_gain
    ),
    interpolate = TRUE
  ) +
  
  geom_contour(
    aes(
      z = marginal_gain
    ),
    breaks = c(
      10,
      15,
      25
    ),
    colour = "white",
    linewidth = 0.25
  ) +
  
  geom_text_contour(
    aes(
      z = marginal_gain
    ),
    breaks = c(
      10,
      15,
      25
    ),
    size = 2.6,
    stroke = 0.2,
    colour = "black"
  ) +
  
  # Cyan line indicates the usage level at which the
  # marginal community protection gain is greatest.
  geom_line(
    data = figure4_peak_curve,
    aes(
      usage * 100,
      eir
    ),
    colour = "#00D0FF",
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  
  scale_fill_viridis_c(
    option = "magma",
    name = "Additional reduction per\n10-percentage-point increase (pp)"
  ) +
  
  figure4_base_layers +
  
  labs(
    x = "ITN usage",
    y = NULL
  )

figure4B
# ============================================================
# COMBINE FIGURE 4
# ============================================================

figure4 <- (
  figure4A +
    figure4B
) +
  plot_annotation(
    tag_levels = "A"
  )


# ============================================================
# DISPLAY FIGURE 4
# ============================================================

figure4


# ============================================================
# SAVE FIGURE 4
# ============================================================

ggsave(
  filename = here::here(
    "04_Figures",
    "Figure_4_community_protection.tiff"
  ),
  plot = figure4,
  device = ragg::agg_tiff,
  width = 22,
  height = 10,
  units = "cm",
  res = 600,
  compression = "lzw",
  background = "white"
)
