#################################
# Figure 3: Population-level impact of ITNs
#
# Created: September 2026
#
# Purpose:
#   Quantify the population-level impact of ITN use,
#   including:
#     A. Infections averted among users and non-users
#     B. Incremental infections averted per 10 percentage-point
#        increase in ITN usage
#     C. Infections averted per ITN
#
# Starting data:
#   figure2_seed_data in figure2_prevalence_reduction.R script
#################################

library(tidyverse)
library(ggsci)
library(patchwork)


# ============================================================
# GENERAL SETTINGS
# ============================================================

# Population denominator for reporting
population_denominator <- 1000

# ITN colours
population_impact_colours <- c(
  "Remaining infections" = "white",
  "Infections averted among ITN users" = "#432CA1",
  "Infections averted among non-users" = "#D55E00"
)

population_impact_order <- c(
  "Infections averted among ITN users",
  "Infections averted among non-users",
  "Remaining infections"
)


# ============================================================
# 1. PREPARE FIGURE 2 DATA
# ============================================================

figure2_seed_data <- figure2_seed_data %>%
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
  )


# ============================================================
# 2. NON-USER PREVALENCE
# ============================================================

figure3_nonuser_prevalence <- figure2_seed_data %>%
  filter(
    year >= 2023,
    year <= 2031,
    population == "ITN non-users"
  ) %>%
  group_by(
    seed,
    EIR,
    EIR_cat,
    futNetcovstart2023
  ) %>%
  summarise(
    nonuser_prevalence = mean(
      prevalenceRate,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# ============================================================
# 3. BASELINE NON-USER PREVALENCE
# ============================================================

figure3_baseline_nonuser <- figure3_nonuser_prevalence %>%
  filter(
    futNetcovstart2023 == 0
  ) %>%
  select(
    seed,
    EIR,
    EIR_cat,
    baseline_nonuser_prevalence = nonuser_prevalence
  )


# ============================================================
# 4. INFECTIONS AVERTED AMONG NON-USERS
#
# Reported per 1,000 total population.
# ============================================================

figure3_nonuser_impact <- figure3_nonuser_prevalence %>%
  left_join(
    figure3_baseline_nonuser,
    by = c(
      "seed",
      "EIR",
      "EIR_cat"
    )
  ) %>%
  mutate(
    
    # Number of non-users per 1,000 total population
    nonusers_per_1000 =
      population_denominator *
      (1 - futNetcovstart2023),
    
    # Infections averted among non-users
    # per 1,000 total population
    infections_averted_nonusers =
      (
        baseline_nonuser_prevalence -
          nonuser_prevalence
      ) *
      nonusers_per_1000
    
  )


# ============================================================
# 5. USER PREVALENCE
# ============================================================

figure3_user_prevalence <- figure2_seed_data %>%
  filter(
    year >= 2023,
    year <= 2031,
    population == "ITN users"
  ) %>%
  group_by(
    seed,
    EIR,
    EIR_cat,
    futNetcovstart2023
  ) %>%
  summarise(
    user_prevalence = mean(
      prevalenceRate,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# ============================================================
# 6. INFECTIONS AVERTED AMONG USERS
#
# Reported per 1,000 total population.
# ============================================================

figure3_user_impact <- figure3_user_prevalence %>%
  left_join(
    figure3_baseline_nonuser,
    by = c(
      "seed",
      "EIR",
      "EIR_cat"
    )
  ) %>%
  mutate(
    
    # Number of ITN users per 1,000 total population
    users_per_1000 =
      population_denominator *
      futNetcovstart2023,
    
    # Infections averted among users
    # per 1,000 total population
    infections_averted_users =
      (
        baseline_nonuser_prevalence -
          user_prevalence
      ) *
      users_per_1000
    
  )


# ============================================================
# 7. TOTAL POPULATION IMPACT
# ============================================================

figure3_population_impact <- figure3_nonuser_impact %>%
  select(
    seed,
    EIR,
    EIR_cat,
    futNetcovstart2023,
    infections_averted_nonusers,
    baseline_nonuser_prevalence
  ) %>%
  left_join(
    figure3_user_impact %>%
      select(
        seed,
        EIR,
        EIR_cat,
        futNetcovstart2023,
        infections_averted_users
      ),
    by = c(
      "seed",
      "EIR",
      "EIR_cat",
      "futNetcovstart2023"
    )
  ) %>%
  mutate(
    
    # Total infections averted per 1,000 population
    total_infections_averted =
      infections_averted_users +
      infections_averted_nonusers,
    
    # Baseline infections per 1,000 population
    baseline_infections_per_1000 =
      baseline_nonuser_prevalence *
      population_denominator,
    
    # Remaining infections per 1,000 population
    remaining_infections =
      baseline_infections_per_1000 -
      total_infections_averted
    
  )


# ============================================================
# 8. SUMMARISE ACROSS SEEDS
# ============================================================

figure3_population_impact_summary <- figure3_population_impact %>%
  group_by(
    EIR_cat,
    futNetcovstart2023
  ) %>%
  summarise(
    
    # --------------------------------------------------------
    # Total population impact
    # --------------------------------------------------------
    
    mean_total_infections_averted =
      mean(
        total_infections_averted,
        na.rm = TRUE
      ),
    
    q25_total_infections_averted =
      quantile(
        total_infections_averted,
        0.25,
        na.rm = TRUE
      ),
    
    q75_total_infections_averted =
      quantile(
        total_infections_averted,
        0.75,
        na.rm = TRUE
      ),
    
    # --------------------------------------------------------
    # Direct impact among users
    # --------------------------------------------------------
    
    mean_infections_averted_users =
      mean(
        infections_averted_users,
        na.rm = TRUE
      ),
    
    q25_infections_averted_users =
      quantile(
        infections_averted_users,
        0.25,
        na.rm = TRUE
      ),
    
    q75_infections_averted_users =
      quantile(
        infections_averted_users,
        0.75,
        na.rm = TRUE
      ),
    
    # --------------------------------------------------------
    # community impact among non-users
    # --------------------------------------------------------
    
    mean_infections_averted_nonusers =
      mean(
        infections_averted_nonusers,
        na.rm = TRUE
      ),
    
    q25_infections_averted_nonusers =
      quantile(
        infections_averted_nonusers,
        0.25,
        na.rm = TRUE
      ),
    
    q75_infections_averted_nonusers =
      quantile(
        infections_averted_nonusers,
        0.75,
        na.rm = TRUE
      ),
    
    # --------------------------------------------------------
    # user and non-user percentage contribution 
    # --------------------------------------------------------
    
    percent_user_contribution = (mean_infections_averted_users/
      mean_total_infections_averted)*100,
    
    percent_nonuser_contribution = (mean_infections_averted_nonusers/
                                   mean_total_infections_averted)*100,
    
    # --------------------------------------------------------
    # Remaining infections
    # --------------------------------------------------------
    
    mean_remaining_infections =
      mean(
        remaining_infections,
        na.rm = TRUE
      ),
    
    q25_remaining_infections =
      quantile(
        remaining_infections,
        0.25,
        na.rm = TRUE
      ),
    
    q75_remaining_infections =
      quantile(
        remaining_infections,
        0.75,
        na.rm = TRUE
      ),
    
    counterfactual_infections = (mean_infections_averted_nonusers +
                                        mean_infections_averted_users +
                                        mean_remaining_infections),
    
    percent_counterfactual_averted = (mean_total_infections_averted /
                                        counterfactual_infections)*100,
    
    .groups = "drop"
  )


# ============================================================
# FIGURE 3A
# Population-level infections:
# infections averted among users +
# infections averted among non-users +
# remaining infections
# ============================================================

figure3A_data <- figure3_population_impact_summary %>%
  select(
    EIR_cat,
    futNetcovstart2023,
    mean_infections_averted_users,
    mean_infections_averted_nonusers,
    mean_remaining_infections
  ) %>%
  filter(
    futNetcovstart2023 > 0
  ) %>%
  pivot_longer(
    cols = c(
      mean_infections_averted_users,
      mean_infections_averted_nonusers,
      mean_remaining_infections
    ),
    names_to = "impact_component",
    values_to = "infections"
  ) %>%
  mutate(
    
    impact_component = case_when(
      
      impact_component ==
        "mean_infections_averted_users" ~
        "Infections averted among ITN users",
      
      impact_component ==
        "mean_infections_averted_nonusers" ~
        "Infections averted among non-users",
      
      impact_component ==
        "mean_remaining_infections" ~
        "Remaining infections"
      
    ),
    
    impact_component = factor(
      impact_component,
      levels = population_impact_order
    )
    
  )


figure3A <- ggplot(
  figure3A_data,
  aes(
    x = futNetcovstart2023,
    y = infections,
    fill = impact_component,
    group = impact_component
  )
) +
  
  geom_col(
    position = position_stack(reverse = TRUE),
    width = 0.08,
    colour = "grey8",
    linewidth = 0.2
  ) +
  
  facet_wrap(
    ~ EIR_cat,
    nrow = 1
  ) +
  
  scale_fill_manual(
    values = population_impact_colours,
    breaks = population_impact_order
  ) +
  
  scale_x_continuous(
    breaks = seq(0, 0.8, 0.2),
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  
  labs(
    x = "ITN usage",
    y = "Infections per 1,000 population",
    fill = NULL
  ) +
  
  theme_pub() +
  
  theme(
    legend.position = "bottom"
  )


figure3A


# ============================================================
# FIGURE 3B
# Incremental infections averted per
# 10-percentage-point increase in ITN usage
# ============================================================

figure3_incremental_impact <- figure3_population_impact %>%
  group_by(
    seed,
    EIR,
    EIR_cat
  ) %>%
  arrange(
    futNetcovstart2023,
    .by_group = TRUE
  ) %>%
  mutate(
    
    coverage_difference =
      futNetcovstart2023 -
      lag(futNetcovstart2023),
    
    infections_averted_difference =
      total_infections_averted -
      lag(total_infections_averted),
    
    # Convert the observed change between coverage
    # levels to the equivalent impact of a
    # 10-percentage-point increase.
    incremental_infections_averted_per_10pp =
      (
        infections_averted_difference /
          coverage_difference
      ) * 0.10
    
  ) %>%
  ungroup() %>%
  
  group_by(
    EIR_cat,
    futNetcovstart2023
  ) %>%
  summarise(
    
    inc_benefit_mean =
      mean(
        incremental_infections_averted_per_10pp,
        na.rm = TRUE
      ),
    
    inc_benefit_q25 =
      quantile(
        incremental_infections_averted_per_10pp,
        0.25,
        na.rm = TRUE
      ),
    
    inc_benefit_q75 =
      quantile(
        incremental_infections_averted_per_10pp,
        0.75,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


figure3B <- ggplot(
  figure3_incremental_impact %>%
    filter(
      futNetcovstart2023 > 0.1
    ),
  aes(
    x = futNetcovstart2023,
    y = inc_benefit_mean,
    colour = EIR_cat,
    fill = EIR_cat,
    group = EIR_cat
  )
) +
  
  geom_line(
    linewidth = 0.8
  ) +
  
  geom_ribbon(
    aes(
      ymin = inc_benefit_q25,
      ymax = inc_benefit_q75
    ),
    alpha = 0.18,
    colour = NA,
    show.legend = FALSE
  ) +
  
  geom_point(
    size = 1.5
  ) +
  
  scale_color_bmj(
    palette = "default",
    alpha = 0.7
  ) +
  
  scale_fill_bmj(
    palette = "default",
    alpha = 0.7
  ) +
  
  scale_x_continuous(
    breaks = seq(0, 0.8, by = 0.2),
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  
  labs(
    x = "ITN usage reached",
    y = expression(
      atop(
        "Incremental infections averted per",
        "10-percentage-point increase in usage"
      )
    ),
    colour = NULL,
    fill = NULL
  ) +
  
  theme_pub() +
  
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


figure3B


# ============================================================
# FIGURE 3C
# Infections averted per ITN
# ============================================================

figure3_per_net <- figure3_population_impact %>%
  mutate(
    
    # Number of ITNs required per 1,000 population,
    # assuming 1 ITN per 1.8 people.
    nets_per_1000 =
      (1000 * futNetcovstart2023) / 1.8,
    
    infections_averted_per_itn =
      total_infections_averted /
      nets_per_1000
    
  ) %>%
  filter(
    futNetcovstart2023 > 0
  ) %>%
  group_by(
    EIR_cat,
    futNetcovstart2023
  ) %>%
  summarise(
    
    per_itn_mean =
      mean(
        infections_averted_per_itn,
        na.rm = TRUE
      ),
    
    per_itn_q25 =
      quantile(
        infections_averted_per_itn,
        0.25,
        na.rm = TRUE
      ),
    
    per_itn_q75 =
      quantile(
        infections_averted_per_itn,
        0.75,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


figure3C <- ggplot(
  figure3_per_net,
  aes(
    x = futNetcovstart2023,
    y = per_itn_mean,
    colour = EIR_cat,
    fill = EIR_cat,
    group = EIR_cat
  )
) +
  
  geom_line(
    linewidth = 0.8
  ) +
  
  geom_ribbon(
    aes(
      ymin = per_itn_q25,
      ymax = per_itn_q75
    ),
    alpha = 0.18,
    colour = NA,
    show.legend = FALSE
  ) +
  
  geom_point(
    size = 1.5
  ) +
  
  scale_color_bmj(
    palette = "default",
    alpha = 0.7
  ) +
  
  scale_fill_bmj(
    palette = "default",
    alpha = 0.7
  ) +
  
  scale_x_continuous(
    breaks = seq(0, 0.8, by = 0.2),
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  
  labs(
    x = "ITN usage",
    y = "Infections averted per ITN",
    colour = NULL,
    fill = NULL
  ) +
  
  theme_pub() +
  
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


figure3C


# ============================================================
# COMBINE FIGURE 3
# ============================================================

figure3A <- figure3A +
  labs(tag = "A") +
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 14
    ),
    plot.tag.position = c(0, 1)
  )


figure3B <- figure3B +
  labs(tag = "B") +
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 14
    ),
    plot.tag.position = c(0, 1)
  )


figure3C <- figure3C +
  labs(tag = "C") +
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 14
    ),
    plot.tag.position = c(0, 1)
  )


pbc <- (
  figure3B +
    plot_spacer() +
    figure3C
) +
  plot_layout(
    widths = c(1, 0.08, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom"
  )


figure3 <- figure3A / pbc +
  plot_layout(
    heights = c(1, 1)
  )


# Display Figure 3
figure3


# ============================================================
# SAVE FIGURE
# ============================================================

ggsave(
  filename = here::here(
    "04_Figures",
    "Figure_3_population_impact.tiff"
  ),
  plot = figure3,
  width = 18,
  height = 17,
  units = "cm",
  dpi = 600,
  compression = "lzw"
)

#=======================================
# filter for reporting and supplementary
#=======================================
contrib <- figure3_population_impact_summary%>%
  select(EIR_cat, futNetcovstart2023, mean_total_infections_averted,
         q25_total_infections_averted, q75_total_infections_averted,
         percent_user_contribution, percent_nonuser_contribution,
         percent_counterfactual_averted)%>%
  #filter(futNetcovstart2023 == 0.5)
  filter(futNetcovstart2023 %in% c(0.1, 0.8))

