# start from uncertainty script (with EIR_cat) or prev_vs_eir_plot script (with continuous EIR)
library(dplyr)
library(scales)

# from main simulations
dat_main <- figure2_seed_data %>%
  filter(
    futNetcovstart2023==0.5,
    year>2022 & year<2032,
    EIR %in% c(1.5,4, 48), 
    population %in% c("ITN non-users", "ITN users")#, "LLIN users"
  ) %>%
  group_by(EIR, population) %>%
  summarise(
    prev_red = median(prev_reduction, na.rm = TRUE),
    q25  = quantile(prev_reduction, 0.25, na.rm = TRUE),
    q75  = quantile(prev_reduction, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# for alternative simulation
# from main simulations
dat_robus <- fig2_seed_EIR_cat_sen %>%
  filter(
    futNetcovstart2023==0.5,
    year>2022 & year<2032,
    population %in% c("LLIN non-users", "LLIN users")#, "LLIN users"
  ) %>%
  group_by(EIR, population) %>%
  summarise(
    prev_red_robus = median(prev_reduction, na.rm = TRUE),
    q25_robus = quantile(prev_reduction, 0.25, na.rm = TRUE),
    q75_robus  = quantile(prev_reduction, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

#bind columns
joint_dat <- right_join(dat_main, dat_robus, by = c("EIR", "population"))%>%
  filter(population == "LLIN non-users")%>%
  mutate(diff_red = prev_red - prev_red_robus)


library(tidyr); library(dplyr); library(ggplot2)

df_long <- joint_dat %>%
  transmute(EIR, population,
            orig_mid = prev_red, orig_lo = q25, orig_hi = q75,
            robus_mid = prev_red_robus, robus_lo = q25_robus, robus_hi = q75_robus) %>%
  pivot_longer(-c(EIR, population),
               names_to = c("type", ".value"),
               names_pattern = "(orig|robus)_(.*)") %>%
  mutate(type = recode(type, orig = "Main simulation", robus = "Robustness check"))



# flag whether IQRs overlap, for the annotation
overlap_flag <- joint_dat %>%
  mutate(overlap = !(q25 > q75_robus | q25_robus > q75)) %>%
  select(EIR, overlap)

# small multiplicative dodge so points don't sit exactly on top of each other on a log axis
dodge_factor <- 1.05
df_long <- df_long %>%
  mutate(EIR_dodge = ifelse(type == "Main simulation", EIR / dodge_factor, EIR * dodge_factor))

# annotation position: just above the EIR=4 "Original" point
ann <- df_long %>% filter(EIR == 4, type == "Main simulation")

robus <- ggplot(df_long, aes(x = EIR_dodge, y = mid, ymin = lo, ymax = hi, color = type)) +
  geom_pointrange(size = 0.6, linewidth = 0.9) +
  # annotate("text", x = 4, y = ann$hi + 6,
  #          label = "Non-overlapping IQRs", size = 3.3, fontface = "italic") +
  scale_x_log10(breaks = c(1.5, 4, 48), labels = c("1.5", "4", "48")) +
  scale_color_manual(values = c("Main simulation" = "#F8766D", "Robustness check" = "#00BFC4")) +
  labs(x = "EIR (log scale)",
       y = "Relative prevalence reduction\namong non-users (%)",
       color = NULL) +
  theme_minimal(base_size = 13) +
  theme_pub()

ggsave(
  "Fig S2.tiff",
  plot = robus,
  device = ragg::agg_tiff,
  width = 19,
  height = 17,
  units = "cm",
  res = 600,
  compression = "lzw",
  background = "white"
)



