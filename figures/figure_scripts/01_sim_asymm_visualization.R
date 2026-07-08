library(ggplot2)
library(dplyr)

# 1. Load Symmetrical Results (Extract BOTH AncCond 1->0 and phyloglm)
df_sym <- read.csv("01_sim_power_results.csv") %>%
  filter(n.tips == 500, q.rate %in% c(0.5, 1.0, 2.0)) %>%
  rename(q_val = q.rate) %>%
  filter(method %in% c("anccond.10", "phyloglm")) %>%
  mutate(line_group = case_when(
    method == "anccond.10" ~ "AncCond q10 = q01",
    method == "phyloglm"   ~ "Phyloglm q10 = q01"
  ))

# 2. Load Asymmetrical Results (Extract BOTH AncCond 1->0 and phyloglm)
df_asym <- read.csv("01_sim_asymm_power_results.csv") %>%
  filter(n.tips == 500, q01 %in% c(0.5, 1.0, 2.0)) %>%
  rename(q_val = q01) %>%
  filter(method %in% c("anccond.10", "phyloglm")) %>%
  mutate(line_group = case_when(
    method == "anccond.10" ~ "AncCond q10 = 0.1",
    method == "phyloglm"   ~ "Phyloglm q10 = 0.1"
  ))

# 3. Combine and Clean
df_final <- bind_rows(df_sym, df_asym) %>%
  mutate(
    q_label = paste0("q01 = ", sprintf("%.1f", q_val)),
    # Factor levels grouped by regime for a clean legend
    line_group = factor(line_group, levels = c(
      "AncCond q10 = q01", 
      "Phyloglm q10 = q01",
      "AncCond q10 = 0.1", 
      "Phyloglm q10 = 0.1"
    ))
  )

# 4. Create the Plot
p <- ggplot(df_final, aes(x = s, y = rate, color = line_group, group = line_group)) +
  # Significance baseline
  geom_hline(yintercept = 0.05, linetype = "dotted", color = "black", alpha = 0.8) +
  
  geom_line(linewidth = 1.4) +
  geom_point(size = 3) +
  
  # HORIZONTAL FACETS
  facet_grid(. ~ q_label) +
  
  # COLORS: Deep colors for AncCond, lighter tints for phyloglm
  scale_color_manual(
    name = "Simulation Method",
    values = c(
      "AncCond q10 = q01"  = "#21908CFF", # Deep Teal
      "Phyloglm q10 = q01"      = "#7EC7C4",   # Light Pastel Teal
      "AncCond q10 = 0.1" = "#440154FF", # Deep Purple
      "Phyloglm q10 = 0.1"     = "#9B6AA6"    # Light Pastel Purple
    )
  ) +
  
  # Visual refinements
  scale_y_continuous(
    limits = c(0, 1.05), 
    breaks = seq(0, 1, 0.25),
    expand = expansion(mult = c(0, 0.05)) # Y-axis starts exactly at 0
  ) +
  scale_x_continuous(breaks = c(1, 2, 3, 5, 10)) +
  
  labs(
    title = "Power Comparison of Aysmmetric Transition Rates",
    x = "Signal Strength (s)",
    y = "Power",
    color = "Method"
  ) +
  
  theme_bw() + 
  theme(
    plot.title = element_text(face = "bold", size = 16),
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11)
  ) +
  # Stack legend in 2 columns for a clean 2x2 block
  guides(color = guide_legend(override.aes = list(linewidth = 1.5, size = 4), ncol = 2))

# 5. Display and Save
print(p)
#ggsave("anccond_h2h_final_titled.png", p, width = 14, height = 6, dpi = 300)