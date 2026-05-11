library(ggplot2)
library(dplyr)

# 1. Load Symmetrical Results (includes AncCond 1->0 and phyloglm)
df_sym <- read.csv("01_sim_power_results.csv") %>%
  filter(n.tips == 500, q.rate %in% c(0.5, 1.0, 2.0)) %>%
  rename(q_val = q.rate) %>%
  filter(method %in% c("anccond.10", "phyloglm")) %>%
  mutate(line_group = case_when(
    method == "anccond.10" ~ "AncCond 1->0 (Symmetrical)",
    method == "phyloglm"   ~ "Phyloglm"
  ))

# 2. Load Asymmetrical Results (AncCond 1->0)
df_asym <- read.csv("01_sim_asymm_power_results.csv") %>%
  filter(n.tips == 500, q01 %in% c(0.5, 1.0, 2.0)) %>%
  rename(q_val = q01) %>%
  filter(method == "anccond.10") %>%
  mutate(line_group = "AncCond 1->0 (Asymmetrical)")

# 3. Combine and Clean
df_final <- bind_rows(df_sym, df_asym) %>%
  mutate(
    q_label = paste0("q01 = ", sprintf("%.1f", q_val)),
    # Ensure factor levels match the color manual exactly to show in legend
    line_group = factor(line_group, levels = c(
      "AncCond 1->0 (Symmetrical)", 
      "AncCond 1->0 (Asymmetrical)", 
      "Phyloglm"
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
  
  # COLORS: Teal (Sym), Purple (Asym), Golden Yellow (phyloglm)
  scale_color_manual(
    name = "Simulation Method",
    values = c(
      "AncCond 1->0 (Symmetrical)"  = "#21908CFF", # Teal
      "AncCond 1->0 (Asymmetrical)" = "#440154FF", # Purple
      "Phyloglm"                    = "#D4A017"    # Golden Yellow
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
    title = "Power Comparison for 1->0 Transitions",
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
  # Boost legend icon size for visibility
  guides(color = guide_legend(override.aes = list(linewidth = 1.5, size = 4)))

# 5. Display and Save
print(p)
ggsave("anccond_h2h_final_titled.png", p, width = 14, height = 6, dpi = 300)