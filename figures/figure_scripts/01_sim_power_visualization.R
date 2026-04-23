library(ggplot2)
library(dplyr)

# 1. Load the results
df <- read.csv("01_sim_power_results.csv")

# 2. Preparation: Clean up names and REMOVE PGLS
df_plot <- df %>%
  mutate(method = case_when(
    method == "anccond.01" ~ "AncCond (0->1)",
    method == "anccond.10" ~ "AncCond (1->0)",
    method == "phyloglm"   ~ "phyloglm",
    TRUE ~ method
  )) %>%
  # Filter to keep only the three methods you want
  filter(method %in% c("AncCond (0->1)", "AncCond (1->0)", "phyloglm")) %>%
  filter(!is.na(rate))

# Factorizing to ensure the order in the legend matches your preference
df_plot$method <- factor(df_plot$method, levels = c(
  "AncCond (0->1)", 
  "AncCond (1->0)", 
  "phyloglm"
))

# 3. Labeling for the facets
n_labs <- setNames(paste0("n = ", unique(df_plot$n.tips)), unique(df_plot$n.tips))
q_labs <- setNames(paste0("q = ", unique(df_plot$q.rate)), unique(df_plot$q.rate))

# 4. Create the Plot
p <- ggplot(df_plot, aes(x = s, y = rate, color = method, group = method)) +
  # Significance baseline (Alpha = 0.05)
  geom_hline(yintercept = 0.05, linetype = "dotted", color = "black", alpha = 0.6) +
  
  # Lines and Points
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  
  # THE GRID
  facet_grid(q.rate ~ n.tips, 
             labeller = labeller(q.rate = as_labeller(q_labs), 
                                 n.tips = as_labeller(n_labs))) +
  
  # MATCHING FPR COLORS: Using specific Viridis Hex Codes
  scale_color_manual(values = c(
    "AncCond (0->1)" = "#440154FF", # Deep Purple
    "AncCond (1->0)" = "#21908CFF", # Teal
    "phyloglm"       = "#FDE725FF"  # Bright Yellow
  )) +
  
  # Visual refinements
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25)) +
  scale_x_continuous(breaks = unique(df_plot$s)) +
  labs(
    title = "Method Performance: Power Comparison",
    subtitle = "Matched colors to FPR Plot (Purple/Teal = AncCond, Yellow = phyloglm)",
    x = "Signal Strength (s)",
    y = "Rate (n.sig / n.valid)",
    color = "Method"
  ) +
  theme_bw() + 
  theme(
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

# 5. Display and Save
print(p)
ggsave("anccond_power.png", p, width = 12, height = 10, dpi = 300)