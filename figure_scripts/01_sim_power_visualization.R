library(ggplot2)
library(dplyr)

# 1. Load the results
df <- read.csv("01_sim_power_results.csv")

# 2. Preparation: Clean up method names and ensure consistent factor levels
df_plot <- df %>%
  mutate(method = case_when(
    method == "anccond.01" ~ "AncCond (0->1)",
    method == "anccond.10" ~ "AncCond (1->0)",
    method == "phyloglm"   ~ "phyloglm",
    method == "pgls"       ~ "PGLS",
    TRUE ~ method
  )) %>%
  filter(!is.na(rate))

# Factorizing ensures the colors are assigned to the correct methods every time
df_plot$method <- factor(df_plot$method, levels = c(
  "AncCond (0->1)", 
  "AncCond (1->0)", 
  "phyloglm", 
  "PGLS"
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
  
  # THE GRID: Rows = Transition Rate (q), Columns = Taxa Size (n)
  facet_grid(q.rate ~ n.tips, 
             labeller = labeller(q.rate = as_labeller(q_labs), 
                                 n.tips = as_labeller(n_labs))) +
  
  # CALL DEFAULT COLORS MANUALLY
  # These are the exact hex codes ggplot2 uses for a 4-color palette
  scale_color_manual(values = c(
    "AncCond (0->1)" = "#F8766D", # Default Red
    "AncCond (1->0)" = "#7CAE00", # Default Green
    "phyloglm"       = "#00BFC4", # Default Cyan/Blue
    "PGLS"           = "#C77CFF"  # Default Purple
  )) +
  
  # Visual refinements
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25)) +
  scale_x_continuous(breaks = unique(df_plot$s)) +
  labs(
    title = "Method Performance: Directional Rejection Rates",
    subtitle = "Y-axis: n.sig / n.valid | Rows: Transition Rate (q) | Cols: Taxa Size (n)",
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
ggsave("anccond_dots_grid_default_colors.png", p, width = 12, height = 10, dpi = 300)