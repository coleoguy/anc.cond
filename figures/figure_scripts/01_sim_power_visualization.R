library(ggplot2)
library(dplyr)

# 1. Load the results
df <- read.csv("01_sim_power_results.csv")

# 2. Preparation: Clean up names, remove PGLS, and filter specific parameters
df_plot <- df %>%
  mutate(method = case_when(
    method == "anccond.01" ~ "AncCond (0->1)",
    method == "anccond.10" ~ "AncCond (1->0)",
    method == "phyloglm"   ~ "phyloglm",
    TRUE ~ method
  )) %>%
  # Filter to keep only the three methods you want
  filter(method %in% c("AncCond (0->1)", "AncCond (1->0)", "phyloglm")) %>%
  filter(!is.na(rate)) %>%
  # REMOVE n=100 and q=0.5 for visualization
  filter(n.tips != 100) %>%
  filter(q.rate != 0.5)

# Factorizing to ensure the order in the legend matches your preference
df_plot$method <- factor(df_plot$method, levels = c(
  "AncCond (0->1)", 
  "AncCond (1->0)", 
  "phyloglm"
))

# 3. Labeling for the facets (These will now automatically exclude 100 and 0.5)
n_labs <- setNames(paste0("n = ", unique(df_plot$n.tips)), unique(df_plot$n.tips))
q_labs <- setNames(paste0("q = ", unique(df_plot$q.rate)), unique(df_plot$q.rate))

# 4. Create the Plot
p <- ggplot(df_plot, aes(x = s, y = rate, color = method, group = method)) +
  geom_hline(yintercept = 0.05, linetype = "dotted", color = "black", alpha = 0.6) +
  geom_line(linewidth = 1) + # Increased line thickness slightly for visibility
  geom_point(size = 3) +    # Increased point size
  
  facet_grid(q.rate ~ n.tips, 
             labeller = labeller(q.rate = as_labeller(q_labs), 
                                 n.tips = as_labeller(n_labs))) +
  
  # UPDATED COLORS: Darkened phyloglm to a deeper gold/yellow
  scale_color_manual(values = c(
    "AncCond (0->1)" = "#440154FF", 
    "AncCond (1->0)" = "#21908CFF", 
    "phyloglm"       = "#D4A017"    # Darker Gold/Mustard for better contrast
  )) +
  
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25)) +
  scale_x_continuous(breaks = unique(df_plot$s)) +
  labs(
    title = "Method Performance: Power Comparison",
    x = "Signal Strength (s)",
    y = "Power",
    color = "Method"
  ) +
  theme_bw() + 
  theme(
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    
    # BIGGER LEGEND
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.width = unit(1.5, "cm") # Makes the legend lines longer/easier to see
  )

print(p)
ggsave("anccond_power_filtered.png", p, width = 12, height = 10, dpi = 300)