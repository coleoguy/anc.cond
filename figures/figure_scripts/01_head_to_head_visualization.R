library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Load results
df <- read.csv("01_sim_power_results.csv")

# 2. Reshape for Head-to-Head
df_h2h <- df %>%
  filter(s > 1) %>% 
  select(n.tips, q.rate, s, method, rate) %>%
  pivot_wider(names_from = method, values_from = rate) %>%
  mutate(
    AncCond_power = (anccond.01 + anccond.10) / 2,
    phyloglm_power = phyloglm
  ) %>%
  filter(!is.na(AncCond_power) & !is.na(phyloglm_power))

# 3. Create the Plot
p <- ggplot(df_h2h, aes(x = phyloglm_power, y = AncCond_power)) +
  # Darker diagonal reference line
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.6) +
  
  # Point Layer
  geom_point(aes(size = n.tips, 
                 color = as.factor(s), 
                 shape = as.factor(q.rate)), 
             alpha = 0.75) + # Slightly higher alpha for visibility
  
  # Numerical Corner Annotations
  annotate("text", x = 0, y = 1, label = "AncCond Better", 
           fontface = "italic", color = "gray40", 
           hjust = -0.1, vjust = 0.4) + 
  
  annotate("text", x = 1, y = 0, label = "Phyloglm Better", 
           fontface = "italic", color = "gray40", 
           hjust = 1.1, vjust = 0.3) + 
  
  # MATCHING THE THEME: Using the same palette logic (Purple -> Teal -> Dark Gold)
  # Adjust the number of colors in 'values' if you have more/fewer signal strengths (s)
  scale_color_manual(
    name = "Signal Strength (s)",
    values = c(
      "2" = "#440154FF", # Deep Purple
      "3" = "#31688EFF", # Steel Blue
      "5" = "#40917A", # Teal
      "10" = "#D4A017"    # Dark Gold (The color we just set for phyloglm)
    ),
    guide = guide_legend(override.aes = list(size = 4))
  ) +
  
  scale_shape_manual(values = c(18, 17, 15, 16), 
                     name = "Trans. Rate (q)", 
                     guide = guide_legend(override.aes = list(size = 4))) +
  
  scale_size_continuous(range = c(4, 11), 
                        breaks = c(50, 100, 200, 500), 
                        name = "Tree Size (n)") +
  
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  
  labs(
    title = "Head-to-Head Power Comparison",
    x = "Phyloglm Power",
    y = "AncCond Power"
  ) +
  
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.box = "vertical",
    legend.position = "right",
    
    # ENHANCED LEGEND
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    legend.spacing.y = unit(0.2, "cm")
  )


# 4. Display and Save
print(p)
ggsave("h2h_power_matched.png", p, width = 10, height = 8, dpi = 300)