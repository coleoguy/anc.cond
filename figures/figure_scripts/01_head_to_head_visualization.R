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
             alpha = 0.7) +
  
  # UPDATED: Numerical Corner Annotations
  # Top-left: Set to 0 and 1. Adjust 'hjust' and 'vjust' to nudge the text.
  annotate("text", x = 0, y = 1, label = "AncCond better", 
           fontface = "italic", color = "gray40", 
           hjust = -0.1,  # Negative nudges it right
           vjust = 0.4) + # Above 1 nudges it down
  
  # Bottom-right: Set to 1 and 0. 
  annotate("text", x = 1, y = 0, label = "phyloglm better", 
           fontface = "italic", color = "gray40", 
           hjust = 1.1,   # Above 1 nudges it left
           vjust = 0.3) + # Negative nudges it up
  
  # Aesthetics
  scale_color_viridis_d(name = "Signal Strength (s)") +
  scale_shape_manual(values = c(18, 17, 15, 16), name = "Trans. Rate (q)") +
  scale_size_continuous(range = c(3, 10), 
                        breaks = c(100, 200, 300, 400, 500), 
                        name = "Tree Size (n)") +
  
  # Fixed 1:1 Aspect Ratio
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  
  labs(
    title = "Head-to-Head Power Comparison",
    x = "phyloglm power",
    y = "AncCond power"
  ) +
  
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.box = "vertical",
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

# 4. Display and Save
print(p)
ggsave("h2h_power.png", p, width = 10, height = 7, dpi = 300)