library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Prepare the data (calculating mean and min/max across taxa)
df_fpr_final <- df_asym %>%
  filter(s == 1, q10 == 0.1) %>%
  select(n.tips, q01, anccond.01, anccond.10, phyloglm) %>%
  pivot_longer(cols = c(anccond.01, anccond.10, phyloglm), 
               names_to = "method", 
               values_to = "fpr") %>%
  mutate(
    method_label = case_when(
      method == "anccond.01" ~ "AncCond (0 -> 1)",
      method == "anccond.10" ~ "AncCond (1 -> 0)",
      method == "phyloglm"   ~ "phyloglm"
    ),
    # Ensure 1.0 and 2.0 stay formatted correctly
    symmetry_status = case_when(
      q01 == 0.1 ~ "Symmetrical (0.1/0.1)",
      TRUE ~ paste0("Asym (", sprintf("%.1f", q01), "/0.1)")
    )
  ) %>%
  # Group by method and symmetry to find the taxa-driven range
  group_by(method_label, symmetry_status) %>%
  summarize(
    mean_fpr = mean(fpr, na.rm = TRUE),
    min_fpr  = min(fpr, na.rm = TRUE),
    max_fpr  = max(fpr, na.rm = TRUE),
    .groups = "drop"
  )

# Set factor levels for X-axis
df_fpr_final$symmetry_status <- factor(df_fpr_final$symmetry_status, 
                                       levels = c("Symmetrical (0.1/0.1)", 
                                                  "Asym (0.5/0.1)", 
                                                  "Asym (1.0/0.1)", 
                                                  "Asym (2.0/0.1)"))

# 2. Create the Figure
# We use 'position_dodge' so the lines/bars for different methods don't sit on top of each other
pd <- position_dodge(width = 0.3)

p_fpr_bars <- ggplot(df_fpr_final, aes(x = symmetry_status, y = mean_fpr, 
                                       color = method_label, group = method_label)) +
  # Alpha = 0.05 baseline
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  
  # THE BARS: Error bars showing the range across taxa (50 to 500)
  geom_errorbar(aes(ymin = min_fpr, ymax = max_fpr), 
                width = 0.2, position = pd, linewidth = 0.8) +
  
  # THE LINES: Connecting the means
  geom_line(position = pd, linewidth = 1) +
  
  # THE POINTS: Mean values
  geom_point(position = pd, size = 3) +
  
  # Aesthetics
  scale_color_viridis_d(name = "Method / Direction", option = "viridis") +
  scale_y_continuous(limits = c(0, 0.25), breaks = seq(0, 0.25, 0.05)) +
  
  labs(
    title = "False Positive Rate: Impact of Asymmetry",
    subtitle = "Bars represent the range of FPR across tip sizes (50-500)",
    x = "Transition Rate Symmetry (q01 / q10)",
    y = "False Positive Rate (FPR)"
  ) +
  
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

# 3. Save and Display
print(p_fpr_bars)
ggsave("fpr_adversarial_errorbars.png", p_fpr_bars, width = 10, height = 7, dpi = 300)