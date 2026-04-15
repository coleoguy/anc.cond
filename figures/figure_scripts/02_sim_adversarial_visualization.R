# ============================================================================
# Plotting Script: AncCond Adversarial Results Visualization
# ============================================================================

library(ggplot2)
library(tidyr)
library(dplyr)

# 1. Load the results
# Ensure the file path matches your working directory
res_df <- read.csv("02_sim_OU_adversarial_results.csv")

# 2. Reshape data for plotting (Wide to Long)
# This moves 'FP.standard' and 'FP.pruned' into a single 'Method' column
plot_data <- res_df %>%
  pivot_longer(
    cols = c(FP.standard, FP.pruned),
    names_to = "Method",
    values_to = "FP_Rate"
  ) %>%
  mutate(
    Method = recode(Method, "FP.standard" = "Standard", "FP.pruned" = "Pruned"),
    # Create a cleaner label for the facets
    alpha_label = paste0("Alpha (pull): ", alpha)
  )

# 3. Create the Visualization
p <- ggplot(plot_data, aes(x = delta.theta, y = FP_Rate, color = Method, group = Method)) +
  # Add a reference line for the nominal significance level (0.05)
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", alpha = 0.6) +
  # Add points and lines
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  # Facet by alpha to see how pull strength changes the behavior
  facet_wrap(~ alpha_label) +
  # Aesthetics and Labels
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_color_manual(values = c("Standard" = "#0072B2", "Pruned" = "#D55E00")) +
  labs(
    title = "AncCond False Positive Rates: OU Adversarial Scenario",
    subtitle = "Directionality Test: Discrete trait drives Continuous evolution (Y -> X)",
    x = "Delta Theta (Difference between optima)",
    y = "False Positive Rate (p < 0.05)",
    caption = "Red dashed line indicates the expected 5% error rate."
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "gray95", color = NA),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# 4. Save the plot
print(p)
ggsave("anccond_adversarial_plot.png", plot = p, width = 10, height = 5, dpi = 300)