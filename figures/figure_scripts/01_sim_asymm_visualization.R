library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Load results
df_sym  <- read.csv("01_sim_power_results.csv")
df_asym <- read.csv("01_sim_asymm_power_results.csv")

# 2. Process Symmetrical AncCond Data
df_sym_plot <- df_sym %>%
  filter(method %in% c("anccond.01", "anccond.10")) %>%
  mutate(
    method_label = case_when(
      method == "anccond.01" ~ "AncCond (0->1) Sym",
      method == "anccond.10" ~ "AncCond (1->0) Sym"
    )
  ) %>%
  select(n.tips, q.rate, s, method = method_label, rate)

# 3. Process Asymmetrical AncCond Data (Rare loss direction only)
df_asym_plot <- df_asym %>%
  filter(q10 == 0.1) %>% 
  mutate(
    q.rate = q01, 
    method = "AncCond (1->0) Asym (q10=0.1)",
    rate = anccond.10 
  ) %>%
  select(n.tips, q.rate, s, method, rate)

# 4. Combine and define factor levels for consistent legend ordering
df_final <- rbind(df_sym_plot, df_asym_plot) %>%
  filter(!is.na(rate))

df_final$method <- factor(df_final$method, levels = c(
  "AncCond (0->1) Sym",
  "AncCond (1->0) Sym",
  "AncCond (1->0) Asym (q10=0.1)"
))

# 5. Facet labels
n_labs <- setNames(paste0("n = ", unique(df_final$n.tips)), unique(df_final$n.tips))
q_labs <- setNames(paste0("q(forward) = ", unique(df_final$q.rate)), unique(df_final$q.rate))

# 6. Create the Plot
p <- ggplot(df_final, aes(x = s, y = rate, color = method, group = method)) +
  # Alpha = 0.05 baseline
  geom_hline(yintercept = 0.05, linetype = "dotted", color = "black", alpha = 0.6) +
  
  # Lines and Points
  # 'linetype = "solid"' is default, but we'll ensure all are solid here
  geom_line(linewidth = 0.9, linetype = "solid") +
  geom_point(size = 2.5) +
  
  # THE GRID: Rows = Transition Rate (q), Columns = Taxa Size (n)
  facet_grid(q.rate ~ n.tips, 
             labeller = labeller(q.rate = as_labeller(q_labs), 
                                 n.tips = as_labeller(n_labs))) +
  
  # CUSTOM COLOR MAPPING
  # Using ggplot2 default hex codes for red/green and a standard blue for the asym line
  scale_color_manual(values = c(
    "AncCond (0->1) Sym"            = "#F8766D", # Default ggplot Red
    "AncCond (1->0) Sym"            = "#7CAE00", # Default ggplot Green
    "AncCond (1->0) Asym (q10=0.1)" = "#619CFF"  # Standard ggplot Blue
  )) +
  
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25)) +
  scale_x_continuous(breaks = unique(df_final$s)) +
  labs(
    title = "AncCond Power: Symmetrical vs. Asymmetrical Evolution",
    subtitle = "Red/Green = Symmetrical (Mk) | Blue = Asymmetrical (Source-Sink: q10=0.1)",
    x = "Signal Strength (s)",
    y = "Rate (n.sig / n.valid)",
    color = "Model & Direction"
  ) +
  theme_bw() + 
  theme(
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

# 7. Finalize
print(p)
ggsave("anccond_asym_blue_comparison.png", p, width = 12, height = 10, dpi = 300)