library(ggplot2)
par(mfrow = c(1,1), mar = c(5,4,4,2) + .1)

## --- Unidirectional ---
load('../../results/scaling_unidirectional_results.RData')
if (exists("results") && !exists("uni_scaling_results")) {
  uni_scaling_results <- results
}

# Parse the comma-separated 'value' column: take the 1st number per row (unidirectional)
uni_mat <- do.call(
  cbind,
  lapply(uni_scaling_results, function(df) {
    vals <- strsplit(as.character(df$value), ",", fixed = TRUE)
    p <- sapply(vals, function(v) suppressWarnings(as.numeric(trimws(v[1]))))
    matrix(p, ncol = 1, dimnames = list(NULL, "pval"))
  })
)

# Proportion <= 0.05 per scale
probs <- colMeans(uni_mat <= 0.05, na.rm = TRUE)


## --- Bidirectional ---
load('../../results/scaling_bidirectional_results.RData')
if (exists("results") && !exists("bi_scaling_results")) {
  bi_scaling_results <- results
}

# Parse the comma-separated 'value' column: take first two numbers per row as (p12, p21)
bi_mat <- do.call(
  cbind,
  lapply(bi_scaling_results, function(df) {
    vals <- strsplit(as.character(df$value), ",", fixed = TRUE)
    m <- t(sapply(vals, function(v) {
      a <- suppressWarnings(as.numeric(trimws(v)))
      if (length(a) == 0) a <- c(NA_real_, NA_real_)
      if (length(a) == 1) a <- c(a[1], NA_real_)
      if (length(a) >= 2) a <- a[1:2]
      a
    }))
    colnames(m) <- c("p12","p21")
    m
  })
)

# Compute proportions <= 0.025 per scale (two cols per scale: p12, p21)
n_scales <- ncol(bi_mat) / 2
probs2 <- sapply(seq_len(n_scales), function(i) {
  cols <- (2*i - 1):(2*i)
  vals <- bi_mat[, cols, drop = FALSE]
  round(sum(vals <= 0.025, na.rm = TRUE) / sum(!is.na(vals)), 2)
})


## --- Build results data frame + plot (your original style) ---
n_scales <- length(probs)  # number of unidirectional scaling factors

results <- as.data.frame(matrix(NA_real_, 2 * n_scales, 3))
colnames(results) <- c("value", "scaling.factor", "type")

results[, 1] <- c(probs2, probs)
results[, 2] <- rep(seq_len(n_scales), 2)

x <- factor(c("bidirectional false positive",
              rep("bidirectional power", n_scales - 1),
              "unidirectional false positive",
              rep("unidirectional power", n_scales - 1)))
x <- factor(x, levels(x)[c(2,4,1,3)])
results[, 3] <- x

ggplot(results, aes(y = value, x = scaling.factor, group = type, colour = type, shape = type)) +
  geom_hline(yintercept = .05, linetype = "dashed") +
  geom_point(stat = "identity", size = 4, alpha = .8) +
  geom_line() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text = element_text(family = "sans", face = "plain", color = "#000000", size = 15),
        legend.text = element_text(size = rel(0.6))) +
  scale_size(range = c(1, 4)) +
  xlab("scaling factor") +
  ylab("percent significant") +
  scale_colour_manual(labels = c("bidirectional power",
                                 "unidirectional power",
                                 "bidirectional false positive",
                                 "unidirectional false positive"),
                      values = c("#377eb8","#e41a1c","#0d0101","#0d0101")) +
  scale_shape_manual(labels = c("bidirectional power",
                                "unidirectional power",
                                "bidirectional false positive",
                                "unidirectional false positive"),
                     values = c(16,16,0,4))
