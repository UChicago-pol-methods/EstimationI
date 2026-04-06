## Coverage simulation: classical vs robust SEs under heteroskedasticity
## Varies: sample size (n) and regressor kurtosis (light vs heavy tails)
## DGP: y_i = beta * x_i + e_i, with sigma_i^2 = x_i^2

library(ggplot2)
library(sandwich)
library(lmtest)

set.seed(42)
B <- 5000
beta_true <- 1

# Three regressor distributions with increasing kurtosis
# t(df) kurtosis = 3 + 6/(df-4) for df > 4
# Normal: kurtosis = 3
# t(5):   kurtosis = 9
# t(4.22): kurtosis ≈ 30  (like CPS wages in Hansen's example)
x_dists <- list(
  "Normal (kurtosis = 3)"            = function(n) rnorm(n),
  "t(5) (kurtosis = 9)"              = function(n) rt(n, df = 5),
  "t(4.2) (kurtosis ~ 30, like wages)" = function(n) rt(n, df = 4.2)
)

sample_sizes <- c(30, 50, 100, 250, 500, 1000)

results <- data.frame()

for (dist_name in names(x_dists)) {
  draw_x <- x_dists[[dist_name]]
  for (n in sample_sizes) {
    cover_classical <- 0
    cover_hc0 <- 0
    cover_hc2 <- 0

    for (b in 1:B) {
      x <- draw_x(n)
      e <- rnorm(n, 0, abs(x))   # sigma_i = |x_i|, so sigma_i^2 = x_i^2
      y <- beta_true * x + e

      fit <- lm(y ~ x - 1)       # no intercept to match the theory

      # Classical CI
      se_classical <- summary(fit)$coefficients[1, 2]
      ci_class <- coef(fit) + c(-1, 1) * qnorm(0.975) * se_classical
      cover_classical <- cover_classical + (ci_class[1] <= beta_true & beta_true <= ci_class[2])

      # HC0
      se_hc0 <- sqrt(vcovHC(fit, type = "HC0")[1, 1])
      ci_hc0 <- coef(fit) + c(-1, 1) * qnorm(0.975) * se_hc0
      cover_hc0 <- cover_hc0 + (ci_hc0[1] <= beta_true & beta_true <= ci_hc0[2])

      # HC2
      se_hc2 <- sqrt(vcovHC(fit, type = "HC2")[1, 1])
      ci_hc2 <- coef(fit) + c(-1, 1) * qnorm(0.975) * se_hc2
      cover_hc2 <- cover_hc2 + (ci_hc2[1] <= beta_true & beta_true <= ci_hc2[2])
    }

    results <- rbind(results, data.frame(
      Distribution = dist_name,
      n = n,
      SE_Type = c("Classical", "HC0 (White)", "HC2"),
      Coverage = c(cover_classical, cover_hc0, cover_hc2) / B
    ))
  }
}

# Preserve facet order
results$Distribution <- factor(results$Distribution, levels = names(x_dists))
results$SE_Type <- factor(results$SE_Type, levels = c("Classical", "HC0 (White)", "HC2"))

p <- ggplot(results, aes(x = n, y = Coverage, color = SE_Type, shape = SE_Type)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black", linewidth = 0.5) +
  facet_wrap(~ Distribution) +
  scale_y_continuous(limits = c(0.45, 1), labels = scales::percent_format()) +
  scale_x_log10(breaks = c(30, 100, 300, 1000)) +
  labs(
    x = "Sample size (n)",
    y = "Coverage of nominal 95% CI",
    color = "SE estimator",
    shape = "SE estimator"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  )

ggsave("coverage_sim.pdf", p, width = 10, height = 4.2)
ggsave("coverage_sim.png", p, width = 10, height = 4.2, dpi = 200)
cat("Done.\n")
