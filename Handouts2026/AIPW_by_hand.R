# AIPW by hand on Project STAR (Tennessee class-size experiment)
# Companion to Handouts2026/AIPW_review_handout.pdf
#
# Project STAR randomized Tennessee students (within schools) to small or
# regular-sized kindergarten classes. We treat the small-vs-regular
# comparison as a clean RCT (dropping the "regular + aide" arm).
#   D = 1 for small class, 0 for regular class
#   Y = kindergarten reading score (continuous)
# Treatment effect of interest: tau = E[Y(1) - Y(0)].
#
# The script: (1) difference in means, (2) outcome regressions, (3) the
# AIPW estimator computed by hand, (4) the influence-function SE.
#
# Two versions of the outcome model are reported:
#   - in-sample: mu_d fitted and predicted on the same data (transparent
#     but uses each unit's outcome in its own mu_d prediction);
#   - cross-fit (K = 5 folds): mu_d for each unit is built from data
#     that excludes it, satisfying the sample-splitting condition
#     behind the o_p(n^{-1/2}) remainder bound.

library(AER)
data("STAR")

covars <- c("gender", "ethnicity", "birth", "lunchk", "schoolk",
            "experiencek", "tethnicityk")
dat <- subset(STAR,
              !is.na(stark) & !is.na(readk) &
              stark %in% c("small", "regular") &
              ethnicity %in% c("cauc", "afam"))    # 99% of the sample
dat <- dat[complete.cases(dat[, covars]), ]
dat$ethnicity <- droplevels(dat$ethnicity)
dat$D <- as.numeric(dat$stark == "small")
dat$Y <- dat$readk

Y    <- dat$Y
D    <- dat$D
X_df <- dat[, covars]
n    <- length(Y)
p    <- mean(D)                         # within-school randomization

## 1. Difference in means and its SE
tau_dm <- mean(Y[D == 1]) - mean(Y[D == 0])
se_dm  <- sqrt(var(Y[D == 1]) / sum(D == 1) +
               var(Y[D == 0]) / sum(D == 0))

## 2. AIPW with in-sample outcome models
##    Correct specification of mu_d is not required for consistency under
##    randomization, but a better fit lowers the AIPW asymptotic variance.
dat_fit <- cbind(Y = Y, X_df)
fit1    <- lm(Y ~ ., data = dat_fit, subset = D == 1)
fit0    <- lm(Y ~ ., data = dat_fit, subset = D == 0)
mu1     <- predict(fit1, newdata = X_df)
mu0     <- predict(fit0, newdata = X_df)

## helper: AIPW point estimate and IF-based SE given mu1, mu0 vectors
aipw_summary <- function(mu1, mu0) {
  psi <- mu1 - mu0 +
         D       * (Y - mu1) / p -
         (1 - D) * (Y - mu0) / (1 - p)
  tau <- mean(psi)
  list(tau = tau,
       se  = sqrt(mean((psi - tau)^2) / n))
}

aipw_in <- aipw_summary(mu1, mu0)

## 3. AIPW with cross-fit outcome models (K = 5 folds)
##    For each unit i, mu_d(X_i) is constructed from data in the other
##    K - 1 folds, so mu_hat is independent of unit i -- this is what
##    Neyman orthogonality plus L^2-consistency of mu_d needs to make
##    the feasible-AIPW remainder o_p(n^{-1/2}).
set.seed(20260511)
K <- 5
folds <- sample(rep(1:K, length.out = n))

mu1_cf <- rep(NA_real_, n)
mu0_cf <- rep(NA_real_, n)
for (k in 1:K) {
  out <- folds == k
  tr  <- !out
  fit1_k <- lm(Y ~ ., data = dat_fit, subset = tr & D == 1)
  fit0_k <- lm(Y ~ ., data = dat_fit, subset = tr & D == 0)
  mu1_cf[out] <- predict(fit1_k, newdata = X_df[out, ])
  mu0_cf[out] <- predict(fit0_k, newdata = X_df[out, ])
}

aipw_cf <- aipw_summary(mu1_cf, mu0_cf)

## 4. 95% CIs and Wald test of H0: tau = 0
z       <- qnorm(0.975)
ci_dm   <- tau_dm        + c(-1, 1) * z * se_dm
ci_in   <- aipw_in$tau   + c(-1, 1) * z * aipw_in$se
ci_cf   <- aipw_cf$tau   + c(-1, 1) * z * aipw_cf$se

## 5. Report
results <- data.frame(
  estimator = c("Difference in means",
                "AIPW (in-sample mu_d)",
                "AIPW (5-fold cross-fit mu_d)"),
  estimate  = c(tau_dm, aipw_in$tau, aipw_cf$tau),
  se        = c(se_dm,  aipw_in$se,  aipw_cf$se),
  ci_lo     = c(ci_dm[1], ci_in[1], ci_cf[1]),
  ci_hi     = c(ci_dm[2], ci_in[2], ci_cf[2])
)
print(results, row.names = FALSE, digits = 4)

R2_mu1 <- summary(fit1)$r.squared
R2_mu0 <- summary(fit0)$r.squared
cat(sprintf("\nIn-sample outcome R^2:  arm 1 = %.3f,  arm 0 = %.3f\n",
            R2_mu1, R2_mu0))
cat(sprintf("Wald z (cross-fit AIPW, H0: tau = 0):  %.2f\n",
            aipw_cf$tau / aipw_cf$se))
cat(sprintf("SE ratio AIPW/DM:        in-sample = %.3f,  cross-fit = %.3f\n",
            aipw_in$se / se_dm,
            aipw_cf$se / se_dm))
