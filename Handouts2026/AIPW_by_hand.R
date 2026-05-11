# AIPW by hand on the Bertrand-Mullainathan audit study
# Companion to Handouts2026/AIPW_review_handout.pdf
#
# Bertrand & Mullainathan (2004) sent ~4870 resumes to job ads in Boston
# and Chicago with names randomly assigned to signal applicant race.
#   D = 1 if African-American-sounding name, 0 otherwise (randomized)
#   Y = 1 if the employer called back
# Treatment effect of interest: tau = E[Y(1) - Y(0)].
#
# The script: (1) difference in means, (2) outcome regressions, (3) the
# AIPW estimator computed by hand, (4) the influence-function SE.
#
# Caveat: for honest inference, the outcome models should be fit on
# sample-split folds. A single in-sample fit is used here for clarity.

library(AER)
data("ResumeNames")

Y    <- as.numeric(ResumeNames$call == "yes")
D    <- as.numeric(ResumeNames$ethnicity == "afam")
X_df <- ResumeNames[, c("gender", "quality", "city", "jobs", "experience",
                        "honors", "computer", "college")]
n    <- length(Y)
p    <- mean(D)                         # ~ 0.5 by construction

## 1. Difference in means and its SE
tau_dm <- mean(Y[D == 1]) - mean(Y[D == 0])
se_dm  <- sqrt(var(Y[D == 1]) / sum(D == 1) +
               var(Y[D == 0]) / sum(D == 0))

## 2. Outcome models, separately by arm
##    Correct specification of mu_d is not required for consistency under
##    randomization, but a better fit lowers the AIPW asymptotic variance.
dat  <- cbind(Y = Y, X_df)
fit1 <- lm(Y ~ ., data = dat, subset = D == 1)
fit0 <- lm(Y ~ ., data = dat, subset = D == 0)
mu1  <- predict(fit1, newdata = X_df)
mu0  <- predict(fit0, newdata = X_df)

## 3. AIPW estimator (the slide-3 formula, by hand)
psi <- mu1 - mu0 +
       D       * (Y - mu1) / p -
       (1 - D) * (Y - mu0) / (1 - p)
tau_aipw <- mean(psi)

## 4. Variance from the empirical influence function
psi_hat <- psi - tau_aipw               # centered influence-function values
V_hat   <- mean(psi_hat^2)
se_aipw <- sqrt(V_hat / n)

## 5. 95% CIs and Wald test of H0: tau = 0
z       <- qnorm(0.975)
ci_dm   <- tau_dm   + c(-1, 1) * z * se_dm
ci_aipw <- tau_aipw + c(-1, 1) * z * se_aipw
wald    <- tau_aipw / se_aipw

## 6. Report
results <- data.frame(
  estimator = c("Difference in means", "AIPW (by hand)"),
  estimate  = c(tau_dm, tau_aipw),
  se        = c(se_dm, se_aipw),
  ci_lo     = c(ci_dm[1], ci_aipw[1]),
  ci_hi     = c(ci_dm[2], ci_aipw[2])
)
print(results, row.names = FALSE, digits = 4)

cat(sprintf("\nWald statistic, AIPW, H0: tau = 0:  z = %.2f\n", wald))
cat(sprintf("SE ratio AIPW/DM:                    %.3f\n", se_aipw / se_dm))
cat(sprintf("Implied variance reduction:          %.0f%%\n",
            100 * (1 - (se_aipw / se_dm)^2)))
