# %% [markdown]
# # Unit 4. Additional single-equation topics
#
# Econometrics II (Henry Redondo, 2026-27). R companion to **Unit 4** of the slides: generated regressors (A1, A2), testing for endogeneity (B1, B2), overidentifying restrictions (C1), functional form and heteroskedasticity (D1, D2).
#
# Each section has the idea with the slide page (page of the PDF), then the code, then how to read the output. Numbers quoted as "slides" are the ones printed in the slides.

# %%
library(tidyverse)
library(broom)
library(estimatr)        # lm_robust(), iv_robust()
library(ivreg)           # ivreg()
library(lmtest)          # resettest(), bptest()
library(modelsummary)

options(digits = 4, readr.show_col_types = FALSE)
theme_set(theme_minimal())

# Henry's data folder. To read local copies, set the environment variable ECMTII_DATA to their folder.
data_dir <- Sys.getenv("ECMTII_DATA", "https://henryredondo.github.io/assets/files/teaching/econometrics-ii")

# %% [markdown]
# ## The data: 3,000 workers (p. 4)
#
# Simulated labour-force survey: log wage, schooling, experience, instruments for schooling (parents' schooling `meduc` and `feduc`, distance to the nearest university `dist`), the number of siblings `sibs`, and the traits of each worker's firm. Ability is never observed. The true return to schooling is 0.08. A separate file has 100 firms where the employer wage premium `q` is measured.
#
# > **Data pending.** On 25 November 2026 the CSV files were not on the course website. The slides give `workers.csv`. They give neither the name of the 100-firm file nor the names of the firm traits: the cell below assumes `firms.csv` and the columns `size`, `exporter`, `capital` (in both files) and `q`. Change `firm_file` and `traits` when the real names are known.

# %%
firm_file <- "firms.csv"                       # assumption
traits    <- c("size", "exporter", "capital")  # assumption: firm size, exporter, capital intensity

d     <- read_csv(file.path(data_dir, "workers.csv")) |> mutate(expersq = exper^2)
firms <- read_csv(file.path(data_dir, firm_file))

glimpse(d)

# %% [markdown]
# ## A1. A premium estimated on 100 firms (p. 12-15)
#
# Wages depend on the employer premium `q`, which is observed only in the firm file. Stage 1 learns there how `q` depends on firm traits; stage 2 uses the prediction $\hat q$ for the workers' firms as a regressor. OLS with a **generated regressor** is consistent, and its reported standard error ignores that stage 1 was estimated.

# %%
# Stage 1: the 100 firms
stage1 <- lm(reformulate(traits, response = "q"), data = firms)
tidy(stage1)

# Stage 2: the 3,000 workers, with the predicted premium
d <- d |> mutate(qhat = predict(stage1, newdata = d))
stage2 <- lm(lwage ~ exper + expersq + qhat, data = d)
tidy(stage2) |> filter(term == "qhat")

# %% [markdown]
# ### Two-stage bootstrap
#
# Resample **both** samples (the 100 firms and the 3,000 workers) and repeat both stages, so the uncertainty of stage 1 is included.

# %%
boot_once <- function() {
  f_b <- slice_sample(firms, prop = 1, replace = TRUE)
  w_b <- slice_sample(d, prop = 1, replace = TRUE)
  s1  <- lm(reformulate(traits, response = "q"), data = f_b)
  w_b$qhat <- predict(s1, newdata = w_b)
  coef(lm(lwage ~ exper + expersq + qhat, data = w_b))[["qhat"]]
}

set.seed(2026)
gamma_boot <- map_dbl(1:500, \(b) boot_once())

tibble(method   = c("OLS with qhat (reported)", "Two-stage bootstrap"),
       gamma    = coef(stage2)[["qhat"]],
       std.error = c(tidy(stage2)$std.error[4], sd(gamma_boot))) |>
  mutate(conf.low = gamma - 1.96 * std.error, conf.high = gamma + 1.96 * std.error)

# %% [markdown]
# **Reading.** Slides: stage 1 gives 0.051, 0.134 and 0.056 (truth 0.060, 0.120, 0.050; s.e. 0.014, 0.051, 0.020). Stage 2 gives a pass-through of 0.517 with a reported s.e. of 0.059; OLS with the true `q` gives 0.519 (0.057); the two-stage bootstrap s.e. is 0.114. The three estimates sit near the truth (0.50): consistency holds. The reported standard error is about half the bootstrap one, because one set of stage 1 estimates built $\hat q$ for all 3,000 workers and more workers cannot average out an error they all share.
#
# The row "OLS with the true `q`" needs the workers' true premium, which is not in the files described in the slides.

# %% [markdown]
# ## A2. What many samples say (p. 16)
#
# > The slides repeat the exercise over 2,000 samples with Henry's data-generating process, which is not published. The process below is a simplified one written for this notebook (500 samples). Its numbers will be of the same order as the slides, not equal.

# %%
delta <- c(0.06, 0.12, 0.05)   # true stage 1 coefficients (p. 12)

draw_traits <- function(n) cbind(size = rnorm(n, 4, 1.2), exporter = rbinom(n, 1, 0.3), capital = rnorm(n, 2, 1))

one_sample <- function(gamma, n_firms = 100, n_workers = 3000) {
  Wf <- draw_traits(n_firms)                                                # firm file
  qf <- drop(Wf %*% delta) + rnorm(n_firms, sd = 0.17)                      # premium, measured with error
  Ww <- draw_traits(n_workers)                                              # workers' firms
  exper <- pmax(0, rnorm(n_workers, 16, 8))
  lwage <- 1 + 0.035 * exper - 0.0006 * exper^2 + gamma * (drop(Ww %*% delta) + rnorm(n_workers, sd = 0.08)) +
    rnorm(n_workers, sd = 0.3)
  qhat <- drop(cbind(1, Ww) %*% coef(lm.fit(cbind(1, Wf), qf)))
  est  <- coef(summary(lm(lwage ~ exper + I(exper^2) + qhat)))["qhat", ]
  c(gamma_hat = est[["Estimate"]], se = est[["Std. Error"]])
}

set.seed(2026)
mc <- bind_rows(
  map(1:500, \(i) one_sample(gamma = 0.5)) |> bind_rows() |> mutate(gamma = 0.5),
  map(1:500, \(i) one_sample(gamma = 0))   |> bind_rows() |> mutate(gamma = 0)
)

mc |> summarise(mean_estimate  = mean(gamma_hat),
                true_spread    = sd(gamma_hat),
                mean_reported_se = mean(se),
                rejection_rate = mean(abs((gamma_hat - gamma) / se) > 1.96),   # naive 5% test of the true null
                .by = gamma)

# %% [markdown]
# **Reading.** Slides, with a true pass-through of 0.50: mean estimate 0.486, true spread 0.105, average reported s.e. 0.058 (45% too small), and the naive 5% test of the true null rejects 28.8% of the time (9.0% with the bootstrap). With a true value of 0 the naive test rejects 5.0%: the estimation error of stage 1 enters multiplied by the pass-through, so it vanishes under that null.

# %% [markdown]
# ## B1. Is schooling endogenous? (p. 23-25)
#
# The regression-based test in three steps: (1) reduced form of `educ` on all the exogenous variables, keep the residuals $\hat v$; (2) add $\hat v$ to the wage equation (the **control-function** regression); (3) t test of its coefficient. The coefficient on `educ` in step 2 is numerically the 2SLS estimate.

# %%
fs <- lm(educ ~ exper + expersq + meduc + feduc + dist, data = d)
d  <- d |> mutate(vhat = resid(fs))

cf  <- lm(lwage ~ educ + exper + expersq + vhat, data = d)
ols <- lm(lwage ~ educ + exper + expersq, data = d)
iv  <- ivreg(lwage ~ educ + exper + expersq | exper + expersq + meduc + feduc + dist, data = d)

modelsummary(list("OLS" = ols, "2SLS" = iv, "Control function" = cf),
             coef_map = c("educ", "vhat"), fmt = 4, gof_map = "nobs", output = "markdown")

# t test of vhat: classical and robust
t_vhat <- tidy(cf) |> filter(term == "vhat") |> pull(statistic)
c(classical = t_vhat,
  robust    = tidy(lm_robust(lwage ~ educ + exper + expersq + vhat, data = d, se_type = "HC1"))$statistic[5])

# %%
# The built-in versions after 2SLS: weak instruments, Wu-Hausman, Sargan
diag_iv <- summary(iv, diagnostics = TRUE)$diagnostics
diag_iv

# The control-function coefficient on educ is 2SLS, and t^2 is the Wu-Hausman F
stopifnot(all.equal(coef(cf)[["educ"]], coef(iv)[["educ"]]),
          all.equal(t_vhat^2, diag_iv["Wu-Hausman", "statistic"]))

# %% [markdown]
# The same in one line with `iv_robust`, whose diagnostics use the robust variance.

# %%
iv_r <- iv_robust(lwage ~ educ + exper + expersq | exper + expersq + meduc + feduc + dist,
                  data = d, se_type = "HC1", diagnostics = TRUE)
tidy(iv_r) |> filter(term == "educ")
iv_r$diagnostic_first_stage_fstatistic
iv_r$diagnostic_endogeneity_test

# %% [markdown]
# **Reading.** Slides: OLS 0.1103 (s.e. 0.0028), 2SLS 0.0826 (0.0057), coefficient on $\hat v$ 0.0370 (0.0065), t = 5.72 (robust 5.86), so exogeneity of schooling is rejected; t² = 32.75 is the Wu-Hausman F of `ivreg`. First stage F = 334. The truth is 0.08: OLS overstates it by 0.030 because ability is in the error and raises schooling. The price of IV is a standard error twice the OLS one.

# %% [markdown]
# ## B2. Large difference, no rejection? (p. 26)
#
# The power of the endogeneity test comes from the strength of the instruments. We simulate surveys with weaker and weaker instruments and count how often the test rejects.
#
# > Simplified data-generating process written for this notebook (one instrument, 300 samples per strength); the slides use 600 samples of Henry's process.

# %%
one_survey <- function(strength, n = 3000) {
  ability <- rnorm(n); z <- rnorm(n)
  educ  <- 12 + strength * z + ability + rnorm(n)
  lwage <- 1 + 0.08 * educ + 0.06 * ability + rnorm(n, sd = 0.15)   # OLS bias: 0.03
  vhat  <- resid(lm(educ ~ z))
  cf    <- coef(summary(lm(lwage ~ educ + vhat)))
  c(F = summary(lm(educ ~ z))$fstatistic[["value"]],
    gap = abs(cf["educ", "Estimate"] - coef(lm(lwage ~ educ))[["educ"]]),   # |IV - OLS|
    rejects = abs(cf["vhat", "t value"]) > 1.96)
}

set.seed(2026)
map(c(0.02, 0.05, 0.1, 0.5), \(s) map(1:300, \(i) one_survey(s)) |> bind_rows() |> mutate(strength = s)) |>
  list_rbind() |>
  summarise(mean_F = mean(F), median_gap_iv_ols = median(gap), share_rejects = mean(rejects), .by = strength)

# %% [markdown]
# **Reading.** Slides: with a mean F of 1.2, 4.7, 15.5 and 361 the test rejects 4.8%, 12.5%, 32.5% and 100% of the time, while the average distance between IV and OLS is 0.092, 0.051, 0.041 and 0.030. With weak instruments IV and OLS can be far apart and the test still fails to reject, because the IV standard error is too large. (The table above reports the median distance, which is less sensitive to the extreme estimates that weak instruments produce.)

# %% [markdown]
# ## C1. Sargan, and what it cannot see (p. 29-33)
#
# With more instruments than endogenous regressors the exogeneity of the instruments can be partly tested. The **Sargan** statistic is N·R² of the regression of the 2SLS residuals on all the exogenous variables, χ² with (number of instruments − number of endogenous regressors) degrees of freedom. It assumes homoskedasticity; the robust version is Hansen's J.

# %%
# By hand: N x R2 of the 2SLS residuals on z
aux    <- lm(resid(iv) ~ exper + expersq + meduc + feduc + dist, data = d)
sargan <- nobs(aux) * summary(aux)$r.squared
c(R2 = summary(aux)$r.squared, sargan = sargan, p = pchisq(sargan, df = 2, lower.tail = FALSE))

stopifnot(all.equal(sargan, diag_iv["Sargan", "statistic"]))

# Robust version (Hansen J), from iv_robust
iv_r$diagnostic_overid_test

# %%
# The same test for other instrument lists
overid <- function(instruments) {
  f <- as.formula(str_glue("lwage ~ educ + exper + expersq | exper + expersq + {instruments}"))
  m <- ivreg(f, data = d)
  s <- summary(m, diagnostics = TRUE)$diagnostics["Sargan", ]
  tibble(instruments, beta_educ = coef(m)[["educ"]], sargan = s[["statistic"]], df = s[["df1"]], p = s[["p-value"]])
}

map(c("meduc + feduc + dist", "meduc + feduc + dist + sibs", "meduc + feduc"), overid) |> list_rbind()

# %% [markdown]
# **Reading.** Slides: with the valid instruments R² = 0.00121, Sargan = 3,000 × 0.00121 = 3.63 (p = 0.163) and robust J = 4.02 (p = 0.134): no rejection. Adding `sibs`, which is correlated with ability, gives 88.5 with 3 d.f. and 2SLS jumps to 0.106.
#
# What Sargan cannot see: the test compares instruments with each other. If `meduc` and `feduc` were invalid in the same way (family background raising wages directly), the test with those two alone would not reject (the slides simulate that world: 2SLS 0.136, p = 0.328); only adding a valid instrument of a different kind (`dist`) reveals it.

# %% [markdown]
# ## D1. RESET on the wage equation (p. 38)
#
# RESET adds powers of the fitted values and tests them jointly. It assumes exogenous regressors, so we apply it to the control-function equation, with and without the square of experience.

# %%
cf_linear <- lm(lwage ~ educ + exper + vhat, data = d)   # forgot exper^2

resettest(cf_linear, power = 2:3, type = "fitted")
resettest(cf, power = 2:3, type = "fitted")

# %%
# The three versions by hand for the linear equation: F, Wald and LM (N x R2)
aug   <- update(cf_linear, . ~ . + I(fitted(cf_linear)^2) + I(fitted(cf_linear)^3))
F_by_hand <- anova(cf_linear, aug)$F[2]
lm_aux <- lm(resid(cf_linear) ~ educ + exper + vhat + I(fitted(cf_linear)^2) + I(fitted(cf_linear)^3), data = d)

c(F = F_by_hand, Wald = 2 * F_by_hand, LM = nobs(lm_aux) * summary(lm_aux)$r.squared)

stopifnot(all.equal(F_by_hand, resettest(cf_linear, power = 2:3, type = "fitted")$statistic[["RESET"]]))

# %% [markdown]
# **Reading.** Slides: linear in experience, F = 14.25, Wald 28.50, N·R² 28.28, p < 0.001; with the square, F = 1.61, Wald 3.22, N·R² 3.22, p = 0.200. RESET rejects the equation that is too simple and does not reject the one with the quadratic. The classical F is not robust to heteroskedasticity, which is the next point.

# %% [markdown]
# ## D2. Is the variance constant? (p. 37, 39)
#
# Breusch-Pagan: N·R² of the squared residuals on the regressors. White (special case): on the fitted values and their squares. After 2SLS: squared 2SLS residuals on the exogenous variables.

# %%
bptest(cf)                                        # Breusch-Pagan, 4 regressors
bptest(cf, ~ fitted(cf) + I(fitted(cf)^2))        # White, special case

# Breusch-Pagan by hand
bp_aux <- lm(resid(cf)^2 ~ educ + exper + expersq + vhat, data = d)
stopifnot(all.equal(nobs(bp_aux) * summary(bp_aux)$r.squared, bptest(cf)$statistic[["BP"]]))

# After 2SLS: squared residuals on the exogenous z (5 regressors)
nr2_iv <- nobs(iv) * summary(lm(resid(iv)^2 ~ exper + expersq + meduc + feduc + dist, data = d))$r.squared
c(NR2_after_2SLS = nr2_iv, p = pchisq(nr2_iv, df = 5, lower.tail = FALSE))

# %%
# What it costs: the standard error of educ
modelsummary(list("CF, classical" = cf, "CF, robust" = cf, "2SLS, robust" = iv),
             vcov = list("classical", "HC1", "HC1"),
             coef_map = "educ", fmt = 4, gof_map = "nobs", output = "markdown")

# %% [markdown]
# **Reading.** Slides: Breusch-Pagan 290.9 and White 288.5 (2 d.f.), both p < 0.001; after 2SLS, N·R² = 68.1. In the simulation the s.d. of the wage noise rises with schooling. The coefficient on `educ` stays at 0.0826; in the augmented regression its standard error moves from 0.0056 (classical) to 0.0059 (robust), and 2SLS gives 0.0057. Report the robust one.

# %% [markdown]
# ## Try it yourself (p. 27, 40)
#
# 1. Compare the Wu-Hausman row of `diag_iv` with the squared t statistic on `vhat`.
# 2. Add `sibs` to the instruments. What does the Sargan row do?
# 3. Which specification does RESET reject, and why?
# 4. How many degrees of freedom does the Breusch-Pagan statistic use here, and why?
