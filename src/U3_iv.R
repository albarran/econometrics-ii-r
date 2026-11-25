# %% [markdown]
# # Unit 3. Instrumental variables
#
# Econometrics II (Henry Redondo, 2026-27). R companion to **Unit 3** of the slides, one section per block:
#
# - **Block 3A.** The simple IV estimator (A0 to A4)
# - **Block 3B.** Two-stage least squares and weak instruments (B1 to B4, two Monte Carlo experiments)
# - **Blocks 3C and 3D.** GMM, the Hansen J test and multiple indicators (C1 to C3, D1)
#
# Each analysis has the idea in two or three lines with the slide page (page of the PDF of its block), then the code, then how to read the output. Numbers quoted as "slides" are the ones printed in the slides.
#
# Run the notebook from the top: each block uses objects created in the previous ones.

# %%
library(tidyverse)
library(broom)
library(estimatr)        # lm_robust(), iv_robust()
library(ivreg)           # ivreg()
library(gmm)             # gmm(), specTest()
library(marginaleffects) # hypotheses()
library(modelsummary)
library(fixest)          # feols(), one optional cell

options(digits = 4, readr.show_col_types = FALSE)
theme_set(theme_minimal())

# Henry's data folder. To read local copies, set the environment variable ECMTII_DATA to their folder.
data_dir <- Sys.getenv("ECMTII_DATA", "https://henryredondo.github.io/assets/files/teaching/econometrics-ii")

# %% [markdown]
# ## Block 3A. The simple IV estimator
#
# Page numbers in this block refer to `Unit3A_slides.pdf`.

# %% [markdown]
# ### Warm-up: IV in a simulated world
#
# Before the class data, a simulation where we know the truth (script `iv_regression.R` of the previous course). Earnings depend on education (true coefficient 2), experience and ability. Ability is unobserved and raises education, so OLS is inconsistent. The instrument `z` moves education and is independent of everything else.

# %%
set.seed(12345)
n <- 5000
sim <- tibble(
  ability  = rnorm(n, mean = 1, sd = 2),
  z        = rnorm(n),
  exper    = abs(rnorm(n, mean = 10, sd = 4)),
  educ     = 12 + 1.2 * ability + 0.8 * z + rnorm(n),
  u        = rnorm(n),
  earnings = 1 + 2 * educ + 1.5 * exper + 1.5 * ability + u
)

# The two conditions. Exogeneity can be checked here only because we simulated u.
sim |> summarise(cor_z_u = cor(z, u), cor_z_ability = cor(z, ability))
lm_robust(educ ~ exper + z, data = sim, se_type = "HC1") |>
  tidy() |>
  filter(term == "z")

# %%
models <- list(
  "OLS"         = lm_robust(earnings ~ educ + exper, data = sim, se_type = "HC1"),
  "IV"          = ivreg(earnings ~ educ + exper | z + exper, data = sim),
  "OLS, oracle" = lm_robust(earnings ~ educ + exper + ability, data = sim, se_type = "HC1")
)
modelsummary(models, gof_map = "nobs", output = "markdown")

# %% [markdown]
# **Reading.** The truth is 2. OLS is well above it, because `educ` picks up the effect of ability. IV recovers a value close to 2 without observing ability, with a larger standard error. The last column needs ability, which no real analyst has.
#
# In `ivreg` the formula is `y ~ regressors | instruments`: the exogenous regressors (`exper`) appear on both sides of the bar.

# %% [markdown]
# ### The data: Programa Transparenta (p. 4-5)
#
# 5,000 simulated firms. `transp` = 1 if the firm adopted pay transparency by 2025. A lottery assigned firms to three arms: control (`arm` = 0), letter (1), letter and software (2). The outcome is the change in the log gender wage gap between 2021 and 2025.
#
# > **Data pending.** On 25 November 2026 the CSV files were not on the course website. The cell below reads them from the address where they should be. The oracle file name (`transparency_oracle.csv`) and the firm identifier (`id`) are assumptions to confirm.

# %%
d <- read_csv(file.path(data_dir, "transparency_firms.csv")) |>
  mutate(
    sector      = factor(sector),
    province    = factor(province),
    gap2021     = log(wage_men_2021 / wage_women_2021),
    dgap        = log(wage_men_2025 / wage_women_2025) - gap2021,
    growth_men  = log(wage_men_2025 / wage_men_2021),
    growth_wom  = log(wage_women_2025 / wage_women_2021),
    lsize       = log(size),
    share_women = n_women / size,
    invited     = as.integer(arm >= 1),   # either arm
    letter      = as.integer(arm == 1),   # letter only
    software    = as.integer(arm == 2)    # letter + software
  ) |>
  # The oracle: the firm's egalitarian culture, which an analyst never observes
  left_join(read_csv(file.path(data_dir, "transparency_oracle.csv")), by = "id")

glimpse(d)

# %%
# Controls used in every regression, the two instruments of Block 3B,
# and a helper that builds formulas with them
X <- "lsize + share_women + sector + province"
Z <- "letter + software"
fml <- function(x) as.formula(str_glue(x, .envir = parent.frame()), env = globalenv())

deparse1(fml("dgap ~ transp + {X}"))

# %% [markdown]
# ### A0. Was the lottery fair? (p. 24)
#
# An instrument that comes from a lottery should be unrelated to anything that was fixed before the programme. We compare pre-programme means by arm and test that the 2021 gap is equal across arms.

# %%
d |>
  summarise(
    firms = n(),
    gap_2021 = mean(gap2021),
    employees = mean(size),
    share_women = mean(share_women),
    parent_abroad = mean(parent_abroad),
    fair_score = mean(fair_score),
    adopted_2025 = mean(transp),
    .by = arm
  ) |>
  arrange(arm)

# %%
# Joint test: the two arm coefficients are zero (robust)
lm_robust(gap2021 ~ factor(arm), data = d, se_type = "HC1") |> hypotheses(joint = "arm")
lm_robust(lsize ~ factor(arm), data = d, se_type = "HC1") |> hypotheses(joint = "arm")

# %% [markdown]
# **Reading.** Slides: the three arms look alike before the programme (gap 0.164, 0.163, 0.159; F = 0.74, p = 0.48; firm size p = 0.57), and adoption differs a lot afterwards (0.253, 0.305, 0.487). The lottery moved adoption and nothing else.

# %% [markdown]
# ### A1. The naive answer, and what the oracle says (p. 25)
#
# OLS of `dgap` on `transp` and the controls, with robust (HC1) standard errors. Culture is in the error term: egalitarian firms adopt more often and close their gap faster.

# %%
ols <- lm_robust(fml("dgap ~ transp + {X}"), data = d, se_type = "HC1")
ols_oracle <- lm_robust(fml("dgap ~ transp + culture + {X}"), data = d, se_type = "HC1")

modelsummary(list("OLS" = ols, "OLS + culture" = ols_oracle),
  coef_map = c("transp", "culture"), fmt = 4, gof_map = "nobs", output = "markdown"
)

d |> summarise(mean_culture = mean(culture), .by = transp)

# %% [markdown]
# **Reading.** Slides: OLS gives −0.0561 (s.e. 0.0019). With culture in the regression the coefficient falls to −0.0310, so almost half of the OLS number was culture. Adopters have mean culture +0.54 and non-adopters −0.37. The target for the rest of the unit is −0.031, without using `culture`.

# %% [markdown]
# ### A2. The IV estimator, by hand (p. 28)
#
# Instrument: `invited` (either arm). With one endogenous regressor and one instrument, the IV estimator is the **Wald ratio**: reduced form divided by first stage.

# %%
first_stage <- lm_robust(fml("transp ~ invited + {X}"), data = d, se_type = "HC1")
reduced_form <- lm_robust(fml("dgap ~ invited + {X}"), data = d, se_type = "HC1")
iv <- ivreg(fml("dgap ~ transp + {X} | invited + {X}"), data = d)

wald <- coef(reduced_form)[["invited"]] / coef(first_stage)[["invited"]]
c(
  first_stage = coef(first_stage)[["invited"]], reduced_form = coef(reduced_form)[["invited"]],
  wald = wald, ivreg = coef(iv)[["transp"]]
)

# The ratio is the ivreg estimate, exactly
stopifnot(all.equal(wald, coef(iv)[["transp"]]))

# %%
modelsummary(
  list("First stage (transp)" = first_stage, "Reduced form (dgap)" = reduced_form, "IV (dgap)" = iv),
  vcov = list(NULL, NULL, "HC0"), # the lm_robust models keep their HC1; IV with HC0 as in the slides
  coef_map = c("invited", "transp"), fmt = 4, gof_map = "nobs", output = "markdown"
)

# %% [markdown]
# The same in one line with `estimatr::iv_robust`, which computes the robust standard errors itself. Its default is HC2; the slides use HC0 for IV.

# %%
iv_r <- iv_robust(fml("dgap ~ transp + {X} | invited + {X}"), data = d, se_type = "HC0", diagnostics = TRUE)
tidy(iv_r) |> filter(term == "transp")

# Same coefficients and same standard errors as ivreg + sandwich HC0
stopifnot(
  all.equal(coef(iv_r), coef(iv)),
  all.equal(iv_r$std.error, sqrt(diag(sandwich::vcovHC(iv, type = "HC0"))))
)

# First-stage F: classical (ivreg) and robust (iv_robust)
summary(iv, diagnostics = TRUE)$diagnostics["Weak instruments", ]
iv_r$diagnostic_first_stage_fstatistic

# %% [markdown]
# **Reading.** Slides: first stage 0.1395 (s.e. 0.0129), reduced form −0.0052 (0.0020), IV −0.0375 (0.0132), F = 118. Invited firms were 14 pp more likely to adopt and their gap fell 0.5 pp more (the intention-to-treat effect). If the invitation works only through adoption, the effect of adopting is 0.52/0.14 ≈ 3.7 pp. IV rescales the effect of the lottery by how much the lottery moved the regressor.

# %% [markdown]
# ### A3. What bad instruments do (p. 30)
#
# Two other candidates are strongly related to adoption: `parent_abroad` and `peer_share`. Both are also related to culture, so they fail exogeneity. We run the same IV with each instrument.

# %%
iv_with <- function(z) {
  m <- iv_robust(fml("dgap ~ transp + {X} | {z} + {X}"), data = d, se_type = "HC0", diagnostics = TRUE)
  tibble(
    instrument = z,
    first_stage_F = m$diagnostic_first_stage_fstatistic[["value"]],
    beta = coef(m)[["transp"]],
    se = m$std.error[["transp"]],
    cor_culture = cor(d[[z]], d$culture)
  )
}

map(c("invited", "parent_abroad", "peer_share"), iv_with) |>
  list_rbind() |>
  add_row(instrument = "OLS", beta = coef(ols)[["transp"]], se = ols$std.error[["transp"]])

# %% [markdown]
# **Reading.** Slides: F = 118, 72 and 93, so all three are relevant. The estimates are −0.038 (invited), −0.111 (parent_abroad) and −0.188 (peer_share), against −0.056 for OLS and −0.031 for the oracle. The two invalid instruments land further from the truth than OLS, with small standard errors. Their correlation with culture is 0.14 and 0.42 (0.001 for `invited`). A strong first stage says nothing about validity; here we can see it only because we have the oracle.

# %% [markdown]
# ### A4. The price of identification (p. 40)
#
# IV uses only the variation in `transp` that comes from the lottery, so it is much less precise than OLS.

# %%
bind_rows("OLS" = tidy(ols), "IV (invited)" = tidy(iv_r), .id = "model") |>
  filter(term == "transp") |>
  select(model, estimate, std.error, conf.low, conf.high)

iv_r$std.error[["transp"]] / ols$std.error[["transp"]]

# %%
modelplot(list("OLS" = ols, "IV (invited)" = iv_r), coef_map = "transp") +
  geom_vline(xintercept = coef(ols_oracle)[["transp"]], linetype = "dashed") +
  labs(x = "Effect of transparency on the change in the gap (dashed: oracle)", y = NULL)

# %% [markdown]
# **Reading.** Slides: OLS −0.056, CI [−0.060, −0.052]; IV −0.038, CI [−0.063, −0.012]. The IV standard error is 6.8 times the OLS one. The OLS interval is narrow and leaves the oracle value out; the IV interval is wide and contains it.

# %% [markdown]
# ### Try it yourself (p. 42)
#
# 1. Divide the reduced-form coefficient by the first-stage one with more decimals. Do you get the `ivreg` estimate exactly?
# 2. Drop the controls from the three regressions (`fml("dgap ~ transp | invited")`). How much does the IV estimate move? Why so little?
# 3. Compare the IV standard error from `vcov(iv)` (classical), `sandwich::vcovHC(iv, type = "HC0")` and `iv_robust(..., se_type = "HC2")`.

# %% [markdown]
# ## Block 3B. Two-stage least squares and weak instruments
#
# Page numbers in this block refer to `Unit3B_slides.pdf`.

# %% [markdown]
# ### B1. The first stage with both arms (p. 32)
#
# Two instruments for one endogenous regressor: the equation is overidentified. The first stage regresses `transp` on both instruments and the controls; relevance is the joint test that both instrument coefficients are zero.

# %%
fs <- lm_robust(fml("transp ~ {Z} + {X}"), data = d, se_type = "HC1")
tidy(fs) |> filter(term %in% c("letter", "software")) |> select(term, estimate, std.error, statistic)

# Joint robust F test of H0: letter = software = 0
hypotheses(fs, joint = c("letter", "software"))

# %% [markdown]
# **Reading.** Slides: letter 0.0482 (s.e. 0.0153, t = 3.2), software 0.2366 (0.0165, t = 14.3), joint F = 103.6. A letter raises adoption by 5 pp and the software by 24 pp. Keep in mind the letter's t of 3.2 for B4.

# %% [markdown]
# ### B2. 2SLS with both arms, and the two-step trap (p. 33)
#
# 2SLS combines the instruments through the first-stage fitted values. Doing the two stages by hand gives the right coefficient and the **wrong standard errors**, because stage 2 forms residuals with the fitted `transp` where it should use the observed one.

# %%
iv2 <- ivreg(fml("dgap ~ transp + {X} | {Z} + {X}"), data = d)

# By hand: stage 1 fitted values, then OLS in stage 2
d <- d |> mutate(transp_hat = fitted(lm(fml("transp ~ {Z} + {X}"), data = d)))
stage2 <- lm(fml("dgap ~ transp_hat + {X}"), data = d)

# Same coefficient
stopifnot(all.equal(coef(stage2)[["transp_hat"]], coef(iv2)[["transp"]]))

modelsummary(
  list("ivreg, robust" = iv2, "ivreg, classical" = iv2,
       "by hand, robust" = stage2, "by hand, classical" = stage2),
  vcov = list("HC0", "classical", "HC1", "classical"),
  coef_map = c("transp" = "transp", "transp_hat" = "transp"),
  fmt = 4, gof_map = "nobs", output = "markdown"
)

# %%
# Why: the residual standard deviation
X_obs <- model.matrix(fml("~ transp + {X}"), data = d)             # observed regressors
u_ok  <- d$dgap - drop(X_obs %*% coef(iv2))                         # correct 2SLS residuals
stopifnot(all.equal(unname(u_ok), unname(resid(iv2))))
c(correct = sqrt(sum(u_ok^2) / df.residual(iv2)), by_hand = sigma(stage2))

# %% [markdown]
# **Reading.** Slides: −0.0408 with both methods; robust s.e. 0.0098 (`ivreg`) against 0.0105 (by hand), classical 0.0095 against 0.0102; residual s.d. 0.0647 (correct) against 0.0692 (by hand). The 95% CI is [−6.0, −2.2] pp. With `invited` alone the s.e. was 0.0132: using both arms buys back precision. Use a 2SLS command.

# %% [markdown]
# #### 2SLS in matrices, once (p. 6)
#
# $\hat\beta_{2SLS} = (\hat X'X)^{-1}\hat X'y$ with $\hat X = Z(Z'Z)^{-1}Z'X$.

# %%
y     <- d$dgap
Z_all <- model.matrix(fml("~ {Z} + {X}"), data = d)                 # instruments, exogenous regressors included
X_hat <- Z_all %*% solve(crossprod(Z_all), crossprod(Z_all, X_obs))
b_2sls <- solve(crossprod(X_hat, X_obs), crossprod(X_hat, y))

stopifnot(all.equal(drop(b_2sls), coef(iv2), check.attributes = FALSE))
b_2sls["transp", ]

# %% [markdown]
# #### The same in one line
#
# `estimatr::iv_robust` gives 2SLS with robust standard errors and the diagnostics. `fixest::feols` uses the syntax `y ~ exogenous | endogenous ~ instruments`.

# %%
iv2_r <- iv_robust(fml("dgap ~ transp + {X} | {Z} + {X}"), data = d, se_type = "HC0", diagnostics = TRUE)
stopifnot(all.equal(coef(iv2_r), coef(iv2)),
          all.equal(iv2_r$std.error, sqrt(diag(sandwich::vcovHC(iv2, type = "HC0")))))
tidy(iv2_r) |> filter(term == "transp")

iv2_fe <- feols(fml("dgap ~ {X} | transp ~ {Z}"), data = d, vcov = "hetero")
stopifnot(all.equal(coef(iv2_fe)[["fit_transp"]], coef(iv2)[["transp"]]))
coeftable(iv2_fe)["fit_transp", ]
fitstat(iv2_fe, ~ ivf + ivwald)   # first-stage F: classical and robust

# %% [markdown]
# `feols` with `vcov = "hetero"` applies the small-sample factor N/(N−k) (HC1), so its standard error differs from HC0 in the fourth significant digit.

# %% [markdown]
# ### B3. Who moves, men or women? (p. 37)
#
# The gap can close because women's wages rise or because men's wages grow more slowly. Same 2SLS, with wage growth of each group as the outcome.

# %%
iv_men <- iv_robust(fml("growth_men ~ transp + {X} | {Z} + {X}"), data = d, se_type = "HC0")
iv_wom <- iv_robust(fml("growth_wom ~ transp + {X} | {Z} + {X}"), data = d, se_type = "HC0")

modelsummary(list("Men" = iv_men, "Women" = iv_wom, "Gap (men - women)" = iv2_r),
             coef_map = "transp", fmt = 4, gof_map = "nobs", output = "markdown")

# 2SLS is linear in y: the two effects add up to the effect on the gap
stopifnot(all.equal(coef(iv_men)[["transp"]] - coef(iv_wom)[["transp"]], coef(iv2)[["transp"]]))

# %% [markdown]
# **Reading.** Slides: men −0.0294 (s.e. 0.0064), women +0.0113 (0.0066), difference −0.0408. Most of the effect is slower wage growth for men. OLS gives −0.036 and +0.020: it overstates both, because culture moves both.

# %% [markdown]
# ### B4. What if we only had the letter? (p. 38)
#
# Drop the software arm and use the letter as the only instrument. Its first stage is weak.

# %%
iv_arm <- function(data, z) {
  m <- iv_robust(fml("dgap ~ transp + {X} | {z} + {X}"), data = data, se_type = "HC0", diagnostics = TRUE)
  tidy(m) |> filter(term == "transp") |>
    transmute(instrument = z, n = nobs(m), first_stage_F = m$diagnostic_first_stage_fstatistic[["value"]],
              estimate, std.error, conf.low, conf.high)
}

bind_rows(
  iv_arm(filter(d, arm != 2), "letter"),     # control and letter arms
  iv_arm(filter(d, arm != 1), "software")    # control and software arms
)

# %% [markdown]
# **Reading.** Slides, letter only (N = 3,813): first stage 0.048 (0.015), F = 9.8, 2SLS −0.027 with s.e. 0.045 and CI [−0.116, 0.062]. Software only: F = 208, −0.041 (0.010). With the letter alone the interval runs from "closes the gap by 12 pp" to "widens it by 6 pp". One sample cannot tell whether −0.027 is close to the truth; a simulation can.

# %% [markdown]
# ### One weak instrument: Monte Carlo (p. 39)
#
# We draw many samples from a simulated world with a known effect and look at the distribution of the 2SLS estimate.
#
# > The slides use Henry's script `weak_iv_mc.R`, which is not published. The data-generating process below is a simplified one written for this notebook (no controls, 500 samples), calibrated to the first stages of B1. The percentages will be close to the slides, not equal.

# %%
truth <- -0.032

# One sample: lottery, culture (unobserved), adoption and change in the gap
draw <- function(n, arms = c(0, 1, 2), cells = 18) {
  tibble(
    arm     = sample(arms, n, replace = TRUE, prob = c(0.5, 0.25, 0.25)[arms + 1]),
    cell    = sample(cells, n, replace = TRUE),
    culture = rnorm(n),
    transp  = as.integer(runif(n) < 0.25 + 0.05 * (arm == 1) + 0.24 * (arm == 2) + 0.17 * culture),
    dgap    = truth * transp - 0.030 * culture + rnorm(n, sd = 0.065)
  )
}

# 2SLS and first-stage F in matrices (fast): returns the coefficient on x
tsls <- function(y, x, Zmat) {
  Zc   <- cbind(1, Zmat)
  xhat <- drop(Zc %*% solve(crossprod(Zc), crossprod(Zc, x)))
  q    <- ncol(Zmat)
  f    <- ((sum((x - mean(x))^2) - sum((x - xhat)^2)) / q) / (sum((x - xhat)^2) / (length(x) - q - 1))
  c(beta = cov(xhat, y) / cov(xhat, x), F = f, ols = cov(x, y) / var(x))
}

one_rep <- function(design) {
  s <- switch(design,
    "both arms, 5,000 firms"   = draw(5000),
    "letter only, 5,000 firms" = draw(5000, arms = c(0, 1)),
    "letter only, 1,500 firms" = draw(1500, arms = c(0, 1)))
  Zmat <- if (design == "both arms, 5,000 firms") cbind(s$arm == 1, s$arm == 2) else cbind(s$arm == 1)
  tsls(s$dgap, s$transp, Zmat)
}

set.seed(2026)
designs <- c("both arms, 5,000 firms", "letter only, 5,000 firms", "letter only, 1,500 firms")
mc <- map(designs, \(dsg) map(1:500, \(i) one_rep(dsg)) |> bind_rows() |> mutate(design = dsg)) |>
  list_rbind() |>
  mutate(design = fct_inorder(design)) |>
  filter(is.finite(beta))   # drops the rare sample whose first stage is exactly zero

mc |> summarise(median_F = median(F), median_beta = median(beta), mean_ols = mean(ols),
                share_5pp_off = mean(abs(beta - truth) > 0.05), .by = design)

# %%
ggplot(mc, aes(beta, colour = design)) +
  geom_density() +
  geom_vline(xintercept = truth, linetype = "dashed") +
  coord_cartesian(xlim = c(-0.2, 0.1)) +
  labs(x = "Estimated effect of transparency (dashed: truth)", y = NULL, colour = NULL) +
  theme(legend.position = "bottom")

# %% [markdown]
# **Reading.** Slides (1,000 samples): median F of 117, 12 and 4; share of estimates more than 5 pp from the truth 0%, 27% and 54%. With one weak instrument the distribution stays roughly centred on the truth and its tails explode: any single estimate can land almost anywhere.

# %% [markdown]
# ### Many weak instruments: Monte Carlo (p. 40)
#
# "Use more information": interact the letter with the 18 sector-by-province cells and use the 18 interactions as instruments. Each one carries almost no signal.

# %%
one_rep_many <- function() {
  s <- draw(5000, arms = c(0, 1))
  Zmat <- model.matrix(~ 0 + factor(cell), data = s) * (s$arm == 1)   # letter x 18 cells
  tsls(s$dgap, s$transp, Zmat)
}

set.seed(2026)
mc_many <- map(1:500, \(i) one_rep_many()) |> bind_rows()

mc_many |> summarise(median_F = median(F), median_beta = median(beta), mean_ols = mean(ols))

# %% [markdown]
# **Reading.** Slides: joint F ≈ 1.6; the median estimate moves from −0.032 (truth) to −0.046, most of the way to OLS (−0.057). The first stage overfits, the fitted values start to reproduce `transp` itself, endogeneity included, and 2SLS drifts toward OLS with standard errors that look reasonable.

# %% [markdown]
# ### Try it yourself (p. 41)
#
# 1. Do the coefficients for `growth_men` and `growth_wom` add up to the one for `dgap`? (B3 checks it with `stopifnot`.)
# 2. Use the software arm only (`filter(d, arm != 1)`). Compare F and the standard error with B4.
# 3. Why are the by-hand standard errors of B2 wrong?

# %% [markdown]
# ## Blocks 3C and 3D. GMM, overidentification and multiple indicators
#
# Page numbers in this block refer to `Unit3C_slides.pdf` (C1 to C3) and `Unit3D_slides.pdf` (D1 and the final comparison). `iv2` and `iv2_r` are the 2SLS estimates of Block 3B.

# %% [markdown]
# ### C1. Small firms, noisy gaps (p. 24)
#
# The gap is an average over a firm's workers, so it is noisier in small firms. If the variance of the error depends on firm size, 2SLS is no longer the efficient GMM estimator and robust standard errors are needed.

# %%
d |>
  mutate(uhat = resid(iv2),
         size_class = cut(size, c(0, 25, 50, 100, Inf), labels = c("10-25", "26-50", "51-100", "more than 100"))) |>
  summarise(firms = n(), sd_uhat = sd(uhat), .by = size_class) |>
  arrange(size_class)

# %% [markdown]
# **Reading.** Slides: the residual s.d. is 0.082, 0.067, 0.055 and 0.047 from the smallest to the largest firms (965, 1,901, 1,428 and 706 firms). It is 1.8 times larger in the smallest firms: heteroskedasticity.

# %% [markdown]
# ### C2. Two-step efficient GMM (p. 25)
#
# With two instruments there are two moment conditions for one parameter. GMM minimises a weighted distance of the sample moments from zero; the efficient weight is the inverse of the variance of the moments. The minimised criterion is the **Hansen J** statistic, which tests the overidentifying restriction.

# %%
g <- gmm(fml("dgap ~ transp + {X}"), fml("~ {Z} + {X}"), data = as.data.frame(d),   # gmm() needs a plain data frame
         vcov = "MDS", centeredVcov = FALSE)

bind_rows("2SLS, robust (iv_robust)" = tidy(iv2_r), "GMM two-step (gmm)" = tidy(g), .id = "estimator") |>
  filter(term == "transp") |>
  select(estimator, estimate, std.error)

specTest(g)   # Hansen J

# J and its p-value as numbers, to reuse below
hansen <- function(m) { t <- specTest(m)$test; c(J = as.numeric(t[1]), p = as.numeric(t[2])) }

# %% [markdown]
# #### The two steps in matrices
#
# Step 1: 2SLS residuals $\hat u_i$. Step 2: $\hat W = \left(\frac{1}{N}\sum_i \hat u_i^2 z_i z_i'\right)^{-1}$, $\hat\beta_{GMM} = (X'Z\hat W Z'X)^{-1}X'Z\hat W Z'y$ and $J = N\,\bar g'\hat W \bar g$ with $\bar g = \frac{1}{N}Z'(y - X\hat\beta_{GMM})$.

# %%
y  <- d$dgap
Xm <- model.matrix(fml("~ transp + {X}"), data = d)
Zm <- model.matrix(fml("~ {Z} + {X}"), data = d)
N  <- nrow(d)

u1 <- resid(iv2)                                   # step 1
W  <- solve(crossprod(Zm * u1) / N)                # inverse of (1/N) sum u^2 z z'
A  <- crossprod(Xm, Zm) %*% W
b_gmm <- drop(solve(A %*% crossprod(Zm, Xm), A %*% crossprod(Zm, y)))   # step 2
gbar  <- crossprod(Zm, y - Xm %*% b_gmm) / N
J     <- N * drop(t(gbar) %*% W %*% gbar)

c(by_hand = b_gmm[["transp"]], gmm = coef(g)[["transp"]])
c(J_by_hand = J, J_gmm = hansen(g)[["J"]], p_by_hand = pchisq(J, df = ncol(Zm) - ncol(Xm), lower.tail = FALSE))

stopifnot(all.equal(b_gmm[["transp"]], coef(g)[["transp"]], tolerance = 1e-4),
          all.equal(J, hansen(g)[["J"]], tolerance = 0.02))

# %% [markdown]
# **Reading.** Slides: `ivreg` robust −0.0408 (s.e. 0.0098), two-step GMM −0.0408 (0.0098), Hansen J = 0.14 with 1 d.f., p = 0.71. The estimate is the same to the fourth decimal: the lottery makes the instruments unrelated to firm size, hence to the variance of the error, and reweighting the moments changes almost nothing. J asks whether the two arms tell the same story; here they do.
#
# The `gmm` package may rebuild the weight matrix at the final estimate, so its J and standard errors can differ from the by-hand ones in the last digits (the checks above allow for it).

# %% [markdown]
# ### C3. Adding a bad instrument (p. 26)
#
# Add `peer_share`, which is related to culture, to the instrument list. Now there are three moments and one parameter.

# %%
g3 <- gmm(fml("dgap ~ transp + {X}"), fml("~ {Z} + peer_share + {X}"), data = as.data.frame(d),
          vcov = "MDS", centeredVcov = FALSE)

tibble(instruments = c("letter software", "letter software peer_share"),
       beta        = c(coef(g)[["transp"]], coef(g3)[["transp"]]),
       hansen_J    = c(hansen(g)[["J"]], hansen(g3)[["J"]]),
       p_value     = c(hansen(g)[["p"]], hansen(g3)[["p"]]))

# The overidentification test of iv_robust (robust, from 2SLS) for the same instrument list
iv_robust(fml("dgap ~ transp + {X} | {Z} + peer_share + {X}"), data = d,
          se_type = "HC0", diagnostics = TRUE)$diagnostic_overid_test

# %% [markdown]
# **Reading.** Slides: with `peer_share` the estimate doubles (−0.083 against −0.041) and J jumps from 0.14 (p = 0.71) to 61.3 (p < 10⁻¹³; 61.7 in the `gmm` package). The moments contradict each other: no value of the coefficient makes the three close to zero. The test has power only if at least one instrument is valid, and it cannot say which one fails.

# %% [markdown]
# ### D1. Proxy against multiple indicators (p. 8-9)
#
# Without a lottery: `fair_score` and `women_mgmt` are two noisy measures of culture. Used as a control (a **proxy**), each one removes part of the bias. Used as **multiple indicators**, one enters the equation and the other is its instrument, which removes the noise problem. `transp` is its own instrument here: the assumption is that culture is the only reason why adoption is endogenous.

# %%
d |> summarise(cor_fair_culture   = cor(fair_score, culture),
               cor_wmgmt_culture  = cor(women_mgmt, culture),
               cor_fair_wmgmt     = cor(fair_score, women_mgmt))

# %%
ols     <- lm_robust(fml("dgap ~ transp + {X}"), data = d, se_type = "HC1")
proxy1  <- lm_robust(fml("dgap ~ transp + fair_score + {X}"), data = d, se_type = "HC1")
proxy2  <- lm_robust(fml("dgap ~ transp + fair_score + women_mgmt + {X}"), data = d, se_type = "HC1")
ind1    <- iv_robust(fml("dgap ~ transp + fair_score + {X} | transp + women_mgmt + {X}"), data = d,
                     se_type = "HC0", diagnostics = TRUE)
ind2    <- iv_robust(fml("dgap ~ transp + women_mgmt + {X} | transp + fair_score + {X}"), data = d,
                     se_type = "HC0", diagnostics = TRUE)
oracle  <- lm_robust(fml("dgap ~ transp + culture + {X}"), data = d, se_type = "HC1")

d1 <- list("OLS" = ols, "Proxy: fair_score" = proxy1, "Proxy: both" = proxy2,
           "fair_score IV women_mgmt" = ind1, "women_mgmt IV fair_score" = ind2, "Oracle" = oracle)

map(d1, tidy) |> list_rbind(names_to = "specification") |>
  filter(term == "transp") |> select(specification, estimate, std.error)

# First-stage F of the indicator that works as instrument
c(ind1 = ind1$diagnostic_first_stage_fstatistic[["value"]], ind2 = ind2$diagnostic_first_stage_fstatistic[["value"]])

# %% [markdown]
# **Reading.** Slides: correlations with culture 0.73 and 0.54, between indicators 0.41. OLS −0.0561; proxy `fair_score` −0.0444; both proxies −0.0415; multiple indicators −0.0316 and −0.0301 (s.e. 0.0026); oracle −0.0310. The proxies close part of the distance between OLS and the oracle; instrumenting one indicator with the other matches the oracle, with a standard error of 0.0026 against 0.0098 for 2SLS with the lottery. It works because the two measurement errors are unrelated and culture is the only source of endogeneity.

# %% [markdown]
# ### Every estimator against the truth (3D p. 12)
#
# Points and 95% robust confidence intervals for the coefficient on `transp`.

# %%
iv1_r <- iv_robust(fml("dgap ~ transp + {X} | invited + {X}"), data = d, se_type = "HC0")

all_models <- list(
  "OLS"                       = ols,
  "Proxy (fair_score)"        = proxy1,
  "IV (invited)"              = iv1_r,
  "2SLS (both arms)"          = iv2_r,
  "GMM (both arms)"           = g,
  "Multiple indicators"       = ind1,
  "Oracle (culture observed)" = oracle
)

estimates <- map(all_models, \(m) tidy(m, conf.int = TRUE)) |>
  list_rbind(names_to = "estimator") |>
  filter(term == "transp") |>
  mutate(estimator = fct_rev(fct_inorder(estimator)))

ggplot(estimates, aes(estimate, estimator)) +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high)) +
  geom_vline(xintercept = coef(oracle)[["transp"]], linetype = "dashed") +
  labs(x = "Effect of transparency on the change in the gap (dashed: oracle)", y = NULL)

# %% [markdown]
# **Reading.** OLS and the proxy are precise and biased. The lottery-based estimators cover the truth (−0.032 in the simulation); they are noisier, and 2SLS and GMM with both arms are the most precise of them. The indicators match the oracle if their assumptions hold.

# %% [markdown]
# ### Try it yourself (3C p. 27, 3D p. 10)
#
# 1. Compare `coef(g)` with `coef(iv2)`. They agree to four decimals: why is the efficiency gain zero here?
# 2. Replace `peer_share` by `parent_abroad` in C3. Does J catch it too? Why is it harder?
# 3. Why is `transp` not instrumented in D1? What are we assuming about it?
# 4. Could you add `letter` and `software` as extra instruments in D1? For which regressor?
