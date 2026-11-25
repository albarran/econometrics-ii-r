# %% [markdown]
# # Unit 3B. Two-stage least squares and weak instruments
#
# Econometrics II (Henry Redondo, 2026-27). R companion to **Unit 3, Block B** of the slides: analyses B1 to B4 and the two Monte Carlo experiments on weak instruments.
#
# Each section has the idea with the slide page (page of the PDF), then the code, then how to read the output. Numbers quoted as "slides" are the ones printed in the slides.

# %%
library(tidyverse)
library(broom)
library(estimatr)        # lm_robust(), iv_robust()
library(ivreg)           # ivreg()
library(marginaleffects) # hypotheses()
library(modelsummary)
library(fixest)          # feols(), one optional cell

options(digits = 4, readr.show_col_types = FALSE)
theme_set(theme_minimal())

# Henry's data folder. To read local copies, set the environment variable ECMTII_DATA to their folder.
data_dir <- Sys.getenv("ECMTII_DATA", "https://henryredondo.github.io/assets/files/teaching/econometrics-ii")

# %% [markdown]
# ## The data (same as Block A)
#
# > **Data pending.** On 25 November 2026 the CSV files were not on the course website. The cell reads them from the address where they should be.

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
    invited     = as.integer(arm >= 1),
    letter      = as.integer(arm == 1),   # letter only
    software    = as.integer(arm == 2)    # letter + software
  )

X <- "lsize + share_women + sector + province"   # controls
Z <- "letter + software"                         # instruments
fml <- function(x) as.formula(str_glue(x, .envir = parent.frame()), env = globalenv())

# %% [markdown]
# ## B1. The first stage with both arms (p. 32)
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
# ## B2. 2SLS with both arms, and the two-step trap (p. 33)
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
# ### 2SLS in matrices, once (p. 6)
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
# ### The same in one line
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
# ## B3. Who moves, men or women? (p. 37)
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
# ## B4. What if we only had the letter? (p. 38)
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
# ## One weak instrument: Monte Carlo (p. 39)
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
# ## Many weak instruments: Monte Carlo (p. 40)
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
# ## Try it yourself (p. 41)
#
# 1. Do the coefficients for `growth_men` and `growth_wom` add up to the one for `dgap`? (B3 checks it with `stopifnot`.)
# 2. Use the software arm only (`filter(d, arm != 1)`). Compare F and the standard error with B4.
# 3. Why are the by-hand standard errors of B2 wrong?
