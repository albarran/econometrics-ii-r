# %% [markdown]
# # Unit 2D and 2E. Inference, heteroskedasticity and weighted least squares
#
# Econometrics II (Henry Redondo, 2026-27). R companion to **Unit 2, Blocks D and E** of the slides: analyses D2 to D8 and E1.
#
# Each section has the idea with the slide page (page of the PDF), then the code, then how to read the output. Numbers quoted as "slides" are the ones printed in the slides.

# %%
library(tidyverse)
library(broom)
library(estimatr)        # lm_robust()
library(marginaleffects) # hypotheses()
library(lmtest)          # bptest()
library(modelsummary)

options(digits = 4, readr.show_col_types = FALSE)
theme_set(theme_minimal())

# Henry's data folder. To read local copies, set the environment variable ECMTII_DATA to their folder.
data_dir <- Sys.getenv("ECMTII_DATA", "https://henryredondo.github.io/assets/files/teaching/econometrics-ii")

# %% [markdown]
# ## The data and the two hedonic models
#
# > **Data pending.** On 25 November 2026 the CSV files were not on the course website. The cell reads them from the address where they should be.

# %%
d <- read_csv(file.path(data_dir, "idealista_alicante_clase.csv")) |>
  mutate(zona    = fct_relevel(factor(zona), "Benalúa"),
         anio    = factor(anio),
         lprecio = log(precio),
         lm2     = log(m2)) |>
  left_join(read_csv(file.path(data_dir, "idealista_oraculo_clase.csv")), by = "id")

f_niv <- precio  ~ m2  + m2_terraza + terraza + antig + I(antig^2) + bajo + zona + anio
f_log <- lprecio ~ lm2 + terraza + m2_terraza + antig + I(antig^2) + bajo + zona + anio

# Robust (HC1) by default; the plain lm versions are for the classical tests
mn <- lm_robust(f_niv, data = d, se_type = "HC1"); mn_c <- lm(f_niv, data = d)
ml <- lm_robust(f_log, data = d, se_type = "HC1"); ml_c <- lm(f_log, data = d)

# %% [markdown]
# ## D2. Is a terrace worth anything? (2D p. 25-27)
#
# t test and 95% confidence interval for the terrace coefficient in the log model, with robust standard errors.

# %%
tidy(ml) |> filter(term == "terraza")

# The interval by hand: estimate +- t(0.975, N - k) x s.e.
with(tidy(ml) |> filter(term == "terraza"), estimate + c(-1, 1) * qt(0.975, df) * std.error)

# %% [markdown]
# **Reading.** Slides: 0.1041, s.e. 0.0189, t = 5.50, 95% CI [0.0669, 0.1413] with N − k = 976. The hypothesis of no premium is rejected at any conventional level and the interval is narrow. From Block C: holding quality fixed the premium is 0.033, outside this interval. A confidence interval measures sampling noise only; a precise estimate can be precisely wrong.

# %% [markdown]
# ## Joint tests: the four zone coefficients (2D p. 43)
#
# H0: the four zone coefficients are zero. `hypotheses(joint = ...)` takes a regular expression that selects the coefficients and uses the model's own variance: classical for `lm`, robust for `lm_robust`.

# %%
hypotheses(ml_c, joint = "^zona")                          # classical F
hypotheses(ml, joint = "^zona")                            # robust F
hypotheses(ml, joint = "^zona", joint_test = "chisq")      # robust Wald (chi2) = 4 x F

# The classical F is the comparison of restricted and unrestricted sums of squares
mr <- lm(update(f_log, . ~ . - zona), data = d)
anova(mr, ml_c)

# %% [markdown]
# **Reading.** Slides: the classical F is about 175 (4F = 699): the zones matter jointly.

# %% [markdown]
# ### Individual t tests against a joint F test
#
# From the previous course (`joint_test.R`): two highly correlated regressors, each with a true effect. The data cannot tell which of the two matters, and can tell that at least one does.

# %%
set.seed(123)
sim <- tibble(education  = rnorm(100),
              experience = 0.9 * education + rnorm(100, 0, 0.3),
              wage       = 1 + 0.5 * education + 0.5 * experience + rnorm(100))
m_sim <- lm(wage ~ education + experience, data = sim)

tidy(m_sim)
hypotheses(m_sim, joint = c("education", "experience"))
cor(sim$education, sim$experience)

# %% [markdown]
# **Reading.** Both t statistics are below 1.96 in absolute value and the joint F is 31.27: the two coefficients are individually imprecise because the regressors move together, and jointly they are clearly different from zero.

# %% [markdown]
# ## D4. The appraisers' rule (2D p. 30-31)
#
# Appraisers count a square metre of terrace as half a square metre of interior. In the levels model, H0: $\theta = \beta_{m2\_terraza} - 0.5\,\beta_{m2} = 0$. The variance of θ̂ needs the covariance of the two estimates; `hypotheses()` takes care of it.

# %%
bind_rows("Robust (HC1)" = hypotheses(mn,   "m2_terraza - 0.5 * m2 = 0"),
          "Classical"    = hypotheses(mn_c, "m2_terraza - 0.5 * m2 = 0"), .id = "se") |>
  as_tibble() |> select(se, estimate, std.error, statistic, p.value)

# By hand, with the covariance
V <- vcov(mn)
r <- c(m2_terraza = 1, m2 = -0.5)
theta    <- sum(r * coef(mn)[names(r)])
se_theta <- sqrt(drop(t(r) %*% V[names(r), names(r)] %*% r))
stopifnot(all.equal(se_theta, hypotheses(mn, "m2_terraza - 0.5 * m2 = 0")$std.error))
c(theta = theta, se = se_theta, t = theta / se_theta)

# %% [markdown]
# **Reading.** Slides: θ̂ = 667.55, robust s.e. 431.92, t = 1.55 (p = 0.12); classical s.e. 271.40, t = 2.46 (p = 0.014). With robust standard errors the 50% rule cannot be rejected; with classical ones it would be. The levels model is strongly heteroskedastic, so the robust verdict is the valid one. The point estimate says a terrace square metre is worth about 90% of an interior one: the data cannot tell 50% from 90%. (`hypotheses` reports a normal p-value, the slides a t one; with 976 degrees of freedom they coincide to the third decimal.)

# %% [markdown]
# ## D6. The LM test for the zones (2D p. 42-43)
#
# Recipe: (1) estimate the restricted model, without the zones, and keep its residuals; (2) regress them on all the regressors of the unrestricted model; (3) LM = N·R², χ² with 4 degrees of freedom under H0.

# %%
aux <- lm(update(f_log, resid(mr) ~ .), data = d)
LM  <- nobs(aux) * summary(aux)$r.squared

c(R2_aux = summary(aux)$r.squared, LM = LM, p = pchisq(LM, df = 4, lower.tail = FALSE),
  four_F = 4 * anova(mr, ml_c)$F[2])

# %% [markdown]
# **Reading.** Slides: R² of the auxiliary regression 0.4172, LM = 1,000 × 0.4172 = 417.2, against 4F = 699. Same conclusion as the F test, estimating only the restricted model. Under H0 the two are asymptotically equivalent; far from H0 they differ, and LM ≤ W.

# %% [markdown]
# ## D7. Is there heteroskedasticity? (2E p. 4-5)
#
# Breusch-Pagan: N·R² of the regression of the squared residuals on the fitted values.

# %%
bptest(mn_c, ~ fitted(mn_c))
bptest(ml_c, ~ fitted(ml_c))

# By hand for the levels model
bp_aux <- lm(resid(mn_c)^2 ~ fitted(mn_c))
stopifnot(all.equal(nobs(bp_aux) * summary(bp_aux)$r.squared, bptest(mn_c, ~ fitted(mn_c))$statistic[["BP"]]))

# Residual s.d. of the log model, beach against the rest
d |> mutate(u = resid(ml_c), beach = zona == "Playa de San Juan") |> summarise(sd_resid = sd(u), .by = beach)

# %% [markdown]
# **Reading.** Slides: χ²(1) = 122.80 in levels (p ≈ 10⁻²⁸) and 17.37 in logs (p = 3 × 10⁻⁵). Logs remove most of the heteroskedasticity: the residual s.d. is 0.148 away from the beach and 0.201 at the beach. Two responses: keep OLS with robust standard errors, or change the estimator and give noisy observations less weight (GLS).

# %% [markdown]
# ## D8. Weighted least squares with aggregated data (2E p. 12-13)
#
# Suppose we only receive averages: the mean log price per m² in each zone-by-year cell. A cell mean from 25 flats is more precise than one from 5, Var(ȳ_c) = σ²/n_c, so GLS weights each cell by its number of listings.

# %%
cells <- d |>
  mutate(lpm2 = log(precio / m2)) |>
  summarise(lpm2 = mean(lpm2), n = n(), .by = c(zona, anio))

cells |> summarise(cells = n(), min_n = min(n), mean_n = mean(n), max_n = max(n))

ols_cells <- lm(lpm2 ~ zona + anio, data = cells)
wls_cells <- lm(lpm2 ~ zona + anio, data = cells, weights = n)

modelsummary(list("OLS (equal weights)" = ols_cells, "WLS (weights n)" = wls_cells),
             coef_map = c("zonaCentro" = "Centro", "zonaPlaya de San Juan" = "Playa de San Juan", "anio2025" = "2025"),
             fmt = 4, gof_map = "nobs", output = "markdown")

# %% [markdown]
# **Reading.** Slides: 70 cells with between 5 and 25 listings (mean 14.3). Centro 0.0347 (0.0207) with OLS and 0.0419 (0.0181) with WLS; Playa de San Juan 0.5034 (0.0207) and 0.5034 (0.0196); 2025 0.5376 (0.0347) and 0.5447 (0.0317). The coefficients hardly move and the standard errors fall a little, because cell sizes are not very unequal. With cells of 5 and 5,000 observations weighting matters a lot.

# %% [markdown]
# ## E1. How much is a terrace worth? Every version (2E p. 19-20)
#
# The premium of a typical terrace (mean size among flats with one), in per cent, from each specification of the unit: $100\,(e^{\hat\beta_{terraza} + \hat\beta_{m2\_terraza}\,\bar s} - 1)$. The first row has no terrace size in the model.

# %%
m2t <- mean(d$m2_terraza[d$terraza == 1])

premium <- function(m) {
  b <- coef(m)
  100 * (exp(b[["terraza"]] + (if ("m2_terraza" %in% names(b)) b[["m2_terraza"]] * m2t else 0)) - 1)
}

specs <- list(
  "No controls"                 = lm(lprecio ~ terraza, data = d),
  "Hedonic model"               = ml_c,
  "+ proxy (number of photos)"  = lm(update(f_log, . ~ . + num_fotos), data = d),
  "+ quality (the oracle)"      = lm(update(f_log, . ~ . + calidad), data = d)
)

tibble(specification = names(specs),
       premium_pct   = unname(map_dbl(specs, premium)),
       on_130000_eur = 130000 * premium_pct / 100,
       R2            = unname(map_dbl(specs, \(m) summary(m)$r.squared)))

# %% [markdown]
# **Reading.** Slides: 46.0% with no controls (59,819 € on a 130,000 € flat), 19.2% with the hedonic model, 15.2% with the proxy, 11.4% with the oracle; the truth in the simulation is 13.0% (16,845 €). The naive number is three and a half times the truth and the hedonic model still overstates it by almost half. Each layer of econometrics moved the answer towards the truth, and none of them showed up in R².

# %% [markdown]
# ## Try it yourself (2D p. 47, 2E p. 22)
#
# 1. Test the 50% rule in the log model. What exactly is the hypothesis there?
# 2. Run `bptest(ml_c)`, with all the regressors in place of the fitted values. Same verdict? How many degrees of freedom?
# 3. Feasible GLS: regress `log(resid(ml_c)^2)` on the regressors, set `w = exp(fitted)` and refit with `weights = 1/w`. Does the terrace coefficient move? Do the standard errors fall?
# 4. In D8, why would robust standard errors on the unweighted cell regression not make WLS unnecessary?
