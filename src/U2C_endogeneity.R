# %% [markdown]
# # Unit 2C. Robust standard errors, omitted variables, proxies and measurement error
#
# Econometrics II (Henry Redondo, 2026-27). R companion to **Unit 2, Block C** of the slides: analyses D1 (classical against robust standard errors) and C1 to C4, plus two simulations from the previous course (proxy variables and measurement error).
#
# Each section has the idea with the slide page (page of the PDF), then the code, then how to read the output. Numbers quoted as "slides" are the ones printed in the slides.

# %%
library(tidyverse)
library(broom)
library(estimatr)        # lm_robust()
library(marginaleffects) # hypotheses()
library(modelsummary)

options(digits = 4, readr.show_col_types = FALSE)
theme_set(theme_minimal())

# Henry's data folder. To read local copies, set the environment variable ECMTII_DATA to their folder.
data_dir <- Sys.getenv("ECMTII_DATA", "https://henryredondo.github.io/assets/files/teaching/econometrics-ii")

# %% [markdown]
# ## The data, with the oracle
#
# The unobserved quality of each flat (`calidad`) is kept in a separate file. An analyst never has it; we open it to see what we are missing.
#
# > **Data pending.** On 25 November 2026 the CSV files were not on the course website. The cell reads them from the address where they should be.

# %%
d <- read_csv(file.path(data_dir, "idealista_alicante_clase.csv")) |>
  mutate(zona        = fct_relevel(factor(zona), "Benalúa"),
         anio        = factor(anio),
         lprecio     = log(precio),
         lm2         = log(m2),
         lm2_anuncio = log(m2_anuncio)) |>
  left_join(read_csv(file.path(data_dir, "idealista_oraculo_clase.csv")), by = "id")

f_niv <- precio  ~ m2  + m2_terraza + terraza + antig + I(antig^2) + bajo + zona + anio
f_log <- lprecio ~ lm2 + terraza + m2_terraza + antig + I(antig^2) + bajo + zona + anio

# %% [markdown]
# ## D1. Classical against robust standard errors (p. 17-18)
#
# The classical formula assumes that the error variance is the same for every flat. The robust (White, here HC1) formula does not. `modelsummary` can show the same model with both.

# %%
mn <- lm(f_niv, data = d)
ml <- lm(f_log, data = d)

modelsummary(list("Levels, classical" = mn, "Levels, robust" = mn, "Logs, classical" = ml, "Logs, robust" = ml),
             vcov = c("classical", "HC1", "classical", "HC1"),
             coef_map = c("m2" = "m2 / log(m2)", "lm2" = "m2 / log(m2)", "m2_terraza" = "m2_terraza",
                          "terraza" = "terraza", "bajo" = "bajo"),
             fmt = fmt_significant(4), gof_map = "nobs", output = "markdown")

# %% [markdown]
# ### The sandwich, once (p. 15)
#
# $\widehat{Var}(\hat\beta) = (X'X)^{-1}\left(\sum_i \hat u_i^2 x_i x_i'\right)(X'X)^{-1}$, times N/(N−k) for HC1. `lm_robust(se_type = "HC1")` does this in one step and is what the rest of the notebooks use.

# %%
X  <- model.matrix(ml)
u  <- resid(ml)
N  <- nrow(X); k <- ncol(X)
bread <- solve(crossprod(X))
meat  <- crossprod(X * u)                       # sum of u_i^2 x_i x_i'
V_hc1 <- N / (N - k) * bread %*% meat %*% bread

ml_r <- lm_robust(f_log, data = d, se_type = "HC1")
stopifnot(all.equal(sqrt(diag(V_hc1)), ml_r$std.error, check.attributes = FALSE))

c(by_hand = sqrt(diag(V_hc1))[["terraza"]], lm_robust = ml_r$std.error[["terraza"]],
  classical = tidy(ml)$std.error[tidy(ml)$term == "terraza"])

# %% [markdown]
# **Reading.** Slides, ratio robust/classical: in levels `m2` 1.56 (46.12 → 71.93), `m2_terraza` 1.60, `terraza` 1.24, `bajo` 0.84; in logs 1.05, 1.26, 1.11 and 0.99. In euros a 400,000 € beach flat misses by far more euros than a 70,000 € flat: strong heteroskedasticity, and the classical s.e. of `m2` is 36% too small. In logs the two are much closer. Robust errors can also be smaller (`bajo`). Rule: report robust standard errors by default.

# %% [markdown]
# ## C1. Opening the oracle on omitted quality (p. 43-44)
#
# Short regression (the analyst's), long regression (with `calidad`) and auxiliary regression (`calidad` on the regressors of the short one). In the sample the omitted-variable formula is an identity: short − long = γ̂ × δ̂.

# %%
short <- lm_robust(f_log, data = d, se_type = "HC1")
long  <- lm_robust(update(f_log, . ~ . + calidad), data = d, se_type = "HC1")
aux   <- lm(update(f_log, calidad ~ .), data = d)

modelsummary(list("Short (analyst)" = short, "Long (+ calidad)" = long),
             coef_map = c("terraza", "calidad"), fmt = 4, gof_map = c("nobs", "r.squared"), output = "markdown")

gamma_hat <- coef(long)[["calidad"]]
delta_hat <- coef(aux)[["terraza"]]
c(short_minus_long = coef(short)[["terraza"]] - coef(long)[["terraza"]],
  gamma_x_delta    = gamma_hat * delta_hat,
  delta            = delta_hat,
  cor_terraza_calidad = cor(d$terraza, d$calidad))

stopifnot(all.equal(coef(short)[["terraza"]] - coef(long)[["terraza"]], gamma_hat * delta_hat))

# %% [markdown]
# **Reading.** Slides: `terraza` 0.1041 (short) and 0.0329 (long); `calidad` 0.0975; δ̂ = 0.7301; 0.1041 − 0.0329 = 0.0975 × 0.7301 = 0.0712; R² 0.900 and 0.933; correlation 0.41. Better flats have terraces (δ > 0) and quality raises prices (γ > 0), so the analyst's terrace premium is about three times the premium holding quality fixed.

# %% [markdown]
# ## C2. The number of photos as a proxy for quality (p. 49-50)
#
# A proxy is an observed variable related to the omitted one. Adding it as a control removes part of the bias.

# %%
proxy <- lm_robust(update(f_log, . ~ . + num_fotos), data = d, se_type = "HC1")
tidy(proxy) |> filter(term %in% c("terraza", "num_fotos")) |> select(term, estimate, std.error)

c(cor_fotos_calidad = cor(d$num_fotos, d$calidad),
  share_of_bias_removed = (coef(short)[["terraza"]] - coef(proxy)[["terraza"]]) /
                          (coef(short)[["terraza"]] - coef(long)[["terraza"]]))

# %% [markdown]
# **Reading.** Slides: `terraza` 0.0705 (0.0169), `num_fotos` 0.0115 (0.0008), correlation with quality 0.74. The proxy removes (0.1041 − 0.0705)/(0.1041 − 0.0329) = 47% of the bias. The coefficient on the photos has no causal meaning: the photos are there to absorb quality.

# %% [markdown]
# ## C3. A bad control: the energy certificate (p. 69)
#
# The energy certificate is partly a consequence of the age of the building. Controlling for it removes part of the effect of age that we want to measure. We compare the price difference between a flat of 50 years and one of 10, which is $40\beta_{antig} + 2400\beta_{antig^2}$ (2400 = 50² − 10²), with and without the control.

# %%
mc <- lm_robust(update(f_log, . ~ . + cert_energ), data = d, se_type = "HC1")

bind_rows("Hedonic model"  = hypotheses(short, "40 * antig + 2400 * `I(antig^2)` = 0"),
          "+ cert_energ"   = hypotheses(mc,    "40 * antig + 2400 * `I(antig^2)` = 0"),
          .id = "model") |>
  as_tibble() |> select(model, estimate, std.error, conf.low, conf.high)

# %% [markdown]
# **Reading.** With the certificate in the regression, the age effect is measured "holding the certificate fixed": the part of the effect of age that goes through energy efficiency is taken out. If the question is how much value a flat loses with age, the certificate is a bad control. (The slides give the code for this analysis and no table.)

# %% [markdown]
# ## C4. The square metres in the ad (p. 62-63)
#
# The ad declares an area (`m2_anuncio`) that measures the land-registry area (`m2`) with error. Classical measurement error in a regressor shrinks its coefficient towards zero by the factor $\sigma^2_{r^*}/(\sigma^2_{r^*} + \sigma^2_e)$, where $r^*$ is the part of the true regressor not explained by the other regressors.

# %%
m_true  <- lm(f_log, data = d)
m_noisy <- lm(update(f_log, . ~ . - lm2 + lm2_anuncio), data = d)

# Theoretical attenuation factor
r_star <- resid(lm(update(f_log, lm2 ~ . - lm2), data = d))   # log(m2) net of the other regressors
e      <- d$lm2_anuncio - d$lm2                                # measurement error (in logs)

tibble(sd_error_m2        = sd(d$m2_anuncio - d$m2),
       cor_m2_m2_anuncio  = cor(d$m2, d$m2_anuncio),
       elasticity_true    = coef(m_true)[["lm2"]],
       elasticity_noisy   = coef(m_noisy)[["lm2_anuncio"]],
       ratio              = elasticity_noisy / elasticity_true,
       attenuation_theory = var(r_star) / (var(r_star) + var(e)))

# %% [markdown]
# **Reading.** Slides: s.d. of the error 11.6 m², correlation 0.93; elasticity 0.900 with the land-registry area and 0.712 with the declared one; ratio 0.791 against a theoretical factor of 0.779. A correlation of 0.93 sounds like the same variable and still costs 21% of the elasticity.

# %% [markdown]
# ## Simulation: omitted ability and proxies
#
# From the previous course (`Proxy_variable.do`). True model: `log_earn = 1 + 0.08 educ + 0.03 exper + 0.5 ability + u`, with ability unobserved and correlated with education. Two test scores serve as proxies, one strongly and one weakly related to ability.

# %%
set.seed(12345)
n <- 10000
sim <- tibble(
  ability    = rnorm(n),
  educ       = 12 + 2 * ability + rnorm(n),
  exper      = abs(rnorm(n, 10, 5)),
  log_earn   = 1 + 0.08 * educ + 0.03 * exper + 0.5 * ability + rnorm(n),
  sat_strong = 1000 + 200 * ability + rnorm(n, 0, 50),
  sat_weak   = 1000 + 50 * ability + rnorm(n, 0, 200)
)

modelsummary(
  list("True model"   = lm(log_earn ~ educ + exper + ability, data = sim),
       "Omit ability" = lm(log_earn ~ educ + exper, data = sim),
       "Strong proxy" = lm(log_earn ~ educ + exper + sat_strong, data = sim),
       "Weak proxy"   = lm(log_earn ~ educ + exper + sat_weak, data = sim)),
  coef_map = c("educ", "exper"), fmt = 4, gof_map = "nobs", output = "markdown")

sim |> summarise(across(c(sat_strong, sat_weak, educ), \(x) cor(x, ability)))

# %% [markdown]
# **Reading.** The truth for `educ` is 0.08. Omitting ability pushes the estimate up (the bias is 0.5 × the coefficient of ability on education, about 0.2). The strong proxy removes most of the bias and the weak one very little.

# %%
# The average over many samples: the bias stays, however large the sample
set.seed(1)
b_educ <- map_dbl(1:200, \(i) {
  ability <- rnorm(n); educ <- 12 + 2 * ability + rnorm(n); exper <- abs(rnorm(n, 10, 5))
  log_earn <- 1 + 0.08 * educ + 0.03 * exper + 0.5 * ability + rnorm(n)
  coef(lm(log_earn ~ educ + exper))[["educ"]]
})
c(mean = mean(b_educ), sd = sd(b_educ))

# %% [markdown]
# ## Simulation: measurement error in y and in x
#
# From the previous course (`Measurement_error.do`). Error in the dependent variable leaves the coefficients consistent and raises the standard errors. Classical error in a regressor biases its coefficient towards zero.

# %%
set.seed(12345)
sim2 <- tibble(
  ability      = rnorm(n),
  educ         = 12 + 2 * ability + rnorm(n),
  exper        = abs(rnorm(n, 10, 5)),
  log_earn     = 1 + 0.08 * educ + 0.03 * exper + 0.5 * ability + rnorm(n),
  log_earn_err = log_earn + rnorm(n, 0, 0.3),   # error in y
  educ_err     = educ + rnorm(n, 0, 2)          # classical error in x
)

modelsummary(
  list("No error"   = lm(log_earn ~ educ + exper + ability, data = sim2),
       "Error in y" = lm(log_earn_err ~ educ + exper + ability, data = sim2),
       "Error in x" = lm(log_earn ~ educ_err + exper + ability, data = sim2)),
  coef_map = c("educ" = "educ", "educ_err" = "educ", "ability" = "ability"),
  fmt = 4, gof_map = "nobs", output = "markdown")

# %% [markdown]
# **Reading.** With error in `y` the coefficient on education stays near 0.08 with a slightly larger standard error. With error in education the coefficient shrinks towards zero, and the coefficient on ability moves too: the error in one regressor contaminates the others.

# %% [markdown]
# ## Try it yourself (p. 69)
#
# 1. Verify γ̂ × δ̂ = short − long to the fourth decimal for `m2_terraza`. Which way is the bias? Explain it with the formula.
# 2. Use both `num_fotos` and `cert_energ` as controls. Which question does that model answer?
# 3. What would measurement error in the price do?
