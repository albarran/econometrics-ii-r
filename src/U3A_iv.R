# %% [markdown]
# # Unit 3A. The simple IV estimator
#
# Econometrics II (Henry Redondo, 2026-27). R companion to **Unit 3, Block A** of the slides: the same analyses (A0 to A4) with `tidyverse`, `estimatr`, `ivreg` and `modelsummary`.
#
# Each section has the idea in two or three lines with the slide page (page of the PDF), then the code, then how to read the output. Numbers quoted as "slides" are the ones printed in the slides.

# %%
library(tidyverse)
library(broom)
library(estimatr)        # lm_robust(), iv_robust()
library(ivreg)           # ivreg()
library(marginaleffects) # hypotheses()
library(modelsummary)

options(digits = 4, readr.show_col_types = FALSE)
theme_set(theme_minimal())

# Henry's data folder. To read local copies, set the environment variable ECMTII_DATA to their folder.
data_dir <- Sys.getenv("ECMTII_DATA", "https://henryredondo.github.io/assets/files/teaching/econometrics-ii")

# %% [markdown]
# ## Warm-up: IV in a simulated world
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
lm_robust(educ ~ exper + z, data = sim, se_type = "HC1") |> tidy() |> filter(term == "z")

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
# ## The data: Programa Transparenta (slides p. 4-5)
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
    lsize       = log(size),
    share_women = n_women / size,
    invited     = as.integer(arm >= 1)
  )

# The oracle: the firm's egalitarian culture, which an analyst never observes.
oracle <- read_csv(file.path(data_dir, "transparency_oracle.csv"))
d <- left_join(d, oracle, by = "id")

glimpse(d)

# %%
# Controls used in every regression, and a helper that builds formulas with them
X <- "lsize + share_women + sector + province"
fml <- function(x) as.formula(str_glue(x, .envir = parent.frame()), env = globalenv())

deparse1(fml("dgap ~ transp + {X}"))

# %% [markdown]
# ## A0. Was the lottery fair? (p. 24)
#
# An instrument that comes from a lottery should be unrelated to anything that was fixed before the programme. We compare pre-programme means by arm and test that the 2021 gap is equal across arms.

# %%
d |>
  summarise(
    firms         = n(),
    gap_2021      = mean(gap2021),
    employees     = mean(size),
    share_women   = mean(share_women),
    parent_abroad = mean(parent_abroad),
    fair_score    = mean(fair_score),
    adopted_2025  = mean(transp),
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
# ## A1. The naive answer, and what the oracle says (p. 25)
#
# OLS of `dgap` on `transp` and the controls, with robust (HC1) standard errors. Culture is in the error term: egalitarian firms adopt more often and close their gap faster.

# %%
ols        <- lm_robust(fml("dgap ~ transp + {X}"), data = d, se_type = "HC1")
ols_oracle <- lm_robust(fml("dgap ~ transp + culture + {X}"), data = d, se_type = "HC1")

modelsummary(list("OLS" = ols, "OLS + culture" = ols_oracle),
             coef_map = c("transp", "culture"), fmt = 4, gof_map = "nobs", output = "markdown")

d |> summarise(mean_culture = mean(culture), .by = transp)

# %% [markdown]
# **Reading.** Slides: OLS gives −0.0561 (s.e. 0.0019). With culture in the regression the coefficient falls to −0.0310, so almost half of the OLS number was culture. Adopters have mean culture +0.54 and non-adopters −0.37. The target for the rest of the unit is −0.031, without using `culture`.

# %% [markdown]
# ## A2. The IV estimator, by hand (p. 28)
#
# Instrument: `invited` (either arm). With one endogenous regressor and one instrument, the IV estimator is the **Wald ratio**: reduced form divided by first stage.

# %%
first_stage  <- lm_robust(fml("transp ~ invited + {X}"), data = d, se_type = "HC1")
reduced_form <- lm_robust(fml("dgap ~ invited + {X}"), data = d, se_type = "HC1")
iv           <- ivreg(fml("dgap ~ transp + {X} | invited + {X}"), data = d)

wald <- coef(reduced_form)[["invited"]] / coef(first_stage)[["invited"]]
c(first_stage = coef(first_stage)[["invited"]], reduced_form = coef(reduced_form)[["invited"]],
  wald = wald, ivreg = coef(iv)[["transp"]])

# The ratio is the ivreg estimate, exactly
stopifnot(all.equal(wald, coef(iv)[["transp"]]))

# %%
modelsummary(
  list("First stage (transp)" = first_stage, "Reduced form (dgap)" = reduced_form, "IV (dgap)" = iv),
  vcov = list(NULL, NULL, "HC0"),   # the lm_robust models keep their HC1; IV with HC0 as in the slides
  coef_map = c("invited", "transp"), fmt = 4, gof_map = "nobs", output = "markdown"
)

# %% [markdown]
# The same in one line with `estimatr::iv_robust`, which computes the robust standard errors itself. Its default is HC2; the slides use HC0 for IV.

# %%
iv_r <- iv_robust(fml("dgap ~ transp + {X} | invited + {X}"), data = d, se_type = "HC0", diagnostics = TRUE)
tidy(iv_r) |> filter(term == "transp")

# Same coefficients and same standard errors as ivreg + sandwich HC0
stopifnot(all.equal(coef(iv_r), coef(iv)),
          all.equal(iv_r$std.error, sqrt(diag(sandwich::vcovHC(iv, type = "HC0")))))

# First-stage F: classical (ivreg) and robust (iv_robust)
summary(iv, diagnostics = TRUE)$diagnostics["Weak instruments", ]
iv_r$diagnostic_first_stage_fstatistic

# %% [markdown]
# **Reading.** Slides: first stage 0.1395 (s.e. 0.0129), reduced form −0.0052 (0.0020), IV −0.0375 (0.0132), F = 118. Invited firms were 14 pp more likely to adopt and their gap fell 0.5 pp more (the intention-to-treat effect). If the invitation works only through adoption, the effect of adopting is 0.52/0.14 ≈ 3.7 pp. IV rescales the effect of the lottery by how much the lottery moved the regressor.

# %% [markdown]
# ## A3. What bad instruments do (p. 30)
#
# Two other candidates are strongly related to adoption: `parent_abroad` and `peer_share`. Both are also related to culture, so they fail exogeneity. We run the same IV with each instrument.

# %%
iv_with <- function(z) {
  m <- iv_robust(fml("dgap ~ transp + {X} | {z} + {X}"), data = d, se_type = "HC0", diagnostics = TRUE)
  tibble(instrument    = z,
         first_stage_F = m$diagnostic_first_stage_fstatistic[["value"]],
         beta          = coef(m)[["transp"]],
         se            = m$std.error[["transp"]],
         cor_culture   = cor(d[[z]], d$culture))
}

map(c("invited", "parent_abroad", "peer_share"), iv_with) |>
  list_rbind() |>
  add_row(instrument = "OLS", beta = coef(ols)[["transp"]], se = ols$std.error[["transp"]])

# %% [markdown]
# **Reading.** Slides: F = 118, 72 and 93, so all three are relevant. The estimates are −0.038 (invited), −0.111 (parent_abroad) and −0.188 (peer_share), against −0.056 for OLS and −0.031 for the oracle. The two invalid instruments land further from the truth than OLS, with small standard errors. Their correlation with culture is 0.14 and 0.42 (0.001 for `invited`). A strong first stage says nothing about validity; here we can see it only because we have the oracle.

# %% [markdown]
# ## A4. The price of identification (p. 40)
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
# ## Try it yourself (p. 42)
#
# 1. Divide the reduced-form coefficient by the first-stage one with more decimals. Do you get the `ivreg` estimate exactly?
# 2. Drop the controls from the three regressions (`fml("dgap ~ transp | invited")`). How much does the IV estimate move? Why so little?
# 3. Compare the IV standard error from `vcov(iv)` (classical), `sandwich::vcovHC(iv, type = "HC0")` and `iv_robust(..., se_type = "HC2")`.
