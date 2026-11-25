# %% [markdown]
# # Unit 2B. Interpreting coefficients
#
# Econometrics II (Henry Redondo, 2026-27). R companion to **Unit 2, Block B** of the slides: analyses B1 to B9 (levels and logs, quadratics, dummies, interactions, the hedonic price index).
#
# Each section has the idea with the slide page (page of the PDF), then the code, then how to read the output. Numbers quoted as "slides" are the ones printed in the slides.

# %%
library(tidyverse)
library(broom)
library(estimatr)        # lm_robust()
library(marginaleffects) # hypotheses()
library(modelsummary)
library(fixest)          # feols(), one optional cell

options(digits = 4, readr.show_col_types = FALSE)
theme_set(theme_minimal())

# Henry's data folder. To read local copies, set the environment variable ECMTII_DATA to their folder.
data_dir <- Sys.getenv("ECMTII_DATA", "https://henryredondo.github.io/assets/files/teaching/econometrics-ii")

# %% [markdown]
# ## The data and the two hedonic models
#
# > **Data pending.** On 25 November 2026 the CSV files were not on the course website. The cell reads them from the address where they should be. The slides use a beach dummy `costa` that is not in their variable list: if the file does not bring it, it is built as "zone is Playa de San Juan" (assumption to confirm).

# %%
d <- read_csv(file.path(data_dir, "idealista_alicante_clase.csv")) |>
  mutate(zona    = fct_relevel(factor(zona), "Benalúa"),   # base group
         anio    = factor(anio),
         lprecio = log(precio),
         lm2     = log(m2))
if (!"costa" %in% names(d)) d <- d |> mutate(costa = as.integer(zona == "Playa de San Juan"))

f_niv <- precio  ~ m2  + m2_terraza + terraza + antig + I(antig^2) + bajo + zona + anio
f_log <- lprecio ~ lm2 + terraza + m2_terraza + antig + I(antig^2) + bajo + zona + anio

# Robust (HC1) standard errors are the default in the whole unit
mn <- lm_robust(f_niv, data = d, se_type = "HC1")   # levels
ml <- lm_robust(f_log, data = d, se_type = "HC1")   # logs

# %% [markdown]
# ## B1. Coefficients in euros and in percentages (p. 7, 16)
#
# In levels a coefficient is euros per unit of the regressor. With log price on the left: a log regressor gives an elasticity, and a regressor in levels gives (approximately) 100·β per cent per unit.

# %%
modelsummary(list("Levels (€)" = mn, "Logs" = ml),
             coef_map = c("m2" = "m2 / log(m2)", "lm2" = "m2 / log(m2)", "m2_terraza" = "m2_terraza"),
             fmt = fmt_significant(4), gof_map = c("nobs", "r.squared"), output = "markdown")

# %% [markdown]
# **Reading.** Slides, levels: `m2` 1,615.76 (robust s.e. 71.93), `m2_terraza` 1,475.43 (431.08). Logs: `lm2` 0.8997 (0.0186), `m2_terraza` 0.0053 (0.0014). One more interior square metre goes with 1,616 € more, holding the rest fixed. In logs, 1% more area goes with 0.90% more price, and one more square metre of terrace with about 0.53% more.

# %% [markdown]
# ## B2. Does an old building lose value forever? (p. 13)
#
# With `antig` and its square, the effect of age is $\beta_1 + 2\beta_2\,antig$ and the profile turns at $-\beta_1/(2\beta_2)$. That is a non-linear function of the coefficients: `hypotheses()` computes its standard error by the delta method.

# %%
tidy(ml) |> filter(str_detect(term, "antig")) |> select(term, estimate, std.error)

hypotheses(ml, "-antig / (2 * `I(antig^2)`) = 0")

# %%
b <- coef(ml)
ggplot() +
  geom_function(fun = \(a) 100 * (exp(b[["antig"]] * a + b[["I(antig^2)"]] * a^2) - 1)) +
  xlim(0, 90) +
  labs(x = "Age when listed (years)", y = "Price vs. a new flat (%)")

# %% [markdown]
# **Reading.** Slides: coefficients −0.0142 and 0.0001297; minimum at 54.9 years. A flat loses value with age up to about 55 years and recovers afterwards (old buildings in the centre, renovated). The test of `= 0` printed by `hypotheses` is of no interest here; the estimate and its confidence interval are.

# %% [markdown]
# ## B3. The ground-floor discount (p. 27)

# %%
modelsummary(list("Levels (€)" = mn, "Logs" = ml), coef_map = "bajo",
             fmt = fmt_significant(4), gof_map = "nobs", output = "markdown")

# Exact percentage effect of a dummy in a log model: 100 * (exp(beta) - 1)
hypotheses(ml, "100 * (exp(bajo) - 1) = 0")

# %% [markdown]
# **Reading.** Slides: −22,357.73 € (s.e. 3,317.22) in levels and −0.1210 (0.0157) in logs. A ground-floor flat is about 22,000 € cheaper; in per cent, 100·(e^−0.121 − 1) = −11.4%. With a dummy `alto` = 1 − `bajo` the coefficient would change sign and the intercept would absorb the difference.

# %% [markdown]
# ## B4. Five zones, base group Benalúa (p. 31)
#
# A factor with five levels enters as four dummies; each coefficient compares its zone with the base group.

# %%
tidy(ml) |>
  filter(str_detect(term, "^zona")) |>
  transmute(zone = str_remove(term, "^zona"), estimate, std.error,
            exact_effect_pct = 100 * (exp(estimate) - 1))

# %% [markdown]
# **Reading.** Slides: Centro 0.0559 (+5.8%), Carolinas −0.1193 (−11.2%), San Blas −0.1747 (−16.0%), Playa de San Juan 0.3172 (+37.3%). Beach flats are listed at more than double the price of Benalúa on average, and the estimated premium is 37%: beach flats are also larger, newer and have more terraces, and the regression holds those fixed.

# %% [markdown]
# ## B5. Changing the base group (p. 35)

# %%
ml_centro <- lm_robust(f_log, data = mutate(d, zona = fct_relevel(zona, "Centro")), se_type = "HC1")
tidy(ml_centro) |> filter(str_detect(term, "^zona")) |> select(term, estimate, std.error)

# Same model, different parametrisation: the new coefficient is a difference of the old ones
stopifnot(all.equal(coef(ml_centro)[["zonaPlaya de San Juan"]],
                    coef(ml)[["zonaPlaya de San Juan"]] - coef(ml)[["zonaCentro"]]),
          all.equal(ml_centro$r.squared, ml$r.squared))

# %% [markdown]
# **Reading.** Slides: with Centro as base, the coefficient of Playa de San Juan is 0.2612 = 0.3172 − 0.0559. The fit and all the other coefficients are unchanged.

# %% [markdown]
# ## B6. The terrace premium (p. 39, 42)
#
# The model has a dummy `terraza` and the terrace size `m2_terraza` (zero without a terrace). The premium of having a terrace of size s is $\beta_{terraza} + \beta_{m2\_terraza}\,s$.

# %%
modelsummary(list("Levels (€)" = mn, "Logs" = ml), coef_map = c("terraza", "m2_terraza"),
             fmt = fmt_significant(4), gof_map = "nobs", output = "markdown")

# Premium of a typical terrace (mean size among flats with a terrace), in % and in euros
m2t <- round(mean(d$m2_terraza[d$terraza == 1]), 1)
hypotheses(ml, str_glue("100 * (exp(terraza + m2_terraza * {m2t}) - 1) = 0"))
hypotheses(mn, str_glue("terraza + m2_terraza * {m2t} = 0"))

# Two dummies, no interaction: terrace and ground floor together
hypotheses(ml, "terraza + bajo = 0")

# %% [markdown]
# **Reading.** Slides: `terraza` 9,008.58 (5,273.38) in levels and 0.1041 (0.0189) in logs; `m2_terraza` 1,475.43 and 0.0053. A typical terrace of 13.4 m² is worth 9,009 + 1,475 × 13.4 ≈ 28,800 € in the levels model and 19.2% in the log model. Without an interaction, a ground-floor flat with a terrace differs from the reference by 0.1041 − 0.1210 = −0.0169.

# %% [markdown]
# ## B7. Is a terrace worth more at the beach? (p. 46)
#
# An interaction between two dummies lets the terrace premium differ at the beach.

# %%
m7 <- lm_robust(update(f_log, . ~ . + terraza:costa), data = d, se_type = "HC1")
tidy(m7) |> filter(term %in% c("terraza", "terraza:costa")) |> select(term, estimate, std.error, statistic)

# Terrace coefficient at the beach
hypotheses(m7, "terraza + `terraza:costa` = 0")

# %% [markdown]
# **Reading.** Slides: `terraza` 0.1026 (0.0188), interaction 0.0502 (0.0418). The point estimate says 5 points more at the beach, with a standard error that does not allow us to distinguish it from zero.

# %% [markdown]
# ## B8. Does the area elasticity differ at the beach? (p. 50)

# %%
m8 <- lm_robust(update(f_log, . ~ . + lm2:costa), data = d, se_type = "HC1")
tidy(m8) |> filter(term %in% c("lm2", "lm2:costa")) |> select(term, estimate, std.error, statistic)

# %% [markdown]
# **Reading.** Slides: elasticity 0.9081 (0.0191) away from the beach, interaction −0.0438 (0.0563), t = −0.78. No evidence of a different elasticity at the beach.

# %% [markdown]
# ## B9. Year dummies are a house-price index (p. 52)
#
# Each year coefficient is the price change of an identical flat relative to the base year, 2013. The index is $100\,e^{\hat\gamma_t}$.

# %%
index <- tidy(ml) |>
  filter(str_detect(term, "^anio")) |>
  transmute(year = parse_number(term), hedonic = 100 * exp(estimate)) |>
  add_row(year = 2013, hedonic = 100, .before = 1)

raw <- d |>
  summarise(eur_m2 = mean(precio / m2), .by = anio) |>
  mutate(year = parse_number(as.character(anio)), raw_mean = 100 * eur_m2 / eur_m2[year == 2013])

left_join(index, select(raw, year, raw_mean), by = "year") |>
  pivot_longer(-year, names_to = "index") |>
  ggplot(aes(year, value, colour = index)) +
  geom_line() + geom_point() +
  labs(x = "Listing year", y = "Index, 2013 = 100", colour = NULL)

# %% [markdown]
# **Reading.** The hedonic index follows the price of a flat of constant characteristics. The raw mean of €/m² also moves when the mix of flats listed each year changes.

# %% [markdown]
# ## The same in fixest
#
# `fixest::feols` with `vcov = "hetero"` gives the same coefficients and robust standard errors; `i(zona, ref = ...)` sets the base group inside the formula.

# %%
ml_fe <- feols(lprecio ~ lm2 + terraza + m2_terraza + antig + I(antig^2) + bajo +
                 i(zona, ref = "Benalúa") + anio, data = d, vcov = "hetero")

stopifnot(all.equal(coef(ml_fe)[["terraza"]], coef(ml)[["terraza"]]))
c(lm_robust_HC1 = ml$std.error[["terraza"]], feols_hetero = se(ml_fe)[["terraza"]])

# %% [markdown]
# ## Try it yourself (p. 54)
#
# 1. Use `fct_relevel(zona, "Playa de San Juan")`. Which coefficients change? Does R²?
# 2. Add `terraza:bajo`. Is a ground-floor terrace worth less?
# 3. Build the index by month with `factor(fecha)`. Is it too noisy?
# 4. `lm_robust` uses HC2 by default and `sandwich::vcovHC` uses HC3; the slides use HC1. Compare the three for `terraza`.
