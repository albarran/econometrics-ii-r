# %% [markdown]
# # Unit 2. The linear regression model
#
# Econometrics II (Henry Redondo, 2026-27). R companion to **Unit 2** of the slides, one section per block:
#
# - **Block 2A.** OLS mechanics and goodness of fit (A1 to A4)
# - **Block 2B.** Interpreting coefficients (B1 to B9)
# - **Block 2C.** Robust standard errors, omitted variables, proxies and measurement error (D1, C1 to C4)
# - **Blocks 2D and 2E.** Inference, heteroskedasticity and weighted least squares (D2 to D8, E1)
#
# Each analysis has the idea in two or three lines with the slide page (page of the PDF of its block), then the code, then how to read the output. Numbers quoted as "slides" are the ones printed in the slides.
#
# Run the notebook from the top: each block uses objects created in the previous ones.

# %%
library(tidyverse)
library(broom)           # tidy(), glance(), augment()
library(estimatr)        # lm_robust()
library(marginaleffects) # hypotheses()
library(lmtest)          # bptest()
library(modelsummary)
library(fixest)          # feols(), one optional cell

options(digits = 4, readr.show_col_types = FALSE)
theme_set(theme_minimal())

# Henry's data folder. To read local copies, set the environment variable ECMTII_DATA to their folder.
data_dir <- Sys.getenv("ECMTII_DATA", "https://henryredondo.github.io/assets/files/teaching/econometrics-ii")

# %% [markdown]
# ## The data: 1,000 listings in Alicante, 2013-2026 (2A p. 3-4)
#
# Simulated listings: asking price, floor area, terrace, age, ground floor, zone and date. Because they are simulated, the truth is known and each estimate can be checked against it.
#
# > **Data pending.** On 25 November 2026 the CSV files were not on the course website. The cell reads them from the address where they should be. The slides use a beach dummy `costa` that is not in their variable list: if the file does not bring it, it is built as "zone is Playa de San Juan" (assumption to confirm).

# %%
d <- read_csv(file.path(data_dir, "idealista_alicante_clase.csv")) |>
  mutate(zona        = fct_relevel(factor(zona), "Benalúa"),   # base group
         anio        = factor(anio),
         lprecio     = log(precio),
         lm2         = log(m2),
         lm2_anuncio = log(m2_anuncio))
if (!"costa" %in% names(d)) d <- d |> mutate(costa = as.integer(zona == "Playa de San Juan"))

# The two hedonic models used in the whole unit
f_niv <- precio  ~ m2  + m2_terraza + terraza + antig + I(antig^2) + bajo + zona + anio
f_log <- lprecio ~ lm2 + terraza + m2_terraza + antig + I(antig^2) + bajo + zona + anio

d |> summarise(n = n(), median_price = median(precio), mean_price = mean(precio), mean_m2 = mean(m2),
               share_terrace = mean(terraza), mean_m2_terrace = mean(m2_terraza[terraza == 1]))
d |> summarise(mean_price = mean(precio), .by = zona)

# %% [markdown]
# **Reading.** Slides: median price 131,750 €, mean 158,666 €; mean area 91.5 m²; 43.5% have a terrace (13.4 m² on average); mean price from 108,589 € in San Blas to 295,367 € in Playa de San Juan.

# %% [markdown]
# ## Block 2A. OLS: mechanics and goodness of fit
#
# Page numbers in this block refer to `Unit2A_slides.pdf`.

# %% [markdown]
# ### A1. Euros per square metre (p. 19)
#
# The simple regression of price on floor area. `tidy()` returns the coefficient table and `glance()` the fit statistics, both as data frames.

# %%
m1 <- lm(precio ~ m2, data = d)
tidy(m1)
glance(m1) |> select(r.squared, nobs)

# %%
ggplot(d, aes(m2, precio)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  labs(x = "Floor area (m2)", y = "Asking price (€)")

# %% [markdown]
# **Reading.** Slides: slope 2,206.21 (s.e. 74.97), intercept −43,151.50, R² = 0.465. One more square metre goes with 2,206 € more on the asking price. The intercept only positions the line: nobody lists a flat of 0 m². Floor area alone accounts for 46% of the variation in prices.

# %% [markdown]
# #### The OLS formula, once (p. 9-12)
#
# $\hat\beta = (X'X)^{-1}X'y$. This is what `lm` computes (with a more stable algorithm).

# %%
X <- cbind(1, d$m2)
y <- d$precio
beta_hat <- solve(t(X) %*% X, t(X) %*% y)

stopifnot(all.equal(drop(beta_hat), coef(m1), check.attributes = FALSE))
drop(beta_hat)

# %% [markdown]
# ### A2. The first-order conditions hold in the data (p. 20)
#
# OLS chooses the coefficients so that the residuals sum to zero and are orthogonal to every regressor. `augment()` adds fitted values (`.fitted`) and residuals (`.resid`) to the data.

# %%
augment(m1) |>
  summarise(sum_resid     = sum(.resid),
            sum_m2_resid  = sum(m2 * .resid),
            cor_resid_m2  = cor(.resid, m2),
            mean_fitted   = mean(.fitted),
            mean_y        = mean(precio))

# %% [markdown]
# **Reading.** Slides: 2.9 × 10⁻⁸, 2.2 × 10⁻⁶ and a correlation of about 10⁻¹⁶: zeros up to rounding error; the mean of the fitted values equals the mean of the price (158,666.30). These equalities hold by construction in every sample, true model or not. So "the residuals are uncorrelated with `m2`" is an algebraic fact and says nothing about whether `m2` is exogenous.

# %% [markdown]
# ### A3. SST = SSE + SSR with the full hedonic model (p. 34)
#
# With a constant in the model, the total sum of squares splits into the explained part (SSE) and the residual part (SSR), and R² = SSE/SST.

# %%
m3 <- lm(f_niv, data = d)

ss <- augment(m3) |>
  summarise(SSE = sum((.fitted - mean(precio))^2),
            SSR = sum(.resid^2),
            SST = sum((precio - mean(precio))^2)) |>
  mutate(R2 = SSE / SST)
ss

stopifnot(all.equal(ss$SSE + ss$SSR, ss$SST),
          all.equal(ss$R2, glance(m3)$r.squared))

glance(m3) |> transmute(r.squared, adj.r.squared, parameters = df + 1)

# %% [markdown]
# **Reading.** Slides: SSE = 7.6798 × 10¹², SSR = 1.5507 × 10¹², SST = 9.2305 × 10¹², R² = 0.8320, adjusted R² = 0.8280, k = 24. Adding terrace, age, floor, zone and year raises R² from 0.46 to 0.83. R² cannot tell whether the terrace coefficient is right; that is Block C.

# %% [markdown]
# ### A4. More parameters, worse predictions (p. 37-38)
#
# R² measures the fit in the estimation sample. To see how the model predicts flats it has not seen: estimate on the odd ids, predict the even ids, and compute the same ratio there.

# %%
train <- filter(d, id %% 2 == 1)
test  <- filter(d, id %% 2 == 0)

# The hedonic model, and a much larger one: zone x year, and area slopes by zone and by year
f_big <- update(f_log, . ~ . + zona:anio + lm2:zona + lm2:anio)

r2_out <- function(model, newdata) {
  pred <- suppressWarnings(predict(model, newdata = newdata))   # the big model can be rank deficient
  1 - sum((newdata$lprecio - pred)^2) / sum((newdata$lprecio - mean(newdata$lprecio))^2)
}

map(list("Hedonic model" = f_log, "+ zone x year, slopes by zone and year" = f_big), \(f) {
  m <- lm(f, data = train)
  tibble(k = length(coef(m)), R2 = glance(m)$r.squared, adj_R2 = glance(m)$adj.r.squared,
         R2_out_of_sample = r2_out(m, test))
}) |> list_rbind(names_to = "model")

# %% [markdown]
# **Reading.** Slides: with 24 parameters R² = 0.904, adjusted 0.899, out of sample 0.890; with 93 parameters R² = 0.917, adjusted 0.899, out of sample 0.865. Sixty-nine more parameters raise R² by 1.3 points and lower it by 2.5 points on new flats. Adjusted R², AIC, BIC and cross-validation are versions of one idea: penalise complexity.

# %% [markdown]
# ### Try it yourself (p. 39)
#
# 1. Regress `precio` on `m2` without a constant, `lm(precio ~ 0 + m2, data = d)`. Do the residuals still sum to zero? Why?
# 2. Add `factor(mes)` to the model. What happens to R² in and out of sample?
# 3. Add a column of pure noise, `mutate(d, ruido = rnorm(n()))`, to `lm(precio ~ m2)`. Does R² rise? Does adjusted R²?

# %% [markdown]
# ## Block 2B. Interpreting coefficients
#
# Page numbers in this block refer to `Unit2B_slides.pdf`. From here on, robust (HC1) standard errors are the default: `lm_robust(se_type = "HC1")`. The plain `lm` versions are kept for the classical formulas and tests.

# %%
mn <- lm_robust(f_niv, data = d, se_type = "HC1"); mn_c <- lm(f_niv, data = d)   # levels
ml <- lm_robust(f_log, data = d, se_type = "HC1"); ml_c <- lm(f_log, data = d)   # logs

# %% [markdown]
# ### B1. Coefficients in euros and in percentages (p. 7, 16)
#
# In levels a coefficient is euros per unit of the regressor. With log price on the left: a log regressor gives an elasticity, and a regressor in levels gives (approximately) 100·β per cent per unit.

# %%
modelsummary(list("Levels (€)" = mn, "Logs" = ml),
             coef_map = c("m2" = "m2 / log(m2)", "lm2" = "m2 / log(m2)", "m2_terraza" = "m2_terraza"),
             fmt = fmt_significant(4), gof_map = c("nobs", "r.squared"), output = "markdown")

# %% [markdown]
# **Reading.** Slides, levels: `m2` 1,615.76 (robust s.e. 71.93), `m2_terraza` 1,475.43 (431.08). Logs: `lm2` 0.8997 (0.0186), `m2_terraza` 0.0053 (0.0014). One more interior square metre goes with 1,616 € more, holding the rest fixed. In logs, 1% more area goes with 0.90% more price, and one more square metre of terrace with about 0.53% more.

# %% [markdown]
# ### B2. Does an old building lose value forever? (p. 13)
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
# ### B3. The ground-floor discount (p. 27)

# %%
modelsummary(list("Levels (€)" = mn, "Logs" = ml), coef_map = "bajo",
             fmt = fmt_significant(4), gof_map = "nobs", output = "markdown")

# Exact percentage effect of a dummy in a log model: 100 * (exp(beta) - 1)
hypotheses(ml, "100 * (exp(bajo) - 1) = 0")

# %% [markdown]
# **Reading.** Slides: −22,357.73 € (s.e. 3,317.22) in levels and −0.1210 (0.0157) in logs. A ground-floor flat is about 22,000 € cheaper; in per cent, 100·(e^−0.121 − 1) = −11.4%. With a dummy `alto` = 1 − `bajo` the coefficient would change sign and the intercept would absorb the difference.

# %% [markdown]
# ### B4. Five zones, base group Benalúa (p. 31)
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
# ### B5. Changing the base group (p. 35)

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
# ### B6. The terrace premium (p. 39, 42)
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
# ### B7. Is a terrace worth more at the beach? (p. 46)
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
# ### B8. Does the area elasticity differ at the beach? (p. 50)

# %%
m8 <- lm_robust(update(f_log, . ~ . + lm2:costa), data = d, se_type = "HC1")
tidy(m8) |> filter(term %in% c("lm2", "lm2:costa")) |> select(term, estimate, std.error, statistic)

# %% [markdown]
# **Reading.** Slides: elasticity 0.9081 (0.0191) away from the beach, interaction −0.0438 (0.0563), t = −0.78. No evidence of a different elasticity at the beach.

# %% [markdown]
# ### B9. Year dummies are a house-price index (p. 52)
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
# ### The same in fixest
#
# `fixest::feols` with `vcov = "hetero"` gives the same coefficients and robust standard errors; `i(zona, ref = ...)` sets the base group inside the formula.

# %%
ml_fe <- feols(lprecio ~ lm2 + terraza + m2_terraza + antig + I(antig^2) + bajo +
                 i(zona, ref = "Benalúa") + anio, data = d, vcov = "hetero")

stopifnot(all.equal(coef(ml_fe)[["terraza"]], coef(ml)[["terraza"]]))
c(lm_robust_HC1 = ml$std.error[["terraza"]], feols_hetero = se(ml_fe)[["terraza"]])

# %% [markdown]
# ### Try it yourself (p. 54)
#
# 1. Use `fct_relevel(zona, "Playa de San Juan")`. Which coefficients change? Does R²?
# 2. Add `terraza:bajo`. Is a ground-floor terrace worth less?
# 3. Build the index by month with `factor(fecha)`. Is it too noisy?
# 4. `lm_robust` uses HC2 by default and `sandwich::vcovHC` uses HC3; the slides use HC1. Compare the three for `terraza`.

# %% [markdown]
# ## Block 2C. Robust standard errors, omitted variables, proxies and measurement error
#
# Page numbers in this block refer to `Unit2C_slides.pdf`. It closes with two simulations from the previous course (proxy variables and measurement error).

# %% [markdown]
# ### D1. Classical against robust standard errors (p. 17-18)
#
# The classical formula assumes that the error variance is the same for every flat. The robust (White, here HC1) formula does not. `modelsummary` can show the same model with both.

# %%
modelsummary(list("Levels, classical" = mn_c, "Levels, robust" = mn_c, "Logs, classical" = ml_c, "Logs, robust" = ml_c),
             vcov = c("classical", "HC1", "classical", "HC1"),
             coef_map = c("m2" = "m2 / log(m2)", "lm2" = "m2 / log(m2)", "m2_terraza" = "m2_terraza",
                          "terraza" = "terraza", "bajo" = "bajo"),
             fmt = fmt_significant(4), gof_map = "nobs", output = "markdown")

# %% [markdown]
# #### The sandwich, once (p. 15)
#
# $\widehat{Var}(\hat\beta) = (X'X)^{-1}\left(\sum_i \hat u_i^2 x_i x_i'\right)(X'X)^{-1}$, times N/(N−k) for HC1. `lm_robust(se_type = "HC1")` does this in one step.

# %%
X  <- model.matrix(ml_c)
u  <- resid(ml_c)
N  <- nrow(X); k <- ncol(X)
bread <- solve(crossprod(X))
meat  <- crossprod(X * u)                       # sum of u_i^2 x_i x_i'
V_hc1 <- N / (N - k) * bread %*% meat %*% bread

stopifnot(all.equal(sqrt(diag(V_hc1)), ml$std.error, check.attributes = FALSE))

c(by_hand = sqrt(diag(V_hc1))[["terraza"]], lm_robust = ml$std.error[["terraza"]],
  classical = tidy(ml_c)$std.error[tidy(ml_c)$term == "terraza"])

# %% [markdown]
# **Reading.** Slides, ratio robust/classical: in levels `m2` 1.56 (46.12 → 71.93), `m2_terraza` 1.60, `terraza` 1.24, `bajo` 0.84; in logs 1.05, 1.26, 1.11 and 0.99. In euros a 400,000 € beach flat misses by far more euros than a 70,000 € flat: strong heteroskedasticity, and the classical s.e. of `m2` is 36% too small. In logs the two are much closer. Robust errors can also be smaller (`bajo`). Rule: report robust standard errors by default.

# %% [markdown]
# ### The oracle
#
# The unobserved quality of each flat (`calidad`) is kept in a separate file. An analyst never has it; we open it to see what we are missing.

# %%
d <- d |> left_join(read_csv(file.path(data_dir, "idealista_oraculo_clase.csv")), by = "id")

# %% [markdown]
# ### C1. Opening the oracle on omitted quality (p. 43-44)
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
# ### C2. The number of photos as a proxy for quality (p. 49-50)
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
# ### C3. A bad control: the energy certificate (p. 69)
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
# ### C4. The square metres in the ad (p. 62-63)
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
# ### Simulation: omitted ability and proxies
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
# ### Simulation: measurement error in y and in x
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
# ### Try it yourself (p. 69)
#
# 1. Verify γ̂ × δ̂ = short − long to the fourth decimal for `m2_terraza`. Which way is the bias? Explain it with the formula.
# 2. Use both `num_fotos` and `cert_energ` as controls. Which question does that model answer?
# 3. What would measurement error in the price do?

# %% [markdown]
# ## Blocks 2D and 2E. Inference, heteroskedasticity and weighted least squares
#
# Page numbers in this block refer to `Unit2D_slides.pdf` and `Unit2E_slides.pdf`, as marked in each heading.

# %% [markdown]
# ### D2. Is a terrace worth anything? (2D p. 25-27)
#
# t test and 95% confidence interval for the terrace coefficient in the log model, with robust standard errors.

# %%
tidy(ml) |> filter(term == "terraza")

# The interval by hand: estimate +- t(0.975, N - k) x s.e.
with(tidy(ml) |> filter(term == "terraza"), estimate + c(-1, 1) * qt(0.975, df) * std.error)

# %% [markdown]
# **Reading.** Slides: 0.1041, s.e. 0.0189, t = 5.50, 95% CI [0.0669, 0.1413] with N − k = 976. The hypothesis of no premium is rejected at any conventional level and the interval is narrow. From Block C: holding quality fixed the premium is 0.033, outside this interval. A confidence interval measures sampling noise only; a precise estimate can be precisely wrong.

# %% [markdown]
# ### Joint tests: the four zone coefficients (2D p. 43)
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
# #### Individual t tests against a joint F test
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
# ### D4. The appraisers' rule (2D p. 30-31)
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
# ### D6. The LM test for the zones (2D p. 42-43)
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
# ### D7. Is there heteroskedasticity? (2E p. 4-5)
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
# ### D8. Weighted least squares with aggregated data (2E p. 12-13)
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
# ### E1. How much is a terrace worth? Every version (2E p. 19-20)
#
# The premium of a typical terrace (mean size among flats with one), in per cent, from each specification of the unit: $100\,(e^{\hat\beta_{terraza} + \hat\beta_{m2\_terraza}\,\bar s} - 1)$. The first row has no terrace size in the model.

# %%
# m2t: mean terrace size among flats with a terrace (Block 2B)
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
# ### Try it yourself (2D p. 47, 2E p. 22)
#
# 1. Test the 50% rule in the log model. What exactly is the hypothesis there?
# 2. Run `bptest(ml_c)`, with all the regressors in place of the fitted values. Same verdict? How many degrees of freedom?
# 3. Feasible GLS: regress `log(resid(ml_c)^2)` on the regressors, set `w = exp(fitted)` and refit with `weights = 1/w`. Does the terrace coefficient move? Do the standard errors fall?
# 4. In D8, why would robust standard errors on the unweighted cell regression not make WLS unnecessary?
