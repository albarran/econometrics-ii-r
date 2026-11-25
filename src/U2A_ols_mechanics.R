# %% [markdown]
# # Unit 2A. OLS: mechanics and goodness of fit
#
# Econometrics II (Henry Redondo, 2026-27). R companion to **Unit 2, Block A** of the slides: analyses A1 to A4.
#
# Each section has the idea with the slide page (page of the PDF), then the code, then how to read the output. Numbers quoted as "slides" are the ones printed in the slides.

# %%
library(tidyverse)
library(broom)   # tidy(), glance(), augment()

options(digits = 4, readr.show_col_types = FALSE)
theme_set(theme_minimal())

# Henry's data folder. To read local copies, set the environment variable ECMTII_DATA to their folder.
data_dir <- Sys.getenv("ECMTII_DATA", "https://henryredondo.github.io/assets/files/teaching/econometrics-ii")

# %% [markdown]
# ## The data: 1,000 listings in Alicante, 2013-2026 (p. 3-4)
#
# Simulated listings: asking price, floor area, terrace, age, ground floor, zone and date. Because they are simulated, the truth is known and each estimate can be checked against it.
#
# > **Data pending.** On 25 November 2026 the CSV files were not on the course website. The cell reads them from the address where they should be.

# %%
d <- read_csv(file.path(data_dir, "idealista_alicante_clase.csv")) |>
  mutate(zona    = fct_relevel(factor(zona), "Benalúa"),   # base group
         anio    = factor(anio),
         lprecio = log(precio),
         lm2     = log(m2))

# The two hedonic models used in the whole unit
f_niv <- precio  ~ m2  + m2_terraza + terraza + antig + I(antig^2) + bajo + zona + anio
f_log <- lprecio ~ lm2 + terraza + m2_terraza + antig + I(antig^2) + bajo + zona + anio

d |> summarise(n = n(), median_price = median(precio), mean_price = mean(precio), mean_m2 = mean(m2),
               share_terrace = mean(terraza), mean_m2_terrace = mean(m2_terraza[terraza == 1]))
d |> summarise(mean_price = mean(precio), .by = zona)

# %% [markdown]
# **Reading.** Slides: median price 131,750 €, mean 158,666 €; mean area 91.5 m²; 43.5% have a terrace (13.4 m² on average); mean price from 108,589 € in San Blas to 295,367 € in Playa de San Juan.

# %% [markdown]
# ## A1. Euros per square metre (p. 19)
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
# ### The OLS formula, once (p. 9-12)
#
# $\hat\beta = (X'X)^{-1}X'y$. This is what `lm` computes (with a more stable algorithm).

# %%
X <- cbind(1, d$m2)
y <- d$precio
beta_hat <- solve(t(X) %*% X, t(X) %*% y)

stopifnot(all.equal(drop(beta_hat), coef(m1), check.attributes = FALSE))
drop(beta_hat)

# %% [markdown]
# ## A2. The first-order conditions hold in the data (p. 20)
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
# ## A3. SST = SSE + SSR with the full hedonic model (p. 34)
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
# ## A4. More parameters, worse predictions (p. 37-38)
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
# ## Try it yourself (p. 39)
#
# 1. Regress `precio` on `m2` without a constant, `lm(precio ~ 0 + m2, data = d)`. Do the residuals still sum to zero? Why?
# 2. Add `factor(mes)` to the model. What happens to R² in and out of sample?
# 3. Add a column of pure noise, `mutate(d, ruido = rnorm(n()))`, to `lm(precio ~ m2)`. Does R² rise? Does adjusted R²?
