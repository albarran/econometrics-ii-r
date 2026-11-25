# %% [markdown]
# # Unit 3C and 3D. GMM, overidentification and multiple indicators
#
# Econometrics II (Henry Redondo, 2026-27). R companion to **Unit 3, Blocks C and D** of the slides: analyses C1 to C3 (efficient GMM and the Hansen J test), D1 (proxy against multiple indicators) and the final comparison of every estimator.
#
# Each section has the idea with the slide page (page of the PDF), then the code, then how to read the output. Numbers quoted as "slides" are the ones printed in the slides.

# %%
library(tidyverse)
library(broom)
library(estimatr)        # lm_robust(), iv_robust()
library(ivreg)           # ivreg()
library(gmm)             # gmm(), specTest()
library(modelsummary)

options(digits = 4, readr.show_col_types = FALSE)
theme_set(theme_minimal())

# Henry's data folder. To read local copies, set the environment variable ECMTII_DATA to their folder.
data_dir <- Sys.getenv("ECMTII_DATA", "https://henryredondo.github.io/assets/files/teaching/econometrics-ii")

# %% [markdown]
# ## The data (same as Blocks A and B)
#
# > **Data pending.** On 25 November 2026 the CSV files were not on the course website. The cell reads them from the address where they should be. The oracle file name (`transparency_oracle.csv`) and the firm identifier (`id`) are assumptions to confirm.

# %%
d <- read_csv(file.path(data_dir, "transparency_firms.csv")) |>
  mutate(
    sector      = factor(sector),
    province    = factor(province),
    gap2021     = log(wage_men_2021 / wage_women_2021),
    dgap        = log(wage_men_2025 / wage_women_2025) - gap2021,
    lsize       = log(size),
    share_women = n_women / size,
    invited     = as.integer(arm >= 1),
    letter      = as.integer(arm == 1),
    software    = as.integer(arm == 2)
  ) |>
  left_join(read_csv(file.path(data_dir, "transparency_oracle.csv")), by = "id") |>
  as.data.frame()   # gmm() needs a plain data frame

X <- "lsize + share_women + sector + province"   # controls
Z <- "letter + software"                         # instruments
fml <- function(x) as.formula(str_glue(x, .envir = parent.frame()), env = globalenv())

iv2 <- ivreg(fml("dgap ~ transp + {X} | {Z} + {X}"), data = d)   # 2SLS of Block B

# %% [markdown]
# ## C1. Small firms, noisy gaps (p. 24)
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
# ## C2. Two-step efficient GMM (p. 25)
#
# With two instruments there are two moment conditions for one parameter. GMM minimises a weighted distance of the sample moments from zero; the efficient weight is the inverse of the variance of the moments. The minimised criterion is the **Hansen J** statistic, which tests the overidentifying restriction.

# %%
g <- gmm(fml("dgap ~ transp + {X}"), fml("~ {Z} + {X}"), data = d,
         vcov = "MDS", centeredVcov = FALSE)

iv2_r <- iv_robust(fml("dgap ~ transp + {X} | {Z} + {X}"), data = d, se_type = "HC0", diagnostics = TRUE)

bind_rows("2SLS, robust (iv_robust)" = tidy(iv2_r), "GMM two-step (gmm)" = tidy(g), .id = "estimator") |>
  filter(term == "transp") |>
  select(estimator, estimate, std.error)

specTest(g)   # Hansen J

# J and its p-value as numbers, to reuse below
hansen <- function(m) { t <- specTest(m)$test; c(J = as.numeric(t[1]), p = as.numeric(t[2])) }

# %% [markdown]
# ### The two steps in matrices
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
# ## C3. Adding a bad instrument (p. 26)
#
# Add `peer_share`, which is related to culture, to the instrument list. Now there are three moments and one parameter.

# %%
g3 <- gmm(fml("dgap ~ transp + {X}"), fml("~ {Z} + peer_share + {X}"), data = d,
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
# ## D1. Proxy against multiple indicators (p. 8-9)
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
# ## Every estimator against the truth (3D p. 12)
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
# ## Try it yourself (3C p. 27, 3D p. 10)
#
# 1. Compare `coef(g)` with `coef(iv2)`. They agree to four decimals: why is the efficiency gain zero here?
# 2. Replace `peer_share` by `parent_abroad` in C3. Does J catch it too? Why is it harder?
# 3. Why is `transp` not instrumented in D1? What are we assuming about it?
# 4. Could you add `letter` and `software` as extra instruments in D1? For which regressor?
