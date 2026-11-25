# %% [markdown]
# # Unit 1. Conditional expectation and linear projection
#
# Econometrics II (Henry Redondo, 2026-27). R companion to **Unit 1** of the slides. The unit is theory and has no code in the slides: this notebook checks its main results in a simulated population (the code is written for this notebook; the figure follows slides p. 34-37). No data files are needed.

# %%
library(tidyverse)
library(broom)

options(digits = 4)
theme_set(theme_minimal())

# %% [markdown]
# ## A simulated population
#
# `x` is years of education (0 to 20) and the conditional expectation of `y` is concave in `x`: $E(y \mid x) = 1 + 3\,(1 - e^{-0.15x})$. With a very large sample, sample means are close to population expectations.

# %%
set.seed(2026)
cef <- function(x) 1 + 3 * (1 - exp(-0.15 * x))

pop <- tibble(x = sample(0:20, 200000, replace = TRUE, prob = dbinom(0:20, 20, 0.55)),
              y = cef(x) + rnorm(200000, sd = 0.4))

# %% [markdown]
# ## The CEF and the linear projection (p. 34-37, 43)
#
# The CEF is estimated by the mean of `y` at each value of `x`. The linear projection $L(y \mid 1, x) = \beta_0 + \beta_1 x$ has $\beta_1 = Cov(x, y)/Var(x)$ and $\beta_0 = E(y) - \beta_1 E(x)$: what `lm` computes.

# %%
cef_hat <- pop |> summarise(cef_hat = mean(y), n = n(), .by = x) |> arrange(x)

lp <- lm(y ~ x, data = pop)
b1 <- cov(pop$x, pop$y) / var(pop$x)
b0 <- mean(pop$y) - b1 * mean(pop$x)
stopifnot(all.equal(c(b0, b1), unname(coef(lp))))
coef(lp)

# %%
ggplot(slice_sample(pop, n = 1500), aes(x, y)) +
  geom_point(alpha = 0.15) +
  geom_line(aes(y = cef_hat, colour = "E(y | x)"), data = cef_hat, linewidth = 1) +
  geom_abline(aes(intercept = b0, slope = b1, colour = "L(y | 1, x)"), linewidth = 1) +
  labs(x = "x (e.g. years of education)", colour = NULL)

# %% [markdown]
# **Reading.** The CEF is a curve; the linear projection is the best straight-line approximation to it. It fits well where most of the population is and misses at the extremes.

# %% [markdown]
# ## Law of iterated expectations (p. 15-17)
#
# $E(y) = E[E(y \mid x)]$: the average of the conditional means, weighted by how frequent each `x` is, is the unconditional mean.

# %%
c(mean_y = mean(pop$y), mean_of_cef = with(cef_hat, weighted.mean(cef_hat, n)))

# %% [markdown]
# ## The best predictor (p. 20, 37)
#
# Among all functions of `x`, the CEF has the smallest mean squared error; among linear functions, the linear projection does. Any other line does worse.

# %%
pop <- pop |> left_join(select(cef_hat, x, cef_hat), by = "x") |> mutate(lp_hat = fitted(lp))

pop |> summarise(mse_cef        = mean((y - cef_hat)^2),
                 mse_lp         = mean((y - lp_hat)^2),
                 mse_other_line = mean((y - (b0 + 0.1 + (b1 - 0.01) * x))^2))

# %% [markdown]
# ## Variance decomposition (p. 28)
#
# $Var(y) = Var[E(y \mid x)] + E[Var(y \mid x)]$: explained plus unexplained variance.

# %%
pop |> summarise(var_y = mean((y - mean(y))^2),
                 var_of_cef  = mean((cef_hat - mean(y))^2),
                 mean_cond_var = mean((y - cef_hat)^2)) |>
  mutate(sum = var_of_cef + mean_cond_var)

# %% [markdown]
# ## Projecting the CEF (p. 56)
#
# The linear projection of `y` on `x` and the linear projection of the CEF on `x` are the same line.

# %%
lp_of_cef <- lm(cef_hat ~ x, data = pop)
stopifnot(all.equal(coef(lp_of_cef), coef(lp)))
rbind(lp_of_y = coef(lp), lp_of_cef = coef(lp_of_cef))

# %% [markdown]
# ## CEF error and LP error (p. 26, 57-58)
#
# The CEF error $u = y - E(y \mid x)$ has mean zero at **every** value of `x`. The LP error $e = y - L(y \mid 1, x)$ is uncorrelated with `x` by construction, and its mean still depends on `x` when the CEF is not linear.

# %%
pop |> summarise(cor_x_cef_error = cor(x, y - cef_hat), cor_x_lp_error = cor(x, y - lp_hat))

pop |>
  summarise(`CEF error` = mean(y - cef_hat), `LP error` = mean(y - lp_hat), .by = x) |>
  pivot_longer(-x, names_to = "error", values_to = "mean") |>
  ggplot(aes(x, mean, colour = error)) +
  geom_line() + geom_point() +
  labs(y = "Mean of the error at each x", colour = NULL)

# %% [markdown]
# **Reading.** Both correlations are zero. The mean of the LP error is positive in the middle and negative at the extremes: zero correlation with `x` is weaker than mean independence. This difference is what separates "the linear projection always exists" from "the linear model is the CEF", and it comes back in Unit 2 with the assumptions of OLS.
