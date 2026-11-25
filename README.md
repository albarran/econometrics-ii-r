# Econometrics II in R

R notebooks for Econometrics II (Master in Economics, Universidad de Alicante, 2026-27), taught by Henry Redondo. They follow his [slides](https://henryredondo.github.io/teaching/econometrics-ii/) analysis by analysis (A1, B2, ...) with `tidyverse`, `estimatr`, `ivreg`, `marginaleffects` and `modelsummary`. Each section gives the idea with the slide page, the code, and how to read the output.

| Notebook                                                                       | Topic                                                           | Slides                       | Colab                                                                                                                                                                                                 | Status                |
| ------------------------------------------------------------------------------ | --------------------------------------------------------------- | ---------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------- |
| [U1_cef_lp](notebooks/U1_cef_lp.ipynb)                                         | Conditional expectation and linear projection                   | Unit 1                       | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/albarran/econometrics-ii-r/blob/main/colab/U1_cef_lp.ipynb)                     | runs (no data needed) |
| [U2A_ols_mechanics](notebooks/U2A_ols_mechanics.ipynb)                         | OLS mechanics, R², out-of-sample fit                            | Unit 2A (A1-A4)              | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/albarran/econometrics-ii-r/blob/main/colab/U2A_ols_mechanics.ipynb)             | data pending          |
| [U2B_interpreting_coefficients](notebooks/U2B_interpreting_coefficients.ipynb) | Logs, quadratics, dummies, interactions, price index            | Unit 2B (B1-B9)              | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/albarran/econometrics-ii-r/blob/main/colab/U2B_interpreting_coefficients.ipynb) | data pending          |
| [U2C_endogeneity](notebooks/U2C_endogeneity.ipynb)                             | Robust s.e., omitted variables, proxies, measurement error      | Unit 2C (D1, C1-C4)          | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/albarran/econometrics-ii-r/blob/main/colab/U2C_endogeneity.ipynb)               | data pending          |
| [U2DE_inference_gls](notebooks/U2DE_inference_gls.ipynb)                       | t, F, Wald and LM tests, heteroskedasticity, WLS                | Unit 2D-2E (D2-D8, E1)       | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/albarran/econometrics-ii-r/blob/main/colab/U2DE_inference_gls.ipynb)            | data pending          |
| [U3A_iv](notebooks/U3A_iv.ipynb)                                               | The simple IV estimator                                         | Unit 3A (A0-A4)              | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/albarran/econometrics-ii-r/blob/main/colab/U3A_iv.ipynb)                        | data pending          |
| [U3B_2sls_weak](notebooks/U3B_2sls_weak.ipynb)                                 | 2SLS and weak instruments                                       | Unit 3B (B1-B4, Monte Carlo) | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/albarran/econometrics-ii-r/blob/main/colab/U3B_2sls_weak.ipynb)                 | data pending          |
| [U3CD_gmm_indicators](notebooks/U3CD_gmm_indicators.ipynb)                     | GMM, Hansen J, multiple indicators                              | Unit 3C-3D (C1-C3, D1)       | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/albarran/econometrics-ii-r/blob/main/colab/U3CD_gmm_indicators.ipynb)           | data pending          |
| [U4_additional_topics](notebooks/U4_additional_topics.ipynb)                   | Generated regressors, endogeneity, Sargan, RESET, Breusch-Pagan | Unit 4                       | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/albarran/econometrics-ii-r/blob/main/colab/U4_additional_topics.ipynb)          | data pending          |

**Data pending** (25 November 2026): the notebooks read Henry's CSV files from `https://henryredondo.github.io/assets/files/teaching/econometrics-ii/`, where they are not published yet. Until they are, those notebooks stop with an error at the cell that reads the data; the cells before it run (for example, the simulation that opens `U3A_iv`). The code has been run from top to bottom on stand-in data with the same variable names (`tools/mock_data.R`), so the numbers have not been compared with the slides yet. Names that the slides do not give and the notebooks assume: `transparency_oracle.csv` with columns `id` and `culture`; a firm identifier `id` in `transparency_firms.csv`; `firms.csv` with `size`, `exporter`, `capital` and `q` for Unit 4; the dummy `costa` in Unit 2B.

If you have the files in a folder, set the environment variable `ECMTII_DATA` to that folder before starting Jupyter and the notebooks read them from there.

## Run in Colab

Click the badge. The first code cell installs the packages that the Colab R runtime lacks (binaries, one or two minutes); then *Runtime > Run all*.

## Run on your computer

Follow the [DMA setup](https://albarran.github.io/DMA/setup.html), open a notebook of the folder `notebooks/` in Jupyter and choose the kernel **R (DMA)**.

## How the notebooks are built

One source per notebook in `src/` (an R script in percent format). `python3 tools/build.py` writes the two versions without outputs: `notebooks/` (kernel `ir-dma`, "R (DMA)") and `colab/` (kernel `ir`, with the installation cell). Edit the source, never the `.ipynb`.

```bash
python3 tools/build.py
Rscript tools/mock_data.R /tmp/ecmtii-mock        # stand-in data, for testing only
ECMTII_DATA=/tmp/ecmtii-mock tools/test.sh        # runs every notebook in a temporary copy
```
