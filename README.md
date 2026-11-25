# Econometrics II in R

R notebooks for Econometrics II (Master in Economics, Universidad de Alicante, 2026-27), taught by Henry Redondo. They follow his [slides](https://henryredondo.github.io/teaching/econometrics-ii/) analysis by analysis (A1, B2, ...) with `tidyverse`, `estimatr`, `ivreg`, `marginaleffects` and `modelsummary`. Each section gives the idea with the slide page, the code, and how to read the output.

| Notebook                                                     | Unit                                                  | Sections                                                                                                                                                                                            | Colab                                                                                                                                                                                        | Status                |
| ------------------------------------------------------------ | ----------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------- |
| [U1_cef_lp](notebooks/U1_cef_lp.ipynb)                       | Unit 1. Conditional expectation and linear projection | CEF and linear projection, law of iterated expectations, best predictor, variance decomposition                                                                                                     | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/albarran/econometrics-ii-r/blob/main/colab/U1_cef_lp.ipynb)            | runs (no data needed) |
| [U2_ols](notebooks/U2_ols.ipynb)                             | Unit 2. The linear regression model                   | 2A OLS mechanics and fit (A1-A4); 2B interpreting coefficients (B1-B9); 2C robust s.e., omitted variables, proxies, measurement error (D1, C1-C4); 2D-2E tests, heteroskedasticity, WLS (D2-D8, E1) | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/albarran/econometrics-ii-r/blob/main/colab/U2_ols.ipynb)               | data pending          |
| [U3_iv](notebooks/U3_iv.ipynb)                               | Unit 3. Instrumental variables                        | 3A simple IV (A0-A4); 3B 2SLS and weak instruments (B1-B4, Monte Carlo); 3C-3D GMM, Hansen J, multiple indicators (C1-C3, D1)                                                                       | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/albarran/econometrics-ii-r/blob/main/colab/U3_iv.ipynb)                | data pending          |
| [U4_additional_topics](notebooks/U4_additional_topics.ipynb) | Unit 4. Additional single-equation topics             | Generated regressors (A1-A2); endogeneity test (B1-B2); Sargan and Hansen J (C1); RESET, Breusch-Pagan, White (D1-D2)                                                                               | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/albarran/econometrics-ii-r/blob/main/colab/U4_additional_topics.ipynb) | data pending          |

One notebook per unit, with one section per block of the slides.

**Download** both versions (local and Colab) in one file: [econometrics-ii-r.zip](https://github.com/albarran/econometrics-ii-r/releases/latest/download/econometrics-ii-r.zip).

**Data pending** (25 November 2026): the notebooks read Henry's CSV files from `https://henryredondo.github.io/assets/files/teaching/econometrics-ii/`, where they are not published yet. Until they are, those notebooks stop with an error at the cell that reads the data; the cells before it run (for example, the simulation that opens `U3_iv`). The code has been run from top to bottom on stand-in data with the same variable names (`tools/mock_data.R`), so the numbers have not been compared with the slides yet. Names that the slides do not give and the notebooks assume: `transparency_oracle.csv` with columns `id` and `culture`; a firm identifier `id` in `transparency_firms.csv`; `firms.csv` with `size`, `exporter`, `capital` and `q` for Unit 4; the dummy `costa` in Unit 2B.

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
