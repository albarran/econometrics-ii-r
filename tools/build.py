#!/usr/bin/env python3
"""Build the notebooks from src/*.R (one source per notebook, percent format).

    python3 tools/build.py

Writes two versions of every source, both without outputs:
  notebooks/<name>.ipynb  kernel "R (DMA)" (ir-dma), for local Jupyter
  colab/<name>.ipynb      kernel "ir", with a first cell that installs in Colab
                          the packages the notebook loads and does nothing elsewhere

Source format: a cell starts at a line "# %%" (code) or "# %% [markdown]"
(markdown, every line prefixed with "# ").
"""
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REPO = "albarran/econometrics-ii-r"
BASE = {"tidyverse"}  # already in the Colab R runtime; listed anyway, installed only if missing

COLAB_CELL = '''# Colab only: installs the packages this notebook needs (1-2 minutes the first time).
# On your own computer this cell does nothing.
need <- c({pkgs})
if (Sys.getenv("COLAB_RELEASE_TAG") != "" || dir.exists("/content")) {{
  os <- readLines("/etc/os-release")
  codename <- sub(".*=", "", grep("^(UBUNTU|VERSION)_CODENAME=", os, value = TRUE)[1])
  if (is.na(codename)) codename <- "jammy"
  # binaries from Posit Package Manager (no compilation); they need R's user agent
  options(repos = c(CRAN = sprintf("https://packagemanager.posit.co/cran/__linux__/%s/latest", codename)),
          HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(),
                                  paste(getRversion(), R.version$platform, R.version$arch, R.version$os)))
  missing <- setdiff(need, rownames(installed.packages()))
  if (length(missing) > 0) install.packages(missing, Ncpus = 2)
}}'''

KERNELS = {
    "notebooks": {"name": "ir-dma", "display_name": "R (DMA)", "language": "R"},
    "colab": {"name": "ir", "display_name": "R", "language": "R"},
}


def parse(text):
    cells, kind, buf = [], None, []

    def flush():
        if kind is None:
            return
        lines = buf[:]
        while lines and not lines[-1].strip():
            lines.pop()
        while lines and not lines[0].strip():
            lines.pop(0)
        if kind == "markdown":
            lines = [re.sub(r"^# ?", "", l) for l in lines]
        if lines:
            cells.append((kind, "\n".join(lines)))

    for line in text.splitlines():
        if line.startswith("# %%"):
            flush()
            kind = "markdown" if "[markdown]" in line else "code"
            buf = []
        elif kind is not None:
            buf.append(line)
    flush()
    return cells


def cell(kind, source):
    c = {"cell_type": kind, "metadata": {}, "source": source.splitlines(keepends=True)}
    if kind == "code":
        c.update(execution_count=None, outputs=[])
    return c


def notebook(cells, kernel, name):
    meta = {"kernelspec": kernel, "language_info": {"name": "R"}}
    if kernel["name"] == "ir":
        meta["colab"] = {"name": name + ".ipynb", "provenance": []}
    return {"cells": [cell(k, s) for k, s in cells], "metadata": meta, "nbformat": 4, "nbformat_minor": 5}


def main():
    for src in sorted((ROOT / "src").glob("*.R")):
        cells = parse(src.read_text(encoding="utf-8"))
        code = "\n".join(s for k, s in cells if k == "code")
        pkgs = sorted(set(re.findall(r"library\((\w+)\)", code)) | set(re.findall(r"\b(\w+)::", code)))
        pkgs = [p for p in pkgs if p not in ("stats", "utils", "base")]
        install = COLAB_CELL.format(pkgs=", ".join(f'"{p}"' for p in pkgs))
        first_code = next(i for i, (k, _) in enumerate(cells) if k == "code")
        versions = {"notebooks": cells, "colab": cells[:first_code] + [("code", install)] + cells[first_code:]}
        for folder, cs in versions.items():
            out = ROOT / folder / (src.stem + ".ipynb")
            out.parent.mkdir(exist_ok=True)
            out.write_text(json.dumps(notebook(cs, KERNELS[folder], src.stem), indent=1, ensure_ascii=False) + "\n",
                           encoding="utf-8")
        print(f"{src.stem}: {len(cells)} cells, packages: {' '.join(pkgs)}")


if __name__ == "__main__":
    main()
