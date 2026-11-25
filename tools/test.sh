#!/usr/bin/env bash
# Run every notebooks/*.ipynb top to bottom in a temporary copy (the repository copies stay
# without outputs). Usage: tools/test.sh [name ...]
#   ECMTII_DATA=<folder>  read the data from a local folder (e.g. the one written by
#                         tools/mock_data.R) instead of Henry's website
#   JUPYTER=<path>        jupyter to use (default: ~/venvs/dma/bin/jupyter)
set -u
cd "$(dirname "$0")/.."
JUPYTER="${JUPYTER:-$HOME/venvs/dma/bin/jupyter}"
out="$(mktemp -d "${TMPDIR:-/tmp}/ecmtii-test.XXXXXX")"
[ $# -gt 0 ] && set -- "${@/%/.ipynb}" && set -- "${@/#/notebooks/}" || set -- notebooks/*.ipynb
fail=0
for nb in "$@"; do
  if "$JUPYTER" nbconvert --to notebook --execute --ExecutePreprocessor.kernel_name=ir-dma \
       --ExecutePreprocessor.timeout=1200 --output-dir "$out" "$nb" > "$out/log.txt" 2>&1; then
    echo "OK    $nb"
  else
    echo "FAIL  $nb"; tail -25 "$out/log.txt"; fail=1
  fi
done
echo "Executed copies: $out"
exit $fail
