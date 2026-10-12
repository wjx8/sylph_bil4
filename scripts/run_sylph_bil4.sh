#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
PYTHON_BIN=${PYTHON_BIN:-python}
exec "$PYTHON_BIN" "$ROOT/scripts/run_sylph_bil4.py" --package-root "$ROOT" "$@"
