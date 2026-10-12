#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
PYTHON_BIN=${PYTHON_BIN:-python}
SYLPH_BIN=${SYLPH_BIN:-sylph}

command -v "$PYTHON_BIN" >/dev/null 2>&1 || {
  echo "ERROR: Python not found: $PYTHON_BIN" >&2
  exit 1
}
command -v "$SYLPH_BIN" >/dev/null 2>&1 || {
  echo "ERROR: Sylph not found: $SYLPH_BIN" >&2
  echo "Activate the environment first: mamba activate sylph-bil4" >&2
  exit 1
}

version=$("$SYLPH_BIN" --version)
[[ "$version" == "sylph 0.5.1" ]] || {
  echo "ERROR: expected sylph 0.5.1, observed: $version" >&2
  exit 1
}

echo "Python: $($PYTHON_BIN --version 2>&1)"
echo "Sylph: $version"
PYTHONPATH="$ROOT/scripts${PYTHONPATH:+:$PYTHONPATH}" \
  "$PYTHON_BIN" "$ROOT/scripts/validate_assets.py" --root "$ROOT"
echo "Installation and bundled reference assets are ready."
