#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
SYLPH_BIN=${SYLPH_BIN:-sylph}
THREADS=${THREADS:-8}
DB_DIR="$ROOT/database/sylph_BIL4"
DB="$DB_DIR/BIL4.syldb"

command -v "$SYLPH_BIN" >/dev/null 2>&1 || {
  echo "ERROR: Sylph not found: $SYLPH_BIN" >&2
  exit 1
}
[[ "$($SYLPH_BIN --version)" == "sylph 0.5.1" ]] || {
  echo "ERROR: Sylph 0.5.1 is required" >&2
  exit 1
}

python "$ROOT/scripts/validate_assets.py" --root "$ROOT" --skip-database

if [[ -s "$DB" ]]; then
  echo "BIL4.syldb already exists and the reference MD5 values passed; keeping it."
  exit 0
fi

tmp_dir=$(mktemp -d "$DB_DIR/.build.XXXXXX")
cleanup() {
  rm -rf "$tmp_dir"
}
trap cleanup EXIT

"$SYLPH_BIN" sketch "$DB_DIR"/genomes/*.fna \
  -o "$tmp_dir/BIL4" \
  -c 200 \
  -t "$THREADS"
[[ -s "$tmp_dir/BIL4.syldb" ]] || {
  echo "ERROR: Sylph database build failed" >&2
  exit 1
}
mv "$tmp_dir/BIL4.syldb" "$DB"
echo "Built $DB"
