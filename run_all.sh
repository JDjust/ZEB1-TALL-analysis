#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
"${ZEB1_RSCRIPT:-Rscript}" --vanilla "$repo_dir/scripts/06_figures/run_all.R" "$@"
