#!/usr/bin/env bash
# Compatibility wrapper for the old section-15 cohort command.
set -e

script_dir="$(dirname "${BASH_SOURCE[0]}")"

echo "WARNING: 15c-ld_run_upload.sh is deprecated."
echo "Running the staged section-15 workflow instead: 15a, 15b, then check_upload.sh 15 upload."

bash "${script_dir}/15a-ld_prepare_stats.sh" "$@"
bash "${script_dir}/15b-ld_compress_data.sh" "$@"
bash "${script_dir}/check_upload.sh" 15 upload "$@"
