#!/usr/bin/env bash
# ld_decrypt_cohort.sh — central decrypt + reassemble section-15 15c uploads
# into the cohort_stats/ tree that central accumulation consumes. No network I/O.
#
# Usage: ld_decrypt_cohort.sh <input_dir> <output_dir> <study_name> <gpg_passphrase_file>
#   input_dir   dir holding <study>_chr<C>_15_*.tgz.aes + .md5sum (left untouched)
#   output_dir  where the merged cohort_stats/ tree is rebuilt
#   study_name  cohort identifier without the 15c chromosome suffix
#   gpg_passphrase_file  readable file containing this cohort's passphrase
#
# Override the gpg binary/wrapper for testing via the GPG env var.
set -euo pipefail

if [ "$#" -ne 4 ]; then
	echo "Usage: $0 <input_dir> <output_dir> <study_name> <gpg_passphrase_file>" >&2
	exit 2
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python "${script_dir}/ld_reassemble_cohort.py" "$1" "$2" "$3" \
	--gpg-passphrase-file "$4"
