#!/usr/bin/env python
import argparse
from pathlib import Path

import ld_aggregate as agg
from ld_hail import init_hail, read_a_block_chunk


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Section-15b central LD aggregation")
    p.add_argument("--mode", required=True, choices=["accumulate", "finalise"])
    p.add_argument("--precursor-dir", required=True)
    p.add_argument("--log-file", required=True)
    p.add_argument("--cohort-dir", help="15a cohort output dir (accumulate mode)")
    p.add_argument("--force", action="store_true")
    p.add_argument("--panel-dir", help="output panel dir (finalise mode)")
    p.add_argument("--maf-threshold", type=float, default=agg.DEFAULT_MAF_THRESHOLD)
    p.add_argument("--min-adj-diag", type=float, default=agg.DEFAULT_MIN_ADJ_DIAG)
    p.add_argument("--min-cohorts", type=int, default=None)
    p.add_argument("--a-block-size", type=int, default=None)
    p.add_argument("--a-max-dense-gb", type=float, default=1.0)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    if args.mode == "accumulate":
        if not args.cohort_dir:
            raise SystemExit("--cohort-dir is required for accumulate mode")
        hail_log = Path(args.log_file).parent / "hail.log"
        init_hail(hail_log)
        agg.accumulate(args.cohort_dir, args.precursor_dir,
                       chunk_reader=read_a_block_chunk, force=args.force)
        print(f"[section15b] accumulated cohort into {args.precursor_dir}", flush=True)
    else:
        # finalise branch added in Task 11
        raise SystemExit("finalise mode not yet implemented")


if __name__ == "__main__":
    main()
