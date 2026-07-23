#!/usr/bin/env python
"""Reassemble section-15 15c per-chromosome uploads for 15b accumulate."""

from __future__ import annotations

import argparse
import copy
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import tarfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import ld_checksums


AUTOSOME_ORDER = {str(i): i for i in range(1, 23)}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Decrypt and reassemble 15c per-chromosome LD uploads"
    )
    parser.add_argument("input_dir", help="directory holding 15c .tgz.aes uploads")
    parser.add_argument("output_dir", help="merged cohort directory for 15b accumulate")
    parser.add_argument("study_name", help="cohort study_name without chromosome suffix")
    parser.add_argument(
        "--gpg",
        default=os.environ.get("GPG", "gpg"),
        help="gpg binary or wrapper to use for decryption",
    )
    return parser.parse_args()


def md5_file(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def verify_plaintext_md5(tgz_path: Path, md5sum_path: Path) -> None:
    if not md5sum_path.is_file():
        raise FileNotFoundError(f"Missing plaintext md5 file: {md5sum_path}")
    expected = md5sum_path.read_text(encoding="utf-8").split()[0]
    actual = md5_file(tgz_path)
    if actual != expected:
        raise ValueError(
            f"Plaintext md5 mismatch for {tgz_path.name}: expected {expected}, got {actual}"
        )


def safe_extract(tgz_path: Path, dest: Path) -> None:
    dest = dest.resolve()
    with tarfile.open(tgz_path, "r:gz") as tar:
        for member in tar.getmembers():
            target = (dest / member.name).resolve()
            if not str(target).startswith(str(dest) + os.sep):
                raise ValueError(f"Unsafe tar member path in {tgz_path}: {member.name}")
        tar.extractall(dest)


def decrypt_archive(input_dir: Path, staging_dir: Path, base: str, gpg: str) -> Path:
    aes_path = input_dir / f"{base}.tgz.aes"
    if not aes_path.is_file():
        raise FileNotFoundError(f"Missing encrypted archive: {aes_path}")
    tgz_path = staging_dir / f"{base}.tgz"
    if tgz_path.exists():
        tgz_path.unlink()
    subprocess.run([gpg, "--output", str(tgz_path), "--decrypt", str(aes_path)], check=True)
    verify_plaintext_md5(tgz_path, input_dir / f"{base}.md5sum")
    return tgz_path


def discover_chromosomes(input_dir: Path, study_name: str) -> list[str]:
    pattern = re.compile(rf"^{re.escape(study_name)}_chr([0-9]+)_15_scaffold\.tgz\.aes$")
    chromosomes = []
    for path in input_dir.iterdir():
        match = pattern.match(path.name)
        if match:
            chrom = match.group(1)
            if chrom not in AUTOSOME_ORDER:
                raise ValueError(f"Unsupported chromosome in 15c upload name: chr{chrom}")
            chromosomes.append(chrom)
    if not chromosomes:
        raise FileNotFoundError(
            f"No 15c scaffold archives found for study '{study_name}' in {input_dir}"
        )
    return sorted(set(chromosomes), key=lambda c: AUTOSOME_ORDER[c])


def chunk_bases(input_dir: Path, study_name: str, chrom: str) -> list[str]:
    prefix = f"{study_name}_chr{chrom}_15_chr{chrom}_chunk_"
    bases = [p.name[:-8] for p in input_dir.glob(f"{prefix}*.tgz.aes")]
    if not bases:
        raise FileNotFoundError(
            f"No chunk archives found for {study_name} chr{chrom} in {input_dir}"
        )
    return sorted(bases)


def restore_chromosome(
    input_dir: Path,
    staging_dir: Path,
    study_name: str,
    chrom: str,
    gpg: str,
) -> Path:
    chrom_dir = staging_dir / f"chr{chrom}"
    chrom_dir.mkdir(parents=True, exist_ok=True)

    scaffold_base = f"{study_name}_chr{chrom}_15_scaffold"
    scaffold_tgz = decrypt_archive(input_dir, staging_dir, scaffold_base, gpg)
    safe_extract(scaffold_tgz, chrom_dir)
    scaffold_tgz.unlink()

    a_blocks_dir = chrom_dir / "A_blocks"
    a_blocks_dir.mkdir(parents=True, exist_ok=True)
    for base in chunk_bases(input_dir, study_name, chrom):
        tgz_path = decrypt_archive(input_dir, staging_dir, base, gpg)
        safe_extract(tgz_path, a_blocks_dir)
        tgz_path.unlink()

    ld_checksums.verify_cohort_checksums(chrom_dir)
    return chrom_dir


def read_manifest(chrom_dir: Path) -> dict:
    return json.loads((chrom_dir / "manifest.json").read_text(encoding="utf-8"))


def check_contract(chrom: str, template: dict, manifest: dict) -> None:
    fields = ["study_name", "genome_build", "autosomes_only", "covariate_schema"]
    for field in fields:
        if manifest.get(field) != template.get(field):
            raise ValueError(f"Manifest field {field!r} differs for chr{chrom}")
    template_a = template["A_blocks"]
    manifest_a = manifest["A_blocks"]
    for field in ["radius_bp", "block_size", "chunk_rows", "max_dense_gb"]:
        if manifest_a.get(field) != template_a.get(field):
            raise ValueError(f"A_blocks field {field!r} differs for chr{chrom}")


def merged_counts(manifests: list[dict], n_variants: int) -> dict:
    keys = manifests[0]["variant_index"]["counts"].keys()
    source = [m["variant_index"]["counts"] for m in manifests]
    counts = {key: 0 for key in keys}
    for key in keys:
        values = [int(c[key]) for c in source]
        if key in {"total_rows", "excluded_non_autosomal"}:
            counts[key] = values[0] if len(set(values)) == 1 else max(values)
        elif key == "excluded_other_chromosome":
            counts[key] = 0
        elif key == "kept_count":
            counts[key] = int(n_variants)
        else:
            counts[key] = int(sum(values))
    return counts


def write_merged_variants(chrom_dirs: list[Path], out_path: Path) -> int:
    n_rows = 0
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out_path, "wt") as out_handle:
        wrote_header = False
        for chrom_dir in chrom_dirs:
            with gzip.open(chrom_dir / "variants.tsv.gz", "rt") as in_handle:
                header = in_handle.readline()
                if not wrote_header:
                    out_handle.write(header)
                    wrote_header = True
                for line in in_handle:
                    out_handle.write(line)
                    n_rows += 1
    return n_rows


def write_merged_b(chrom_dirs: list[Path], out_path: Path, n_variants: int) -> tuple[int, int]:
    matrices = [np.load(d / "B.npy", allow_pickle=False) for d in chrom_dirs]
    merged = np.vstack(matrices)
    if merged.shape[0] != n_variants:
        raise ValueError(
            f"Merged B.npy has {merged.shape[0]} rows but variants.tsv.gz has {n_variants}"
        )
    np.save(out_path, merged, allow_pickle=False)
    return int(merged.shape[0]), int(merged.shape[1])


def write_merged_d(chrom_dirs: list[Path], out_path: Path) -> tuple[int, int]:
    matrices = [np.load(d / "D.npy", allow_pickle=False) for d in chrom_dirs]
    first = matrices[0]
    for idx, matrix in enumerate(matrices[1:], start=2):
        if not np.array_equal(matrix, first):
            raise ValueError(f"D.npy for chromosome input {idx} differs from the first")
    np.save(out_path, first, allow_pickle=False)
    return int(first.shape[0]), int(first.shape[1])


def copy_a_blocks(chrom_dirs: list[Path], output_dir: Path) -> None:
    out_a = output_dir / "A_blocks"
    out_a.mkdir(parents=True, exist_ok=True)
    for chrom_dir in chrom_dirs:
        for source in sorted((chrom_dir / "A_blocks").glob("chr*")):
            target = out_a / source.name
            if target.exists():
                raise FileExistsError(f"Duplicate A_blocks chromosome directory: {target.name}")
            shutil.copytree(source, target)


def build_manifest(
    template: dict,
    manifests: list[dict],
    chromosomes: list[str],
    n_variants: int,
    b_shape: tuple[int, int],
    d_shape: tuple[int, int],
) -> dict:
    merged = copy.deepcopy(template)
    merged["generated_at_utc"] = datetime.now(timezone.utc).isoformat()
    merged["status"] = "reassembled"
    merged["variant_index"]["chromosome_filter"] = "all"
    merged["variant_index"]["filters_applied"][0] = "autosomes 1-22"
    merged["variant_index"]["counts"] = merged_counts(manifests, n_variants)
    merged["variant_index"]["n_samples_used"] = template["variant_index"].get("n_samples_used")
    merged["B_block"]["shape"] = list(b_shape)
    merged["A_blocks"]["chromosomes"] = {}
    for chrom, manifest in zip(chromosomes, manifests):
        chrom_meta = manifest["A_blocks"]["chromosomes"].get(chrom)
        if chrom_meta is None:
            raise ValueError(f"Manifest for chr{chrom} lacks A_blocks.chromosomes.{chrom}")
        merged["A_blocks"]["chromosomes"][chrom] = chrom_meta
    merged["reassembly"] = {
        "source": "15c per-chromosome uploads",
        "reassembled_at_utc": merged["generated_at_utc"],
        "chromosomes": [f"chr{c}" for c in chromosomes],
        "D_shape": list(d_shape),
        "source_manifest_chromosome_filters": {
            f"chr{c}": m["variant_index"].get("chromosome_filter")
            for c, m in zip(chromosomes, manifests)
        },
    }
    return merged


def write_qc_report(output_dir: Path, chromosomes: list[str], n_variants: int) -> None:
    text = [
        "Section 15 central cohort reassembly completed successfully.",
        f"Chromosomes reassembled: {', '.join('chr' + c for c in chromosomes)}",
        f"Merged variants: {n_variants}",
        "Output is suitable for one 15b accumulate call.",
    ]
    (output_dir / "qc_report.txt").write_text("\n".join(text) + "\n", encoding="utf-8")


def ensure_publishable(output_dir: Path) -> bool:
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
        return True
    entries = [p for p in output_dir.iterdir() if p.name != ".reassembly_tmp"]
    if not entries:
        return True
    if (output_dir / "manifest.json").is_file() and (output_dir / "checksums.json").is_file():
        ld_checksums.verify_cohort_checksums(output_dir)
        print(f"[ld_reassemble] existing complete output verified: {output_dir}", flush=True)
        return False
    raise FileExistsError(
        f"Output directory is not empty and is not a verified completed reassembly: {output_dir}"
    )


def reassemble(input_dir: Path, output_dir: Path, study_name: str, gpg: str) -> None:
    if not ensure_publishable(output_dir):
        return

    tmp_root = output_dir / ".reassembly_tmp"
    if tmp_root.exists():
        shutil.rmtree(tmp_root)
    restore_root = tmp_root / "restored"
    merged_root = tmp_root / "merged"
    restore_root.mkdir(parents=True)
    merged_root.mkdir(parents=True)

    chromosomes = discover_chromosomes(input_dir, study_name)
    chrom_dirs = [restore_chromosome(input_dir, restore_root, study_name, c, gpg) for c in chromosomes]
    manifests = [read_manifest(d) for d in chrom_dirs]

    template = manifests[0]
    for chrom, manifest in zip(chromosomes, manifests):
        check_contract(chrom, template, manifest)

    n_variants = write_merged_variants(chrom_dirs, merged_root / "variants.tsv.gz")
    b_shape = write_merged_b(chrom_dirs, merged_root / "B.npy", n_variants)
    d_shape = write_merged_d(chrom_dirs, merged_root / "D.npy")
    copy_a_blocks(chrom_dirs, merged_root)
    manifest = build_manifest(template, manifests, chromosomes, n_variants, b_shape, d_shape)
    (merged_root / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_qc_report(merged_root, chromosomes, n_variants)
    ld_checksums.write_cohort_checksums(merged_root)
    ld_checksums.verify_cohort_checksums(merged_root)

    for path in merged_root.iterdir():
        shutil.move(str(path), output_dir / path.name)
    shutil.rmtree(tmp_root)
    print(
        f"[ld_reassemble] reassembled {study_name} ({len(chromosomes)} chromosome(s)) into {output_dir}",
        flush=True,
    )


def main() -> None:
    args = parse_args()
    reassemble(Path(args.input_dir), Path(args.output_dir), args.study_name, args.gpg)


if __name__ == "__main__":
    main()
