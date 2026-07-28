# resources/genetics/ld_checksums.py
"""Fast file checksums for section-15 federated artefact integrity.

Cohort-side (15a) emits a ``checksums.json`` over the artefacts that section-15
packaging and central aggregation consume (``B.npy``, ``D.npy``,
``variants.tsv.gz``, ``manifest.json``, and every file under ``A_blocks/``).
Central aggregate ``accumulate`` re-verifies before folding a cohort in, so a
corrupted or truncated upload is caught *before* the cohort's raw data is discarded.

Hashing uses ``hashlib.blake2b`` (Python stdlib): meaningfully faster than
sha256 in software, cryptographically strong, and adds no new dependency to the
federated cohort-site install footprint. This is integrity-against-corruption,
not adversarial tamper-proofing.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

CHECKSUM_ALGORITHM = "blake2b"
CHECKSUM_FILENAME = "checksums.json"
# Per-cohort artefacts that section-15 packaging and central aggregation consume.
CORE_ARTIFACTS = ("B.npy", "D.npy", "variants.tsv.gz", "manifest.json")
_READ_CHUNK = 1 << 20  # 1 MiB streaming reads


class ChecksumError(RuntimeError):
    """Raised when a cohort artefact fails integrity verification."""


def hash_file(path: str | Path) -> str:
    """Return the blake2b hex digest of a file, read in streamed blocks."""
    h = hashlib.blake2b()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(_READ_CHUNK), b""):
            h.update(block)
    return h.hexdigest()


def iter_artifact_relpaths(cohort_dir: str | Path) -> list[str]:
    """Sorted POSIX relative paths of every artefact to checksum for a cohort.

    Covers the core files plus every file under ``A_blocks/`` recursively;
    excludes the checksums file itself.
    """
    cohort_dir = Path(cohort_dir)
    rels: list[str] = []
    for name in CORE_ARTIFACTS:
        if (cohort_dir / name).is_file():
            rels.append(name)
    a_blocks = cohort_dir / "A_blocks"
    if a_blocks.is_dir():
        for p in a_blocks.rglob("*"):
            if p.is_file():
                rels.append(p.relative_to(cohort_dir).as_posix())
    return sorted(rels)


def compute_checksums(cohort_dir: str | Path, relpaths: list[str]) -> dict[str, str]:
    """Map each relative path to its blake2b digest."""
    cohort_dir = Path(cohort_dir)
    return {rel: hash_file(cohort_dir / rel) for rel in relpaths}


def write_cohort_checksums(cohort_dir: str | Path) -> dict:
    """Enumerate, hash, and atomically write ``checksums.json``; return the doc.

    Call this as the final cohort-side step, after ``manifest.json`` is written.
    """
    cohort_dir = Path(cohort_dir)
    relpaths = iter_artifact_relpaths(cohort_dir)
    doc = {"algorithm": CHECKSUM_ALGORITHM,
           "files": compute_checksums(cohort_dir, relpaths)}
    path = cohort_dir / CHECKSUM_FILENAME
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    tmp.write_text(json.dumps(doc, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(tmp, path)
    return doc


def read_cohort_checksums(cohort_dir: str | Path) -> dict:
    """Read ``checksums.json``; raise if it is absent."""
    path = Path(cohort_dir) / CHECKSUM_FILENAME
    if not path.is_file():
        raise ChecksumError(
            f"No {CHECKSUM_FILENAME} in {cohort_dir}; cannot verify cohort integrity")
    return json.loads(path.read_text(encoding="utf-8"))


def verify_cohort_checksums(cohort_dir: str | Path) -> None:
    """Recompute every listed digest and hard-fail on any mismatch or absence."""
    cohort_dir = Path(cohort_dir)
    doc = read_cohort_checksums(cohort_dir)
    algo = doc.get("algorithm")
    if algo != CHECKSUM_ALGORITHM:
        raise ChecksumError(
            f"{CHECKSUM_FILENAME} uses algorithm {algo!r}, "
            f"expected {CHECKSUM_ALGORITHM!r}")
    files = doc.get("files") or {}
    if not files:
        raise ChecksumError(f"{CHECKSUM_FILENAME} lists no files to verify")
    for rel, expected in sorted(files.items()):
        target = cohort_dir / rel
        if not target.is_file():
            raise ChecksumError(f"Checksummed artefact missing: {rel}")
        actual = hash_file(target)
        if actual != expected:
            raise ChecksumError(
                f"Checksum mismatch for {rel}: expected {expected}, got {actual}")
