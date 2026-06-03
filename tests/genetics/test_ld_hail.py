import numpy as np
import pytest
from hail.linalg import BlockMatrix

import ld_hail
from synthetic_plink import write_synthetic_plink as _write_synthetic_plink


def test_load_genotype_matrixtable_filters_chromosome_rows_and_orders_columns(
    hail_session, tmp_path
):
    from ld_hail import load_genotype_matrixtable
    from ld_qc import build_variant_index

    bfile = tmp_path / "test"
    _write_synthetic_plink(
        bfile,
        calls=[
            [0, 1, 2, 0],
            [2, None, 1, 0],
            [1, 1, 1, 1],
            [0, 2, 0, 2],
        ],
        variants=[
            ("22", 100, "A", "G"),
            ("22", 200, "C", "T"),
            ("21", 300, "T", "A"),
            ("22", 400, "G", "C"),
        ],
        samples=["IID1", "IID2", "IID3", "IID4"],
    )

    variant_index = build_variant_index(f"{bfile}.bim", chromosome="22")
    mt = load_genotype_matrixtable(
        bfile=bfile,
        chromosome="22",
        final_samples=["IID3", "IID1"],
        variant_index=variant_index,
        n_partitions=2,
    )

    assert mt.count_rows() == 3
    assert mt.count_cols() == 2

    expected_ids = [v["variant_id"] for v in variant_index["variants"]]
    assert mt.variant_id.collect() == expected_ids
    assert [str(s) for s in mt.s.collect()] == ["IID3", "IID1"]


def test_load_genotype_matrixtable_raises_on_missing_sample(hail_session, tmp_path):
    from ld_hail import load_genotype_matrixtable
    from ld_qc import build_variant_index

    bfile = tmp_path / "test"
    _write_synthetic_plink(
        bfile,
        calls=[[0, 1], [2, 1]],
        variants=[("22", 100, "A", "G"), ("22", 200, "C", "T")],
        samples=["IID1", "IID2"],
    )
    variant_index = build_variant_index(f"{bfile}.bim", chromosome="22")

    with pytest.raises(ValueError, match="final_samples missing"):
        load_genotype_matrixtable(
            bfile=bfile,
            chromosome="22",
            final_samples=["IID1", "IID_GHOST"],
            variant_index=variant_index,
            n_partitions=2,
        )


def test_prepare_for_cross_products_diagnostics_and_imputation(hail_session, tmp_path):
    from ld_hail import load_genotype_matrixtable, prepare_for_cross_products
    from ld_qc import build_variant_index

    bfile = tmp_path / "test"
    _write_synthetic_plink(
        bfile,
        calls=[
            [0, 1, 2, 0],
            [1, None, 1, 2],
        ],
        variants=[("22", 100, "A", "G"), ("22", 200, "C", "T")],
        samples=["IID1", "IID2", "IID3", "IID4"],
    )
    variant_index = build_variant_index(f"{bfile}.bim", chromosome="22")
    mt = load_genotype_matrixtable(
        bfile=bfile,
        chromosome="22",
        final_samples=["IID1", "IID2", "IID3", "IID4"],
        variant_index=variant_index,
        n_partitions=2,
    )

    mt_imputed, diagnostics = prepare_for_cross_products(mt)

    assert [d["variant_id"] for d in diagnostics] == [
        v["variant_id"] for v in variant_index["variants"]
    ]
    assert diagnostics[0]["n_nonmissing"] == 4
    assert diagnostics[0]["n_imputed"] == 0
    assert diagnostics[0]["genotype_mean"] == pytest.approx(0.75)
    assert diagnostics[1]["n_nonmissing"] == 3
    assert diagnostics[1]["n_imputed"] == 1
    assert diagnostics[1]["genotype_mean"] == pytest.approx(4.0 / 3.0)

    rows = mt_imputed.entries().select("variant_id", "GT_dosage").collect()
    by_pair = {(r.variant_id, r.s): float(r.GT_dosage) for r in rows}
    assert by_pair[("22:100:A:G", "IID2")] == pytest.approx(1.0)
    assert by_pair[("22:200:C:T", "IID2")] == pytest.approx(4.0 / 3.0)


def test_prepare_for_cross_products_raises_on_all_missing_variant(
    hail_session, tmp_path
):
    from ld_hail import load_genotype_matrixtable, prepare_for_cross_products
    from ld_qc import build_variant_index

    bfile = tmp_path / "test"
    _write_synthetic_plink(
        bfile,
        calls=[[None, None], [0, 1]],
        variants=[("22", 100, "A", "G"), ("22", 200, "C", "T")],
        samples=["IID1", "IID2"],
    )
    variant_index = build_variant_index(f"{bfile}.bim", chromosome="22")
    mt = load_genotype_matrixtable(
        bfile=bfile,
        chromosome="22",
        final_samples=["IID1", "IID2"],
        variant_index=variant_index,
        n_partitions=2,
    )
    with pytest.raises(ValueError, match="entirely missing"):
        prepare_for_cross_products(mt)


def test_compute_b_block_matches_direct_numpy(hail_session, tmp_path):
    import numpy as np

    from ld_hail import (
        compute_b_block,
        load_genotype_matrixtable,
        prepare_for_cross_products,
    )
    from ld_qc import build_variant_index

    bfile = tmp_path / "test"
    calls = [
        [0, 1, 2, 0],
        [1, None, 1, 2],
        [2, 0, 1, 1],
    ]
    _write_synthetic_plink(
        bfile,
        calls=calls,
        variants=[
            ("22", 100, "A", "G"),
            ("22", 200, "C", "T"),
            ("22", 300, "T", "A"),
        ],
        samples=["IID1", "IID2", "IID3", "IID4"],
    )
    variant_index = build_variant_index(f"{bfile}.bim", chromosome="22")
    final_samples = ["IID1", "IID2", "IID3", "IID4"]
    mt = load_genotype_matrixtable(
        bfile=bfile,
        chromosome="22",
        final_samples=final_samples,
        variant_index=variant_index,
        n_partitions=2,
    )
    mt_imputed, _ = prepare_for_cross_products(mt)

    covariate_matrix = np.array(
        [
            [1.0, 30.0, 1.0],
            [1.0, 45.0, 2.0],
            [1.0, 60.0, 1.0],
            [1.0, 55.0, 2.0],
        ],
        dtype=np.float64,
    )

    b = compute_b_block(mt_imputed, covariate_matrix)

    x_imputed = np.array(
        [
            [0.0, 1.0, 2.0, 0.0],
            [1.0, 4.0 / 3.0, 1.0, 2.0],
            [2.0, 0.0, 1.0, 1.0],
        ]
    )
    expected = x_imputed @ covariate_matrix
    assert b.shape == (3, 3)
    assert b.dtype == np.float64
    np.testing.assert_allclose(b, expected, rtol=1e-12)


def test_compute_b_block_rejects_mismatched_sample_count(hail_session, tmp_path):
    import numpy as np

    from ld_hail import (
        compute_b_block,
        load_genotype_matrixtable,
        prepare_for_cross_products,
    )
    from ld_qc import build_variant_index

    bfile = tmp_path / "test"
    _write_synthetic_plink(
        bfile,
        calls=[[0, 1, 2, 0]],
        variants=[("22", 100, "A", "G")],
        samples=["IID1", "IID2", "IID3", "IID4"],
    )
    variant_index = build_variant_index(f"{bfile}.bim", chromosome="22")
    mt = load_genotype_matrixtable(
        bfile=bfile,
        chromosome="22",
        final_samples=["IID1", "IID2", "IID3", "IID4"],
        variant_index=variant_index,
        n_partitions=2,
    )
    mt_imputed, _ = prepare_for_cross_products(mt)

    bad_covariates = np.zeros((3, 3), dtype=np.float64)
    with pytest.raises(ValueError, match="must match exactly"):
        compute_b_block(mt_imputed, bad_covariates)


def test_compute_a_block_banded_drops_out_of_band_blocks(hail_session, tmp_path):
    import random

    import numpy as np
    from hail.linalg import BlockMatrix

    from ld_hail import (
        compute_a_block_banded,
        load_genotype_matrixtable,
        prepare_for_cross_products,
    )
    from ld_qc import build_variant_index

    rng = random.Random(42)
    n_variants = 17
    n_samples = 4
    block_size = 8
    radius_bp = 4_000_000

    alleles_cycle = [("A", "G"), ("C", "T"), ("T", "C"), ("G", "A")]
    variants = [
        ("22", (i + 1) * 1_000_000, *alleles_cycle[i % len(alleles_cycle)])
        for i in range(n_variants)
    ]
    calls = [
        [rng.randint(0, 2) for _ in range(n_samples)] for _ in range(n_variants)
    ]
    samples = [f"IID{i + 1}" for i in range(n_samples)]

    bfile = tmp_path / "test"
    _write_synthetic_plink(bfile, calls=calls, variants=variants, samples=samples)
    variant_index = build_variant_index(f"{bfile}.bim", chromosome="22")

    mt = load_genotype_matrixtable(
        bfile=bfile,
        chromosome="22",
        final_samples=samples,
        variant_index=variant_index,
        n_partitions=2,
    )
    mt_imputed, _ = prepare_for_cross_products(mt)

    out_dir = tmp_path / "A_blocks"
    meta = compute_a_block_banded(
        mt_imputed,
        variant_index,
        out_dir=out_dir,
        radius_bp=radius_bp,
        block_size=block_size,
    )

    chr22_meta = meta["chromosomes"]["22"]
    assert chr22_meta["n_variants"] == n_variants
    assert chr22_meta["max_idx_distance_in_band"] == 4

    a_back = BlockMatrix.read(str(out_dir / "chr22")).to_numpy()
    x_imputed = np.array(calls, dtype=np.float64)
    full_a = x_imputed @ x_imputed.T

    # Build the expected kept/dropped mask matching Hail's sparsify_band(0, 4, blocks_only=True).
    # A block (r, c) is KEPT iff the band [0, 4] overlaps the block's j-i range.
    n_row_blocks = (n_variants + block_size - 1) // block_size
    dropped_mask = np.zeros((n_variants, n_variants), dtype=bool)
    for r in range(n_row_blocks):
        for c in range(n_row_blocks):
            i_lo, i_hi = r * block_size, min((r + 1) * block_size, n_variants)
            j_lo, j_hi = c * block_size, min((c + 1) * block_size, n_variants)
            block_jmi_min = j_lo - (i_hi - 1)
            block_jmi_max = (j_hi - 1) - i_lo
            kept = max(block_jmi_min, 0) <= min(block_jmi_max, 4)
            if not kept:
                dropped_mask[i_lo:i_hi, j_lo:j_hi] = True

    # At least one block must be dropped for the banding test to be meaningful.
    assert dropped_mask.any()

    np.testing.assert_array_equal(a_back[dropped_mask], 0.0)
    np.testing.assert_allclose(a_back[~dropped_mask], full_a[~dropped_mask], rtol=1e-12)


def test_compute_a_block_banded_writes_per_chromosome_directories(
    hail_session, tmp_path
):
    from hail.linalg import BlockMatrix

    from ld_hail import (
        compute_a_block_banded,
        load_genotype_matrixtable,
        prepare_for_cross_products,
    )
    from ld_qc import build_variant_index

    bfile = tmp_path / "test"
    _write_synthetic_plink(
        bfile,
        calls=[
            [0, 1, 2, 0],
            [1, 0, 1, 2],
            [2, 1, 0, 1],
        ],
        variants=[
            ("1", 100, "A", "G"),
            ("1", 200, "C", "T"),
            ("2", 300, "T", "A"),
        ],
        samples=["IID1", "IID2", "IID3", "IID4"],
    )
    variant_index = build_variant_index(f"{bfile}.bim")
    mt = load_genotype_matrixtable(
        bfile=bfile,
        chromosome=None,
        final_samples=["IID1", "IID2", "IID3", "IID4"],
        variant_index=variant_index,
        n_partitions=2,
    )
    mt_imputed, _ = prepare_for_cross_products(mt)

    out_dir = tmp_path / "A_blocks"
    meta = compute_a_block_banded(
        mt_imputed, variant_index, out_dir=out_dir, block_size=8,
    )

    assert set(meta["chromosomes"].keys()) == {"1", "2"}
    assert (out_dir / "chr1").is_dir()
    assert (out_dir / "chr2").is_dir()

    a1 = BlockMatrix.read(str(out_dir / "chr1")).to_numpy()
    a2 = BlockMatrix.read(str(out_dir / "chr2")).to_numpy()
    assert a1.shape == (2, 2)
    assert a2.shape == (1, 1)


def test_load_genotype_matrixtable_genotype_calls_round_trip(hail_session, tmp_path):
    """Confirm Hail decodes PLINK calls so mt.GT.n_alt_alleles() matches input."""
    from ld_hail import load_genotype_matrixtable
    from ld_qc import build_variant_index

    bfile = tmp_path / "test"
    calls = [
        [0, 1, 2, None],
        [2, 0, 1, 1],
    ]
    _write_synthetic_plink(
        bfile,
        calls=calls,
        variants=[("22", 100, "A", "G"), ("22", 200, "C", "T")],
        samples=["IID1", "IID2", "IID3", "IID4"],
    )
    variant_index = build_variant_index(f"{bfile}.bim", chromosome="22")

    mt = load_genotype_matrixtable(
        bfile=bfile,
        chromosome="22",
        final_samples=["IID1", "IID2", "IID3", "IID4"],
        variant_index=variant_index,
        n_partitions=2,
    )
    mt = mt.annotate_entries(n_alt=mt.GT.n_alt_alleles())

    rows = mt.entries().select("variant_id", "n_alt").collect()
    by_variant: dict[str, dict[str, int | None]] = {}
    for row in rows:
        by_variant.setdefault(row.variant_id, {})[row.s] = row.n_alt
    expected_ids = [v["variant_id"] for v in variant_index["variants"]]
    for variant_idx, vid in enumerate(expected_ids):
        for sample_idx, iid in enumerate(["IID1", "IID2", "IID3", "IID4"]):
            assert by_variant[vid][iid] == calls[variant_idx][sample_idx], (
                f"mismatch at variant {vid} sample {iid}"
            )


def test_read_a_block_chunk_round_trip(hail_session, tmp_path):
    arr = np.array([[5.0, 2.0, 0.0], [0.0, 6.0, 3.0], [0.0, 0.0, 7.0]])
    bm = BlockMatrix.from_numpy(arr, block_size=8)
    out = tmp_path / "chunk"
    bm.write(str(out), overwrite=True)
    back = ld_hail.read_a_block_chunk(out)
    np.testing.assert_array_equal(back, arr)


def test_write_r_block_chunk_sparsifies_row_intervals(hail_session, tmp_path):
    dense = np.array([[1.0, 0.4, 0.2], [0.0, 1.0, 0.5], [0.0, 0.0, 1.0]],
                     dtype=np.float64)
    starts = np.array([0, 1, 2], dtype=np.int64)
    stops = np.array([3, 3, 3], dtype=np.int64)   # full upper triangle within window
    out = tmp_path / "r_chunk"
    ld_hail.write_r_block_chunk(dense, starts, stops, out, block_size=8)
    back = BlockMatrix.read(str(out)).to_numpy()
    # Hail BlockMatrix is float64-only (from_numpy upcasts), so the round-trip
    # dtype is float64 — float32 storage is not achievable via Hail BlockMatrix.
    assert back.dtype == np.float64
    np.testing.assert_allclose(np.triu(back), np.triu(dense), rtol=1e-6)
    # below-diagonal entries were sparsified away
    assert back[1, 0] == 0.0
