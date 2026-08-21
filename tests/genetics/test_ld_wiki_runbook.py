from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
WIKI = REPO_ROOT / "godmc_phase2.wiki" / "Run-federated-LD-reference-panel.md"


def test_runbook_documents_staged_workflow_and_is_current():
    text = WIKI.read_text()
    assert "15a-ld_prepare_stats.sh" in text
    assert "15b-ld_compress_data.sh" in text
    assert "check_upload.sh 15 upload" in text
    # documents the disk expectation that motivates the design
    assert "150 GB" in text
    # documents the Imperial endpoint config the operator must set
    assert "ld_upload_password_file" in text
    # documents resumability
    assert ".prepared" in text and ".packaged" in text
    # stale status removed
    assert "not yet implemented" not in text
