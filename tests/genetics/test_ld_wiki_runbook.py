from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
WIKI = REPO_ROOT / "godmc_phase2.wiki" / "Run-federated-LD-reference-panel.md"


def test_runbook_documents_wrapper_and_is_current():
    text = WIKI.read_text()
    # documents the one cohort command
    assert "15c-ld_run_upload.sh" in text
    # documents the disk expectation that motivates the design
    assert "150 GB" in text
    # documents the Imperial endpoint config the operator must set
    assert "imperial_user" in text and "imperial_key" in text
    # documents resumability
    assert ".uploaded" in text
    # stale status removed
    assert "not yet implemented" not in text
