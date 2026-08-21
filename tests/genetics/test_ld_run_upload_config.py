from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def test_config_example_has_upload_password_file():
    text = (REPO_ROOT / "config.example").read_text()
    assert "ld_upload_password_file=" in text
    assert "GODMC_UPLOAD_PASSWORD" not in text


def test_parameters_has_section_15_upload_endpoint():
    text = (REPO_ROOT / "resources" / "parameters").read_text()
    assert "godmc_upload_url=" in text
