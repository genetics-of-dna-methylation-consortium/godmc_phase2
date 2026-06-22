from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def test_config_example_has_imperial_creds():
    text = (REPO_ROOT / "config.example").read_text()
    assert "imperial_user=" in text
    assert "imperial_key=" in text


def test_parameters_has_imperial_endpoint_and_chromosomes():
    text = (REPO_ROOT / "resources" / "parameters").read_text()
    assert "imperial_host=" in text
    assert "imperial_path=" in text
    # default autosome list, env-overridable
    assert 'ld_chromosomes="${ld_chromosomes:-' in text
    # first and last autosome present in the default
    assert "1 2 3" in text and "22" in text
