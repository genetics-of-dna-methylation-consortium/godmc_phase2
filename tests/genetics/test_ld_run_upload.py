import subprocess
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
WRAPPER = REPO_ROOT / "15c-ld_run_upload.sh"


def test_compatibility_wrapper_parses():
    result = subprocess.run(
        ["bash", "-n", str(WRAPPER)], capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr


def test_compatibility_wrapper_runs_staged_commands_in_order(tmp_path):
    wrapper = tmp_path / WRAPPER.name
    wrapper.write_text(WRAPPER.read_text())
    calls = tmp_path / "calls.txt"

    scripts = {
        "15a-ld_prepare_stats.sh": "15a",
        "15b-ld_compress_data.sh": "15b",
        "check_upload.sh": "upload",
    }
    for filename, label in scripts.items():
        script = tmp_path / filename
        script.write_text(
            "#!/usr/bin/env bash\n"
            f'printf "%s\\n" "{label}:$*" >> "{calls}"\n'
        )
        script.chmod(0o755)

    result = subprocess.run(
        ["bash", str(wrapper), "-c", "test.config"],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert calls.read_text().splitlines() == [
        "15a:-c test.config",
        "15b:-c test.config",
        "upload:15 upload -c test.config",
    ]
    assert "deprecated" in result.stdout.lower()
