import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
CHECK_UPLOAD = REPO_ROOT / "check_upload.sh"


def test_check_upload_parses():
    r = subprocess.run(["bash", "-n", str(CHECK_UPLOAD)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


def test_section15_block_present_and_calls_helper():
    text = CHECK_UPLOAD.read_text()
    assert '$1 = "15"' in text
    assert "ld_encrypt_cohort.sh" in text
    # uploads both the encrypted archives and their checksums
    assert "${study_name}_15_*.tgz.aes" in text
    assert "${study_name}_15_*.md5sum" in text
