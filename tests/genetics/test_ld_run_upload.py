import os
import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "15c-ld_run_upload.sh"


def _call(func_call, env=None, extra=""):
    """Source the wrapper (main does not run) and invoke one shell snippet."""
    full_env = dict(os.environ)
    if env:
        full_env.update(env)
    return subprocess.run(
        ["bash", "-c", f'source "{SCRIPT}"; {extra} {func_call}'],
        capture_output=True, text=True, env=full_env,
    )


def test_parses():
    r = subprocess.run(["bash", "-n", str(SCRIPT)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


def test_chr_done_detects_sentinel(tmp_path):
    outdir = tmp_path / "chr22"
    outdir.mkdir()
    # no sentinel yet
    assert _call(f'chr_done "{outdir}"').returncode != 0
    (outdir / ".uploaded").touch()
    assert _call(f'chr_done "{outdir}"').returncode == 0


def test_ship_file_uses_override_and_propagates_exit(tmp_path):
    f = tmp_path / "artifact.aes"
    f.write_text("payload")
    dest = tmp_path / "remote"
    dest.mkdir()
    ok = tmp_path / "ship_ok.sh"
    ok.write_text(f'#!/usr/bin/env bash\ncp "$1" "{dest}/"\n')
    ok.chmod(0o755)
    r = _call(f'ship_file "{f}"', env={"LD_SHIP_CMD": str(ok)})
    assert r.returncode == 0, r.stderr
    assert (dest / "artifact.aes").read_text() == "payload"

    fail = tmp_path / "ship_fail.sh"
    fail.write_text("#!/usr/bin/env bash\nexit 7\n")
    fail.chmod(0o755)
    assert _call(f'ship_file "{f}"', env={"LD_SHIP_CMD": str(fail)}).returncode == 7
