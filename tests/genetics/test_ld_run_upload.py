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


PASSPHRASE = "testpass"


def _gpg_wrapper(tmp_path):
    home = tmp_path / "gnupg"
    home.mkdir(mode=0o700)
    w = tmp_path / "gpg_batch.sh"
    w.write_text(
        "#!/usr/bin/env bash\n"
        f'exec gpg --homedir "{home}" --batch --yes '
        f'--pinentry-mode loopback --passphrase "{PASSPHRASE}" "$@"\n'
    )
    w.chmod(0o755)
    return w


def _ship_ok(tmp_path):
    dest = tmp_path / "remote"
    dest.mkdir(exist_ok=True)
    s = tmp_path / "ship_ok.sh"
    s.write_text(f'#!/usr/bin/env bash\ncp "$1" "{dest}/"\n')
    s.chmod(0o755)
    return s, dest


def _make_chunk(ablocks, chrom="chr22", chunk="chunk_0"):
    d = ablocks / chrom / chunk
    d.mkdir(parents=True)
    (d / "part-00000").write_bytes(b"blockmatrix-bytes")
    (d / "metadata.json").write_text('{"block":1}')
    return d


def test_process_chunk_success_ships_and_deletes(tmp_path):
    ablocks = tmp_path / "A_blocks"
    chunk = _make_chunk(ablocks)
    out = tmp_path / "upload"; out.mkdir()
    gpg = _gpg_wrapper(tmp_path)
    ship, dest = _ship_ok(tmp_path)
    r = _call(
        f'process_chunk "{chunk}" "{ablocks}" "{out}" "study_chr22"',
        env={"GPG": str(gpg), "LD_SHIP_CMD": str(ship)},
    )
    assert r.returncode == 0, r.stderr
    base = "study_chr22_15_chr22_chunk_0"
    # shipped to the fake remote
    assert (dest / f"{base}.tgz.aes").is_file()
    assert (dest / f"{base}.md5sum").is_file()
    # local source + staged artifacts removed
    assert not chunk.exists()
    assert list(out.glob("*")) == []


def test_process_chunk_failed_ship_keeps_everything(tmp_path):
    ablocks = tmp_path / "A_blocks"
    chunk = _make_chunk(ablocks)
    out = tmp_path / "upload"; out.mkdir()
    gpg = _gpg_wrapper(tmp_path)
    fail = tmp_path / "ship_fail.sh"
    fail.write_text("#!/usr/bin/env bash\nexit 1\n")
    fail.chmod(0o755)
    r = _call(
        f'process_chunk "{chunk}" "{ablocks}" "{out}" "study_chr22"',
        env={"GPG": str(gpg), "LD_SHIP_CMD": str(fail)},
    )
    assert r.returncode == 1
    # source chunk NOT deleted on a failed upload
    assert chunk.exists()
    assert (chunk / "part-00000").is_file()
