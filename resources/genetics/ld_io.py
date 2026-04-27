import gzip
import json
from pathlib import Path
from typing import Iterable

type JsonValue = (
    str | int | float | bool | None | list[JsonValue] | dict[str, JsonValue]
)


def ensure_dir(path: str | Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_json(path: str | Path, payload: JsonValue) -> None:
    path = Path(path)
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def write_text(path: str | Path, text: str) -> None:
    path = Path(path)
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8") as handle:
        handle.write(text)


def write_gzip_text(path: str | Path, text: str) -> None:
    path = Path(path)
    ensure_dir(path.parent)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(text)


def write_gzip_lines(path: str | Path, lines: Iterable[str]) -> None:
    path = Path(path)
    ensure_dir(path.parent)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for line in lines:
            handle.write(line)
            handle.write("\n")
