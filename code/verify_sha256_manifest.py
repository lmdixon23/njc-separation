#!/usr/bin/env python3
"""Verify the repository's canonical-file SHA-256 manifest."""

from __future__ import annotations

import hashlib
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "SHA256SUMS.txt"
LINE_RE = re.compile(r"^([0-9a-f]{64})  (.+)$")


def fail(message: str) -> None:
    print(f"VERDICT: SHA-256 MANIFEST FAILED: {message}", file=sys.stderr)
    raise SystemExit(1)


def main() -> None:
    if not MANIFEST.is_file():
        fail("SHA256SUMS.txt is missing")

    seen: set[str] = set()
    checked = 0
    for line_number, raw in enumerate(MANIFEST.read_text(encoding="utf-8").splitlines(), 1):
        if not raw or raw.startswith("#"):
            continue
        match = LINE_RE.fullmatch(raw)
        if match is None:
            fail(f"malformed line {line_number}")
        expected, rel_text = match.groups()
        rel = Path(rel_text)
        if rel.is_absolute() or ".." in rel.parts or rel_text.startswith(("/", "\\")):
            fail(f"unsafe path on line {line_number}: {rel_text}")
        normalized = rel.as_posix()
        if normalized in seen:
            fail(f"duplicate path: {normalized}")
        seen.add(normalized)
        path = ROOT / rel
        if not path.is_file():
            fail(f"missing file: {normalized}")
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != expected:
            fail(f"hash mismatch: {normalized}")
        checked += 1

    if checked == 0:
        fail("manifest contains no files")
    print("VERDICT: SHA-256 MANIFEST VERIFIED")
    print(f"files = {checked}")


if __name__ == "__main__":
    main()
