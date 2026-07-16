#!/usr/bin/env python3
"""Fast fail-closed smoke tests for CHECK.py.

The complete adversarial suite remains in run_negative_suite.py and is intended
for release audits, not every CI push.
"""

from __future__ import annotations

import copy
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WITNESS = json.loads((ROOT / "data/witness.json").read_text(encoding="utf-8"))
CERT = json.loads((ROOT / "data/certificate.json").read_text(encoding="utf-8"))


def run_raw(witness_text: str, certificate_text: str) -> tuple[int, str]:
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        (work / "data").mkdir()
        (work / "data/witness.json").write_text(witness_text, encoding="utf-8")
        (work / "data/certificate.json").write_text(certificate_text, encoding="utf-8")
        shutil.copy2(ROOT / "CHECK.py", work / "CHECK.py")
        result = subprocess.run(
            [sys.executable, "CHECK.py"],
            cwd=work,
            capture_output=True,
            text=True,
            timeout=180,
            check=False,
        )
        return result.returncode, result.stdout + result.stderr


def require(name: str, condition: bool, detail: str) -> None:
    if not condition:
        raise AssertionError(f"{name}: {detail}")
    print(f"{name}: PASS")


def main() -> None:
    code, out = run_raw(json.dumps(WITNESS), "{")
    require("malformed-certificate-json", code == 2 and "EXIT2" in out, f"exit={code}")

    cert = copy.deepcopy(CERT)
    del cert["bnb"]
    code, out = run_raw(json.dumps(WITNESS), json.dumps(cert))
    require("missing-required-certificate-field", code == 2 and "schema" in out, f"exit={code}")

    witness = copy.deepcopy(WITNESS)
    witness["A"][0][0] = witness["A"][1][0] = 0.0
    code, out = run_raw(json.dumps(witness), json.dumps(CERT))
    require("corrupted-minor-structure", code == 1 and "O1" in out, f"exit={code}")

    witness = copy.deepcopy(WITNESS)
    witness["witness_pair"]["q"] = copy.deepcopy(witness["witness_pair"]["p"])
    code, out = run_raw(json.dumps(witness), json.dumps(CERT))
    require("zero-divided-difference", code == 2 and "EXIT2" in out, f"exit={code}")

    print("VERDICT: FAIL-CLOSED CI SMOKE TESTS PASSED")


if __name__ == "__main__":
    main()
