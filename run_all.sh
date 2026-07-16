#!/usr/bin/env bash
# Reproduce the maintained deterministic claims in mathematical order.
set -euo pipefail

PY=python3
"$PY" -c "" 2>/dev/null || PY=python
"$PY" -c "" 2>/dev/null || { echo "no working Python interpreter found"; exit 1; }

printf '%s\n' '== Canonical-file integrity =='
"$PY" code/verify_sha256_manifest.py

printf '%s\n' '== Witness matrices, mixed minors, and separation pair =='
"$PY" code/verify_witness.py

printf '%s\n' '== Independent headline reconstruction =='
"$PY" BLIND_CHECK.py

printf '%s\n' '== Floating-point positivity reconstruction =='
"$PY" code/certify.py

printf '%s\n' '== Proof-bearing positivity certificate =='
CERT_BACKUP="$(mktemp "${TMPDIR:-/tmp}/njc-certificate.XXXXXX")"
cp data/certificate.json "$CERT_BACKUP"

restore_certificate() {
  cp "$CERT_BACKUP" data/certificate.json
  rm -f "$CERT_BACKUP"
}

trap restore_certificate EXIT

"$PY" code/emit_certificate.py
"$PY" CHECK.py

restore_certificate
trap - EXIT

printf '%s\n' '== Certified segment lower bound =='
"$PY" code/segment_bound.py

printf '%s\n' '== Sector exclusion certificates =='
"$PY" code/sector_lemma.py

printf '%s\n' '== Fail-closed CI smoke tests =='
"$PY" tests/run_ci_negative_smoke.py

printf '%s\n' '== Final canonical-file integrity =='
"$PY" code/verify_sha256_manifest.py

printf '%s\n' 'VERDICT: ALL MAINTAINED CHECKS PASSED'
