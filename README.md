# Averaged-Jacobian failure at the first open width of the Neural Jacobian Conjecture

![verify](https://github.com/lmdixon23/njc-separation/actions/workflows/verify.yml/badge.svg)

Code, canonical data, and reproducibility material for the accompanying paper.

## Main result

The repository supports an explicit sigmoid ridge network at `n=2`, `N=4` with an everywhere-positive Jacobian determinant and a segment-averaged Jacobian with negative determinant. The example separates pointwise Jacobian positivity from the averaged-Jacobian univalence mechanism that proves the smaller-width cases. This paper does not claim that the witness is injective; that question remains open here.

## Reproduce the maintained claims

Use Python 3.13 and the locked package versions:

```bash
python -m pip install -r requirements-lock.txt
bash run_all.sh
```

Run from the repository root. The final expected line is:

```text
VERDICT: ALL MAINTAINED CHECKS PASSED
```

The default suite is deterministic. It verifies `SHA256SUMS.txt`, recomputes the witness properties, runs the independent reconstruction, checks the proof-bearing positivity certificate, certifies the segment bound and sector exclusions, and exercises a fast fail-closed mutation suite.

The stochastic frontier search is retained for provenance but is not part of default CI. To regenerate it separately, run `python code/njc_v3.py`. Its search path and stopping time can vary across BLAS and library versions even with pinned seeds; the shipped `data/e3.pkl`, witness, and certificate are the artifacts used by the paper.

## Claim-to-script map

| Paper claim | Primary implementation | Evidence |
|---|---|---|
| Witness matrices, mixed minors, and separation pair | `code/verify_witness.py` | deterministic recomputation from canonical data |
| Global positivity of the witness Jacobian | `code/emit_certificate.py`, `CHECK.py` | structure-only certificate checked with outward-rounded and exact arithmetic |
| Floating-point positivity reconstruction | `code/certify.py` | corroborating numerical reconstruction |
| Certified positivity along the displayed segment | `code/segment_bound.py` | outward-rounded one-dimensional subdivision |
| Sector exclusion for possible collisions | `code/sector_lemma.py` | finite sign-pattern enumeration and explicit dual certificates |
| Recorded frontier-study counts | `code/njc_v3.py`, `data/e3.pkl` | platform-dependent search record, not an enumerative theorem |
| Independent headline recomputation | `BLIND_CHECK.py` | fail-closed implementation using only the canonical JSON payloads |

The detailed public crosswalk is `verification/CLAIM-TO-ARTIFACT-MAP.md`.

## Certificate boundary

`data/certificate.json` is the proof-bearing certificate for the compact and tail positivity argument. `CHECK.py` recomputes the required bounds and partition conditions rather than trusting verdict fields. `tests/run_ci_negative_smoke.py` provides fast fail-closed mutation checks for every push. The complete 27-case adversarial suite remains available as `python tests/run_negative_suite.py --require-complete` for release audits; it is intentionally excluded from default CI because several rounding knife-edge constructions are computationally expensive.

## Independent reconstruction and integrity

`BLIND_CHECK.py` recomputes the headline witness quantities without importing the primary implementation under `code/`. `SHA256SUMS.txt` binds the canonical public source, code, data, tests, and verification records; `code/verify_sha256_manifest.py` fails on missing, altered, duplicated, or unsafe entries.

## Scope

- The separation theorem and positivity certificate are maintained claims.
- The frontier-study counts describe the shipped search run and may not reproduce identically on another platform.
- The sector result gives necessary conditions for a collision; it does not prove or disprove injectivity of the witness.

## Repository map

- `paper/`: canonical manuscript source
- `code/`: search, reconstruction, and certificate scripts
- `data/`: canonical witness, study data, and proof-bearing certificate
- `tests/`: fail-closed negative tests
- `verification/`: concise public status and evidence documentation

## Cite

See `CITATION.cff`.
