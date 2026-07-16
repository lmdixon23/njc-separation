# Claim-to-artifact map

| Paper item | Primary artifact | Secondary check |
|---|---|---|
| Witness matrices and mixed minors | `data/witness.json`, `code/verify_witness.py` | `BLIND_CHECK.py` |
| Negative segment-average determinant | `code/verify_witness.py` | `BLIND_CHECK.py` |
| Global determinant positivity | `data/certificate.json`, `CHECK.py` | `code/certify.py` |
| Certificate construction | `code/emit_certificate.py` | schema and arithmetic checks in `CHECK.py` |
| Fail-closed certificate behavior | `tests/run_ci_negative_smoke.py` | fast malformed-input and corruption checks; full release suite in `tests/run_negative_suite.py` |
| Segment determinant lower bound | `code/segment_bound.py` | sampled values in `code/verify_witness.py` |
| Sector exclusion | `code/sector_lemma.py` | explicit matrix-vector certificates printed by the script |
| Frontier-study record | `data/e3.pkl`, `code/njc_v3.py` | optional exploratory `code/verify_frontier.py` |

The public artifacts support the scoped claims listed in `STATUS.md`; they do not establish injectivity of the witness.
