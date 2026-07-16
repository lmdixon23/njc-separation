# Public claim status

This file records the mathematical and computational scope without review chronology or internal workflow history.

| Claim | Public status | Scope or caveat |
|---|---|---|
| Explicit mixed-minor witness | Established | Canonical matrices and pair are stored in `data/witness.json` |
| `det DF>0` on all of `R^2` | Certified computationally | Compact branch-and-bound plus directional tail certificate checked by `CHECK.py` |
| Negative determinant of a segment-averaged Jacobian | Established | Direct recomputation from the canonical pair |
| Singular intermediate averaged Jacobian | Established | Continuity consequence of endpoint signs |
| Positive determinant along the displayed segment | Certified computationally | Outward-rounded one-dimensional subdivision |
| Sector exclusion for collision directions | Established | Finite sign-pattern certificates |
| Collision magnitude bounds | Established as necessary conditions | Do not close injectivity |
| Injectivity of the witness | Open in this paper | No claim is made either way |
| Frontier-study counts | Recorded computational result | Search regeneration is platform-dependent |

The scripts named in `CLAIM-TO-ARTIFACT-MAP.md` reproduce the maintained finite claims.
