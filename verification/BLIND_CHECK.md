# Independent headline reconstruction

`BLIND_CHECK.py` reads only `data/witness.json` and independently recomputes:

1. the six mixed-minor products;
2. the determinant of the segment-averaged Jacobian at the canonical pair;
3. the minimum sampled determinant along the displayed segment.

Run from the repository root:

```bash
python BLIND_CHECK.py
```

The script exits nonzero if the recomputed values disagree with the canonical expectations beyond the stated tolerance.
