#!/usr/bin/env python3
"""Independent headline reconstruction using only canonical JSON payloads."""

from __future__ import annotations

import itertools
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def fail(message: str) -> None:
    print(f"BLIND CHECK FAILED: {message}", file=sys.stderr)
    raise SystemExit(1)


def require(condition: bool, message: str) -> None:
    if not condition:
        fail(message)


def det2(m):
    return m[0][0] * m[1][1] - m[0][1] * m[1][0]


def sig(t):
    if t >= 0:
        e = math.exp(-t)
        return 1.0 / (1.0 + e)
    e = math.exp(t)
    return e / (1.0 + e)


def sp(t):
    s = sig(t)
    return s * (1.0 - s)


def matvec(B, x, c):
    return [B[i][0] * x[0] + B[i][1] * x[1] + c[i] for i in range(4)]


def det_adb(A, d, B):
    m = [
        [
            sum(A[r][i] * d[i] * B[i][k] for i in range(4))
            for k in range(2)
        ]
        for r in range(2)
    ]
    return det2(m)


def close(actual, expected, atol=2e-12, rtol=2e-10):
    return abs(actual - expected) <= atol + rtol * abs(expected)


def main() -> None:
    with (ROOT / "data/witness.json").open(encoding="utf-8") as handle:
        w = json.load(handle)
    with (ROOT / "data/expected_values.json").open(encoding="utf-8") as handle:
        expected = json.load(handle)

    A, B, c = w["A"], w["B"], w["c"]
    p, q = w["witness_pair"]["p"], w["witness_pair"]["q"]

    minors = {}
    for i, j in itertools.combinations(range(4), 2):
        ai = [[A[0][i], A[0][j]], [A[1][i], A[1][j]]]
        bi = [B[i], B[j]]
        minors[f"{i}{j}"] = det2(ai) * det2(bi)

    negative = [(key, value) for key, value in minors.items() if value < 0]
    expected_negative = expected["item1_negative_minor"]

    require(len(negative) == 1, f"expected one negative minor, found {negative}")
    expected_key = "".join(str(i) for i in expected_negative["I"])
    require(
        negative[0][0] == expected_key,
        f"negative minor key {negative[0][0]} != {expected_key}",
    )
    require(
        close(negative[0][1], expected_negative["value"]),
        f"negative minor value {negative[0][1]} != {expected_negative['value']}",
    )

    lp, lq = matvec(B, p, c), matvec(B, q, c)
    divided = [
        (sig(b) - sig(a)) / (b - a) if a != b else sp(a)
        for a, b in zip(lp, lq)
    ]
    avg_det = det_adb(A, divided, B)
    require(
        close(avg_det, expected["item2_seg_avg_det"]),
        f"segment-average determinant {avg_det} != {expected['item2_seg_avg_det']}",
    )

    values = []
    for k in range(201):
        t = k / 200.0
        x = [(1.0 - t) * p[j] + t * q[j] for j in range(2)]
        values.append(det_adb(A, [sp(z) for z in matvec(B, x, c)], B))
    segment_min = min(values)

    require(
        close(segment_min, expected["item3_min_detDF_on_segment"]),
        f"sampled segment minimum {segment_min} != "
        f"{expected['item3_min_detDF_on_segment']}",
    )
    require(
        avg_det < 0 < segment_min,
        f"required sign separation failed: avg_det={avg_det}, "
        f"segment_min={segment_min}",
    )

    print("minor_products =", minors)
    print("negative_minor =", negative)
    print("segment_average_det =", avg_det)
    print("sampled_segment_min_detDF =", segment_min)
    print("BLIND CHECK PASSED")


if __name__ == "__main__":
    main()
