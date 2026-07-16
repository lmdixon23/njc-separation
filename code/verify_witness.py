#!/usr/bin/env python3
"""Deterministically verify the canonical NJC separation witness."""

from __future__ import annotations

import itertools
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def det2(m):
    return m[0][0] * m[1][1] - m[0][1] * m[1][0]


def sigmoid(t: float) -> float:
    if t >= 0:
        e = math.exp(-t)
        return 1.0 / (1.0 + e)
    e = math.exp(t)
    return e / (1.0 + e)


def sigmoid_prime(t: float) -> float:
    s = sigmoid(t)
    return s * (1.0 - s)


def logits(B, x, c):
    return [B[i][0] * x[0] + B[i][1] * x[1] + c[i] for i in range(4)]


def det_adb(A, d, B):
    m = [
        [sum(A[r][i] * d[i] * B[i][k] for i in range(4)) for k in range(2)]
        for r in range(2)
    ]
    return det2(m)


def close(actual: float, expected: float, *, atol: float = 2e-12, rtol: float = 2e-10) -> bool:
    return abs(actual - expected) <= atol + rtol * abs(expected)


def main() -> None:
    witness = json.loads((ROOT / "data/witness.json").read_text(encoding="utf-8"))
    A, B, c = witness["A"], witness["B"], witness["c"]
    p = witness["witness_pair"]["p"]
    q = witness["witness_pair"]["q"]
    claimed = witness["claimed"]

    minors = {}
    for i, j in itertools.combinations(range(4), 2):
        ai = [[A[0][i], A[0][j]], [A[1][i], A[1][j]]]
        bi = [B[i], B[j]]
        minors[f"{i}{j}"] = det2(ai) * det2(bi)
    for key, expected in claimed["minors"].items():
        if not close(minors[key], expected):
            raise AssertionError((key, minors[key], expected))

    lp, lq = logits(B, p, c), logits(B, q, c)
    divided = []
    for a, b in zip(lp, lq):
        if a == b:
            divided.append(sigmoid_prime(a))
        else:
            divided.append((sigmoid(b) - sigmoid(a)) / (b - a))
    avg_det = det_adb(A, divided, B)
    if not close(avg_det, claimed["seg_avg_det"]):
        raise AssertionError(("segment average", avg_det, claimed["seg_avg_det"]))

    segment_values = []
    for k in range(201):
        t = k / 200.0
        x = [(1.0 - t) * p[j] + t * q[j] for j in range(2)]
        d = [sigmoid_prime(z) for z in logits(B, x, c)]
        segment_values.append(det_adb(A, d, B))
    segment_min = min(segment_values)
    if not close(segment_min, claimed["min_detDF_on_segment"]):
        raise AssertionError(("segment minimum", segment_min, claimed["min_detDF_on_segment"]))

    negative = [key for key, value in minors.items() if value < 0]
    if negative != ["23"] or avg_det >= 0 or segment_min <= 0:
        raise AssertionError({"negative_minors": negative, "avg_det": avg_det, "segment_min": segment_min})

    print("VERDICT: CANONICAL WITNESS VERIFIED")
    print("minor_products =", minors)
    print("segment_average_det =", avg_det)
    print("sampled_segment_min_detDF =", segment_min)


if __name__ == "__main__":
    main()
