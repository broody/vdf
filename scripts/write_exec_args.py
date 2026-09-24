#!/usr/bin/env python3
"""Regenerate exec args from the vector files (repo-local paths).

Input layout (challenge L and r = 2^T mod L are derived in-program):
  [N(k), mu_n(k+1), x(k), y(k), pi(k), T(2), nonce, ctx, ct_len, ct(ct_len)]
with k = 8 (512-bit) or 32 (2048-bit). The executable is passed the flat
array prefixed by its length (scarb's Span<felt252> argument encoding).
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from gen_cairo_tests import barrett_mu  # noqa: E402


def exec_args(v: dict) -> list[int]:
    T = v["T"]
    ct = [int(c, 16) for c in v["ciphertext"]]
    return (
        v["N_limbs"] + barrett_mu(int(v["N"], 16), v["n_limbs"]) + v["x_limbs"]
        + v["y_limbs"] + v["pi_limbs"] + [T & ((1 << 64) - 1), T >> 64]
        + [v["nonce"], int(v["ctx"], 16), len(ct)] + ct
    )


if __name__ == "__main__":
    for bits in (512, 2048):
        v = json.load(open(os.path.join(ROOT, "vectors", f"vdf_vectors_{bits}.json")))[0]
        flat = exec_args(v)
        out = os.path.join(ROOT, "vectors", f"exec_args_{bits}.json")
        json.dump(["0x%x" % len(flat)] + ["0x%x" % x for x in flat], open(out, "w"))
        expect = [v["instance_hash"], v["ctx"], v["ct_hash"], v["pt_hash"]]
        print(f"{bits} exec args: {len(flat)} felts -> {out}; expected output {expect}")
