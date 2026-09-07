#!/usr/bin/env python3
"""Shared helpers + generators for Cairo tests and executable args.

`barrett_mu` is imported by the other scripts; `limbs_cairo` formats u64 limb
arrays for Cairo source. Running this file prints a full-verify Cairo test per
vector (the canonical regenerator for the checked-in tests.cairo is
regen_clean_tests.py).
"""

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
VECTORS = json.load(open(os.path.join(ROOT, "vectors", "vdf_vectors_512.json")))


def limbs_cairo(limbs: list[int]) -> str:
    inner = ", ".join(f"{x}" for x in limbs)
    return f"array![{inner}].span()"


def barrett_mu(N_int: int, k: int) -> list[int]:
    b = 1 << 64
    n = k + 1
    mask = (1 << 64) - 1
    mu = (b ** (2 * k)) // N_int
    return [(mu >> (64 * i)) & mask for i in range(n)]


def gen_tests():
    print("use vdf::{verify_vdf, N_LIMBS, L_LIMBS};")
    print()
    print("#[cfg(test)]")
    print("mod tests {")
    print("    use super::*;")
    print()
    for v in VECTORS:
        name = v["name"]
        N_int = int(v["N"], 16)
        L_int = int(v["L"], 16)
        muN = barrett_mu(N_int, v["n_limbs"])
        muL = barrett_mu(L_int, 4)
        T = v["T"]
        T_limbs = [T & ((1 << 64) - 1), (T >> 64) & ((1 << 64) - 1)]

        print("    #[test]")
        print(f"    fn {name}_wesolowski_verifies() {{")
        print(f"        let N = {limbs_cairo(v['N_limbs'])};")
        print(f"        let mu_n = {limbs_cairo(muN)};")
        print(f"        let x = {limbs_cairo(v['x_limbs'])};")
        print(f"        let y = {limbs_cairo(v['y_limbs'])};")
        print(f"        let pi = {limbs_cairo(v['pi_limbs'])};")
        print(f"        let mu_l = {limbs_cairo(muL)};")
        print(f"        let T_limbs = {limbs_cairo(T_limbs)};")
        print(
            "        assert(verify_vdf(N, mu_n, x, y, pi, mu_l, T_limbs, 2), 'VDF must verify');"
        )
        print("    }")
        print()
    print("}")


if __name__ == "__main__":
    gen_tests()
