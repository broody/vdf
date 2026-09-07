#!/usr/bin/env python3
"""Starknet's Poseidon permutation + sponge, pure Python.

Ported from the Cairo corelib semantics (core::poseidon, v2.18.0) and the
cairo-lang reference implementation:
  - field p = 2^251 + 17*2^192 + 1, state m = 3 (rate 2, capacity 1)
  - R_F = 8 full rounds (4 + 4), R_P = 83 partial rounds, cube S-box
  - round constants: sha256(f"Hades{idx}") mod p, row-major 91 x 3
  - MDS = [[3, 1, 1], [1, -1, 1], [1, 1, -2]]
  - partial rounds apply the S-box to the LAST state element only

Sponge (matches core::poseidon::PoseidonTrait / poseidon_hash_span):
  state (0, 0, 0); absorb pairs via hades(s0+x, s1+y, s2); a trailing odd
  element x leaves state (s0+x, s1, s2) and finalize pads s1 += 1; an even
  count (including empty) finalizes with s0 += 1; digest is the first element.

Cross-checked against the corelib doc vector:
  poseidon_hash_span([1, 2]) == 0x0371cb6995ea5e7effcd2e174de264b5b407027a75a231a70c2c8d196107f0e7
"""

import hashlib

P = 2**251 + 17 * 2**192 + 1
R_F = 8
R_P = 83
M = 3


def _round_constant(idx: int) -> int:
    return int(hashlib.sha256(f"Hades{idx}".encode()).hexdigest(), 16) % P


ARK = [[_round_constant(M * i + j) for j in range(M)] for i in range(R_F + R_P)]


def _mds(v: list[int]) -> list[int]:
    return [
        (3 * v[0] + v[1] + v[2]) % P,
        (v[0] - v[1] + v[2]) % P,
        (v[0] + v[1] - 2 * v[2]) % P,
    ]


def hades_permutation(s0: int, s1: int, s2: int) -> tuple[int, int, int]:
    v = [s0, s1, s2]
    n = R_F + R_P
    for r in range(n):
        full = r < R_F // 2 or r >= R_F // 2 + R_P
        v = [(v[i] + ARK[r][i]) % P for i in range(M)]
        if full:
            v = [pow(x, 3, P) for x in v]
        else:
            v[2] = pow(v[2], 3, P)
        v = _mds(v)
    return v[0], v[1], v[2]


def poseidon_hash_span(values: list[int]) -> int:
    """core::poseidon::poseidon_hash_span over an arbitrary-length input."""
    s0, s1, s2 = 0, 0, 0
    i = 0
    while i + 1 < len(values) + 1 and i + 1 <= len(values) - 1:
        s0, s1, s2 = hades_permutation(s0 + values[i], s1 + values[i + 1], s2)
        i += 2
    if i < len(values):  # one trailing element
        s0 += values[i]
        s1 += 1
    else:
        s0 += 1
    return hades_permutation(s0, s1, s2)[0]


if __name__ == "__main__":
    got = poseidon_hash_span([1, 2])
    want = 0x0371CB6995EA5E7EFFCD2E174DE264B5B407027A75A231A70C2C8D196107F0E7
    print(f"poseidon_hash_span([1,2]) = {got:#x}")
    assert got == want, "MISMATCH vs corelib doc vector"
    print("matches corelib doc vector")
