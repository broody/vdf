#!/usr/bin/env python3
"""Model v2: Cairo-ready bigint with u64 limbs + u256 accumulation.

Mirrors the rewrite that avoids felt252 ops Cairo 2.13 does not support
(no PartialOrd, no Rem/Div, no shifts on felt252). All arithmetic lives in
u64 (limbs) and u256 (accumulation/carry/subtraction), which Cairo supports
natively. Validated against the same Wesolowski vectors.
"""

import json

LIMB_BITS = 64
MASK = (1 << LIMB_BITS) - 1
TWO64 = 1 << LIMB_BITS


def to_limbs(v: int, n: int) -> list[int]:
    return [(v >> (LIMB_BITS * i)) & MASK for i in range(n)]


def from_limbs(limbs: list[int]) -> int:
    v = 0
    for i, limb in enumerate(limbs):
        v |= limb << (LIMB_BITS * i)
    return v


def bigint_mul(a: list[int], b: list[int], an: int, bn: int,
               out_limbs: int | None = None) -> list[int]:
    """Schoolbook product, output-limb-at-a-time with chained carry.

    out[t] = (sum_{i+j=t} a[i]*b[j]) + carry[t-1]
    limb   = out[t] mod 2^64
    carry  = out[t] / 2^64
    Only contributing (i, j) pairs are visited: i in [max(0, t-bn+1), min(an-1, t)].
    `out_limbs` truncates the product to the low `out_limbs` columns (mod b^out_limbs);
    the carry out of the last kept column is dropped, which is exact for callers that
    only consume low limbs (Barrett's q3n). The Cairo version mirrors this loop shape
    and accumulates each column in a u128 (hi, lo) pair instead of u256 — same math.
    """
    total = an + bn if out_limbs is None else min(an + bn, out_limbs)
    out = []
    carry = 0
    for t in range(total):
        acc = carry
        lo_i = max(0, t - (bn - 1))
        hi_i = min(an - 1, t)
        for i in range(lo_i, hi_i + 1):
            acc += a[i] * b[t - i]
        out.append(acc & MASK)
        carry = acc >> LIMB_BITS
    return out


def sub_limbs(a: list[int], b: list[int], k: int) -> list[int]:
    """Limb-wise a - b with borrow, u256-style minuend/subtrahend check."""
    out = []
    borrow = 0
    for i in range(k):
        minuend = a[i]
        subtrahend = b[i] + borrow
        if minuend >= subtrahend:
            out.append(minuend - subtrahend)
            borrow = 0
        else:
            out.append(minuend + TWO64 - subtrahend)
            borrow = 1
    return out


def ge_limbs(a: list[int], b: list[int], k: int) -> bool:
    for i in range(k - 1, -1, -1):
        if a[i] != b[i]:
            return a[i] > b[i]
    return True  # equal


def barrett_reduce(p: list[int], N: list[int], mu: list[int], k: int) -> list[int]:
    """Barrett reduction, HAC 14.42: r keeps k+1 limbs through normalization.

    Requires N's top limb nonzero (b^(k-1) <= N < b^k) and p < b^2k. Then the
    quotient estimate satisfies floor(p/N) - 2 <= q3 <= floor(p/N), so
    0 <= p - q3*N < 3N < b^(k+1) and r fits k+1 limbs with NO truncation.
    (The previous version kept only k limbs: whenever p - q3*N >= b^k the
    residue wrapped mod b^k and the two conditional subtractions returned a
    wrong result — e.g. N = 2^512-1, a = 2^448-1, b = 2^448+1 was off by one.)
    """
    q1 = p[k - 1:]  # k+1 limbs
    q2 = bigint_mul(q1, mu, k + 1, k + 1)  # 2k+2 limbs
    q3 = q2[k + 1:]  # k+1 limbs
    q3n = bigint_mul(q3, N, k + 1, k, out_limbs=k + 1)  # low k+1 limbs only

    # r = (p - q3n) mod b^(k+1) with borrow (exact: 0 <= p - q3n < b^(k+1))
    r = []
    borrow = 0
    for i in range(k + 1):
        minuend = p[i]
        subtrahend = q3n[i] + borrow
        if minuend >= subtrahend:
            r.append(minuend - subtrahend)
            borrow = 0
        else:
            r.append(minuend + TWO64 - subtrahend)
            borrow = 1

    # normalize: while r >= N, r -= N (at most 2), over k+1 limbs
    N_ext = N + [0]
    while ge_limbs(r, N_ext, k + 1):
        r = sub_limbs(r, N_ext, k + 1)
    assert r[k] == 0, "normalized residue must fit k limbs"
    return r[:k]


def modmul(a: list[int], b: list[int], N: list[int], mu: list[int], k: int) -> list[int]:
    return barrett_reduce(bigint_mul(a, b, k, k), N, mu, k)


def modpow(base: list[int], exp: list[int], exp_n: int, N: list[int], mu: list[int], k: int) -> list[int]:
    one = [1] + [0] * (k - 1)
    result = barrett_reduce(bigint_mul(one, one, k, k), N, mu, k)
    base_red = modmul(base, one, N, mu, k)
    for i in range(exp_n):
        e = exp[i]
        for bit in range(LIMB_BITS):
            if (e >> bit) & 1:
                result = modmul(result, base_red, N, mu, k)
            base_red = modmul(base_red, base_red, N, mu, k)
    return result


def limbs_eq(a: list[int], b: list[int], k: int) -> bool:
    return a[:k] == b[:k]


def barrett_mu(N_int: int, k: int) -> list[int]:
    return to_limbs((1 << (64 * 2 * k)) // N_int, k + 1)


def mu_valid(m_int: int, mu: list[int], k: int) -> bool:
    """Mirror of Cairo mu_valid: mu == floor(b^(2k)/m)."""
    mu_int = from_limbs(mu)
    b2k = 1 << (64 * 2 * k)
    return mu_int * m_int <= b2k < (mu_int + 1) * m_int


def verify_vector(v: dict) -> bool:
    """Mirror of Cairo verify_vdf: challenge derived in-program, mu validated."""
    from vdf_reference import derive_challenge

    kN = v["n_limbs"]
    N = v["N_limbs"]
    x = v["x_limbs"]
    y = v["y_limbs"]
    pi = v["pi_limbs"]
    T = v["T"]
    N_int = int(v["N"], 16)
    kL = 4

    L_int = derive_challenge(N_int, int(v["x"], 16), int(v["y"], 16), T, kN)
    assert L_int == int(v["L"], 16), "stored L must match the derived challenge"
    L = to_limbs(L_int, kL)

    muN = barrett_mu(N_int, kN)
    muL = barrett_mu(L_int, kL)
    assert mu_valid(N_int, muN, kN)
    assert mu_valid(L_int, muL, kL)

    r_calc = modpow(to_limbs(2, kL), to_limbs(T, 2), 2, L, muL, kL)
    assert from_limbs(r_calc) == int(v["r"], 16), f"r mismatch: {from_limbs(r_calc):x}"

    piL = modpow(pi, L, kL, N, muN, kN)
    xr = modpow(x, r_calc, kL, N, muN, kN)
    lhs = modmul(piL, xr, N, muN, kN)
    return from_limbs(lhs) == int(v["y"], 16)


if __name__ == "__main__":
    import os
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    vectors = json.load(open(os.path.join(root, "vectors", "vdf_vectors_512.json")))
    for v in vectors:
        print(f"{v['name']}: {'PASS' if verify_vector(v) else 'FAIL'}")
