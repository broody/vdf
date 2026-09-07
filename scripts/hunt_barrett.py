#!/usr/bin/env python3
"""Targeted hunt for Barrett boundary failures.

Failure mechanism (k-limb truncation of r = p - q3n):
  true residue path needs r = p - q3*N, which lies in [0, 3N).
  Code keeps only k limbs, so values >= b^k wrap. Since q3 <= floor(p/N),
  r = (p mod N) + e*N with e in {0,1,2}. Wrap needs r >= b^k.
"""
import json
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "ref"))
from cairo_model_v2 import barrett_mu, barrett_reduce, from_limbs, to_limbs

v = json.load(open(os.path.join(ROOT, "vectors", "vdf_vectors_512.json")))[0]
kN = v["n_limbs"]
N_int = int(v["N"], 16)
N = v["N_limbs"]
mu_n = barrett_mu(N_int, kN)
B_K = 1 << (64 * kN)
print(f"N = {N_int:#x}")
print(f"N/b^k = {N_int / B_K:.6f}   (failure region needs N > b^k/3 = {1/3:.4f})")


def check(p_int, k, Nl, N_int, mu):
    got = from_limbs(barrett_reduce(to_limbs(p_int, 2 * k), Nl, mu, k))
    want = p_int % N_int
    return got, want


rng = random.Random(0xBA11E17)
fails = []
trials = 200_000
for _ in range(trials):
    # uniform p over the full 2k-limb range maximizes estimate-error odds
    p_int = rng.randrange(0, 1 << (64 * 2 * kN))
    got, want = check(p_int, kN, N, N_int, mu_n)
    if got != want:
        fails.append((p_int, got, want))
        if len(fails) >= 5:
            break

print(f"uniform-p failures: {len(fails)} (stopped after {trials} max trials)")
for p_int, got, want in fails[:3]:
    q, rem = divmod(p_int, N_int)
    # recompute the estimate error e = floor(p/N) - q3
    q1 = p_int >> (64 * (kN - 1))
    q2 = q1 * from_limbs(mu_n)
    q3 = q2 >> (64 * (kN + 1))
    e = q - q3
    print(f"  p={p_int:#x}")
    print(f"    p mod N = {rem:#x}  (p mod N)/N = {rem / N_int:.4f}")
    print(f"    estimate error e = {e}, r_true = p - q3*N = {(rem + e * N_int) / B_K:.4f} b^k")
    print(f"    got  = {got:#x}")
    print(f"    want = {want:#x}")
    print(f"    got == (r_true - b^k)? {got == rem + e * N_int - B_K}")

# adversarial: maximize remainder AND error odds -> p = q*N + (N-1) for random q
adv = 0
for _ in range(200_000):
    q = rng.randrange(0, (1 << (64 * 2 * kN)) // N_int)
    p_int = q * N_int + (N_int - 1)
    if p_int >= 1 << (64 * 2 * kN):
        continue
    got, want = check(p_int, kN, N, N_int, mu_n)
    if got != want:
        adv += 1
        if adv == 1:
            print(f"\nadversarial fail: q={q:#x}")
            print(f"  p    = {p_int:#x}")
            print(f"  got  = {got:#x}")
            print(f"  want = {want:#x}")
        if adv >= 3:
            break
print(f"adversarial failures: {adv}")
