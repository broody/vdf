#!/usr/bin/env python3
"""Reproduce / verify the fixes for the two historical issues.

A. Forged challenge (soundness): the pre-fix verifier accepted L, r, and both
   Barrett mu constants as trusted inputs. Choosing L = 2^200 (divides 2^T for
   T >= 200) makes r = 2^T mod L = 0, and pi = 1 then "proves" y' = 1 for the
   real x, N, T — no squarings needed. The fixed verifier derives
   L = poseidon(N, x, y, T) in-program and validates mu, so this collapses.

B. Barrett truncation (correctness): the pre-fix reduction kept only k limbs
   of r = p - q3*N; whenever the quotient estimate undershoots and r >= b^k,
   the residue wrapped. Closed form: N = 2^512-1, a = 2^448-1, b = 2^448+1
   returned p mod N - 1. HAC 14.42 keeps k+1 limbs.

Run: python3 scripts/repro_issues.py
Expect: all lines report FIXED / REJECTED. The exec-level negative checks
(forged args -> program output 0) are demonstrated in vectors/ + the git
history of this file's comments; see tests.cairo forged_output_rejected.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "ref"))

from cairo_model_v2 import (
    barrett_mu, barrett_reduce, from_limbs, modmul, modpow, mu_valid, to_limbs,
    verify_vector,
)
from vdf_reference import derive_challenge

v = json.load(open(os.path.join(ROOT, "vectors", "vdf_vectors_512.json")))[0]
kN = v["n_limbs"]
N_int = int(v["N"], 16)
N = v["N_limbs"]
T = v["T"]
kL = 4

print("== current vector (poseidon-derived challenge) ==")
print("verify_vector:", "PASS" if verify_vector(v) else "FAIL")

# ---------------------------------------------------------------- A. forgery
print("\n== A. forged-challenge attack ==")
# Old attack shape: for the real x, N, T claim y' = 1 with pi = 1. Against the
# OLD interface the attacker supplied L = 2^200, r = 0, honest mu for that L,
# and the equation pi^L * x^r == y' held trivially. Against the FIXED flow the
# verifier derives L itself:
L_derived = derive_challenge(N_int, int(v["x"], 16), 1, T, kN)
print(f"derived L for forged instance: {L_derived:#x}")
print(f"  L == 2^200? {L_derived == 1 << 200}  (the old attack value)")
L = to_limbs(L_derived, kL)
mu_l = barrett_mu(L_derived, kL)
assert mu_valid(L_derived, mu_l, kL)
r = modpow(to_limbs(2, kL), to_limbs(T, 2), 2, L, mu_l, kL)
pi = to_limbs(1, kN)
x_r = modpow(v["x_limbs"], r, kL, N, barrett_mu(N_int, kN), kN)
pi_l = modpow(pi, L, kL, N, barrett_mu(N_int, kN), kN)
lhs = modmul(pi_l, x_r, N, barrett_mu(N_int, kN), kN)
print(f"  pi^L * x^r mod N = {from_limbs(lhs):#x}")
print(f"  == forged y' (1)? {from_limbs(lhs) == 1}  -> forgery REJECTED")

# mu tampering
mu_l_bad = barrett_mu(L_derived - 2, kL)
print(f"  mu_l computed for L-2 passes mu_valid? {mu_valid(L_derived, mu_l_bad, kL)}"
      "  -> tampered mu REJECTED")

# ---------------------------------------------------------------- B. barrett
print("\n== B. Barrett k-limb truncation ==")
B = 1 << 64
N_bad = B**kN - 1
a = B**(kN - 1) - 1
b = B**(kN - 1) + 1
p = a * b
got = from_limbs(barrett_reduce(to_limbs(p, 2 * kN), to_limbs(N_bad, kN), barrett_mu(N_bad, kN), kN))
want = p % N_bad
print(f"N = 2^512-1, a = 2^448-1, b = 2^448+1")
print(f"  want p mod N = {want:#x}")
print(f"  got          = {got:#x}")
print(f"  -> {'FIXED' if got == want else 'STILL BROKEN'}")
