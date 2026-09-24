#!/usr/bin/env python3
"""Reproduce / verify the fixes for the historical issues.

A. Forged challenge (soundness): the pre-fix verifier accepted L, r, and both
   Barrett mu constants as trusted inputs. Choosing L = 2^200 (divides 2^T for
   T >= 200) makes r = 2^T mod L = 0, and pi = 1 then "proves" y' = 1 for the
   real x, N, T — no squarings needed. The fixed verifier derives
   L = poseidon(N, x, y, T) in-program and validates mu, so this collapses.

B. Barrett truncation (correctness): the pre-fix reduction kept only k limbs
   of r = p - q3*N; whenever the quotient estimate undershoots and r >= b^k,
   the residue wrapped. Closed form: N = 2^512-1, a = 2^448-1, b = 2^448+1
   returned p mod N - 1. HAC 14.42 keeps k+1 limbs.

C. -1 malleability (soundness): the verifier worked in the full group Z_N^*,
   so from an honest y anyone could "prove" N - y (negate pi for an odd
   challenge). Fixed by the signed group: only y < N - y is accepted and the
   equation is checked up to sign.

D. Composite challenge (soundness): with L a random integer, a forger picks
   alpha = 2^T mod M for a smooth M (cheap by CRT), grinds y = x^(alpha + kM)
   until L = H(..y..) divides M, and gets pi = x^((alpha - r)/L) for a wrong y
   -- no factoring of N. Demonstrated here with the challenge truncated to 48
   bits; fixed by hash-to-prime (Miller-Rabin in-program).

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
from vdf_reference import (
    derive_challenge, find_challenge, is_prime_fs, poseidon_hash_span, wesolowski_verify,
)

v = json.load(open(os.path.join(ROOT, "vectors", "vdf_vectors_512.json")))[0]
kN = v["n_limbs"]
N_int = int(v["N"], 16)
x_int = int(v["x"], 16)
y_int = int(v["y"], 16)
N = v["N_limbs"]
T = v["T"]

print("== current vector (prime challenge, signed group) ==")
print("verify_vector:", "PASS" if verify_vector(v) else "FAIL")

# ---------------------------------------------------------------- A. forgery
print("\n== A. forged-challenge attack ==")
# Old attack shape: for the real x, N, T claim y' = 1 with pi = 1. Against the
# OLD interface the attacker supplied L = 2^200, r = 0, honest mu for that L,
# and the equation pi^L * x^r == y' held trivially. The fixed verifier derives
# a prime L itself:
nonce, L_derived = find_challenge(N_int, x_int, 1, T, kN)
print(f"derived L for forged instance (nonce {nonce}): {L_derived:#x}")
print(f"  L == 2^200? {L_derived == 1 << 200}  (the old attack value)")
r = to_limbs(pow(2, T, L_derived), 4)
pi = to_limbs(1, kN)
x_r = modpow(v["x_limbs"], r, 4, N, barrett_mu(N_int, kN), kN)
pi_l = modpow(pi, to_limbs(L_derived, 4), 4, N, barrett_mu(N_int, kN), kN)
lhs = modmul(pi_l, x_r, N, barrett_mu(N_int, kN), kN)
print(f"  pi^L * x^r mod N = {from_limbs(lhs):#x}")
print(f"  == ±forged y' (1)? {from_limbs(lhs) in (1, N_int - 1)}  -> forgery REJECTED")

# mu tampering
mu_n_bad = barrett_mu(N_int - 2, kN)
print(f"  mu_n computed for N-2 passes mu_valid? {mu_valid(N_int, mu_n_bad, kN)}"
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

# ---------------------------------------------------------------- C. -1
print("\n== C. -1 malleability ==")
y_neg = N_int - y_int
nonce_neg, L_neg = find_challenge(N_int, x_int, y_neg, T, kN)
pi_neg = N_int - pow(x_int, (1 << T) // L_neg, N_int)
holds = pow(pi_neg, L_neg, N_int) * pow(x_int, pow(2, T, L_neg), N_int) % N_int == y_neg
print(f"  pi'^L' * x^r' == N - y (old full-group check): {holds}")
ok = wesolowski_verify(N_int, T, x_int, y_neg, pi_neg, kN, nonce_neg)
print(f"  signed-group verifier accepts N - y? {ok}  -> {'STILL BROKEN' if ok else 'REJECTED'}")

# ---------------------------------------------------------------- D. smooth L
print("\n== D. composite (smooth) challenge ==")
TRUNC, SMOOTH = 48, 1 << 16
primes = [q for q in range(2, SMOOTH) if all(q % d for d in range(2, int(q ** 0.5) + 1))]
M = 1
for q in primes:
    M *= q
alpha = 0  # 2^T mod M by CRT, without T squarings
for q in primes:
    Mq = M // q
    alpha = (alpha + pow(2, T, q) * Mq * pow(Mq, -1, q)) % M


def toy_challenge(y):
    h = poseidon_hash_span([N_int, x_int, y, T])
    return (h % (1 << TRUNC)) | (1 << (TRUNC - 8))


xM, e, y = pow(x_int, M, N_int), alpha, pow(x_int, alpha, N_int)
for k in range(1, 10 ** 6):
    L = toy_challenge(y)
    if M % L == 0:
        r = e % L
        pi = pow(x_int, (e - r) // L, N_int)
        ok = pow(pi, L, N_int) * pow(x_int, pow(2, T, L), N_int) % N_int == y
        honest = pow(x_int, 1 << T, N_int)
        print(f"  after {k} grinds: forged y verifies with composite L={L:#x}: {ok}; "
              f"y is the true output: {y in (honest, N_int - honest)}")
        print(f"  Miller-Rabin accepts that L? {is_prime_fs(L, 0)}  -> "
              f"{'STILL BROKEN' if is_prime_fs(L, 0) else 'REJECTED by hash-to-prime'}")
        break
    e, y = e + M, y * xM % N_int
