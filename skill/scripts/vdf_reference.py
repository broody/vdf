#!/usr/bin/env python3
"""RSW timelock encryption + Wesolowski VDF reference implementation.

Generates JSON test vectors consumed by the Cairo verifier tests.

Schemes:
- RSW timelock: encryptor knows phi(N) and derives the key instantly via
  2^T mod phi(N); everyone else must compute x^(2^T) mod N by T sequential
  squarings.
- Wesolowski VDF: after computing y = x^(2^T) mod N, the prover also produces
  pi = x^floor(2^T / L) mod N where L = hash(x, y, T); the verifier checks
  y == pi^L * x^(2^T mod L) mod N in O(log T) group ops.
"""

import hashlib
import json
import secrets
import sys

from poseidon import poseidon_hash_span


def limbify(value: int, limb_bits: int, n_limbs: int) -> list[int]:
    """Split an integer into little-endian limbs."""
    mask = (1 << limb_bits) - 1
    limbs = []
    for i in range(n_limbs):
        limbs.append((value >> (i * limb_bits)) & mask)
    return limbs


def derive_challenge(N: int, x: int, y: int, T: int, n_limbs: int) -> int:
    """Wesolowski Fiat-Shamir challenge, mirroring Cairo derive_challenge:
    L = poseidon(N || x || y || T u64 limbs) with bit 192 forced, so
    2^192 <= L < 2^256 always satisfies the Barrett precondition."""
    limbs = (
        limbify(N, 64, n_limbs)
        + limbify(x, 64, n_limbs)
        + limbify(y, 64, n_limbs)
        + limbify(T, 64, 2)
    )
    return poseidon_hash_span(limbs) | (1 << 192)


def modpow(base: int, exp: int, mod: int) -> int:
    return pow(base, exp, mod)


def sequential_square(x: int, T: int, N: int) -> int:
    """The actual VDF evaluation: T sequential squarings (no shortcut)."""
    y = x
    for _ in range(T):
        y = (y * y) % N
    return y


def wesolowski_prove(N: int, T: int, x: int, y: int, n_limbs: int) -> tuple[int, int]:
    """Produce (pi, L) for the Wesolowski proof. Requires knowing q = floor(2^T / L),
    which the honest prover computes by actually performing the T squarings."""
    L = derive_challenge(N, x, y, T, n_limbs)
    # q = floor(2^T / L), r = 2^T mod L
    q, r = divmod(1 << T, L)
    pi = modpow(x, q, N)
    return pi, L, r


def wesolowski_verify(N: int, T: int, x: int, y: int, pi: int, n_limbs: int) -> bool:
    L = derive_challenge(N, x, y, T, n_limbs)
    r = modpow(2, T, L)  # 2^T mod L — cheap, O(log T)
    lhs = (modpow(pi, L, N) * modpow(x, r, N)) % N
    return lhs == y


def make_vector(
    n_bits: int,
    T: int,
    limb_bits: int = 64,
    seed: int | None = None,
) -> dict:
    """Build one test vector: safe-prime RSA modulus, RSW key, Wesolowski proof."""
    rng = secrets.SystemRandom()
    if seed is not None:
        rng = __import__("random").Random(seed)

    # 1024-bit demo modulus: N = p*q with p,q ~ n_bits/2. Use safe primes so the
    # group order is well-defined and RSW's phi(N) shortcut works.
    def safe_prime(bits: int):
        while True:
            p = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
            if p % 2 == 0:
                continue
            # crude Miller-Rabin
            if not is_probable_prime(p):
                continue
            q = 2 * p + 1
            if is_probable_prime(q):
                return q
            p += 2

    # Simpler: just generate two probable primes.
    def probable_prime(bits: int):
        while True:
            p = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
            if is_probable_prime(p):
                return p

    p = probable_prime(n_bits // 2)
    q = probable_prime(n_bits // 2)
    N = p * q
    phi = (p - 1) * (q - 1)

    x = rng.getrandbits(n_bits - 8) | 1
    x %= N
    if x == 0:
        x = 1

    # RSW: encryptor uses phi(N) to derive the key instantly.
    exp = modpow(2, T, phi)  # 2^T mod phi(N)
    y_fast = modpow(x, exp, N)  # == x^(2^T) mod N, computed in O(log T)

    # Sequential evaluation (what everyone else must do).
    y_slow = sequential_square(x, T, N)
    assert y_fast == y_slow, "RSW shortcut must match sequential squaring"

    # Wesolowski proof.
    n_limbs = (n_bits + limb_bits - 1) // limb_bits
    pi, L, r = wesolowski_prove(N, T, x, y_fast, n_limbs)
    assert wesolowski_verify(N, T, x, y_fast, pi, n_limbs), "proof must verify"

    # RSW timelock: symmetric key from VDF output, AEAD-ish demo.
    msg = b"reveal-opening"
    key = hashlib.sha256(y_fast.to_bytes((n_bits + 7) // 8, "big")).digest()
    ciphertext = bytes(b ^ key[i % len(key)] for i, b in enumerate(msg))
    # Decrypt after sequential eval (the "deadline" step).
    key2 = hashlib.sha256(y_slow.to_bytes((n_bits + 7) // 8, "big")).digest()
    plain = bytes(c ^ key2[i % len(key2)] for i, c in enumerate(ciphertext))
    assert plain == msg

    return {
        "name": f"rsa{n_bits}_T{T}",
        "n_bits": n_bits,
        "limb_bits": limb_bits,
        "n_limbs": n_limbs,
        "T": T,
        "N": hex(N),
        "N_limbs": limbify(N, limb_bits, n_limbs),
        "x": hex(x),
        "x_limbs": limbify(x, limb_bits, n_limbs),
        "y": hex(y_fast),
        "y_limbs": limbify(y_fast, limb_bits, n_limbs),
        "pi": hex(pi),
        "pi_limbs": limbify(pi, limb_bits, n_limbs),
        "L": hex(L),
        "L_limbs": limbify(L, limb_bits, n_limbs),
        "r": hex(r),
        "r_limbs": limbify(r, limb_bits, n_limbs),
        "mu": None,  # Barrett mu computed in Cairo test setup (N known)
        "ciphertext_hex": ciphertext.hex(),
        "plaintext": msg.decode(),
    }


def is_probable_prime(n: int, rounds: int = 16) -> bool:
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d = n - 1
    s = 0
    while d % 2 == 0:
        s += 1
        d //= 2
    for _ in range(rounds):
        a = secrets.randbelow(n - 3) + 2
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = (x * x) % n
            if x == n - 1:
                break
        else:
            return False
    return True


def main() -> None:
    out_path = sys.argv[1] if len(sys.argv) > 1 else "vdf_vectors.json"
    vectors = []
    # T = 2^20 is a "real" VDF parameter (2^T >> 2^256 challenge L, so the
    # Wesolowski proof is non-trivial) while still being vectorizable fast via
    # the phi(N) shortcut. Realistic T would be ~2^30+ sequential squarings,
    # which the encryptor skips using phi(N); solvers cannot.
    vectors.append(make_vector(1024, T=1 << 20, seed=42))
    # Small-T cross-check that the phi(N) shortcut equals actual squarings.
    vectors.append(make_vector(1024, T=64, seed=1337))
    with open(out_path, "w") as f:
        json.dump(vectors, f, indent=2)
    for v in vectors:
        print(
            f"{v['name']}: L-bits={int(v['L'],16).bit_length()} "
            f"pi-bits={int(v['pi'],16).bit_length()} "
            f"r-bits={int(v['r'],16).bit_length()}"
        )
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
