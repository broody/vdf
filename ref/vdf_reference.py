#!/usr/bin/env python3
"""RSW timelock encryption + Wesolowski VDF reference implementation.

Generates JSON test vectors consumed by the Cairo verifier tests.

Schemes:
- RSW timelock: encryptor knows phi(N) and derives the key instantly via
  2^T mod phi(N); everyone else must compute x^(2^T) mod N by T sequential
  squarings.
- Wesolowski VDF in the signed group Z_N^*/{±1}: the output is the
  sign-canonical |y| = min(y, N - y) of y = x^(2^T) mod N; the prover also
  produces pi = x^floor(2^T / L) mod N where L is a prime derived from
  H(N, x, |y|, T, nonce); the verifier checks pi^L * x^(2^T mod L) == ±|y|
  mod N in O(log T) group ops.
- RSW encryption: key = poseidon(KEY_TAG, ctx, |y|), Poseidon keystream
  ct_i = pt_i + poseidon(key, i) over the Stark field.
"""

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


FELT_PRIME = (1 << 251) + 17 * (1 << 192) + 1
MR_ROUNDS = 20


def short_string(text: str) -> int:
    """Cairo short-string felt encoding."""
    return int.from_bytes(text.encode(), "big")


INSTANCE_TAG = short_string("VDF_WES_INSTANCE_V1")
CHALLENGE_TAG = short_string("VDF_WES_CHALLENGE_V1")
MR_TAG = short_string("VDF_WES_MILLER_RABIN_V1")
KEY_TAG = short_string("VDF_RSW_KEY_V1")


def sign_canonical(y: int, N: int) -> int:
    """Representative of ±y in Z_N^*/{±1}: min(y, N - y)."""
    return min(y, N - y)


def instance_hash(N: int, x: int, y: int, T: int, n_limbs: int) -> int:
    """Mirror of Cairo instance_hash: poseidon(TAG, n_limbs, N || x || y || T limbs)."""
    limbs = (
        limbify(N, 64, n_limbs)
        + limbify(x, 64, n_limbs)
        + limbify(y, 64, n_limbs)
        + limbify(T, 64, 2)
    )
    return poseidon_hash_span([INSTANCE_TAG, n_limbs] + limbs)


def challenge_hash(instance: int, nonce: int) -> int:
    return poseidon_hash_span([CHALLENGE_TAG, instance, nonce])


def challenge_from_hash(h: int) -> int:
    """Low 250 bits with bits 250 and 0 forced: odd, 2^250 <= L < 2^251."""
    return (h & ((1 << 250) - 1)) | (1 << 250) | 1


def mr_bases(seed: int, L: int) -> list[int]:
    """Fiat-Shamir Miller-Rabin bases: two Poseidon outputs as a 504-bit
    integer (low || high), reduced mod L."""
    bases = []
    for i in range(MR_ROUNDS):
        lo = poseidon_hash_span([MR_TAG, seed, 2 * i])
        hi = poseidon_hash_span([MR_TAG, seed, 2 * i + 1])
        bases.append((lo + (hi << 256)) % L)
    return bases


def is_prime_fs(L: int, seed: int) -> bool:
    """Mirror of Cairo is_probable_prime (L odd, bases from `seed`)."""
    d, s = L - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in mr_bases(seed, L):
        z = pow(a, d, L)
        if z in (1, L - 1):
            continue
        for _ in range(s - 1):
            z = z * z % L
            if z == L - 1:
                break
            if z == 1:
                return False
        else:
            return False
    return True


def derive_challenge(N: int, x: int, y: int, T: int, n_limbs: int, nonce: int) -> tuple[int, int]:
    """(L, seed) for a given nonce; L need not be prime (the verifier checks)."""
    seed = challenge_hash(instance_hash(N, x, y, T, n_limbs), nonce)
    return challenge_from_hash(seed), seed


def find_challenge(N: int, x: int, y: int, T: int, n_limbs: int) -> tuple[int, int]:
    """Smallest nonce whose challenge passes the in-program primality test.
    Returns (nonce, L)."""
    for nonce in range(1 << 16):
        L, seed = derive_challenge(N, x, y, T, n_limbs, nonce)
        if is_prime_fs(L, seed):
            return nonce, L
    raise ValueError("no prime challenge below nonce 2^16")


def modpow(base: int, exp: int, mod: int) -> int:
    return pow(base, exp, mod)


def sequential_square(x: int, T: int, N: int) -> int:
    """The actual VDF evaluation: T sequential squarings (no shortcut)."""
    y = x
    for _ in range(T):
        y = (y * y) % N
    return y


def wesolowski_prove(N: int, T: int, x: int, y: int, n_limbs: int) -> tuple[int, int, int, int]:
    """Produce (pi, L, r, nonce) for the sign-canonical output y. Requires
    q = floor(2^T / L), which the honest prover gets from its T squarings."""
    nonce, L = find_challenge(N, x, y, T, n_limbs)
    q, r = divmod(1 << T, L)
    pi = modpow(x, q, N)
    return pi, L, r, nonce


def wesolowski_verify(N: int, T: int, x: int, y: int, pi: int, n_limbs: int, nonce: int) -> bool:
    """Mirror of Cairo verify_vdf."""
    if N % 2 == 0 or not (0 <= y < N - y) or not 0 <= nonce < 1 << 16:
        return False
    L, seed = derive_challenge(N, x, y, T, n_limbs, nonce)
    if not is_prime_fs(L, seed):
        return False
    lhs = (modpow(pi, L, N) * modpow(x, pow(2, T, L), N)) % N
    return lhs in (y, N - y)


def derive_key(ctx: int, y: int, n_limbs: int) -> int:
    return poseidon_hash_span([KEY_TAG, ctx, n_limbs] + limbify(y, 64, n_limbs))


def encrypt(key: int, pt: list[int]) -> list[int]:
    return [(p + poseidon_hash_span([key, i])) % FELT_PRIME for i, p in enumerate(pt)]


def decrypt(key: int, ct: list[int]) -> list[int]:
    return [(c - poseidon_hash_span([key, i])) % FELT_PRIME for i, c in enumerate(ct)]


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

    # RSW: encryptor uses phi(N) to derive the output instantly.
    exp = modpow(2, T, phi)  # 2^T mod phi(N)
    y_fast = modpow(x, exp, N)  # == x^(2^T) mod N, computed in O(log T)

    # Sequential evaluation (what everyone else must do).
    y_slow = sequential_square(x, T, N)
    assert y_fast == y_slow, "RSW shortcut must match sequential squaring"
    y = sign_canonical(y_fast, N)

    # Wesolowski proof over the sign-canonical output.
    n_limbs = (n_bits + limb_bits - 1) // limb_bits
    pi, L, r, nonce = wesolowski_prove(N, T, x, y, n_limbs)
    assert wesolowski_verify(N, T, x, y, pi, n_limbs, nonce), "proof must verify"

    # RSW timelock: Poseidon keystream under a key from the VDF output. The
    # plaintext is an example bid opening (amount, salt, tag).
    ctx = poseidon_hash_span([short_string("auction"), 7, 3])
    pt = [1_250_000, rng.getrandbits(250), short_string("bid")]
    ct = encrypt(derive_key(ctx, y, n_limbs), pt)
    assert decrypt(derive_key(ctx, sign_canonical(y_slow, N), n_limbs), ct) == pt

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
        "y": hex(y),
        "y_limbs": limbify(y, limb_bits, n_limbs),
        "pi": hex(pi),
        "pi_limbs": limbify(pi, limb_bits, n_limbs),
        "nonce": nonce,
        "L": hex(L),
        "L_limbs": limbify(L, limb_bits, 4),
        "r": hex(r),
        "r_limbs": limbify(r, limb_bits, 4),
        "instance_hash": hex(instance_hash(N, x, y, T, n_limbs)),
        "ctx": hex(ctx),
        "plaintext": [hex(v) for v in pt],
        "ciphertext": [hex(v) for v in ct],
        "ct_hash": hex(poseidon_hash_span(ct)),
        "pt_hash": hex(poseidon_hash_span(pt)),
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
    """Regenerate vectors/vdf_vectors_{512,2048}.json (T = 2^16, fixed seeds)."""
    import os

    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    for n_bits, seed in ((512, 512), (2048, 2048)):
        v = make_vector(n_bits, T=1 << 16, seed=seed)
        out_path = os.path.join(root, "vectors", f"vdf_vectors_{n_bits}.json")
        with open(out_path, "w") as f:
            json.dump([v], f, indent=2)
        print(f"{v['name']}: nonce={v['nonce']} L-bits={int(v['L'], 16).bit_length()} -> {out_path}")


if __name__ == "__main__":
    main()
