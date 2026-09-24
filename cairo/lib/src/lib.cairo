/// Wesolowski VDF verifier over an RSA-style modulus, as a Cairo *executable*
/// (not a contract). This is the program you prove with stwo-cairo: run it
/// offchain, produce a CairoProof, and verify the proof onchain via the Cairo
/// verifier (recursive) or a Rust verifier — proving that the T sequential
/// squarings happened without ever running them in a Starknet tx.
///
/// Scheme:
///   prover computes y = x^(2^T) mod N (T sequential squarings, offchain),
///   outputs its sign-canonical form |y| = min(y, N - y), and
///   pi = x^floor(2^T / L) mod N where L is a prime derived from
///   H(N, x, |y|, T, nonce) (hash-to-prime, checked in-program).
///   verifier checks, in O(log T) group ops, in the signed group Z_N^* / {±1}:
///     r = 2^T mod L
///     pi^L * x^r == ±|y|  (mod N)
///
/// Combined with RSW timelock encryption (encryptor knows phi(N) and derives
/// the key instantly; everyone else must square T times), this gives "decryption
/// capability materializes at the deadline" with no standing decryptor, no
/// committee, and no TEE.
///
/// Big-int representation: little-endian limbs of 64 bits. Two limb domains:
///   N (512-bit demo modulus)   -> N_LIMBS = 8 (32 for 2048-bit)
///   L (251-bit prime challenge) -> native u256; 4 limbs as an exponent
/// Barrett reduction mod N uses mu = floor(b^(2k)/N) with k+1 limbs, computed
/// offchain and passed in (validated in-program, see verify_vdf). Arithmetic
/// mod L uses the native u512-by-u256 division.
///
/// Limbs are u64; multiplication columns accumulate exactly in felt252 and
/// split once into u128 halves for the carry chain. ref/cairo_model_v2.py
/// mirrors the limb arithmetic and joint exponentiation.

use core::integer::{u512, u512_safe_div_rem_by_u256};
use core::math::u256_mul_mod_n;

pub const LIMB_BITS: u32 = 64;
pub const N_LIMBS: usize = 8; // 512-bit RSA modulus (demo/CI); 32 for 2048-bit production
pub const L_LIMBS: usize = 4;
/// Miller-Rabin rounds for the challenge prime. For a random odd 251-bit
/// candidate with uniform bases, a composite passes 20 rounds with probability
/// below 2^-107 (Damgard-Landrock-Pomerance average-case bound), so grinding
/// nonces for a composite challenge is infeasible.
pub const MR_ROUNDS: usize = 20;

// Poseidon domain tags (short strings).
pub const INSTANCE_TAG: felt252 = 'VDF_WES_INSTANCE_V1';
pub const CHALLENGE_TAG: felt252 = 'VDF_WES_CHALLENGE_V1';
pub const MR_TAG: felt252 = 'VDF_WES_MILLER_RABIN_V1';
pub const KEY_TAG: felt252 = 'VDF_RSW_KEY_V1';

const TWO64: u128 = 0x10000000000000000;
const MASK64: u128 = 0xffffffffffffffff;

/// Schoolbook column multiplication: one output limb at a time with a chained
/// carry, accumulating each column in felt252 and splitting once into u128 halves.
/// With n = min(an, bn), carry < n*2^64 and acc < n*2^128 by induction:
/// n*(2^64-1)^2 + n*2^64 < n*2^128. Since usize is u32, acc < 2^160,
/// below the Cairo field modulus, and the next carry fits u128. For the
/// 2048-bit configuration n <= 33, giving the tighter bound acc < 2^134.
/// Only contributing pairs are visited: i in [max(0, t+1-bn), min(an-1, t)].
/// The product is truncated to the low `out_limbs` columns (the carry out of
/// the last kept column is dropped), which is exact for callers that only
/// consume low limbs (Barrett's q3n). Mirrors cairo_model_v2.bigint_mul.
pub fn bigint_mul(
    a: Span<u64>, b: Span<u64>, an: usize, bn: usize, out_limbs: usize,
) -> Array<u64> {
    let mut total = an + bn;
    if out_limbs < total {
        total = out_limbs;
    }
    let mut out = ArrayTrait::new();
    let mut carry: u128 = 0;
    let mut t: usize = 0;
    while t < total {
        // The bound above makes felt addition and multiplication exact integers.
        let mut acc: felt252 = carry.into();
        let mut i: usize = if t + 1 > bn {
            t + 1 - bn
        } else {
            0
        };
        let hi_i: usize = if t < an {
            t
        } else {
            an - 1
        };
        while i <= hi_i {
            let ai: felt252 = (*a.at(i)).into();
            let bj: felt252 = (*b.at(t - i)).into();
            acc = acc + ai * bj;
            i += 1;
        }
        let split: u256 = acc.into();
        let lo = split.low;
        let hi = split.high;
        let limb: u64 = (lo & MASK64).try_into().unwrap();
        out.append(limb);
        carry = hi * TWO64 + lo / TWO64;
        t += 1;
    }
    out
}

/// True if a >= b (both k limbs, big-endian comparison).
pub fn ge_limbs(a: Span<u64>, b: Span<u64>, k: usize) -> bool {
    let mut i = k;
    loop {
        if i == 0 {
            break;
        }
        i -= 1;
        let ai = *a.at(i);
        let bi = *b.at(i);
        if ai != bi {
            return ai > bi;
        }
    }
    true // equal
}

/// Limb-wise a - b with borrow, k limbs, u128-safe arithmetic.
pub fn sub_limbs(a: Span<u64>, b: Span<u64>, k: usize) -> Array<u64> {
    let mut out = ArrayTrait::new();
    let mut borrow: u64 = 0;
    let mut i: usize = 0;
    while i < k {
        let minuend: u128 = (*a.at(i)).into();
        let subtrahend: u128 = (*b.at(i)).into() + borrow.into();
        if minuend >= subtrahend {
            out.append((minuend - subtrahend).try_into().unwrap());
            borrow = 0;
        } else {
            out.append((minuend + TWO64 - subtrahend).try_into().unwrap());
            borrow = 1;
        }
        i += 1;
    }
    out
}

/// Barrett reduction of a 2k-limb product p modulo k-limb N (HAC 14.42).
/// mu = floor(b^(2k)/N) with k+1 limbs, precomputed offchain.
///
/// Requires N's top limb nonzero (b^(k-1) <= N < b^k). Then the quotient
/// estimate satisfies floor(p/N) - 2 <= q3 <= floor(p/N), so
/// 0 <= p - q3*N < 3N < b^(k+1): r must keep k+1 limbs through the
/// normalization. (Keeping only k limbs wraps whenever p - q3*N >= b^k and
/// returns a wrong residue — e.g. N = 2^512-1, a = 2^448-1, b = 2^448+1 is
/// off by one. See scripts/repro_issues.py.)
pub fn barrett_reduce(p: Span<u64>, N: Span<u64>, mu: Span<u64>, k: usize) -> Array<u64> {
    // q1 = p >> (k-1)   (k+1 limbs)
    let q1 = p.slice(k - 1, k + 1);
    // q2 = q1 * mu; retain all columns because low carries feed the high half.
    let q2 = bigint_mul(q1, mu, k + 1, k + 1, 2 * (k + 1));
    let q3 = q2.span().slice(k + 1, k + 1); // q2 >> (k+1)

    // q3n = low k+1 limbs of q3 * N (only columns 0..k feed r)
    let q3n = bigint_mul(q3, N, k + 1, k, k + 1);

    // r = (p - q3n) mod b^(k+1) with borrow — exact, no truncation
    let mut r = ArrayTrait::new();
    let mut borrow: u64 = 0;
    let mut i: usize = 0;
    while i < k + 1 {
        let minuend: u128 = (*p.at(i)).into();
        let subtrahend: u128 = (*q3n.at(i)).into() + borrow.into();
        if minuend >= subtrahend {
            r.append((minuend - subtrahend).try_into().unwrap());
            borrow = 0;
        } else {
            r.append((minuend + TWO64 - subtrahend).try_into().unwrap());
            borrow = 1;
        }
        i += 1;
    }

    // N zero-extended to k+1 limbs
    let mut n_ext = ArrayTrait::new();
    i = 0;
    while i < k {
        n_ext.append(*N.at(i));
        i += 1;
    }
    n_ext.append(0);

    // normalize: while r >= N, r -= N (at most 2 iterations since r < 3N)
    while ge_limbs(r.span(), n_ext.span(), k + 1) {
        r = sub_limbs(r.span(), n_ext.span(), k + 1);
    }

    // drop the (now zero) top limb
    let mut out = ArrayTrait::new();
    i = 0;
    while i < k {
        out.append(*r.at(i));
        i += 1;
    }
    out
}

/// Modular multiplication: a*b mod N. a, b < N < b^k.
pub fn modmul(a: Span<u64>, b: Span<u64>, N: Span<u64>, mu: Span<u64>, k: usize) -> Array<u64> {
    let p = bigint_mul(a, b, k, k, 2 * k);
    barrett_reduce(p.span(), N, mu, k)
}

/// Modular exponentiation: base^exp mod N, square-and-multiply over the bits
/// of the (exp_n-limb) exponent. Trailing zero limbs/bits are skipped and the
/// squaring after the top set bit is omitted.
pub fn modpow(
    base: Span<u64>, exp: Span<u64>, exp_n: usize, N: Span<u64>, mu: Span<u64>, k: usize,
) -> Array<u64> {
    // one = 1 in the k-limb domain
    let mut one = ArrayTrait::new();
    one.append(1);
    let mut i: usize = 1;
    while i < k {
        one.append(0);
        i += 1;
    }
    let mut result = one.clone();
    let mut base_red = modmul(base, one.span(), N, mu, k);

    // bit length of exp: skip trailing zero limbs, then count the top limb's bits
    let mut top = exp_n;
    loop {
        if top == 0 {
            break;
        }
        if *exp.at(top - 1) != 0 {
            break;
        }
        top -= 1;
    }
    if top == 0 {
        return result; // exp == 0 -> base^0 = 1
    }
    let mut e_top = *exp.at(top - 1);
    let mut top_bits: u32 = 0;
    while e_top != 0 {
        e_top = e_top / 2;
        top_bits += 1;
    }

    i = 0;
    while i < top {
        let mut e = *exp.at(i);
        let limit: u32 = if i == top - 1 {
            top_bits
        } else {
            LIMB_BITS
        };
        let mut bit: u32 = 0;
        while bit < limit {
            if (e % 2) == 1 {
                result = modmul(result.span(), base_red.span(), N, mu, k);
            }
            e = e / 2;
            bit += 1;
            let last = i == top - 1 && bit == limit;
            if !last {
                base_red = modmul(base_red.span(), base_red.span(), N, mu, k);
            }
        }
        i += 1;
    }
    result
}

/// Compare two k-limb arrays for equality.
pub fn limbs_eq(a: Span<u64>, b: Span<u64>, k: usize) -> bool {
    let mut i: usize = 0;
    while i < k {
        if *a.at(i) != *b.at(i) {
            return false;
        }
        i += 1;
    }
    true
}

/// Validate an offchain Barrett constant: mu == floor(b^(2k)/m), i.e.
/// mu * m <= b^(2k) < (mu + 1) * m. Without this check mu is a trusted input
/// and a prover can make barrett_reduce return attacker-chosen residues.
pub fn mu_valid(m: Span<u64>, mu: Span<u64>, k: usize) -> bool {
    // b^(2k) in 2k+1 limbs
    let mut b2k = ArrayTrait::new();
    let mut i: usize = 0;
    while i < 2 * k {
        b2k.append(0);
        i += 1;
    }
    b2k.append(1);

    // mu * m <= b^(2k)
    let prod = bigint_mul(mu, m, k + 1, k, 2 * k + 1);
    if !ge_limbs(b2k.span(), prod.span(), 2 * k + 1) {
        return false;
    }

    // (mu + 1) * m > b^(2k); if mu + 1 overflows k+1 limbs the product is 0
    let mut mu1 = ArrayTrait::new();
    let mut carry: u128 = 1;
    i = 0;
    while i < k + 1 {
        let s: u128 = (*mu.at(i)).into() + carry;
        mu1.append((s & MASK64).try_into().unwrap());
        carry = s / TWO64;
        i += 1;
    }
    let prod2 = bigint_mul(mu1.span(), m, k + 1, k, 2 * k + 1);
    if !ge_limbs(prod2.span(), b2k.span(), 2 * k + 1) {
        return false;
    }
    !limbs_eq(prod2.span(), b2k.span(), 2 * k + 1)
}

/// Hash of the VDF instance: poseidon(INSTANCE_TAG, len(N), N || x || y || T
/// limbs). Seeds the Fiat-Shamir challenge and is returned by the executable as
/// the public, instance-binding output. `y` is the sign-canonical output |y|.
pub fn instance_hash(N: Span<u64>, x: Span<u64>, y: Span<u64>, T_limbs: Span<u64>) -> felt252 {
    let mut input = array![INSTANCE_TAG, N.len().into()];
    append_limbs(ref input, N);
    append_limbs(ref input, x);
    append_limbs(ref input, y);
    append_limbs(ref input, T_limbs);
    core::poseidon::poseidon_hash_span(input.span())
}

fn append_limbs(ref out: Array<felt252>, limbs: Span<u64>) {
    for limb in limbs {
        out.append((*limb).into());
    }
}

/// Fiat-Shamir seed for the challenge: poseidon(CHALLENGE_TAG, instance, nonce).
/// The prover increments `nonce` offchain until the candidate is prime.
pub fn challenge_hash(instance: felt252, nonce: felt252) -> felt252 {
    core::poseidon::poseidon_hash_span(array![CHALLENGE_TAG, instance, nonce].span())
}

/// Challenge candidate from its seed: the low 250 bits with bits 250 and 0
/// forced, so L is odd and 2^250 <= L < 2^251. Primality is checked separately
/// by `is_probable_prime`.
pub fn challenge_from_hash(h: felt252) -> u256 {
    let h: u256 = h.into();
    u256 {
        low: h.low | 1,
        high: (h.high & 0x3ffffffffffffffffffffffffffffff) | 0x4000000000000000000000000000000,
    }
}

/// base^exp mod n over native u256 (u512 wide product, exact remainder).
pub fn u256_pow_mod(base: u256, exp: u256, n: NonZero<u256>) -> u256 {
    let mut result: u256 = 1;
    let mut b = u256_mul_mod_n(base, 1, n);
    let mut e = exp;
    while e != 0 {
        if e % 2 == 1 {
            result = u256_mul_mod_n(result, b, n);
        }
        e = e / 2;
        if e != 0 {
            b = u256_mul_mod_n(b, b, n);
        }
    }
    result
}

/// Miller-Rabin with MR_ROUNDS Fiat-Shamir bases derived from `seed`. Each
/// base is two Poseidon outputs (a 504-bit integer, low || high) reduced mod
/// L, which is statistically uniform. Requires L odd and L > 3.
pub fn is_probable_prime(L: u256, seed: felt252) -> bool {
    let n: NonZero<u256> = L.try_into().unwrap();
    let l_minus_1 = L - 1;
    let mut d = l_minus_1;
    let mut s: usize = 0;
    while d % 2 == 0 {
        d = d / 2;
        s += 1;
    }

    let mut round: usize = 0;
    while round < MR_ROUNDS {
        let lo: u256 = core::poseidon::poseidon_hash_span(
            array![MR_TAG, seed, (2 * round).into()].span(),
        )
            .into();
        let hi: u256 = core::poseidon::poseidon_hash_span(
            array![MR_TAG, seed, (2 * round + 1).into()].span(),
        )
            .into();
        let wide = u512 { limb0: lo.low, limb1: lo.high, limb2: hi.low, limb3: hi.high };
        let (_, base) = u512_safe_div_rem_by_u256(wide, n);

        let mut z = u256_pow_mod(base, d, n);
        if z != 1 && z != l_minus_1 {
            let mut witness = true;
            let mut j: usize = 1;
            while j < s {
                z = u256_mul_mod_n(z, z, n);
                if z == l_minus_1 {
                    witness = false;
                    break;
                }
                if z == 1 {
                    break;
                }
                j += 1;
            }
            if witness {
                return false;
            }
        }
        round += 1;
    }
    true
}

/// A u256 as 4 little-endian u64 limbs (the exponent format of joint_modpow).
fn u256_limbs(v: u256) -> Array<u64> {
    array![
        (v.low & MASK64).try_into().unwrap(), (v.low / TWO64).try_into().unwrap(),
        (v.high & MASK64).try_into().unwrap(), (v.high / TWO64).try_into().unwrap(),
    ]
}

/// True iff y < N - y, i.e. y is the sign-canonical representative of ±y
/// (N odd, so y != N - y). Requires y < N.
fn is_sign_canonical(y: Span<u64>, N: Span<u64>) -> bool {
    if ge_limbs(y, N, N_LIMBS) {
        return false;
    }
    let neg_y = sub_limbs(N, y, N_LIMBS);
    !ge_limbs(y, neg_y.span(), N_LIMBS)
}

/// Wesolowski VDF verification core, in the signed group Z_N^* / {±1}.
/// `y` must be sign-canonical (y < N - y): accepting only one of ±y removes the
/// -1 malleability that let anyone who knew y "prove" N - y. The challenge L is
/// derived in-program from (N, x, y, T, nonce) and must be prime: a composite
/// challenge admits smooth-challenge forgeries of a wrong y without factoring N.
/// Returns true iff N is odd, mu_n is exact, L is prime and
/// pi^L * x^(2^T mod L) == ±y mod N. `T_limbs` is T as two u64 limbs.
pub fn verify_vdf(
    N: Span<u64>,
    mu_n: Span<u64>,
    x: Span<u64>,
    y: Span<u64>,
    pi: Span<u64>,
    T_limbs: Span<u64>,
    nonce: felt252,
) -> bool {
    // 0. N odd, the offchain-computed Barrett constant matches N, and y is
    //    the sign-canonical output.
    if *N.at(0) % 2 == 0 || !mu_valid(N, mu_n, N_LIMBS) {
        return false;
    }
    if !is_sign_canonical(y, N) {
        return false;
    }
    let seed = challenge_hash(instance_hash(N, x, y, T_limbs), nonce);
    let L = challenge_from_hash(seed);
    if !is_probable_prime(L, seed) {
        return false;
    }

    // 1. r = 2^T mod L (T is two u64 limbs)
    let T: u256 = (*T_limbs.at(0)).into() + (*T_limbs.at(1)).into() * TWO64.into();
    let r = u256_pow_mod(2, T, L.try_into().unwrap());

    // 2. lhs = pi^L * x^r mod N
    let lhs = joint_modpow(pi, u256_limbs(L).span(), x, u256_limbs(r).span(), N, mu_n, N_LIMBS);

    // 3. check lhs == ±y
    if limbs_eq(lhs.span(), y, N_LIMBS) {
        return true;
    }
    let neg_lhs = sub_limbs(N, lhs.span(), N_LIMBS);
    limbs_eq(neg_lhs.span(), y, N_LIMBS)
}

/// RSW key from the sign-canonical VDF output: poseidon(KEY_TAG, ctx, len(y),
/// y limbs). `ctx` is the application context (e.g. H(auction, bid)), so one
/// puzzle output never yields the same key in two contexts.
pub fn derive_key(ctx: felt252, y: Span<u64>) -> felt252 {
    let mut input = array![KEY_TAG, ctx, y.len().into()];
    append_limbs(ref input, y);
    core::poseidon::poseidon_hash_span(input.span())
}

/// Poseidon keystream: ks_i = poseidon(key, i); ct_i = pt_i + ks_i in the field.
/// No integrity tag: the application binds the plaintext by comparing its hash
/// against the commitment posted with the ciphertext.
pub fn decrypt(key: felt252, ct: Span<felt252>) -> Array<felt252> {
    let mut pt = ArrayTrait::new();
    let mut i: usize = 0;
    while i < ct.len() {
        let ks = core::poseidon::poseidon_hash_span(array![key, i.into()].span());
        pt.append(*ct.at(i) - ks);
        i += 1;
    }
    pt
}

pub fn encrypt(key: felt252, pt: Span<felt252>) -> Array<felt252> {
    let mut ct = ArrayTrait::new();
    let mut i: usize = 0;
    while i < pt.len() {
        let ks = core::poseidon::poseidon_hash_span(array![key, i.into()].span());
        ct.append(*pt.at(i) + ks);
        i += 1;
    }
    ct
}

/// a^e * b^f mod N for two L_LIMBS-limb exponents, sharing each squaring.
/// Inputs may be noncanonical k-limb bases; reduce them before precomputation.
/// Both zero exponents return one, matching modpow for the supported moduli.
fn joint_modpow(
    a: Span<u64>, e: Span<u64>, b: Span<u64>, f: Span<u64>, N: Span<u64>, mu: Span<u64>, k: usize,
) -> Array<u64> {
    let mut one = array![1_u64];
    let mut i = 1;
    while i < k {
        one.append(0);
        i += 1;
    }
    let ar = modmul(a, one.span(), N, mu, k);
    let br = modmul(b, one.span(), N, mu, k);
    let ab = modmul(ar.span(), br.span(), N, mu, k);
    let mut result = one;
    // Skip leading zero pairs and copy the first factor instead of multiplying one.
    let mut started = false;
    i = L_LIMBS;
    while i != 0 {
        i -= 1;
        let ei = *e.at(i);
        let fi = *f.at(i);
        let mut mask: u64 = 0x8000000000000000;
        while mask != 0 {
            let eb = (ei / mask) % 2 != 0;
            let fb = (fi / mask) % 2 != 0;
            if started {
                result = modmul(result.span(), result.span(), N, mu, k);
            }
            if eb || fb {
                let factor = if eb {
                    if fb {
                        ab.span()
                    } else {
                        ar.span()
                    }
                } else {
                    br.span()
                };
                if started {
                    result = modmul(result.span(), factor, N, mu, k);
                } else {
                    result = factor.into();
                    started = true;
                }
            }
            mask /= 2;
        }
    }
    result
}
#[cfg(test)]
mod arithmetic_tests;

mod tests;
