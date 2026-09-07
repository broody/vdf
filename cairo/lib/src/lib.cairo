/// Wesolowski VDF verifier over an RSA-style modulus, as a Cairo *executable*
/// (not a contract). This is the program you prove with stwo-cairo: run it
/// offchain, produce a CairoProof, and verify the proof onchain via the Cairo
/// verifier (recursive) or a Rust verifier — proving that the T sequential
/// squarings happened without ever running them in a Starknet tx.
///
/// Scheme:
///   prover computes y = x^(2^T) mod N (T sequential squarings, offchain)
///   and pi = x^floor(2^T / L) mod N where L = H(x, y, T).
///   verifier checks, in O(log T) group ops:
///     r = 2^T mod L
///     pi^L * x^r == y  (mod N)
///
/// Combined with RSW timelock encryption (encryptor knows phi(N) and derives
/// the key instantly; everyone else must square T times), this gives "decryption
/// capability materializes at the deadline" with no standing decryptor, no
/// committee, and no TEE.
///
/// Big-int representation: little-endian limbs of 64 bits. Two limb domains:
///   N (512-bit demo modulus)   -> N_LIMBS = 8 (32 for 2048-bit)
///   L (256-bit challenge)      -> L_LIMBS = 4
/// Barrett reduction uses mu = floor(b^(2k)/N) with k+1 limbs, computed
/// offchain and passed in (validated in-program, see verify_vdf).
///
/// Limbs are u64; multiplication columns accumulate exactly in felt252 and
/// split once into u128 halves for the carry chain. ref/cairo_model_v2.py
/// mirrors the limb arithmetic and joint exponentiation.

pub const LIMB_BITS: u32 = 64;
pub const N_LIMBS: usize = 8; // 512-bit RSA modulus (demo/CI); 32 for 2048-bit production
pub const L_LIMBS: usize = 4;

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

/// Hash of the VDF instance: poseidon(N || x || y || T limbs as felt252s).
/// Seeds the Fiat-Shamir challenge and is returned by the executable as the
/// public, instance-binding output.
pub fn instance_hash(N: Span<u64>, x: Span<u64>, y: Span<u64>, T_limbs: Span<u64>) -> felt252 {
    let mut input = ArrayTrait::new();
    let mut i: usize = 0;
    while i < N.len() {
        input.append((*N.at(i)).into());
        i += 1;
    }
    i = 0;
    while i < x.len() {
        input.append((*x.at(i)).into());
        i += 1;
    }
    i = 0;
    while i < y.len() {
        input.append((*y.at(i)).into());
        i += 1;
    }
    i = 0;
    while i < T_limbs.len() {
        input.append((*T_limbs.at(i)).into());
        i += 1;
    }
    core::poseidon::poseidon_hash_span(input.span())
}

/// Fiat-Shamir challenge for the Wesolowski proof, derived in-program:
/// L = instance_hash(N, x, y, T) as 4 u64 limbs with bit 192 forced, so
/// b^3 <= L < b^4 always satisfies the Barrett precondition. Binding L to the
/// instance is what makes the proof sound: a prover cannot pick L to fit a
/// forged y (e.g. L dividing 2^T, which collapses the check to pi == y).
pub fn derive_challenge(
    N: Span<u64>, x: Span<u64>, y: Span<u64>, T_limbs: Span<u64>,
) -> Array<u64> {
    let h: u256 = instance_hash(N, x, y, T_limbs).into();
    let low: u128 = h.low;
    let high: u128 = h.high;
    let mut out = ArrayTrait::new();
    out.append((low & MASK64).try_into().unwrap());
    out.append((low / TWO64).try_into().unwrap());
    out.append((high & MASK64).try_into().unwrap());
    let l3: u64 = (high / TWO64).try_into().unwrap();
    out.append(l3 | 1); // force bit 192: L >= b^3 (Barrett precondition)
    out
}

/// Wesolowski VDF verification core.
/// Derives the challenge L = H(N, x, y, T) in-program, validates both Barrett
/// constants, then returns true iff pi^L * x^(2^T mod L) == y mod N.
pub fn verify_vdf(
    N: Span<u64>,
    mu_n: Span<u64>,
    x: Span<u64>,
    y: Span<u64>,
    pi: Span<u64>,
    mu_l: Span<u64>,
    T_limbs: Span<u64>,
    T_n: usize,
) -> bool {
    // 0. the offchain-computed Barrett constants must match their moduli
    if !mu_valid(N, mu_n, N_LIMBS) {
        return false;
    }
    let L = derive_challenge(N, x, y, T_limbs);
    if !mu_valid(L.span(), mu_l, L_LIMBS) {
        return false;
    }

    // 1. r = 2^T mod L in the L domain
    let mut two = ArrayTrait::new();
    two.append(2);
    let mut i: usize = 1;
    while i < L_LIMBS {
        two.append(0);
        i += 1;
    }
    let r = modpow(two.span(), T_limbs, T_n, L.span(), mu_l, L_LIMBS);

    // 2. lhs = pi^L * x^r mod N
    let lhs = joint_modpow(pi, L.span(), x, r.span(), N, mu_n, N_LIMBS);

    // 3. check lhs == y
    limbs_eq(lhs.span(), y, N_LIMBS)
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
