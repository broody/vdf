use vdf::{N_LIMBS, decrypt, derive_key, instance_hash, verify_vdf};

/// Slice `n` felts starting at `offset` into u64 limbs, returning the array
/// and the new offset.
fn take_u64s(inputs: Span<felt252>, offset: usize, n: usize) -> (Array<u64>, usize) {
    let mut out = ArrayTrait::new();
    let mut j: usize = 0;
    while j < n {
        let v: u64 = (*inputs.at(offset + j)).try_into().unwrap();
        out.append(v);
        j += 1;
    }
    (out, offset + n)
}

/// Executable entrypoint: flat felt252 inputs
///   [N(k), mu_n(k+1), x(k), y(k), pi(k), T(2), nonce, ctx, ct_len, ct(ct_len)]
/// with k = N_LIMBS (8 for the 512-bit demo, 32 for 2048-bit production).
/// `y` is the sign-canonical VDF output min(y, N - y); `nonce` (< 2^16) selects
/// the prime challenge; the challenge L and r = 2^T mod L are derived in-program.
///
/// On success returns [instance_hash(N, x, y, T), ctx, poseidon(ct), poseidon(pt)]
/// where pt is ct decrypted under derive_key(ctx, y). The onchain side recomputes
/// the instance hash and poseidon(ct) from its stored puzzle and ciphertext, and
/// compares poseidon(pt) against the posted commitment. Returns [0, 0, 0, 0] on
/// verification failure; the input length must match exactly.
#[executable]
pub fn vdf_verify(inputs: Span<felt252>) -> Array<felt252> {
    let (N, off) = take_u64s(inputs, 0, N_LIMBS);
    let (mu_n, off) = take_u64s(inputs, off, N_LIMBS + 1);
    let (x, off) = take_u64s(inputs, off, N_LIMBS);
    let (y, off) = take_u64s(inputs, off, N_LIMBS);
    let (pi, off) = take_u64s(inputs, off, N_LIMBS);
    let (T, off) = take_u64s(inputs, off, 2);
    let nonce: u16 = (*inputs.at(off)).try_into().unwrap();
    let ctx = *inputs.at(off + 1);
    let ct_len: usize = (*inputs.at(off + 2)).try_into().unwrap();
    let ct = inputs.slice(off + 3, ct_len);
    assert(inputs.len() == off + 3 + ct_len, 'input length');

    if !verify_vdf(N.span(), mu_n.span(), x.span(), y.span(), pi.span(), T.span(), nonce.into()) {
        return array![0, 0, 0, 0];
    }
    let pt = decrypt(derive_key(ctx, y.span()), ct);
    array![
        instance_hash(N.span(), x.span(), y.span(), T.span()), ctx,
        core::poseidon::poseidon_hash_span(ct), core::poseidon::poseidon_hash_span(pt.span()),
    ]
}
