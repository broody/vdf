use vdf::{verify_vdf, instance_hash, N_LIMBS, L_LIMBS};

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

/// Executable entrypoint: flat felt252 inputs [N(8), mu_n(9), x(8), y(8),
/// pi(8), mu_l(5), T(2)] for the 512-bit demo (N_LIMBS=8; 16 for 1024-bit
/// production). The challenge L and r = 2^T mod L are derived in-program.
///
/// On success returns instance_hash(N, x, y, T): the onchain verifier
/// recomputes it from the public instance, binding this proof to that exact
/// (N, x, y, T). Returns 0 on verification failure.
#[executable]
pub fn vdf_verify(inputs: Span<felt252>) -> felt252 {
    let (N, off) = take_u64s(inputs, 0, N_LIMBS);
    let (mu_n, off) = take_u64s(inputs, off, N_LIMBS + 1);
    let (x, off) = take_u64s(inputs, off, N_LIMBS);
    let (y, off) = take_u64s(inputs, off, N_LIMBS);
    let (pi, off) = take_u64s(inputs, off, N_LIMBS);
    let (mu_l, off) = take_u64s(inputs, off, L_LIMBS + 1);
    let (T, _) = take_u64s(inputs, off, 2);

    if verify_vdf(
        N.span(), mu_n.span(), x.span(), y.span(), pi.span(), mu_l.span(), T.span(), 2,
    ) {
        instance_hash(N.span(), x.span(), y.span(), T.span())
    } else {
        0
    }
}
