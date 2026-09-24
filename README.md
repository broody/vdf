# VDF — Verifiable Delay Function timelock for Starknet (Cairo)

Wesolowski VDF verification + RSW timelock decryption, implemented as a provable
Cairo **executable** program (not a contract). Run the verification offchain,
prove it with the Stwo prover, and verify the proof onchain cheaply — the
Starknet-native pattern for "decryption capability materializes at a deadline"
with **no standing decryptor, no committee, no TEE**.

## Why

Sealed-bid auctions (and any "reveal at deadline" mechanism) need bids hidden
until a deadline with no trusted party holding the reveal key. This repo
provides the cryptographic primitive:

- **RSW timelock**: the encryptor (bidder) knows φ(N), so they derive the
  decryption key *instantly* (`y = x^(2^T mod φ(N)) mod N`). Everyone else must
  compute `y = x^(2^T) mod N` by **T sequential squarings** — the delay. The
  key materializes only at ~time T. No trusted setup: the bidder generates and
  discards the RSA factors.
- **Wesolowski VDF proof**: once someone has done the squarings, a short proof
  lets a verifier check correctness in **O(log T)** group ops:
  `π^L · x^(2^T mod L) == ±y (mod N)` in the signed group `Z_N^*/{±1}`, with
  the claimed output in sign-canonical form `|y| = min(y, N−y)`. The
  Fiat-Shamir challenge `L` is a 251-bit **prime** hashed from
  `(N, x, |y|, T, nonce)` and **checked inside the verifier** — a prover cannot
  choose `L` (or the Barrett constant) to fit a forged `y`; see *Security
  notes* below.
- **In-program decryption**: the executable also decrypts the RSW ciphertext
  under a key derived from `|y|` and outputs `poseidon(pt)` for comparison
  against the onchain commitment.

## Repository layout

```
cairo/
  lib/   # bigint + Wesolowski verifier (u64 limbs, bounded felt column accumulation, Barrett)
  exec/  # #[executable] vdf_verify entrypoint, depends on lib
ref/     # Python reference + exact Cairo-algorithm model + Poseidon (cross-checked with Cairo)
vectors/ # test vectors (RSW keys + Wesolowski proofs)
scripts/ # generators: vectors, Cairo tests, executable args; bug repros
skill/   # the Hermes skill for this work (SKILL.md + references + scripts)
```

`skill/` is the Hermes Agent skill `vdf-timelock-cairo` — the operational playbook
for this exact pipeline (scheme, Cairo bigint pattern, scarb 2.20 executable
setup, proving gotchas), plus `references/whisper-integration-design.md`
capturing the sealed-bid auction design reasoning. Install it with
`hermes skills install` or copy `skill/SKILL.md` to your skills directory.

The Cairo `lib` and `exec` are **two standalone packages** (a scarb workspace
forces `[cairo]` settings at the workspace level, which conflicts with the
executable target's `enable-gas = false`).

## Documentation

The `docs/` directory is a [Vocs](https://vocs.dev) site (same stack as
Whisper's docs). Run it locally:

```sh
cd docs
pnpm install
pnpm dev          # local dev server
pnpm build        # static site output in docs/dist/
```

## Quick start

Requires scarb ≥ 2.20 (ships `scarb execute` / `scarb prove` / `scarb verify`).

```sh
# 1. Validate the bigint math and the security properties
cd cairo/lib
scarb test                      # 42 passed: 31 arithmetic + 11 verifier (full verify, forgery rejection)

# 2. Run the full Wesolowski verification as an executable (gas off)
cd ../exec
python3 ../../scripts/write_exec_args.py  # writes vectors/exec_args_{512,2048}.json
scarb execute --arguments-file ../../vectors/exec_args_512.json --print-program-output
# Program output: <instance_hash> <ctx> <poseidon(ct)> <poseidon(pt)>; [0, 0, 0, 0] means rejection

# 3. Generate and verify a Stwo proof (optimized trace not yet benchmarked)
scarb execute --output=standard --arguments-file ../../vectors/exec_args_512.json
scarb prove --execution-id 1   # -> target/execute/vdf_exec/execution1/proof/proof.json
scarb verify --proof-file target/execute/vdf_exec/execution1/proof/proof.json
```

The executable input is
`[N(k), mu_n(k+1), x(k), y(k), pi(k), T(2), nonce, ctx, ct_len, ct(ct_len)]`
(k = 8 for 512-bit, 32 for 2048-bit; the length must match exactly). `y` is
the sign-canonical output, `nonce` (< 2^16) selects the prime challenge and
`ctx` is an application context felt (e.g. `H(auction, bid)`). On success it
outputs `[instance_hash, ctx, poseidon(ct), poseidon(pt)]`, not a boolean:
the onchain side recomputes `instance_hash` and `poseidon(ct)` from the stored
puzzle and ciphertext, and compares `poseidon(pt)` with the bid commitment.
This binds the Stwo proof to the exact `(N, x, |y|, T)` and ciphertext it
claims. On verification failure it outputs `[0, 0, 0, 0]`.

Measured execution (Scarb/Cairo 2.18.0, dev profile, `stwo_no_ecop`,
`--output=none --print-resource-usage`, T=65536), current program:

| Modulus | Steps |
|---|---:|
| 512-bit demo | 8,491,007 |
| 2048-bit | 53,897,839 |

The 20 Miller–Rabin rounds on the challenge plus in-program decryption add
about 2.9M (512-bit) and 4.5M (2048-bit) steps over the pre-hash-to-prime
program. Historical optimization pass, measured **before** hash-to-prime and
decryption:

| Modulus | Before | Optimized | Step reduction |
|---|---:|---:|---:|
| 512-bit demo | 14,838,481 | 5,629,093 | 62.1% |
| 2048-bit | 155,970,789 | 49,429,902 | 68.3% |

Joint exponentiation, bounded felt accumulation and Barrett span slices
preserved both instance hashes. These are offchain execution measurements, not transaction fees or measured
proving-time savings. The 512-bit full-verification test's gas estimate is
~1.09B (~0.72B before hash-to-prime, ~1.92B before the optimization pass).

To reproduce the 2048-bit execution, set `N_LIMBS = 32` and run the executable
with `--arguments-file ../../vectors/exec_args_2048.json`. Restore `N_LIMBS = 8`
for the default full-verification tests. Both benchmark vectors are checked in.

A **historical** 512-bit run, before these optimizations (14.8M steps), generated
~660MB of prover input and a ~590MB proof in ~90s on a 64GB box. That run hit
`ECDSA segment is not empty` during `scarb verify`, recorded in
[`stwo-cairo#1733`](https://github.com/starkware-libs/stwo-cairo/issues/1733).
The optimized traces have been executed and tested, but have not yet been
proven or proof-verified; the historical proving figures do not describe them.
For production use a 2048-bit modulus (`N_LIMBS = 32`).

## How it works

### Cairo bigint (the non-obvious part)

Cairo `felt252` lacks `PartialOrd`, `Rem`, and shifts, and `/` is field
division — so bigint math uses:

- **u64 limbs** (little-endian), **bounded felt column accumulation**:
  products are accumulated exactly in `felt252`, then split once per column
  into u128 halves for limb extraction and carry. With at most 33 products
  per column at 2048 bits, the sum including carry is below `2^134`, far below
  the field modulus. Integer division and borrow arithmetic remain in u128.
- **Joint exponentiation**: scan the bits of L and r together, sharing
  squarings and selecting a factor from `1`, pi, x or precomputed `pi*x`.
  Bases are reduced before precomputation, including noncanonical inputs.
- **Span slices** for Barrett's q1 and q3 avoid per-reduction array copies;
  q2 still computes all columns because low carries affect its high half.
- **u128 subtraction with borrow** (borrow applies to the subtrahend)
- **Barrett reduction** with `mu = floor(b^(2k)/N)` computed offchain (k+1
  limbs) and **validated in-program** (`mu·m ≤ b^2k < (mu+1)·m`)
- **Per-modulus limb domains**: N (512-bit → 8 limbs) uses Barrett with its
  own constant `mu_n`; arithmetic mod the 251-bit prime challenge L
  (`2^T mod L`, Miller–Rabin) uses Cairo's native u256 with u512-by-u256
  division, so there is no Barrett constant for L. L enters the joint
  exponentiation as a 4-limb exponent.

`ref/cairo_model_v2.py` is the exact Python mirror of the Cairo algorithm —
validate any change there first (it caught two real bugs: borrow direction and
mixed limb domains).

### Security notes (read before building on this)

- **The challenge is derived in-program.** An earlier version accepted `L`,
  `r`, and both `mu` constants as inputs; that let a prover "verify" any
  claimed output (e.g. `L = 2^200 | 2^T`, `π = 1`, `r = 0` proves `y = 1`).
  `scripts/repro_issues.py` demonstrates the attack; `forged_output_rejected`
  in the test suite pins the fix.
- **The challenge is a prime (hash-to-prime).** A later version derived
  `L = poseidon(N, x, y, T)` with bit 192 forced — a random *composite*
  integer. That is forgeable without factoring `N`: pick `α ≡ 2^T mod M` for
  a smooth `M` (cheap by CRT), grind `y = x^e`, `e = α+kM`, until
  `L = H(..y..)` divides `M`, then `π = x^((e−r)/L)` proves a wrong `y`
  (`scripts/repro_issues.py` section D shows it with a 48-bit truncated
  challenge). Now:
  `instance = poseidon('VDF_WES_INSTANCE_V1', len(N limbs), N ‖ x ‖ |y| ‖ T)`
  (u64 limbs, T as 2 limbs),
  `seed = poseidon('VDF_WES_CHALLENGE_V1', instance, nonce)` with a
  prover-chosen `nonce < 2^16` (any that gives a prime; the reference prover takes the first, ~87 tries
  expected), and `L` = the low 250 bits of `seed` with bits 250 and 0 forced
  (odd, `2^250 ≤ L < 2^251`). `L` must pass 20 Miller–Rabin rounds in-program
  with Fiat-Shamir bases
  `(poseidon('VDF_WES_MILLER_RABIN_V1', seed, 2i) + 2^256·poseidon(…, 2i+1)) mod L`;
  a composite random odd 251-bit candidate passes with probability `< 2^-107`
  (Damgård–Landrock–Pomerance average-case bound), so grinding nonces for a
  composite challenge is infeasible. Tests: `challenge_is_prime_and_r_calc`,
  `composite_candidate_fails_primality`, `composite_challenge_rejected`.
- **The group is signed (`Z_N^*/{±1}`).** In the full group anyone who knew
  the honest `y` could "prove" `N−y` too (negate `π` for an odd challenge) —
  two valid outputs per instance. The verifier requires `y < N−y` (and
  `y < N`, `N` odd) and checks `π^L · x^r == ±y`; `π` is therefore only
  defined up to sign (`−π` also verifies). Tests: `negated_output_rejected`,
  `negated_proof_verifies`; repro: `scripts/repro_issues.py` section C.
- **`mu_n` is validated in-program** (`mu·N ≤ b^2k < (mu+1)·N`). It is the
  only Barrett constant supplied; `mu_valid` remains a generic function.
- **Barrett keeps k+1 limbs through normalization** (HAC 14.42). A k-limb
  `r` wraps whenever `p − q3·N ≥ b^k` (reachable: `N = 2^512−1`,
  `a = 2^448−1`, `b = 2^448+1` was off by one) — and in the RSW setting `N`
  is chosen by the bidder, i.e. adversarially. Regression test:
  `barrett_boundary_top_limb_all_ones`; boundary sweeps:
  `scripts/hunt_barrett.py`.
- **The executable output binds the instance** (see Quick start). An onchain
  verifier must treat `(N, x, |y|, T)` and the ciphertext as public inputs and
  compare the proof's four outputs against its own `instance_hash`, `ctx`,
  `poseidon(ct)` and the posted `poseidon(pt)` commitment — otherwise the
  proof is existential over the instance itself.
- **Decryption has no integrity tag by design.** `key = poseidon('VDF_RSW_KEY_V1',
  ctx, len(y limbs), |y| limbs)`; the Poseidon keystream is
  `ct_i = pt_i + poseidon(key, i)` over the Stark field (messages are felt252
  arrays). The application binds the plaintext by comparing `poseidon(pt)`
  against the commitment posted with the ciphertext. Test:
  `rsw_keystream_roundtrip`.
- **The trapdoor holder can prove anything on their own puzzle.** With a
  bidder-generated `N` the bidder knows φ(N) and can produce a valid proof for
  **any** `y`, so a proof is only sound against parties without the trapdoor.
  Applications must treat two different proven outputs for one instance as
  evidence of the trapdoor holder's misbehaviour (e.g. void and slash), not
  accept the first proof as final.

### Onchain verification (the end state)

The `exec/` package is a Cairo *program*; its execution is proven with Stwo and
verified with the Cairo recursive verifier (stwo-cairo ships one) — so a Starknet
tx would verify a proof of the ~53.9M-step (2048-bit) computation. This onchain
integration is still pending.

## Status

- [x] Python reference + vectors (RSW + Wesolowski)
- [x] Cairo verifier (`lib`) + `#[executable]` entrypoint (`exec`)
- [x] In-program Fiat-Shamir challenge + Barrett constant validation (soundness)
- [x] Hash-to-prime challenge (20 Miller–Rabin rounds) + signed group `Z_N^*/{±1}`
- [x] In-program RSW decryption (Poseidon keystream) + commitment output
- [x] Barrett k+1-limb normalization (HAC 14.42) + boundary regression tests
- [x] `scarb execute` → 4-felt output on 512/2048-bit vectors (8.49M / 53.90M steps)
- [x] Historical `scarb execute --output=standard` → prover input
- [x] Historical `scarb prove` (512-bit, before latest optimizations; ~90s / 590MB)
- [ ] Prove and verify the optimized traces (historical verification hit stwo-cairo#1733)
- [ ] Onchain recursive-verifier integration

## References

- STRK20 sealed-bid auction RFP: https://strk20.starknet.io/rfp/sealed-bid-auctions
- Stwo Cairo: https://github.com/starkware-libs/stwo-cairo (migrated to `starkware-libs/proving`)
- Scarb prove/verify: https://docs.swmansion.com/scarb/docs/extensions/prove-and-verify

License: Apache-2.0
