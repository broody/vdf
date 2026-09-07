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
  `y == π^L · x^(2^T mod L) (mod N)` with the Fiat-Shamir challenge
  `L = poseidon(N, x, y, T)` **derived inside the verifier** — a prover cannot
  choose `L` (or the Barrett constants) to fit a forged `y`; see *Security
  notes* below.

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
scarb test                      # 37 passed: arithmetic + full verify + forgery rejection

# 2. Run the full Wesolowski verification as an executable (gas off)
cd ../exec
python3 ../../scripts/write_exec_512.py   # writes vectors/exec_args_512.json
scarb execute --arguments-file ../../vectors/exec_args_512.json --print-program-output
# Program output: <instance_hash>   <- poseidon(N, x, y, T); 0 would mean rejection

# 3. Generate and verify a Stwo proof (optimized trace not yet benchmarked)
scarb execute --output=standard --arguments-file ../../vectors/exec_args_512.json
scarb prove --execution-id 1   # -> target/execute/vdf_exec/execution1/proof/proof.json
scarb verify --proof-file target/execute/vdf_exec/execution1/proof/proof.json
```

On success the executable outputs `poseidon(N, x, y, T)` (a nonzero felt252),
not a boolean: the onchain verifier recomputes that hash from the public
instance, which binds the Stwo proof to the exact `(N, x, y, T)` it claims.

Measured execution (Scarb/Cairo 2.18.0, dev profile, `stwo_no_ecop`,
`--output=none --print-resource-usage`, T=65536):

| Modulus | Before | Optimized | Step reduction |
|---|---:|---:|---:|
| 512-bit demo | 14,838,481 | 5,629,093 | 62.1% |
| 2048-bit | 155,970,789 | 49,429,902 | 68.3% |

Joint exponentiation, bounded felt accumulation and Barrett span slices preserve
both instance hashes. These are offchain execution measurements, not transaction
fees or measured proving-time savings. The 512-bit full-verification test's gas
estimate is ~0.72B, down from ~1.92B.

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
- **Per-modulus limb domains**: N (512-bit → 8 limbs) and the 256-bit challenge
  L (4 limbs) each get their own Barrett constant

`ref/cairo_model_v2.py` is the exact Python mirror of the Cairo algorithm —
validate any change there first (it caught two real bugs: borrow direction and
mixed limb domains).

### Security notes (read before building on this)

- **The challenge is derived in-program.** `L = poseidon(N, x, y, T)` (bit
  192 forced so `L ≥ 2^192` always satisfies Barrett's precondition). An
  earlier version accepted `L`, `r`, and both `mu` constants as inputs; that
  let a prover "verify" any claimed output (e.g. `L = 2^200 | 2^T`, `π = 1`,
  `r = 0` proves `y = 1`). `scripts/repro_issues.py` demonstrates the attack;
  `forged_output_rejected` in the test suite pins the fix. Like the Python
  reference, `L` is a random 251-bit integer, not a verified prime — see
  `docs/src/pages/concepts/wesolowski-vdf.mdx` for the soundness discussion.
- **Barrett keeps k+1 limbs through normalization** (HAC 14.42). A k-limb
  `r` wraps whenever `p − q3·N ≥ b^k` (reachable: `N = 2^512−1`,
  `a = 2^448−1`, `b = 2^448+1` was off by one) — and in the RSW setting `N`
  is chosen by the bidder, i.e. adversarially. Regression test:
  `barrett_boundary_top_limb_all_ones`; boundary sweeps:
  `scripts/hunt_barrett.py`.
- **The executable output binds the instance** (see Quick start). An onchain
  verifier must treat `(N, x, y, T)` as public inputs and compare the proof's
  public output against `poseidon(N, x, y, T)` — otherwise the proof is
  existential over the instance itself.

### Onchain verification (the end state)

The `exec/` package is a Cairo *program*; its execution is proven with Stwo and
verified with the Cairo recursive verifier (stwo-cairo ships one) — so a Starknet
tx would verify a proof of the ~49.4M-step (2048-bit) computation. This onchain
integration is still pending.

## Status

- [x] Python reference + vectors (RSW + Wesolowski)
- [x] Cairo verifier (`lib`) + `#[executable]` entrypoint (`exec`)
- [x] In-program Fiat-Shamir challenge + Barrett constant validation (soundness)
- [x] Barrett k+1-limb normalization (HAC 14.42) + boundary regression tests
- [x] `scarb execute` → instance hash on 512/2048-bit vectors (5.63M / 49.43M steps)
- [x] Historical `scarb execute --output=standard` → prover input
- [x] Historical `scarb prove` (512-bit, before latest optimizations; ~90s / 590MB)
- [ ] Prove and verify the optimized traces (historical verification hit stwo-cairo#1733)
- [ ] Onchain recursive-verifier integration

## References

- STRK20 sealed-bid auction RFP: https://strk20.starknet.io/rfp/sealed-bid-auctions
- Stwo Cairo: https://github.com/starkware-libs/stwo-cairo (migrated to `starkware-libs/proving`)
- Scarb prove/verify: https://docs.swmansion.com/scarb/docs/extensions/prove-and-verify

License: Apache-2.0
