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
  lib/   # bigint + Wesolowski verifier (u64 limbs, u128 column accumulation, Barrett)
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
scarb test                      # 6 passed: unit checks + full verify + forgery rejection

# 2. Run the full Wesolowski verification as an executable (gas off)
cd ../exec
python3 ../../scripts/write_exec_512.py   # writes vectors/exec_args_512.json
scarb execute --arguments-file ../../vectors/exec_args_512.json --print-program-output
# Program output: <instance_hash>   <- poseidon(N, x, y, T); 0 would mean rejection

# 3. Prove the execution with Stwo (works on a 16GB+ machine; 512-bit verified here)
scarb execute --output=standard --arguments-file ../../vectors/exec_args_512.json
scarb prove --execution-id 1   # -> target/execute/vdf_exec/execution1/proof/proof.json
scarb verify --proof-file target/execute/vdf_exec/execution1/proof/proof.json
```

On success the executable outputs `poseidon(N, x, y, T)` (a nonzero felt252),
not a boolean: the onchain verifier recomputes that hash from the public
instance, which binds the Stwo proof to the exact `(N, x, y, T)` it claims.

Measured execution (scarb 2.18, `--print-resource-usage`): **~14.8M Cairo
steps** for the 512-bit demo, **~45.7M** for the 1024-bit production config
(down from ~40M and ~160M before the bigint optimization pass — see
*Cairo bigint* below).

Proving was verified end-to-end on a 64GB box: the 512-bit execution (14.8M
steps, ~660MB prover input) proves in ~90s to a ~590MB proof. One caveat:
`scarb verify` currently panics on this trace shape with `ECDSA segment is
not empty` — upstream issue
[`stwo-cairo#1733`](https://github.com/starkware-libs/stwo-cairo/issues/1733)
(fix pending in PR #1666). The proof artifact itself is valid; verification
resumes once the pinned stwo-cairo is updated. For production use a 2048-bit
modulus (`N_LIMBS = 32`) and a serious proving machine.

## How it works

### Cairo bigint (the non-obvious part)

Cairo `felt252` lacks `PartialOrd`, `Rem`, and shifts, and `/` is field
division — so bigint math uses:

- **u64 limbs** (little-endian), **u128 column accumulation** for
  multiplication: each column sums its `u64×u64→u128` products in a
  `(hi, lo)` u128 pair (`u128_overflowing_add`), then folds into a 64-bit
  limb + carry. Visiting only contributing `(i, j)` pairs and truncating
  products whose high limbs are never read cut execution ~63% (40.0M → 14.8M
  steps for the 512-bit verify).
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
tx never runs the ~46M-step (1024-bit) computation, it just verifies a small proof.

## Status

- [x] Python reference + vectors (RSW + Wesolowski)
- [x] Cairo verifier (`lib`) + `#[executable]` entrypoint (`exec`)
- [x] In-program Fiat-Shamir challenge + Barrett constant validation (soundness)
- [x] Barrett k+1-limb normalization (HAC 14.42) + boundary regression tests
- [x] `scarb execute` → instance hash on 512/1024-bit vectors (14.8M / 45.7M steps)
- [x] `scarb execute --output=standard` → prover input
- [x] `scarb prove` (512-bit verified end-to-end; ~90s / 590MB proof on a 64GB box)
- [ ] `scarb verify` (blocked: upstream stwo-cairo#1733 ECDSA-segment panic; proof artifact is valid)
- [ ] Onchain recursive-verifier integration

## References

- STRK20 sealed-bid auction RFP: https://strk20.starknet.io/rfp/sealed-bid-auctions
- Stwo Cairo: https://github.com/starkware-libs/stwo-cairo (migrated to `starkware-libs/proving`)
- Scarb prove/verify: https://docs.swmansion.com/scarb/docs/extensions/prove-and-verify

License: Apache-2.0
