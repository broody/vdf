---
name: vdf-timelock-cairo
description: Implement VDF (Verifiable Delay Function) timelock decryption and deadline-decryption mechanisms on Starknet — Wesolowski VDF verification in Cairo, RSW timelock encryption, and the scarb 2.20 executable + Stwo proving pipeline (scarb execute → scarb prove → scarb verify). Use when building "decrypt at deadline" primitives, sealed-bid auction reveals, or any Cairo program that must be proven offchain and verified onchain. Covers the validated u64/u128 bigint pattern for Cairo 2.20 (felt252 lacks PartialOrd/Rem/shifts), HAC 14.42 Barrett reduction with in-program constant validation, in-program Fiat-Shamir challenge derivation, and the exact toolchain gotchas (gas caps, OOM, arg formats).
tags: [cairo, starknet, vdf, timelock, wesolowski, rsw, stwo, proving, cryptography]
related_skills: [strk20-privacy-integration, ocr-and-documents]
---

# VDF Timelock Decryption on Starknet (Cairo executable + Stwo proof)

## When to use

- Sealed-bid auctions where bid amounts must stay hidden **until a deadline** with **no standing decryptor** (no single operator holding the reveal key).
- Any "decryption capability materializes at time T" requirement on Starknet.
- Proving an expensive Cairo computation offchain and verifying it onchain cheaply (the Starknet-native pattern: `#[executable]` program → Stwo proof → Cairo recursive verifier).

The two halves are orthogonal:
- **Custody** (who can move funds) → shadow accounts, contract-gated escrow (see `strk20-privacy-integration` skill).
- **Reveal timing** (when amounts become readable) → this skill: RSW timelock + Wesolowski VDF.

For the full sealed-bid-auction design analysis (why bidder-side reveal griefs, the eBay-UX constraint, repo split decision, GitHub publishing notes), read `references/whisper-integration-design.md`.

## The scheme (validated, do not re-derive)

### RSW timelock encryption (encryptor fast, everyone else slow)
- Encryptor generates RSA modulus N = p·q, knows φ(N).
- Derives the timelock key instantly: y = x^(2^T mod φ(N)) mod N (O(log T)).
- Everyone else must compute y = x^(2^T) mod N by **T sequential squarings** — this is the delay.
- No trusted setup needed per-message: the encryptor is the only one who ever needs the factors, and can discard them after encryption. The bidder keeps their own bid secret from everyone including the auctioneer.

### Wesolowski VDF proof (verifier fast, O(log T) group ops)
- Prover computes y = x^(2^T) mod N (the T squarings) plus:
  - L = poseidon(N, x, y, T) — Fiat-Shamir challenge, **derived inside the verifier program** (bit 192 forced so L ≥ 2^192 always satisfies the Barrett precondition). Never accept L as an input: with a freely chosen L the "proof" verifies any claimed y (L = 2^200 divides 2^T ⇒ r = 0, π = 1 proves y = 1). This was a real, demonstrated forgery — see scripts/repro_issues.py.
  - q = floor(2^T / L), r = 2^T mod L
  - π = x^q mod N
- Verifier checks: **y == π^L · x^(2^T mod L) (mod N)** — needs only two modpows with ~256-bit exponents.
- The Barrett constants (mu for N and for L) are computed offchain but **validated in-program** (mu·m ≤ b^2k < (mu+1)·m) — an unvalidated mu is an attacker-controlled reduction.
- The executable outputs poseidon(N, x, y, T) on success (0 on rejection); the onchain verifier recomputes it from the public instance, binding the Stwo proof to that exact (N, x, y, T).
- L is a random 251-bit integer, not a verified prime (same as the Python reference; cf. Pietrzak 2018/627 §4). The known attacks with a composite challenge require factoring an RSA-sized integer. For theorem-level soundness, generate a prime L offchain and check a primality certificate in-program.

Reference implementation: `ref/vdf_reference.py` (RSW + Wesolowski, generates vectors). Python model of the Cairo algorithm: `ref/cairo_model_v2.py` (validate ANY algorithm change here first — it caught two real bugs before Cairo did). Poseidon mirror: `ref/poseidon.py` (cross-checked against the Cairo builtin via a probe executable).

## Cairo 2.20 bigint pattern (validated)

**Never use felt252 for bigint math.** felt252 has no `PartialOrd`, no `Rem`, no shifts, and `/` is field division. Instead:

- **Limbs**: `u64`, little-endian.
- **Multiplication**: schoolbook columns, each accumulated in a `(hi, lo)` u128 pair — `u64×u64→u128` products are native, column overflows count into `hi` via `u128_overflowing_add` (Cairo 2.18+: `core::num::traits::OverflowingAdd::overflowing_add`; the older `core::integer::u128_overflowing_add` is deprecated and returns a Result, not a tuple). Visit only contributing `(i, j)` pairs and truncate products to the low columns the caller reads (Barrett's q3·N needs only k+1). This took the 512-bit verify from 40.0M to ~14.8M Cairo steps vs generic u256 accumulation.
- **Carry**: `lo & 0xffffffffffffffff` → u64 limb, `hi·2^64 + lo/2^64` → carry — constant-divisor division on u128 is integer division.
- **Subtraction with borrow**: compute in `u128` (`minuend >= subtrahend` branch; add 2^64 when borrow needed). The borrow goes on the **subtrahend**, not the minuend — this was a real bug.
- **Bit scan in modpow**: `e % 2` / `e / 2` loop (u64), not shifts; skip trailing zero limbs/bits of the exponent.
- **Barrett reduction (HAC 14.42)**: mu = floor(b^(2k)/N) with k+1 limbs, computed offchain, **validated in-program**. Keep **k+1 limbs** in r = p − q3·N through normalization — the quotient estimate undershoots by ≤2 so r < 3N < b^(k+1), and truncating to k limbs wraps for adversarially chosen N (e.g. N = 2^512−1, a = 2^448−1, b = 2^448+1 → off by one; the RSW bidder chooses N, so this is reachable). Each modulus gets its own limb domain (N=16 limbs for 1024-bit, L=4 limbs) — mixing domains broke r_calc with an infinite normalize loop.

Core file: `lib/src/lib.cairo` — `bigint_mul`, `ge_limbs`, `sub_limbs`, `barrett_reduce`, `modmul`, `modpow`, `limbs_eq`, `mu_valid`, `instance_hash`, `derive_challenge`, `verify_vdf`. All pub.

## Toolchain: scarb 2.20 executable + Stwo proving

### Setup (gotchas are real)
```
# exec/Scarb.toml — executable package
[[target.executable]]          # REQUIRED
[cairo]
enable-gas = false             # REQUIRED — executable target refuses gas
[dependencies]
cairo_execute = ">=2.20.0"     # the plugin that provides #[executable]
```

- **Workspaces break gas config**: per-package `[cairo]` and `[profile]` are ignored inside a workspace — only the workspace manifest's `[cairo]` applies. Split lib (tests need gas) and exec (gas off) into **two standalone packages** with a path dependency: `vdf = { path = "../lib" }`. No root workspace file at all.
- **Cairo 2.20 prelude**: do NOT `use array::ArrayTrait;` or `use option::OptionTrait;` — they're in the prelude; importing by those paths fails with E0006.
- **No closures** capturing mutable state — use a plain helper returning a tuple.
- **Tests**: `cairo-test` caps gas at 2^32. Post-optimization the full 512-bit verify is ~1.9B gas and fits; 1024-bit (16-limb modpows) does not → validate via `scarb execute` (gas off), not tests. `scarb cairo-test` is deprecated → snforge eventually.

### The proving pipeline
```
scarb execute --output=standard --arguments-file args.json
  → target/execute/<pkg>/executionN/prover_input.json
scarb prove --execution-id N
  → target/execute/<pkg>/executionN/proof/proof.json
scarb verify <path to proof.json>
```

- **Arguments file**: JSON array of **hex strings** (`0x...`), with a leading length felt for `Span<felt252>` params (serde reads length first). Decimal strings rejected; `Span` params need `[len, elem0, ...]`.
- **`--output=standard` is required** to emit prover_input.json; default is None (empty dirs, no error). `--output=cairo-pie` only works for the bootloader target.
- **OOM is the big one**: Stwo needs ~16GB RAM to prove at these step counts (~15M steps for the 512-bit verify, ~46M for 1024-bit). Use a 512-bit test modulus (8 limbs, N_LIMBS=8) on small machines; prove the real thing on a 16GB+ box. The execution dir artifacts are portable — copy to a bigger machine and `scarb prove --execution-id N` there.
- `scarb-prove` binary needs `SCARB_TARGET_DIR` and `SCARB_PROFILE=dev` exported when invoked directly (the `scarb` shell sets them).

### Verification onchain (the end state)
- Rust verifier for offchain/CI.
- **Cairo verifier** (stwo-cairo ships one) — runs inside Cairo VM → recursive proving, onchain verification. This is what makes "prove offchain, verify onchain" work: the heavy computation never runs in a tx; the tx just verifies a small proof.

## Repo layout (standalone, per project decision)
```
vdf/
  lib/   # bigint + verify_vdf (pub, gas-enabled for tests)
  exec/  # #[executable] vdf_verify → depends on lib via path
  ref/   # python reference + cairo_model_v2.py + vectors + generators
```

## Workflow when extending

1. Change the Python model (`cairo_model_v2.py`) FIRST, validate against vectors + random difftests.
2. Port to `lib/src/lib.cairo` mirroring limb-for-limb.
3. Validate with `scarb execute` (not tests) on the smallest vector that fits memory.
4. Regenerate test/args files with the generator scripts.
5. Prove on a big machine; keep execution-dir artifacts portable.

## Pitfalls checklist

- [ ] borrow goes on subtrahend, not minuend
- [ ] per-modulus limb domains, mu has k+1 limbs
- [ ] felt252: no PartialOrd/Rem/shifts — u64/u128 only
- [ ] **challenge L derived in-program** (poseidon of the instance), never an input — a chosen L forges any y
- [ ] **mu constants validated in-program** (mu·m ≤ b^2k < (mu+1)·m) — a bogus mu is an attacker-controlled reduction
- [ ] **Barrett r keeps k+1 limbs** through normalization (HAC 14.42) — k-limb r wraps for adversarial N
- [ ] executable outputs poseidon(N, x, y, T); onchain verifier compares it against the public instance
- [ ] workspace overrides per-package [cairo] — two standalone packages instead
- [ ] no ArrayTrait/OptionTrait imports in 2.20
- [ ] args = hex strings + leading length felt
- [ ] `--output=standard` or no prover_input.json
- [ ] cairo-test gas cap 2^32 — 512-bit full verify fits post-optimization; 1024-bit needs scarb execute
- [ ] OOM on small machines — 512-bit modulus for CI/demo, 1024+ for real
