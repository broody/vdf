---
name: vdf-timelock-cairo
description: Implement VDF (Verifiable Delay Function) timelock decryption and deadline-decryption mechanisms on Starknet — Wesolowski VDF verification in Cairo, RSW timelock encryption, and the scarb 2.20 executable + Stwo proving pipeline (scarb execute → scarb prove → scarb verify). Use when building "decrypt at deadline" primitives, sealed-bid auction reveals, or any Cairo program that must be proven offchain and verified onchain. Covers the validated u64/bounded-felt bigint pattern for Cairo 2.20 (felt252 lacks PartialOrd/Rem/shifts), HAC 14.42 Barrett reduction with in-program constant validation, in-program Fiat-Shamir challenge derivation, and the exact toolchain gotchas (gas caps, OOM, arg formats).
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
- Verifier checks: **y == π^L · x^(2^T mod L) (mod N)** — shares squarings with joint exponentiation over two ~256-bit exponents.
- The Barrett constants (mu for N and for L) are computed offchain but **validated in-program** (mu·m ≤ b^2k < (mu+1)·m) — an unvalidated mu is an attacker-controlled reduction.
- The executable outputs poseidon(N, x, y, T) on success (0 on rejection); the onchain verifier recomputes it from the public instance, binding the Stwo proof to that exact (N, x, y, T).
- L is a random 251-bit integer, not a verified prime (same as the Python reference; cf. Pietrzak 2018/627 §4). The known attacks with a composite challenge require factoring an RSA-sized integer. For theorem-level soundness, generate a prime L offchain and check a primality certificate in-program.

Reference implementation: `ref/vdf_reference.py` (RSW + Wesolowski, generates vectors). Python model of the Cairo algorithm: `ref/cairo_model_v2.py` (validate ANY algorithm change here first — it caught two real bugs before Cairo did). Poseidon mirror: `ref/poseidon.py` (cross-checked against the Cairo builtin via a probe executable).

## Cairo 2.20 bigint pattern (validated)

Use felt252 only for **bounded exact accumulation**: it has no `PartialOrd`, no `Rem`, no shifts, and `/` is field division. Keep comparisons, division and borrow arithmetic in integer types:

- **Limbs**: `u64`, little-endian.
- **Multiplication**: schoolbook columns accumulated in felt252, split once per column into a u256's u128 halves. With n contributing u64 products, carry < n*2^64 and acc < n*2^128. Since usize is u32, acc < 2^160 < the Cairo field modulus and carry fits u128; at 2048 bits n <= 33, so acc < 2^134. Preserve these bounds. Visit only contributing pairs and truncate only when callers consume low columns.
- **Joint exponentiation**: reduce both k-limb bases, precompute pi*x, then scan L and r together from the most significant bit. Share each squaring and select from 1, pi, x, pi*x. Copy the first nonzero factor; both zero exponents return one. The separate modpow computes 2^T mod L.
- **Barrett slices**: q1 and q3 use Span::slice, avoiding copies. Do not omit the low arithmetic columns of q1*mu: their carry feeds the retained high half.
- **Carry**: after splitting the column, `lo & 0xffffffffffffffff` → u64 limb, `hi·2^64 + lo/2^64` → carry — constant-divisor division on u128 is integer division.
- **Subtraction with borrow**: compute in `u128` (`minuend >= subtrahend` branch; add 2^64 when borrow needed). The borrow goes on the **subtrahend**, not the minuend — this was a real bug.
- **Bit scan in modpow**: `e % 2` / `e / 2` loop (u64), not shifts; skip trailing zero limbs/bits of the exponent.
- **Barrett reduction (HAC 14.42)**: mu = floor(b^(2k)/N) with k+1 limbs, computed offchain, **validated in-program**. Keep **k+1 limbs** in r = p − q3·N through normalization — the quotient estimate undershoots by ≤2 so r < 3N < b^(k+1), and truncating to k limbs wraps for adversarially chosen N (e.g. N = 2^512−1, a = 2^448−1, b = 2^448+1 → off by one; the RSW bidder chooses N, so this is reachable). Each modulus gets its own limb domain (N=16 limbs for 1024-bit, L=4 limbs) — mixing domains broke r_calc with an infinite normalize loop.

Core file: `cairo/lib/src/lib.cairo` — `bigint_mul`, `ge_limbs`, `sub_limbs`, `barrett_reduce`, `modmul`, `modpow`, `limbs_eq`, `mu_valid`, `instance_hash`, `derive_challenge`, `verify_vdf`, plus the private `joint_modpow` helper.

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
- **Tests**: `cairo-test` caps gas at 2^32. The full 512-bit verify is ~0.72B estimated gas; all 37 tests pass, including 31 arithmetic checks against Python integers. Use `scarb execute` (gas off) for full larger-modulus verification. `scarb cairo-test` is deprecated → snforge eventually.

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
- **Proving resources**: the optimized executions use 5.63M steps at 512 bits and 49.43M at 2048 bits (Scarb 2.18, T=65536), reductions of 62.1% and 68.3%. Proving these optimized traces has not been benchmarked. Historical 14.8M-step proving used a 64GB box; do not reuse its memory/time figures as measurements of the new trace. Use 512 bits for CI/demo and 2048 bits for production.
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
- [ ] felt column bound documented and preserved; comparisons/division/borrow stay in integers
- [ ] **challenge L derived in-program** (poseidon of the instance), never an input — a chosen L forges any y
- [ ] **mu constants validated in-program** (mu·m ≤ b^2k < (mu+1)·m) — a bogus mu is an attacker-controlled reduction
- [ ] **Barrett r keeps k+1 limbs** through normalization (HAC 14.42) — k-limb r wraps for adversarial N
- [ ] executable outputs poseidon(N, x, y, T); onchain verifier compares it against the public instance
- [ ] workspace overrides per-package [cairo] — two standalone packages instead
- [ ] no ArrayTrait/OptionTrait imports in 2.20
- [ ] args = hex strings + leading length felt
- [ ] `--output=standard` or no prover_input.json
- [ ] cairo-test gas cap 2^32 — 512-bit full verify fits; use scarb execute for larger full-verification runs
- [ ] OOM on small machines — 512-bit modulus for CI/demo, 2048-bit for production
