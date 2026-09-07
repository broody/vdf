# Whisper sealed-bid auction: integration design analysis

Companion to VDF_INTEGRATION_PLAN.md in the whisper repo (drafted 2026-08-28).
Captures the design reasoning so a future session doesn't re-derive it.

## The two orthogonal trust axes

Any private auction has TWO independent problems — never conflate them:

| Axis | Question | Mechanism |
|---|---|---|
| **Custody** | Who can MOVE funds? | Shadow accounts / contract-gated escrow (see strk20-privacy-integration skill) |
| **Reveal timing** | WHEN do amounts become readable? | This skill: RSW timelock + Wesolowski VDF |

- **Shadow accounts fix custody, not visibility.** The shipped `shadow_account_anonymizer` (monorepo `packages/`) runs dapp interactions per-commitment through dedicated contracts; spend is restricted to the configured privacy contract. But it holds PUBLIC open-note balances — amounts visible immediately. It cannot hold an encrypted exact-value note: the pool's note ownership IS possession of the viewing private key (`use_note(owner_private_key, ...)`), and a contract cannot hold a private key.
- **The structural fix for contract-held sealed notes is a pool primitive**: "predicate-escrowed notes" where spend authority = an app-contract callback instead of the viewing-key witness, amount staying encrypted in the note. This is basically what the sealed-bid RFP points at; it's a strong PR to `starkware-libs/starknet-privacy`.

## Why bidder-side reveal is bad UX (and the actual fix)

The "collective maximum then bid beneath it" pattern is NOT how whisper works — whisper escrows exact-value notes (the note amount IS the bid; additive tranches increase it). The real UX problem with bidder-side reveal:

- `settle_auction` requires `revealed_bids.len() == expected_len` — ALL accepted bids revealed or settlement is blocked until `abort_after`, then everyone refunds. One no-show (or a high bidder griefing) stalls the whole auction.
- The RFP's fix: **forfeit + force-reveal fallback** — unrevealed bids are forfeit (costs the griever their own escrow), and force-reveal via threshold auditing handles honest offline bidders. Cooperative bidders reveal themselves with a tiny tx.

## The eBay-UX constraint (one-action bidding)

Requirement: bid once, walk away, auction resolves, money comes back. That means force-reveal at the deadline with ZERO bidder action — which forces the "deadline decryption" problem: *what causes decryption capability to become available at the deadline?*

Options (docs table in `privacy-and-trust.mdx`): bidder reveal (bad UX), threshold committee (M-of-N trust, MPC-over-viewing-key is unsolved), TEE (hardware trust), timelock encryption (the clean one), MPC/private-auction-proof (heavy custom crypto).

**RSW timelock + Wesolowski VDF** gives eBay UX with no standing decryptor:
- Bidder generates RSA modulus N=p·q at bid time, knows φ(N) → derives key instantly via `y = x^(2^T mod φ(N)) mod N`, AEAD-encrypts the opening, discards factors.
- At deadline ANYONE computes `y = x^(2^T) mod N` (T sequential squarings — the delay), decrypts, and reveals permissionlessly. Even the auctioneer learns bids only at the deadline.
- Wesolowski proof lets a verifier check the squarings happened in O(log T) group ops: `y == π^L · x^(2^T mod L) (mod N)`, with the challenge `L = poseidon(N, x, y, T)` derived **inside the verifier program** (an input L was demonstrated forgeable in this repo — see SKILL.md pitfalls).
- **Settlement does not need a VDF proof per bid** (added 2026-09): the VDF guarantees the key is *available* to anyone at the deadline; settlement only needs a valid commitment opening. A solver decrypts offchain and submits the opening with a plain `reveal_bid` — the Stwo proof of the VDF verification is only for the permissionless/disputed path.

## Repo decision (as decided with the user)

**Standalone `vdf` crate, consumed by whisper** — not folded into whisper:
1. **Toolchain separation**: VDF needs Rust (Stwo prover, nightly, heavy build) + scarb 2.20; whisper is Cairo+TS+pnpm. Mixing breaks whisper's AGENTS.md workflow.
2. **It's not auction-specific**: RFP names the whole sealed-bid class (DAO grants, OTC, parameter auctions). Upstreamable to `starkware-libs/proving` (absorbed stwo-cairo).
3. **Whisper stays product-shaped** (contracts/sdk/operator/docs).

Integration split:
- **RSW encryption → whisper `sdk/`** (part of the bid flow; bidders shouldn't need a second package)
- **`reveal_bid` entrypoint → whisper `contracts/`** (additive; reuses existing commitment check; no pool-callback breakage — the canonical pool rejects trailing return values per AGENTS.md)
- **Post-deadline VDF eval + relay → whisper `operator/`** (operator becomes relayer, not standing decryptor)
- **Cairo verifier program + proving pipeline → standalone repo**

## Publishing note (2026-08-28)

Repo created at `broodling-bot/vdf` (public) — the broodling-bot PAT creates repos under the BOT's account, NOT under `broody` (separate user account, no create permission; org endpoint returns Not Found). If the repo should live under broody, the user must transfer it from their GitHub settings or create it under broody and re-point the remote. Local canonical checkout: `/root/vdf-repo`.
