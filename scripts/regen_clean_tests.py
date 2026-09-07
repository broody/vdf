#!/usr/bin/env python3
"""Regenerate a clean tests.cairo: 3 gas-cheap passing tests + 1 ignored full-verify."""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
sys.path.insert(0, os.path.join(ROOT, "ref"))
from gen_cairo_tests import barrett_mu  # noqa: E402
import cairo_model_v2 as m  # noqa: E402

v = json.load(open(os.path.join(ROOT, "vectors", "vdf_vectors_512.json")))[0]
kN = v["n_limbs"]
kL = 4
N_int = int(v["N"], 16)
L_int = int(v["L"], 16)
muN = barrett_mu(N_int, kN)
muL = barrett_mu(L_int, kL)
T = v["T"]
T_limbs = [T & ((1 << 64) - 1), (T >> 64) & ((1 << 64) - 1)]


def arr(limbs):
    return "array![" + ", ".join(str(x) for x in limbs) + "].span()"


xsq = m.modmul(v["x_limbs"], v["x_limbs"], v["N_limbs"], muN, kN)
r_calc = m.modpow(m.to_limbs(2, kL), m.to_limbs(T, 2), 2, m.to_limbs(L_int, kL), muL, kL)

out = []
out.append("use vdf::{modmul, modpow, limbs_eq, verify_vdf, N_LIMBS, L_LIMBS};")
out.append("")
out.append("#[cfg(test)]")
out.append("mod tests {")
out.append("    use super::*;")
out.append("")
out.append("    // Gas-cheap checks that exercise both limb domains, plus the full")
out.append("    // verification (positive + forged-negative) which now fits the gas cap.")
out.append("    #[test]")
out.append("    fn modmul_x_squared() {")
out.append(f"        let N = {arr(v['N_limbs'])};")
out.append(f"        let mu_n = {arr(muN)};")
out.append(f"        let x = {arr(v['x_limbs'])};")
out.append("        let got = modmul(x, x, N, mu_n, N_LIMBS);")
out.append(f"        let exp = {arr(xsq)};")
out.append("        assert(limbs_eq(got.span(), exp, N_LIMBS), 'x^2 mod N');")
out.append("    }")
out.append("")
out.append("    // Challenge derivation (in-program Fiat-Shamir) + r = 2^T mod L.")
out.append("    #[test]")
out.append("    fn r_calc_2_to_T_mod_L() {")
out.append(f"        let N = {arr(v['N_limbs'])};")
out.append(f"        let x = {arr(v['x_limbs'])};")
out.append(f"        let y = {arr(v['y_limbs'])};")
out.append(f"        let mu_l = {arr(muL)};")
out.append(f"        let T_limbs = {arr(T_limbs)};")
out.append("        let L = vdf::derive_challenge(N, x, y, T_limbs);")
out.append(f"        let exp_L = {arr(v['L_limbs'][:4])};")
out.append("        assert(limbs_eq(L.span(), exp_L, L_LIMBS), 'L = H(N,x,y,T)');")
out.append("        let mut two = ArrayTrait::new();")
out.append("        two.append(2);")
out.append("        let mut i: usize = 1;")
out.append("        while i < L_LIMBS {")
out.append("            two.append(0);")
out.append("            i += 1;")
out.append("        }")
out.append("        let got = modpow(two.span(), T_limbs, 2, L.span(), mu_l, L_LIMBS);")
out.append(f"        let exp = {arr(r_calc)};")
out.append("        assert(limbs_eq(got.span(), exp, L_LIMBS), '2^T mod L');")
out.append("    }")
out.append("")
out.append("    // Regression for the Barrett k-limb truncation bug: with N = 2^512-1,")
out.append("    // a = 2^448-1, b = 2^448+1 the quotient estimate undershoots by 1 and")
out.append("    // p - q3*N = b^k + (2^384 - 2) >= b^k; keeping only k limbs wrapped r and")
out.append("    // returned p mod N - 1. HAC 14.42 keeps k+1 limbs and gets 2^384 - 1.")
out.append("    #[test]")
out.append("    fn barrett_boundary_top_limb_all_ones() {")
out.append("        let f = 18446744073709551615;")
out.append("        let N = array![f, f, f, f, f, f, f, f].span();")
out.append("        let mu_n = array![1, 0, 0, 0, 0, 0, 0, 0, 1].span();")
out.append("        let a = array![f, f, f, f, f, f, f, 0].span();")
out.append("        let b = array![1, 0, 0, 0, 0, 0, 0, 1].span();")
out.append("        let got = modmul(a, b, N, mu_n, N_LIMBS);")
out.append("        let exp = array![f, f, f, f, f, f, 0, 0].span(); // 2^384 - 1")
out.append("        assert(limbs_eq(got.span(), exp, N_LIMBS), 'barrett boundary');")
out.append("    }")
out.append("")
out.append("    // The Barrett constants are validated in-program: a tampered mu")
out.append("    // (here: top limb decremented) must be rejected in both domains.")
out.append("    #[test]")
out.append("    fn mu_validation_accepts_and_rejects() {")
out.append(f"        let N = {arr(v['N_limbs'])};")
out.append(f"        let mu_n = {arr(muN)};")
out.append(f"        let mu_l = {arr(muL)};")
out.append(f"        let L = {arr(v['L_limbs'][:4])};")
out.append("        assert(vdf::mu_valid(N, mu_n, N_LIMBS), 'mu_n valid');")
out.append("        assert(vdf::mu_valid(L, mu_l, L_LIMBS), 'mu_l valid');")
mun_bad = muN.copy(); mun_bad[-1] -= 1
mul_bad = muL.copy(); mul_bad[-1] -= 1
out.append(f"        let mu_n_bad = {arr(mun_bad)};")
out.append(f"        let mu_l_bad = {arr(mul_bad)};")
out.append("        assert(!vdf::mu_valid(N, mu_n_bad, N_LIMBS), 'mu_n bad');")
out.append("        assert(!vdf::mu_valid(L, mu_l_bad, L_LIMBS), 'mu_l bad');")
out.append("    }")
out.append("")
# mu_l for the FORGED instance (y' = 1): the forger plays along with the
# in-program challenge derivation, so rejection must come from the Wesolowski
# equation itself, not from mu validation.
from vdf_reference import derive_challenge  # noqa: E402

L_forged = derive_challenge(N_int, int(v["x"], 16), 1, T, kN)
mu_l_forged = barrett_mu(L_forged, kL)

out.append("    // The Bug-A forgery (claim y'=1 with pi=1, formerly accepted via a")
out.append("    // freely chosen challenge L = 2^200, r = 0, pi = 1) must now fail: L is")
out.append("    // derived in-program from (N, x, y', T), and mu_l here is the VALID")
out.append("    // constant for that derived L — rejection comes from the Wesolowski")
out.append("    // equation, not mu validation. The exec-level check is in")
out.append("    // scripts/repro_issues.py.")
out.append("    #[test]")
out.append("    fn forged_output_rejected() {")
out.append(f"        let N = {arr(v['N_limbs'])};")
out.append(f"        let mu_n = {arr(muN)};")
out.append(f"        let x = {arr(v['x_limbs'])};")
out.append(f"        let mu_l = {arr(mu_l_forged)};")
out.append("        let y_fake = array![1, 0, 0, 0, 0, 0, 0, 0].span();")
out.append("        let pi_fake = array![1, 0, 0, 0, 0, 0, 0, 0].span();")
out.append(f"        let T_limbs = {arr(T_limbs)};")
out.append("        assert(!verify_vdf(N, mu_n, x, y_fake, pi_fake, mu_l, T_limbs, 2), 'forgery');")
out.append("    }")
out.append("")
out.append("    // Full Wesolowski check — post-optimization it fits cairo-test's")
out.append("    // 2^32 gas cap (~1.9B gas), so it runs in the default test run.")
out.append("    #[test]")
out.append("    fn full_wesolowski_verifies() {")
out.append(f"        let N = {arr(v['N_limbs'])};")
out.append(f"        let mu_n = {arr(muN)};")
out.append(f"        let x = {arr(v['x_limbs'])};")
out.append(f"        let y = {arr(v['y_limbs'])};")
out.append(f"        let pi = {arr(v['pi_limbs'])};")
out.append(f"        let mu_l = {arr(muL)};")
out.append(f"        let T_limbs = {arr(T_limbs)};")
out.append("        assert(verify_vdf(N, mu_n, x, y, pi, mu_l, T_limbs, 2), 'vdf');")
out.append("    }")
out.append("}")
open(os.path.join(ROOT, "cairo", "lib", "src", "tests.cairo"), "w").write("\n".join(out) + "\n")
print("tests.cairo regenerated (6 passing, none ignored)")
