use vdf::{modmul, modpow, limbs_eq, verify_vdf, N_LIMBS, L_LIMBS};

#[cfg(test)]
mod tests {
    use super::*;

    // Gas-cheap checks that exercise both limb domains, plus the full
    // verification (positive + forged-negative) which now fits the gas cap.
    #[test]
    fn modmul_x_squared() {
        let N = array![5135600859310664923, 11520123331373703671, 12696800888177789570, 14289338622341901805, 16749385228720951234, 7503658199550084444, 13926480577912486189, 10906071976544146032].span();
        let mu_n = array![12915993780622577322, 5643813266975672306, 2459985088210256382, 17813654798058237461, 6275531767571678775, 5819627493585523452, 16871388155403298274, 12754440702329786577, 1].span();
        let x = array![3229459620821208167, 11243138961427088649, 7056696807974197263, 12601608781280138922, 9180340365413516663, 10757350401162291613, 4283120750315134194, 45965871558814651].span();
        let got = modmul(x, x, N, mu_n, N_LIMBS);
        let exp = array![2651537663530800705, 1066792907944258969, 15720133730110394003, 12262143101151987429, 42037960474187539, 7790175396831012372, 10253975274584817012, 1208683166737348252].span();
        assert(limbs_eq(got.span(), exp, N_LIMBS), 'x^2 mod N');
    }

    // Challenge derivation (in-program Fiat-Shamir) + r = 2^T mod L.
    #[test]
    fn r_calc_2_to_T_mod_L() {
        let N = array![5135600859310664923, 11520123331373703671, 12696800888177789570, 14289338622341901805, 16749385228720951234, 7503658199550084444, 13926480577912486189, 10906071976544146032].span();
        let x = array![3229459620821208167, 11243138961427088649, 7056696807974197263, 12601608781280138922, 9180340365413516663, 10757350401162291613, 4283120750315134194, 45965871558814651].span();
        let y = array![4082597687400822324, 8131444373239049487, 2149460391222042572, 3462375584892705488, 12943271587032102907, 15240431225920236865, 10708904839994476056, 10465282867214590134].span();
        let mu_l = array![11362295783320014941, 50198858719615909, 16436873375485845336, 16070280500074908154, 63].span();
        let T_limbs = array![65536, 0].span();
        let L = vdf::derive_challenge(N, x, y, T_limbs);
        let exp_L = array![2567734107976243109, 3769667256659455245, 10505807065023532147, 288811737701498671].span();
        assert(limbs_eq(L.span(), exp_L, L_LIMBS), 'L = H(N,x,y,T)');
        let mut two = ArrayTrait::new();
        two.append(2);
        let mut i: usize = 1;
        while i < L_LIMBS {
            two.append(0);
            i += 1;
        }
        let got = modpow(two.span(), T_limbs, 2, L.span(), mu_l, L_LIMBS);
        let exp = array![4119785783969266068, 12671185394234309213, 6080980755253349237, 232501891742019033].span();
        assert(limbs_eq(got.span(), exp, L_LIMBS), '2^T mod L');
    }

    // Regression for the Barrett k-limb truncation bug: with N = 2^512-1,
    // a = 2^448-1, b = 2^448+1 the quotient estimate undershoots by 1 and
    // p - q3*N = b^k + (2^384 - 2) >= b^k; keeping only k limbs wrapped r and
    // returned p mod N - 1. HAC 14.42 keeps k+1 limbs and gets 2^384 - 1.
    #[test]
    fn barrett_boundary_top_limb_all_ones() {
        let f = 18446744073709551615;
        let N = array![f, f, f, f, f, f, f, f].span();
        let mu_n = array![1, 0, 0, 0, 0, 0, 0, 0, 1].span();
        let a = array![f, f, f, f, f, f, f, 0].span();
        let b = array![1, 0, 0, 0, 0, 0, 0, 1].span();
        let got = modmul(a, b, N, mu_n, N_LIMBS);
        let exp = array![f, f, f, f, f, f, 0, 0].span(); // 2^384 - 1
        assert(limbs_eq(got.span(), exp, N_LIMBS), 'barrett boundary');
    }

    // The Barrett constants are validated in-program: a tampered mu
    // (here: top limb decremented) must be rejected in both domains.
    #[test]
    fn mu_validation_accepts_and_rejects() {
        let N = array![5135600859310664923, 11520123331373703671, 12696800888177789570, 14289338622341901805, 16749385228720951234, 7503658199550084444, 13926480577912486189, 10906071976544146032].span();
        let mu_n = array![12915993780622577322, 5643813266975672306, 2459985088210256382, 17813654798058237461, 6275531767571678775, 5819627493585523452, 16871388155403298274, 12754440702329786577, 1].span();
        let mu_l = array![11362295783320014941, 50198858719615909, 16436873375485845336, 16070280500074908154, 63].span();
        let L = array![2567734107976243109, 3769667256659455245, 10505807065023532147, 288811737701498671].span();
        assert(vdf::mu_valid(N, mu_n, N_LIMBS), 'mu_n valid');
        assert(vdf::mu_valid(L, mu_l, L_LIMBS), 'mu_l valid');
        let mu_n_bad = array![12915993780622577322, 5643813266975672306, 2459985088210256382, 17813654798058237461, 6275531767571678775, 5819627493585523452, 16871388155403298274, 12754440702329786577, 0].span();
        let mu_l_bad = array![11362295783320014941, 50198858719615909, 16436873375485845336, 16070280500074908154, 62].span();
        assert(!vdf::mu_valid(N, mu_n_bad, N_LIMBS), 'mu_n bad');
        assert(!vdf::mu_valid(L, mu_l_bad, L_LIMBS), 'mu_l bad');
    }

    // The Bug-A forgery (claim y'=1 with pi=1, formerly accepted via a
    // freely chosen challenge L = 2^200, r = 0, pi = 1) must now fail: L is
    // derived in-program from (N, x, y', T), and mu_l here is the VALID
    // constant for that derived L — rejection comes from the Wesolowski
    // equation, not mu validation. The exec-level check is in
    // scripts/repro_issues.py.
    #[test]
    fn forged_output_rejected() {
        let N = array![5135600859310664923, 11520123331373703671, 12696800888177789570, 14289338622341901805, 16749385228720951234, 7503658199550084444, 13926480577912486189, 10906071976544146032].span();
        let mu_n = array![12915993780622577322, 5643813266975672306, 2459985088210256382, 17813654798058237461, 6275531767571678775, 5819627493585523452, 16871388155403298274, 12754440702329786577, 1].span();
        let x = array![3229459620821208167, 11243138961427088649, 7056696807974197263, 12601608781280138922, 9180340365413516663, 10757350401162291613, 4283120750315134194, 45965871558814651].span();
        let mu_l = array![14337335109646308387, 5504977480557402834, 13963848616343359628, 1805625020813257192, 43].span();
        let y_fake = array![1, 0, 0, 0, 0, 0, 0, 0].span();
        let pi_fake = array![1, 0, 0, 0, 0, 0, 0, 0].span();
        let T_limbs = array![65536, 0].span();
        assert(!verify_vdf(N, mu_n, x, y_fake, pi_fake, mu_l, T_limbs, 2), 'forgery');
    }

    // Full Wesolowski check — post-optimization it fits cairo-test's
    // 2^32 gas cap (~1.9B gas), so it runs in the default test run.
    #[test]
    fn full_wesolowski_verifies() {
        let N = array![5135600859310664923, 11520123331373703671, 12696800888177789570, 14289338622341901805, 16749385228720951234, 7503658199550084444, 13926480577912486189, 10906071976544146032].span();
        let mu_n = array![12915993780622577322, 5643813266975672306, 2459985088210256382, 17813654798058237461, 6275531767571678775, 5819627493585523452, 16871388155403298274, 12754440702329786577, 1].span();
        let x = array![3229459620821208167, 11243138961427088649, 7056696807974197263, 12601608781280138922, 9180340365413516663, 10757350401162291613, 4283120750315134194, 45965871558814651].span();
        let y = array![4082597687400822324, 8131444373239049487, 2149460391222042572, 3462375584892705488, 12943271587032102907, 15240431225920236865, 10708904839994476056, 10465282867214590134].span();
        let pi = array![2963693354413529202, 13731823323608936747, 16862469021914409217, 13927750145564672412, 3042382273525009142, 16591053712920500940, 17270321545452892979, 263572611769821423].span();
        let mu_l = array![11362295783320014941, 50198858719615909, 16436873375485845336, 16070280500074908154, 63].span();
        let T_limbs = array![65536, 0].span();
        assert(verify_vdf(N, mu_n, x, y, pi, mu_l, T_limbs, 2), 'vdf');
    }
}
