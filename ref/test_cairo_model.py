"""Differential checks against Python integers, including noncanonical bases."""

import random
import unittest

import cairo_model_v2 as m


class ArithmeticTests(unittest.TestCase):
    def test_full_and_truncated_products(self):
        rng = random.Random(34767)
        for k in (4, 8, 16, 32):
            for an, bn in ((k, k), (k + 1, k + 1), (k + 1, k)):
                for case in range(200):
                    a = [m.MASK] * an if case == 0 else [rng.getrandbits(64) for _ in range(an)]
                    b = [m.MASK] * bn if case == 0 else [rng.getrandbits(64) for _ in range(bn)]
                    product = m.from_limbs(a) * m.from_limbs(b)
                    for length in (an + bn, bn + 1, rng.randrange(1, an + bn + 1)):
                        got = m.bigint_mul(a, b, an, bn, length)
                        self.assertEqual(m.from_limbs(got), product % (1 << (64 * length)))

    def test_modular_products(self):
        rng = random.Random(20260907)
        for k in (4, 8, 16, 32):
            bound = 1 << (64 * k)
            for n in (bound - 1, bound // m.TWO64 + 1, rng.randrange(bound // 2, bound)):
                modulus, mu = m.to_limbs(n, k), m.barrett_mu(n, k)
                for a, b in ((0, n - 1), (1, n - 1), (n - 1, n - 1),
                             (bound - 1, bound - 1), (rng.randrange(bound), rng.randrange(bound))):
                    got = m.modmul(m.to_limbs(a, k), m.to_limbs(b, k), modulus, mu, k)
                    self.assertEqual(m.from_limbs(got), a * b % n)

    def test_single_exponentiation(self):
        n = (1 << 255) - 19
        modulus, mu = m.to_limbs(n, 4), m.barrett_mu(n, 4)
        for a in (0, 1, n + 1):
            for e in (0, 1, 1 << 64, (1 << 128) - 1):
                got = m.modpow(m.to_limbs(a, 4), m.to_limbs(e, 2), 2, modulus, mu, 4)
                self.assertEqual(m.from_limbs(got), pow(a, e, n))

    def test_joint_exponentiation(self):
        rng = random.Random(20260907)
        for k in (4, 8, 16, 32):
            n = (1 << (64 * k - 1)) + 1
            modulus, mu = m.to_limbs(n, k), m.barrett_mu(n, k)
            for e, f in ((0, 0), (0, 1), (1, 0), (1, 1), (1 << 192, 1 << 128),
                         ((1 << 256) - 1, (1 << 256) - 1),
                         (rng.getrandbits(256), rng.getrandbits(256))):
                a, b = rng.randrange(1 << (64 * k)), rng.randrange(1 << (64 * k))
                got = m.joint_modpow(m.to_limbs(a, k), m.to_limbs(e, 4),
                                     m.to_limbs(b, k), m.to_limbs(f, 4), modulus, mu, k)
                self.assertEqual(m.from_limbs(got), pow(a, e, n) * pow(b, f, n) % n)


if __name__ == "__main__":
    unittest.main()
