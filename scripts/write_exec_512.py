#!/usr/bin/env python3
"""Regenerate 512-bit exec args from the 512 vector file (repo-local paths).

New input layout (challenge L and r are derived in-program):
  [N(8), mu_n(9), x(8), y(8), pi(8), mu_l(5), T(2)]  = 48 felts
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from gen_cairo_tests import barrett_mu  # noqa: E402

v = json.load(open(os.path.join(ROOT, "vectors", "vdf_vectors_512.json")))[0]
N_int = int(v["N"], 16)
L_int = int(v["L"], 16)
muN = barrett_mu(N_int, v["n_limbs"])
muL = barrett_mu(L_int, 4)
T = v["T"]
T_limbs = [T & ((1 << 64) - 1), (T >> 64) & ((1 << 64) - 1)]
flat = v["N_limbs"] + muN + v["x_limbs"] + v["y_limbs"] + v["pi_limbs"] + muL + T_limbs
args = ["0x%x" % len(flat)] + ["0x%x" % x for x in flat]
out = os.path.join(ROOT, "vectors", "exec_args_512.json")
json.dump(args, open(out, "w"))
print("512 exec args:", len(flat), "felts ->", out)
