"""Fit the binding-energy parameters (mismatch stacks, 5' flank / terminal, 3' end) to mfeprimer's reported fpTm / fpDg.

Model for a primer f (the oligo mfeprimer reports, i.e. the 3'-window variant when one was used) on a site s of the same
length (genome plus strand, so a perfect site has s == f) with l5 the genome base just 5' of the site:

  dH = dH_perfect(f) + sum over stacks i=0..n-2 that contain a mismatch:  P_H[stack] - NN_H[f[i:i+2]]
                     + P_H[five-prime end key] + P_H[three-prime end key]
  dG likewise with P_G (37 C, salt-corrected exactly like the perfect duplex: n-1 stacks)
  Tm from dH and dS = (dH - dG) * 1000 / 310.15

Keys:  stack "f1f2/s1s2" (internal),  "T5:f0f1/s0s1" (stack at the 5' end),  "F5:f0/s0/l5",  "E3:f3".
usage: fit_binding.py <joined.tsv> [<joined.tsv> ...]   -> writes harness/params.json, prints residuals per file
"""
import csv
import math
import sys
from collections import defaultdict

import numpy as np
from scipy.optimize import least_squares

import model as M
import thermo_fit as T

LN_CT4 = M.R * math.log(50e-9 / 4)
RC_STACKS = True
AMB = str.maketrans("RYKMSWBDHV", "NNNNNNNNNN")
COMP = str.maketrans("ACGTRYKMSWBDHVN", "TGCAYRMKSWVHDBN")


def rc(x):
    return x.translate(COMP)[::-1]


def nn_g(d):
    h, s = T.NN[d]
    return h - M.T37 * s / 1000.0


def effective(f, s, l5):
    """The part of the binding that counts: scanning from the 3' end, a run of >= 3 consecutive mismatches ends the
    duplex; the run's 3'-most mismatch is kept as the terminal (5') pair and the genome base before it becomes l5.
    (Known exception, 4 of ~46,000 real bindings: a 3-run followed one base later by a 4-run keeps two — FINDINGS §4e.)"""
    n = len(f)
    run = 0
    for i in range(n - 1, -1, -1):
        if f[i] != s[i]:
            run += 1
            if run == 3:
                j = i + 2
                return f[j:], s[j:], s[j - 1]
        else:
            run = 0
    return f, s, l5


def keys_for(f, s, l5, side="F", r3=""):
    """{param key: count} for one binding. The R primer's ends are scored differently (its 3' flank r3 counts too)."""
    f, s, l5 = effective(f, s, l5)
    k = defaultdict(int)
    n = len(f)
    for i in range(n - 1):
        if f[i:i + 2] != s[i:i + 2]:
            pf, ps = f[i:i + 2], s[i:i + 2]
            if side == "R" and RC_STACKS:          # read the stack from the other strand
                pf, ps = rc(pf), rc(ps)
            k[("T5:" if i == 0 else "") + pf + "/" + ps] += 1
    if side == "F":
        k[f"F5:{f[0]}/{s[0]}/{l5}"] += 1
        k[f"E3:{f[-1]}"] += 1
    else:
        k[f"RF5:{f[0]}/{s[0]}/{l5}"] += 1
        k[f"RE3:{f[-1]}/{r3}"] += 1
    return k


def baseline(f, s, l5="A"):
    """dH, dG of the perfect duplex of f minus the perfect stacks that a mismatch replaces."""
    f, s, l5 = effective(f, s, l5)
    h, sd = M.perfect_hs(f)
    g = h - M.T37 * sd / 1000.0
    for i in range(len(f) - 1):
        if f[i:i + 2] != s[i:i + 2]:
            h -= T.NN[f[i:i + 2]][0]
            g -= nn_g(f[i:i + 2])
    return h, g


def load(paths):
    recs = []
    for p in paths:
        for r in csv.DictReader(open(p), delimiter="\t"):
            if r["found"] != "1" or not r["fpTm"]:
                continue
            f = r["fpSeq"]
            if len(f) != len(r["site"]):
                continue
            # every genome ambiguity code (R, Y, ..., N) is scored like N (verified: harness/docs FINDINGS §4c)
            recs.append({"src": p, "f": f, "s": r["site"].translate(AMB), "l5": r["l5"].translate(AMB), "side": r.get("side", "F"),
                         "r3": r.get("r3", "").translate(AMB),
                         "tm": float(r["fpTm"]), "dg": float(r["fpDg"])})
    return recs


def predict(recs, idx, ph, pg):
    out = []
    for r in recs:
        h0, g0 = baseline(r["f"], r["s"], r["l5"])
        h, g = h0, g0
        for k, c in keys_for(r["f"], r["s"], r["l5"], r["side"], r["r3"]).items():
            j = idx.get(k)
            if j is not None:
                h += c * ph[j]
                g += c * pg[j]
        sdel = (h - g) * 1000.0 / M.T37
        out.append((1000.0 * h / (sdel + LN_CT4) - 273.15, g))
    return out


def main():
    recs = load(sys.argv[1:])
    idx = {}
    rows = []
    for r in recs:
        ks = keys_for(r["f"], r["s"], r["l5"], r["side"], r["r3"])
        for k in ks:
            idx.setdefault(k, len(idx))
        rows.append(ks)
    n = len(idx)
    X = np.zeros((len(recs), n))
    H0 = np.zeros(len(recs))
    G0 = np.zeros(len(recs))
    for i, (r, ks) in enumerate(zip(recs, rows)):
        for k, c in ks.items():
            X[i, idx[k]] = c
        H0[i], G0[i] = baseline(r["f"], r["s"], r["l5"])
    tm = np.array([r["tm"] for r in recs])
    dg = np.array([r["dg"] for r in recs])
    pg, *_ = np.linalg.lstsq(X, dg - G0, rcond=None)          # dG is linear: exact least squares

    def res(ph):
        h = H0 + X @ ph
        g = G0 + X @ pg
        s = (h - g) * 1000.0 / M.T37
        return 1000.0 * h / (s + LN_CT4) - 273.15 - tm

    ph = least_squares(res, np.zeros(n), method="lm" if len(recs) > n else "trf").x
    etm = np.abs(res(ph))
    edg = np.abs(G0 + X @ pg - dg)
    print(f"records {len(recs)}  params {n}  rank {np.linalg.matrix_rank(X)}")
    for src in sorted({r['src'] for r in recs}):
        m = np.array([r["src"] == src for r in recs])
        print(f"  {src}: n={m.sum():5}  Tm max {etm[m].max():.3f}  dG max {edg[m].max():.3f}")
    inv = {j: k for k, j in idx.items()}
    M.save_params({"H": {inv[j]: round(float(ph[j]), 4) for j in range(n)},
                   "G": {inv[j]: round(float(pg[j]), 4) for j in range(n)}})


if __name__ == "__main__":
    main()
