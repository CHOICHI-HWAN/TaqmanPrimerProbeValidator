"""Fit only the parameter keys params.json does not have yet, holding every existing key fixed.

usage: fit_partial.py <joined.tsv> [...]   -> adds the new keys to params.json, prints residuals per file
Much faster than a full refit and cannot disturb parameters that are already validated.
"""
import sys

import numpy as np
from scipy.optimize import least_squares

import model as M
from fit_binding import LN_CT4, baseline, keys_for, load


def main():
    p = M.load_params()
    recs = load(sys.argv[1:])
    new = {}
    rows = []
    for r in recs:
        ks = keys_for(r["f"], r["s"], r["l5"], r["side"], r["r3"])
        h0, g0 = baseline(r["f"], r["s"], r["l5"])
        unk = {}
        for k, c in ks.items():
            if k in p["G"]:
                h0 += c * p["H"][k]
                g0 += c * p["G"][k]
            else:
                unk[k] = c
                new.setdefault(k, len(new))
        rows.append((h0, g0, unk, r))
    use = [x for x in rows if x[2]]
    print(f"records {len(recs)}  with new keys {len(use)}  new keys {len(new)}")
    if not new:
        return
    X = np.zeros((len(use), len(new)))
    for i, (_, _, unk, _) in enumerate(use):
        for k, c in unk.items():
            X[i, new[k]] = c
    H0 = np.array([x[0] for x in use])
    G0 = np.array([x[1] for x in use])
    tm = np.array([x[3]["tm"] for x in use])
    dg = np.array([x[3]["dg"] for x in use])
    pg, *_ = np.linalg.lstsq(X, dg - G0, rcond=None)

    def res(ph):
        h = H0 + X @ ph
        g = G0 + X @ pg
        s = (h - g) * 1000.0 / M.T37
        return 1000.0 * h / (s + LN_CT4) - 273.15 - tm

    ph = least_squares(res, np.zeros(len(new))).x
    print(f"  new-key records: Tm max {np.abs(res(ph)).max():.3f}  dG max {np.abs(G0 + X @ pg - dg).max():.3f}  "
          f"rank {np.linalg.matrix_rank(X)}/{len(new)}")
    for k, j in new.items():
        p["H"][k] = round(float(ph[j]), 4)
        p["G"][k] = round(float(pg[j]), 4)
    M.save_params(p)


if __name__ == "__main__":
    main()
