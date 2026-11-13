"""Discovery step for the hairpin energy (kept for the record; hp_energy.py is the final, table-only model):
fit the loop-size and loop-end terms from `mfeprimer hairpin` JSON runs (synth_hairpin.py, --mode dg).

Model (FINDINGS §8): stem NN stacks with the salt term on every stack; one end term per (3' closing base Y, the loop
base 5' of it) carrying dH and dS (enters Tm and dG); one loop-size term carrying dG only (dG, not Tm).
Tm = 1000*dH/(dS + R ln(Ct/4)) - 273.15.  Writes hairpin_params.json.
"""
import os
import glob
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

import model as M
import thermo_fit as T

OUT = Path(__file__).with_name("hairpin_params.json")
P3 = Path(os.environ.get("PRIMER3_CONFIG", "primer3_config"))
SPECIAL = {l.split()[0] for f in ("triloop.dh", "tetraloop.dh") for l in (P3 / f).read_text().splitlines() if l.strip()}
LN = 1.987 * math.log(50e-9 / 4)


def dec(st):
    op = [i for i, c in enumerate(st) if c == "/"]
    cl = [i for i, c in enumerate(st) if c == "\\"]
    return list(zip(op, reversed(cl)))


def stem_hs(stem):
    h = sum(T.NN[stem[k:k + 2]][0] for k in range(len(stem) - 1))
    s = sum(T.NN[stem[k:k + 2]][1] for k in range(len(stem) - 1)) + 0.368 * (len(stem) - 1) * math.log(M.salt_na())
    return h, s


def perfect_rows(paths):
    rows = []
    for p in paths:
        for r in json.load(open(p)):
            if "score" not in r:
                continue
            s, st = r["seq"], r["struct"]
            pr = dec(st)
            n = len(s)
            P = len(pr)
            if pr[0] != (0, n - 1) or pr[-1] != (P - 1, n - P):
                continue
            if s[P - 1:n - P + 1] in SPECIAL:
                continue
            rows.append({"s": s, "P": P, "L": n - 2 * P, "tm": r["tm"], "dg": r["dg"]})
    return rows


def main():
    rows = perfect_rows(sys.argv[1:] or sorted(glob.glob("../work/hpx/dg*/rows.json")))
    ek = sorted({r["s"][len(r["s"]) - r["P"]] + r["s"][len(r["s"]) - r["P"] - 1] for r in rows})
    lk = sorted({r["L"] for r in rows})
    ei = {k: i for i, k in enumerate(ek)}
    li = {k: i for i, k in enumerate(lk)}
    H0, S0, E, L, tm, dg = [], [], [], [], [], []
    for r in rows:
        s, P = r["s"], r["P"]
        h, sd = stem_hs(s[:P])
        H0.append(h)
        S0.append(sd)
        E.append(ei[s[len(s) - P] + s[len(s) - P - 1]])
        L.append(li[r["L"]])
        tm.append(r["tm"])
        dg.append(r["dg"])
    H0, S0, E, L, tm, dg = map(np.array, (H0, S0, E, L, tm, dg))
    G0 = H0 - M.T37 * S0 / 1000
    ne, nl = len(ek), len(lk)
    # dG: E_G (per key) + L_G (per size); anchor L=3 at the SantaLucia value 3.5
    X = np.zeros((len(rows), ne + nl))
    X[np.arange(len(rows)), E] = 1
    X[np.arange(len(rows)), ne + L] = 1
    anchor = np.zeros(ne + nl)
    anchor[ne + li[3]] = 1
    A = np.vstack([X, anchor * 1e3])
    y = np.concatenate([dg - G0, [3.5 * 1e3]])
    cg, *_ = np.linalg.lstsq(A, y, rcond=None)
    gres = X @ cg - (dg - G0)
    eg = cg[:ne]

    def res(hE):
        h = H0 + hE[E]
        sd = S0 + (hE[E] - eg[E]) * 1000 / M.T37
        return 1000 * h / (sd + LN) - 273.15 - tm
    hE = least_squares(res, np.zeros(ne)).x
    tres = res(hE)
    print(f"{len(rows)} hairpins  dG max {abs(gres).max():.5f}  Tm max {abs(tres).max():.5f}")
    out = {"E": {k: {"H": round(float(hE[i]), 4), "S": round(float((hE[i] - eg[i]) * 1000 / M.T37), 4)}
                 for k, i in ei.items()},
           "L": {str(k): round(float(cg[ne + i]), 4) for k, i in li.items()}}
    OUT.write_text(json.dumps(out, indent=1, sort_keys=True) + "\n")
    print(json.dumps(out["L"]))
    for k, v in sorted(out["E"].items()):
        print(k, v, "G", round(v["H"] - M.T37 * v["S"] / 1000, 4))


if __name__ == "__main__":
    main()
