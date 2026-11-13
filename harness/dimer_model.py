"""Superseded first dimer dG model (FINDINGS §7 now uses the exact table model, Go thermo.DimerDG): perfect duplex
of the paired segment (with the symmetry term only for a mirror-symmetric duplex), mismatch stacks shared with the
primer search, and a fitted term per end keyed by (end pair, the two overhanging neighbours) from fit_dimer_ends()."""
import glob
import json
import math
from pathlib import Path

import numpy as np

import dh_parse as P
import model as M
import thermo_fit as T

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
END_FILE = Path(__file__).with_name("dimer_ends.json")


def rc(x):
    return x.translate(str.maketrans("ACGT-", "TGCA-"))[::-1]


def segment(a, br, sh, i, j):
    """A[i..j] against reversed B (B-rev column c-sh faces A column c): f, s (site in A's frame), end keys, symmetry."""
    f = a[i:j + 1]
    s = "".join(COMP[br[x - sh]] for x in range(i, j + 1))
    aL = a[i - 1] if i > 0 else "-"
    bL = br[i - 1 - sh] if 0 <= i - 1 - sh < len(br) else "-"
    aR = a[j + 1] if j + 1 < len(a) else "-"
    bR = br[j + 1 - sh] if 0 <= j + 1 - sh < len(br) else "-"
    sym = f == rc(s) and (aL, bL) == (bR, aR)
    return f, s, f"L:{a[i]}{br[i - sh]}/{aL}{bL}", f"R:{a[j]}{br[j - sh]}/{aR}{bR}", sym


def core_dg(f, s, sym, p):
    na = M.salt_na()
    h = sd = 0.0
    for x in range(len(f) - 1):
        dh, ds = T.NN[f[x:x + 2]]
        h += dh
        sd += ds
    h += 0.2
    sd += -5.7
    for e in (f[0], f[-1]):
        if e in "AT":
            h += 2.2
            sd += 6.9
    if sym:
        sd += -1.4
    sd += 0.368 * (len(f) - 1) * math.log(na)
    g = h - M.T37 * sd / 1000
    for x in range(len(f) - 1):
        if f[x:x + 2] != s[x:x + 2]:
            dh, ds = T.NN[f[x:x + 2]]
            g += p["G"][f[x:x + 2] + "/" + s[x:x + 2]] - (dh - M.T37 * ds / 1000)
    return g


def fit_dimer_ends(files):
    p = M.load_params()
    rows = []
    for fn in files:
        for d in P.dimers(open(fn).read()):
            mk = d["marks"]
            cols = [c for c, ch in enumerate(mk) if ch in ":."]
            f, s, kL, kR, sym = segment(d["a"], d["brev"], d["shift"], cols[0], cols[-1])
            rows.append((core_dg(f, s, sym, p), kL, kR, d["dg"]))
    keys = sorted({r[1] for r in rows} | {r[2] for r in rows})
    idx = {k: i for i, k in enumerate(keys)}
    X = np.zeros((len(rows), len(keys)))
    y = np.zeros(len(rows))
    for n, r in enumerate(rows):
        X[n, idx[r[1]]] += 1
        X[n, idx[r[2]]] += 1
        y[n] = r[3] - r[0]
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = X @ c - y
    END_FILE.write_text(json.dumps({k: round(float(c[i]), 4) for k, i in idx.items()}, indent=1, sort_keys=True))
    return len(rows), len(keys), float(abs(res).max())


def dimer_dg(a, br, sh, i, j, p, ends):
    f, s, kL, kR, sym = segment(a, br, sh, i, j)
    if kL not in ends or kR not in ends:
        return None
    return core_dg(f, s, sym, p) + ends[kL] + ends[kR]


if __name__ == "__main__":
    print(fit_dimer_ends(sorted(glob.glob("../work/dh/db*.txt"))))
