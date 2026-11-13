"""Hairpin dG/Tm of a given structure, reconstructed from `mfeprimer hairpin` (FINDINGS §8).

Terms (dH kcal/mol, dS cal/K/mol), all from SantaLucia & Hicks 2004 tables (perfect stacks as in the binding model,
mismatch stacks / dangles / tri- and tetraloop bonuses as distributed with primer3):
  init 0.2 / -5.7; one stack per pair of adjacent paired columns; at the loop end one more stack read as
  5'-X X-3' / 3'-Y b-5' (X·Y the closing pair, b the loop base 5' of Y; the loop base after X is not used);
  outer end: mismatch stack l X / r Y when both flanks exist, a dangling end when one does, nothing otherwise;
  salt 0.368 ln[Na+] per stack (half per dangle); loop size and tri/tetraloop bonus enter dG only.
Tm = 1000 dH / (dS + R ln(Ct/4)) - 273.15.
"""
import os
import math
from pathlib import Path

import model as M
import thermo_fit as T

P3 = Path(os.environ.get("PRIMER3_CONFIG", "primer3_config"))
B = "ACGT"
COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}
LOOP_SIZES = [3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30]
LOOP_G = [3.5, 3.5, 3.3, 4.0, 4.2, 4.3, 4.5, 4.6, 5.0, 5.1, 5.3, 5.5, 5.7, 6.1, 6.3]
LN = 1.9872 * math.log(50e-9 / 4)


def _tab(name):
    return [float(x) if x != "inf" else math.inf for x in (P3 / name).read_text().split()]


def _4d(h, s):
    th, ts = _tab(h), _tab(s)
    return {B[a] + B[b] + "/" + B[c] + B[d]: (th[i] / 1000, ts[i])
            for a in range(4) for b in range(4) for c in range(4) for d in range(4)
            for i in [((a * 4 + b) * 4 + c) * 4 + d]}


MM = _4d("stackmm.dh", "stackmm.ds")          # 5'ab3'/3'cd5'
# where mfeprimer's table differs from primer3's file (measured): a missing dH, and a few tandem mismatches;
# every other tandem mismatch (two mismatched columns in a row) contributes nothing
MM.update({"TA/AA": (4.7, 12.9), "TT/TG": (5.8, 16.3), "GT/TT": (5.8, 16.3), "TA/TG": (-0.1, -1.7),
           "GT/AT": (-0.1, -1.7), "CT/CA": (0.0, -4.4), "AC/TC": (0.0, -4.4)})
MM = {k: ((0.0, 0.0) if math.isinf(v[0]) or math.isinf(v[1]) else v) for k, v in MM.items()}
_dh, _ds = _tab("dangle.dh"), _tab("dangle.ds")
DANGLE = {(half, B[a] + B[b] + B[c]): (_dh[half * 64 + (a * 4 + b) * 4 + c] / 1000, _ds[half * 64 + (a * 4 + b) * 4 + c])
          for half in (0, 1) for a in range(4) for b in range(4) for c in range(4)}


def _bonus():
    out = {}
    for f in ("triloop", "tetraloop"):
        h = {l.split()[0]: float(l.split()[1]) for l in (P3 / f"{f}.dh").read_text().splitlines() if l.strip()}
        s = {l.split()[0]: float(l.split()[1]) for l in (P3 / f"{f}.ds").read_text().splitlines() if l.strip()}
        for k in h:
            out[k] = round(h[k] / 1000 - M.T37 * s[k] / 1e6, 1)
    return out


BONUS = _bonus()
BONUS["GTTTTC"] = 0.5          # primer3's file has -500 cal; mfeprimer uses +0.5 kcal/mol


def loop_g(L):
    if L <= 10:
        return LOOP_G[LOOP_SIZES.index(L)] if L in LOOP_SIZES else 0.0
    for size, g in zip(LOOP_SIZES, LOOP_G):
        if size > L:
            return g
    return 0.0


def stack(top, bot):
    """5'-top-3' / 3'-bot-5'."""
    if COMP[top[0]] == bot[0] and COMP[top[1]] == bot[1]:
        return T.NN[top]
    return MM[top + "/" + bot]


def energy(s, pairs):
    """pairs: outermost first, as decoded from the structure. Returns (dG, Tm)."""
    n = len(s)
    lnna = math.log(M.salt_na())
    h, sd, units = 0.2, -5.7, 0.0
    for (i, j), (i2, j2) in zip(pairs, pairs[1:]):
        for c in range(i2 - i):          # mismatched columns between two pairs are read as mismatch stacks
            dh, ds = stack(s[i + c] + s[i + c + 1], s[j - c] + s[j - c - 1])
            h, sd, units = h + dh, sd + ds, units + 1
    i, j = pairs[0]
    l = s[i - 1] if i > 0 else ""
    r = s[j + 1] if j + 1 < n else ""
    if l and r:
        dh, ds = stack(l + s[i], r + s[j])
        h, sd, units = h + dh, sd + ds, units + 1
    elif r:
        dh, ds = DANGLE[(0, s[j] + s[i] + r)]
        h, sd, units = h + dh, sd + ds, units + 0.5
    elif l:
        dh, ds = DANGLE[(1, s[j] + s[i] + l)]
        h, sd, units = h + dh, sd + ds, units + 0.5
    i, j = pairs[-1]
    dh, ds = stack(s[i] + s[i], s[j] + s[j - 1])
    h, sd, units = h + dh, sd + ds, units + 1
    sd += 0.368 * units * lnna
    g = h - M.T37 * sd / 1000 + loop_g(j - i - 1) + BONUS.get(s[i:j + 1], 0.0)
    tm = 1000 * h / (sd + LN) - 273.15
    return g, tm
