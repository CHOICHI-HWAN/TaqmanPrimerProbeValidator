"""Fit candidate nearest-neighbour models against `mfeprimer thermo` output (black-box reference).

usage: python3 thermo_fit.py <file.thermo.tsv> --mono 50 --diva 1.5 --dntp 0.25 --oligo 50
Prints the max / mean absolute error of Tm and dG for every model variant so the closest one can be picked.
"""
import argparse
import csv
import math

# SantaLucia & Hicks (2004) nearest-neighbour parameters (AA/TT differs from the 1998 table): dH kcal/mol, dS cal/(K mol)
NN = {"AA": (-7.6, -21.3), "TT": (-7.6, -21.3), "AT": (-7.2, -20.4), "TA": (-7.2, -21.3),
      "CA": (-8.5, -22.7), "TG": (-8.5, -22.7), "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
      "CT": (-7.8, -21.0), "AG": (-7.8, -21.0), "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
      "CG": (-10.6, -27.2), "GC": (-9.8, -24.4), "GG": (-8.0, -19.9), "CC": (-8.0, -19.9)}
COMP = str.maketrans("ACGT", "TGCA")
R = 1.987


def rc(s):
    return s.translate(COMP)[::-1]


def base_hs(s, init):
    h = s_ = 0.0
    for i in range(len(s) - 1):
        dh, ds = NN[s[i:i + 2]]
        h += dh
        s_ += ds
    if init == "1998":
        for end in (s[0], s[-1]):
            if end in "GC":
                h += 0.1; s_ += -2.8
            else:
                h += 2.3; s_ += 4.1
    elif init == "2004":
        h += 0.2; s_ += -5.7
        for end in (s[0], s[-1]):
            if end in "AT":
                h += 2.2; s_ += 6.9
    if s == rc(s):
        s_ += -1.4
    return h, s_


def na_eq(mono, diva, dntp, mode):
    if mode == "vonahsen":
        free = max(diva - dntp, 0.0)
        return mono + 120.0 * math.sqrt(free)
    if mode == "mono":
        return mono
    raise ValueError(mode)


def model(s, mono, diva, dntp, oligo, init, salt_mode, salt_on):
    h, ds = base_hs(s, init)
    na = na_eq(mono, diva, dntp, salt_mode) / 1000.0
    n = len(s)
    if salt_on == "dS":
        ds += 0.368 * (n - 1) * math.log(na)
    ct = oligo * 1e-9
    x = 1.0 if s == rc(s) else 4.0
    tm = 1000.0 * h / (ds + R * math.log(ct / x)) - 273.15
    if salt_on == "Tm":
        tm += 16.6 * math.log10(na)
    dg = h - 310.15 * ds / 1000.0
    return tm, dg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tsv")
    ap.add_argument("--mono", type=float, default=50)
    ap.add_argument("--diva", type=float, default=1.5)
    ap.add_argument("--dntp", type=float, default=0.25)
    ap.add_argument("--oligo", type=float, default=50)
    a = ap.parse_args()
    rows = [r for r in csv.reader(open(a.tsv), delimiter="\t") if r and not r[0].startswith("#")]
    for init in ("1998", "2004"):
        for salt_mode in ("vonahsen", "mono"):
            for salt_on in ("dS", "Tm"):
                etm, edg = [], []
                for r in rows:
                    seq, tm, dg = r[1], float(r[4]), float(r[5])
                    mtm, mdg = model(seq, a.mono, a.diva, a.dntp, a.oligo, init, salt_mode, salt_on)
                    etm.append(abs(mtm - tm)); edg.append(abs(mdg - dg))
                print(f"{init:5} {salt_mode:9} salt@{salt_on:3}  Tm max {max(etm):6.2f} mean {sum(etm)/len(etm):6.3f} | "
                      f"dG max {max(edg):6.2f} mean {sum(edg)/len(edg):6.3f}")


if __name__ == "__main__":
    main()
