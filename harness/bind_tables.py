"""Binding (primer on template) energy as `mfeprimer` computes it for the primer search — from the SantaLucia
tables, no fitted values (FINDINGS §3). Verified on 135,851 full-precision bindings (`mfeprimer -j`).

f primer 5'->3'; s site in primer sense (the bases a perfect primer would have); l5 / r3 the template-side bases
beyond the primer's 5' / 3' end, in primer sense ('-' at the genome end).
"""
import math

import hp_energy as E
import model as M

C = {"A": "T", "T": "A", "G": "C", "C": "G"}
AMB = str.maketrans("RYKMSWBDHV", "NNNNNNNNNN")   # every ambiguity code is scored as N
R_GAS = 1.9872


def stack(top, bot):
    if "N" in top or "N" in bot:               # any N: no dH/dS, but the stack still counts for salt
        return 0.0, 0.0
    return E.stack(top, bot)


def cut_point(f, s):
    """Columns dropped from the 5' end: a mismatch run (N counts) touching the 5' end, or the 3'-most run of three
    or more; the duplex then starts at the last column of that run."""
    mm = [a != b for a, b in zip(f, s)]
    n, d, i = len(f), 0, 0
    while i < n:
        if mm[i]:
            j = i
            while j < n and mm[j]:
                j += 1
            if j - i >= 3 or i == 0:
                d = j
            i = j
        else:
            i += 1
    return d


def energy(f, s, l5, r3, side):
    """dG37, Tm (and dH, dS) of the binding."""
    s, l5, r3 = s.upper().translate(AMB), l5.upper().translate(AMB), r3.upper().translate(AMB)
    d = cut_point(f, s)
    if d:
        f, s = f[d - 1:], s[d - 1:]
    br = "".join(C.get(x, "N") for x in s)      # template written 3'->5' under the primer
    if side == "R":
        # R: read on the plus strand as rc(primer) over the plus-strand site, so each mismatch is scored as the
        # complementary mismatch (A.G -> T.C); perfect stacks are unchanged
        top, bot = "".join(C[x] for x in reversed(f)), s[::-1]
    else:
        top, bot = f, br
    h, sd, units = 0.2, -5.7, 0.0
    for x in range(len(f) - 1):
        dh, ds = stack(top[x:x + 2], bot[x:x + 2])
        h, sd, units = h + dh, sd + ds, units + 1
    if d:
        pass                                    # cut 5' end: blunt, no terminal-AT term
    elif l5 in C:                               # dangling template base (F reads it uncomplemented)
        dh, ds = E.DANGLE[(0, C[f[0]] + f[0] + l5) if side == "F" else (1, C[f[0]] + f[0] + C[l5])]
        h, sd, units = h + dh, sd + ds, units + 0.5
    elif l5 == "N":
        units += 0.5
    elif f[0] in "AT":                          # genome end: terminal AT penalty
        h, sd = h + 2.2, sd + 6.9
    if side == "F":                             # F: the 3' dangling base is read from the primer itself
        dh, ds = E.DANGLE[(1, f[-1] + C[f[-1]] + f[-1])]
        h, sd, units = h + dh, sd + ds, units + 0.5
    elif r3 in C:
        dh, ds = E.DANGLE[(0, f[-1] + C[f[-1]] + C[r3])]
        h, sd, units = h + dh, sd + ds, units + 0.5
    elif r3 == "N":
        units += 0.5
    elif f[-1] in "AT":
        h, sd = h + 2.2, sd + 6.9
    sd += 0.368 * units * math.log(M.salt_na())
    return h - M.T37 * sd / 1000, 1000 * h / (sd + R_GAS * math.log(50e-9 / 4)) - 273.15, h, sd
