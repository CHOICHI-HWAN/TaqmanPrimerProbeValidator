"""Check hp_energy.energy against mfeprimer's dG/Tm for the structures it reported (synth_hairpin.py rows.json)."""
import json
import sys

import numpy as np

import fit_hairpin as F
import hp_energy as E


def main():
    out = []
    for fn in sys.argv[1:]:
        for r in json.load(open(fn)):
            if "score" not in r:
                continue
            g, tm = E.energy(r["seq"], F.dec(r["struct"]))
            out.append((r["dg"] - g, r["tm"] - tm, r["struct"], r["seq"]))
    dg = np.array([o[0] for o in out])
    dt = np.array([o[1] for o in out])
    print(f"{len(out)} hairpins  dG res [{dg.min():.4f}, {dg.max():.4f}]  Tm res [{dt.min():.4f}, {dt.max():.4f}]  "
          f"|dG|>0.005: {(abs(dg) > 0.005).sum()}  |Tm|>0.01: {(abs(dt) > 0.01).sum()}")
    for o in sorted(out, key=lambda o: -abs(o[0]) - abs(o[1]) / 10)[:6]:
        print(f"  {o[0]:+.4f} {o[1]:+.4f} {o[2]} {o[3]}")


if __name__ == "__main__":
    main()
