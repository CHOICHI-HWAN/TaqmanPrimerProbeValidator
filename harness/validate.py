"""Predict fpTm / fpDg of held-out experiment records with the fitted params.json (no refitting) and report errors."""
import sys
from collections import defaultdict

import model as M
from fit_binding import keys_for, load, predict


def main():
    p = M.load_params()
    keys = sorted(p["G"])
    idx = {k: i for i, k in enumerate(keys)}
    ph = [p["H"][k] for k in keys]
    pg = [p["G"][k] for k in keys]
    for path in sys.argv[1:]:
        recs = load([path])
        pred = predict(recs, idx, ph, pg)
        ok = defaultdict(list)
        unseen = 0
        for r, (tm, dg) in zip(recs, pred):
            ks = keys_for(r["f"], r["s"], r["l5"], r["side"], r["r3"])
            if any(k not in idx for k in ks):
                unseen += 1
                continue
            nm = sum(a != b for a, b in zip(r["f"], r["s"]))
            ok[nm].append((abs(tm - r["tm"]), abs(dg - r["dg"])))
        print(f"{path}: {len(recs)} found records, {unseen} use unseen parameters")
        for nm in sorted(ok):
            v = ok[nm]
            print(f"   {nm} mismatches: n={len(v):5}  Tm max {max(x[0] for x in v):.3f}  dG max {max(x[1] for x in v):.3f}")


if __name__ == "__main__":
    main()
