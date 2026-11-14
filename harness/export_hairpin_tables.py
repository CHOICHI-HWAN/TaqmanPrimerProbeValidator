"""Write internal/thermo/hairpin_tables.json from hp_energy (the tables the Go engine embeds)."""
import json
import math
from pathlib import Path

import hp_energy as E

out = {
    "mm": {k: [v[0], v[1]] for k, v in sorted(E.MM.items())},
    "dangle": {f"{half}:{k}": [v[0], v[1]] for (half, k), v in sorted(E.DANGLE.items())
               if not math.isinf(v[0]) and not math.isinf(v[1])},
    "bonus": dict(sorted(E.BONUS.items())),
    "loop_sizes": E.LOOP_SIZES,
    "loop_g": E.LOOP_G,
}
p = Path(__file__).resolve().parent.parent / "internal/thermo/hairpin_tables.json"
p.write_text(json.dumps(out, indent=0) + "\n")
print(p, len(out["mm"]), len(out["dangle"]), len(out["bonus"]))
