"""Working Python model of mfeprimer's binding thermodynamics, grown step by step from the black-box experiments.

Parameters fitted from experiments live in params.json next to this file (written by the fit_* scripts) so the Go
port can load exactly the same numbers.
"""
import json
import math
from pathlib import Path

import thermo_fit as T

R = 1.987
T37 = 310.15
B = "ACGT"
PARAMS = Path(__file__).with_name("params.json")


def load_params():
    return json.loads(PARAMS.read_text()) if PARAMS.exists() else {}


def save_params(p):
    PARAMS.write_text(json.dumps(p, indent=1, sort_keys=True) + "\n")


def salt_na(mono=50, diva=1.5, dntp=0.25):
    return (mono + 120.0 * math.sqrt(max(diva - dntp, 0.0))) / 1000.0


def perfect_hs(f, mono=50, diva=1.5, dntp=0.25):
    """dH, dS (salt-corrected) of the primer against its exact complement (== `mfeprimer thermo`)."""
    h, s = T.base_hs(f, "2004")
    return h, s + 0.368 * (len(f) - 1) * math.log(salt_na(mono, diva, dntp))


def tm_dg(h, s, oligo=50, self_comp=False):
    ct = oligo * 1e-9
    x = 1.0 if self_comp else 4.0
    return 1000.0 * h / (s + R * math.log(ct / x)) - 273.15, h - T37 * s / 1000.0
