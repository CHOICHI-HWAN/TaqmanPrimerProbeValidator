#!/bin/bash
# textdiff.sh <mfeprimer report> <engine report>: byte diff of everything between the timestamp line and the end of the
# parameter block (the footer differs on purpose: mfeprimer's citation vs the engine's own line).
norm() { sed -n '2,$p' "$1" | sed '/Annealing oligos/q'; }
diff <(norm "$1") <(norm "$2") >/dev/null && echo same || echo DIFF
