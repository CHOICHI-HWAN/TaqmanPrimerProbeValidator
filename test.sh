#!/bin/bash
# Checks the core of the pipeline.
#
#   ./test.sh            core checks (about 10 s): tools, then synthetic genomes with a KNOWN answer run through the
#                        real mfeprimer + recount + genome scan + probe + orientation, and the unit tests
#   ./test.sh --full     + every test, `snakemake -n` (the workflow builds), and the audits of the REAL results:
#                        audit_recount (mismatch counts), audit_pairs (no valid F/R pair missed), audit_strands
#                        (reverse-complement handling), and the mfeprimer/genome-scan concordance
#
# Exit status 0 only if every check passed. Run ./setup.sh first on a new machine.
set -uo pipefail
cd "$(dirname "$(readlink -f "$0")")"

FULL=0
for a in "$@"; do
    case "$a" in
        --full) FULL=1 ;;
        -h|--help) sed -n '2,12p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "unknown option: $a (see --help)" >&2; exit 2 ;;
    esac
done

# the conda env that setup.sh recorded (snakemake, mafft, pixy, numpy ... on PATH), if any
ENV_PREFIX="$(sed -n 's/^env_prefix:[[:space:]]*"\{0,1\}\([^"#]*\)"\{0,1\}.*/\1/p' config/config.yaml 2>/dev/null | head -n 1)"
if [ -n "$ENV_PREFIX" ] && [ -d "$ENV_PREFIX/bin" ]; then export PATH="$ENV_PREFIX/bin:$PATH"; fi
PY="$(command -v python3)"
[ -x "/BI/conda/envs/lasv_pipeline/bin/python" ] && ! "$PY" -c 'import pytest, numpy' 2>/dev/null && PY=/BI/conda/envs/lasv_pipeline/bin/python
SNAKEMAKE="$(command -v snakemake || true)"
[ -z "$SNAKEMAKE" ] && [ -x /BI/conda/envs/lasv_pipeline/bin/snakemake ] && SNAKEMAKE=/BI/conda/envs/lasv_pipeline/bin/snakemake

[ -n "$SNAKEMAKE" ] && export PATH="$(dirname "$SNAKEMAKE"):$PATH"     # mafft / pixy / bgzip live next to it in the same env

PASS=0; FAIL=0; SKIP=0; FAILED=()
ok()   { printf '  \033[32m✔\033[0m %s\n' "$1"; PASS=$((PASS+1)); }
bad()  { printf '  \033[31m✘\033[0m %s\n' "$1"; FAIL=$((FAIL+1)); FAILED+=("$1"); }
skip() { printf '  \033[33m-\033[0m %s (skipped: %s)\n' "$1" "$2"; SKIP=$((SKIP+1)); }
section() { printf '\n\033[1m%s\033[0m\n' "$1"; }
run() {   # label, command...   (shows the last lines of output only when it fails)
    local label="$1"; shift
    local log; log="$(mktemp)"
    if "$@" >"$log" 2>&1; then ok "$label"; else bad "$label"; tail -n 12 "$log" | sed 's/^/      /'; fi
    rm -f "$log"
}

section "1. Tools"
"$PY" -c 'import sys; assert sys.version_info >= (3, 9)' 2>/dev/null && ok "python $("$PY" -V 2>&1 | cut -d' ' -f2)" || bad "python >= 3.9"
for mod in numpy pytest; do "$PY" -c "import $mod" 2>/dev/null && ok "python module $mod" || bad "python module $mod (run ./setup.sh)"; done
if [ -n "$SNAKEMAKE" ]; then ok "snakemake $("$SNAKEMAKE" --version)"; else bad "snakemake (run ./setup.sh)"; fi
for t in mafft pixy bgzip tabix; do command -v "$t" >/dev/null 2>&1 && ok "$t" || bad "$t (run ./setup.sh)"; done
MFE_STATUS="$("$PY" - <<'PY' 2>&1
import sys
sys.path.insert(0, "scripts")
from pathlib import Path
from mfe_common import check_mfeprimer, find_mfeprimer
try:
    r = check_mfeprimer(find_mfeprimer("", Path(".").resolve()))
except SystemExit:
    r = {"ok": False, "error": "not found -- set mfeprimer_bin (setup.sh or the Analysis Parameters page)"}
print(("OK " + r["version"]) if r["ok"] else ("NO " + r["error"]))
PY
)"
case "$MFE_STATUS" in OK*) ok "mfeprimer (${MFE_STATUS#OK })" ;; *) bad "mfeprimer: ${MFE_STATUS#NO }" ;; esac

section "2. Core: synthetic genomes with a known answer (real mfeprimer, recount, genome scan, probe, orientation)"
run "known-answer end-to-end" "$PY" -m pytest -q tests/test_core_e2e.py
run "mismatch counted on the genome, scan_fr / choose_fr / probe, orientation" "$PY" -m pytest -q tests/test_parse_spec_best.py tests/test_orientation.py
run "primer/probe helpers, manifests, standards, profile" "$PY" -m pytest -q tests/test_mfe_common.py tests/test_primers_manifest.py tests/test_standards_manifest.py tests/test_virus_profile.py tests/test_params_mfeprimer.py

if [ "$FULL" = 1 ]; then
    section "3. Whole test suite and the workflow"
    run "all tests" "$PY" -m pytest -q
    if [ -n "$SNAKEMAKE" ]; then run "snakemake -n (the workflow builds)" "$SNAKEMAKE" -n; else skip "snakemake -n" "snakemake missing"; fi

    section "4. Audits of the real results"
    if ls output/primers/P*/alignment/mismatch_table.tsv >/dev/null 2>&1; then
        run "audit_recount: no mismatch under-counted (F, R, probe)" "$PY" scripts/audit_recount.py
        run "audit_pairs: no valid F/R pair missed" bash -c "'$PY' scripts/audit_pairs.py | grep -v 'missed by pipeline:   0' | grep -q . && exit 1 || exit 0"
        run "audit_live: the web app's live analysis agrees with the pipeline" "$PY" scripts/audit_live.py
        run "audit_strands: reverse-complement handling consistent" bash -c "'$PY' scripts/audit_strands.py | tail -1 | grep -q \"'bad': 0\""
        if [ -f output/summary/engine_concordance.tsv ]; then
            run "engines: genome scan never worse than mfeprimer" "$PY" - <<'PY'
import csv, sys
rows = list(csv.DictReader(open("output/summary/engine_concordance.tsv"), delimiter="\t"))
bad = [r["set_id"] for r in rows if r["mfeprimer_only_pass"] != "0" or r["probe_mfeprimer_only_perfect"] != "0"]
sys.exit(f"mfeprimer passes genomes the scan fails in: {bad}" if bad else 0)
PY
        else skip "engine concordance" "output/summary/engine_concordance.tsv missing (run snakemake)"; fi
    else
        skip "audits" "no results yet (run snakemake first)"
    fi
fi

section "Summary"
printf '  passed %d   failed %d   skipped %d\n' "$PASS" "$FAIL" "$SKIP"
if [ "$FAIL" -gt 0 ]; then
    printf '  \033[31mFAILED:\033[0m\n'; for f in "${FAILED[@]}"; do printf '    - %s\n' "$f"; done
    exit 1
fi
printf '  \033[32mAll checks passed.\033[0m%s\n' "$([ "$FULL" = 0 ] && echo '  (./test.sh --full also audits the real results)')"
