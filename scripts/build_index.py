#!/usr/bin/env python3
"""Build the mfeprimer index for one segment's filtered database once, so that the many
parallel primer_match/probe_match jobs for that segment don't race to build it themselves."""
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mfe_common import find_mfeprimer  # noqa: E402


def main():
    db_fasta = Path(snakemake.input[0])
    pipeline_root = Path(snakemake.scriptdir).parent
    mfeprimer_bin = find_mfeprimer(snakemake.config.get("mfeprimer_bin", ""), pipeline_root)
    cpu = int(snakemake.config.get("mfeprimer_cpu", 0) or 0) or None
    cmd = [mfeprimer_bin, "index", "-i", str(db_fasta)] + (["--cpu", str(cpu)] if cpu else [])
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(f"[error] mfeprimer index failed for {db_fasta}:\n{result.stderr}")


if __name__ == "__main__":
    main()
