"""Nothing that is (or looks like) a primer, probe, standard-material or genome sequence may be in git: they are confidential and the
repository goes to GitHub as CODE ONLY. Checks every tracked file and every file `git add .` would add (untracked, not ignored).

 1. no sequence / data file types (FASTA, GenBank, spreadsheets, archives, slides);
 2. no long nucleotide literal in a text file (>= 16 plain bases, >= 18 with IUPAC codes). Tests generate what they need (tests/seqgen.py);
    short toy strings, IUPAC *alphabets* (every letter different) and trivially periodic strings ("ACGTACGT...") are fine;
 3. if the real configuration is present locally (config/primers.tsv, config/standards.tsv), no 12-base piece of any of its sequences
    (either strand) may appear in a tracked text file, and no FASTA record at all.
"""
import csv
import re
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
FORBIDDEN_EXT = {".fa", ".fasta", ".fna", ".fastq", ".fq", ".fai", ".gb", ".gbk", ".numbers", ".pptx", ".xlsx", ".gz", ".zst", ".tar"}
SKIP = ("webapp/static/vendor/", "node_modules/", "tests/test_no_sequences_in_repo.py")
PLAIN = re.compile(r"(?<![A-Za-z0-9_])[ACGTUacgtu]{16,}(?![A-Za-z0-9_])")
IUPAC = re.compile(r"(?<![A-Za-z0-9_])[ACGTURYKMSWBDHVNacgturykmswbdhvn]{18,}(?![A-Za-z0-9_])")
RUN = re.compile(r"[ACGTURYKMSWBDHVN]{12,}", re.I)
K = 12
COMP = str.maketrans("ACGTURYKMSWBDHVN", "TGCAAYRMKSWVHDBN")


def candidate_files():
    out = subprocess.run(["git", "ls-files", "--cached", "--others", "--exclude-standard"], cwd=ROOT, capture_output=True, text=True)
    if out.returncode != 0:
        pytest.skip("not a git checkout")
    return [f for f in out.stdout.split("\n") if f and not f.startswith(SKIP)]


def read_text(f):
    try:
        return (ROOT / f).read_text(encoding="utf-8")
    except (UnicodeDecodeError, FileNotFoundError, IsADirectoryError):
        return None


def is_alphabet_or_pattern(run):
    s = run.upper()
    if len(set(s)) == len(s) or len(set(s)) == 1:               # an alphabet, or AAAA... / NNNN...
        return True
    return any(s == (s[:p] * (len(s) // p + 1))[:len(s)] for p in range(1, 5))     # ACGTACGT..., ATATAT...


def confidential_kmers():
    seqs = []
    for name in ("primers.tsv", "standards.tsv"):
        p = ROOT / "config" / name
        if p.exists():
            with open(p, newline="", encoding="utf-8") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    seqs += [v.strip().upper() for v in row.values() if v and re.fullmatch(r"[ACGTURYKMSWBDHVN]{%d,}" % K, v.strip(), re.I)]
    out = set()
    for s in seqs:
        for t in (s, s.translate(COMP)[::-1]):
            out |= {t[i:i + K] for i in range(len(t) - K + 1)}
    return out


def test_no_sequence_files_are_tracked():
    bad = [f for f in candidate_files() if Path(f).suffix.lower() in FORBIDDEN_EXT]
    assert not bad, f"sequence/data files must not be in git: {bad}"


def test_no_long_nucleotide_literals_in_tracked_text():
    found = []
    for f in candidate_files():
        text = read_text(f)
        if text is None:
            continue
        for n, line in enumerate(text.split("\n"), 1):
            for pat in (PLAIN, IUPAC):
                for m in pat.finditer(line):
                    if not is_alphabet_or_pattern(m.group(0)):
                        found.append(f"{f}:{n} ({len(m.group(0))} nt)")
        if re.search(r"^>[^\s>]+[^\n]*\n[ACGTUNacgtun]{20,}", text, re.M):
            found.append(f"{f}: FASTA record")
    assert not found, "sequence-like strings in git (never commit primers/standards/genomes; generate test data):\n  " + "\n  ".join(found[:20])


def test_no_piece_of_the_real_primers_or_standards_in_tracked_text():
    kmers = confidential_kmers()
    if not kmers:
        pytest.skip("no local config/primers.tsv or standards.tsv to compare with")
    found = []
    for f in candidate_files():
        text = read_text(f)
        if text is None:
            continue
        for n, line in enumerate(text.split("\n"), 1):
            for m in RUN.finditer(line):
                s = m.group(0).upper()
                if any(s[i:i + K] in kmers for i in range(len(s) - K + 1)):
                    found.append(f"{f}:{n}")
                    break
    assert not found, "a real primer/standard sequence (or a piece of one) is in git:\n  " + "\n  ".join(found[:20])
