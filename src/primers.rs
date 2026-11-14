//! Reads oligos and expands them the way MFEprimer names them: degenerate bases into concrete sequences (NAME.1,
//! NAME.2, ...) and each concrete sequence into its 3'-window mismatch variants (NAME.1.1, ...).

use std::collections::HashSet;

/// Expansion order of each IUPAC code, as `mfeprimer degen` writes it (docs/FINDINGS.md).
fn degen(b: u8) -> Option<&'static [u8]> {
    Some(match b {
        b'A' => b"A",
        b'C' => b"C",
        b'G' => b"G",
        b'T' => b"T",
        b'R' => b"GA",
        b'Y' => b"TC",
        b'M' => b"AC",
        b'K' => b"GT",
        b'S' => b"GC",
        b'W' => b"TA",
        b'B' => b"GTC",
        b'D' => b"GTA",
        b'H' => b"TAC",
        b'V' => b"GAC",
        b'N' => b"GATC",
        b'I' => b"I", // deoxyinosine pairs with every base: kept as is, never expanded
        _ => return None,
    })
}

/// Substitution order for a mismatch at one position, per original base (read off real reports).
fn sub_order(b: u8) -> &'static [u8] {
    match b {
        b'A' => b"TCG",
        b'C' => b"TAG",
        b'G' => b"TCA",
        b'T' => b"ACG",
        _ => b"",
    }
}

/// One input sequence.
#[derive(Clone, Debug)]
pub struct Oligo {
    pub name: String,
    pub seq: String,
}

/// Reads oligos from a FASTA file (name = first word of the header), or from "name seq" lines.
pub fn read_fasta(path: &str) -> Result<Vec<Oligo>, String> {
    let data = std::fs::read(path).map_err(|e| format!("open {path}: {e}"))?;
    let text = String::from_utf8_lossy(&data);
    let mut out: Vec<Oligo> = Vec::new();
    for line in text.split('\n') {
        let l = line.trim();
        if l.is_empty() {
            continue;
        }
        if let Some(h) = l.strip_prefix('>') {
            match h.split_whitespace().next() {
                Some(n) => out.push(Oligo { name: n.to_string(), seq: String::new() }),
                None => return Err(format!("{path}: empty FASTA header")),
            }
            continue;
        }
        let f: Vec<&str> = l.split_whitespace().collect();
        if !out.is_empty() && f.len() == 1 {
            out.last_mut().unwrap().seq.push_str(&l.to_uppercase());
            continue;
        }
        if f.len() >= 2 {
            out.push(Oligo { name: f[0].to_string(), seq: f[1].to_uppercase() });
            continue;
        }
        return Err(format!("{path}: cannot parse line {l:?}"));
    }
    Ok(out)
}

/// The concrete sequences of a (possibly degenerate) oligo in MFEprimer's order.
pub fn expand(seq: &str) -> Result<Vec<String>, String> {
    let mut out: Vec<Vec<u8>> = vec![Vec::new()];
    for &b in seq.as_bytes() {
        let alts = degen(b).ok_or_else(|| format!("unsupported base '{}' in {}", b as char, seq))?;
        let mut next = Vec::with_capacity(out.len() * alts.len());
        for p in &out {
            for &a in alts {
                let mut x = p.clone();
                x.push(a);
                next.push(x);
            }
        }
        out = next;
    }
    Ok(out.into_iter().map(|v| String::from_utf8(v).unwrap()).collect())
}

/// One 3'-window variant of a concrete primer: index 0 is the primer itself (name NAME.k), index i > 0 is NAME.k.i.
#[derive(Clone, Debug)]
pub struct Variant {
    pub index: usize,
    pub seq: String,
}

/// The concrete primer and every sequence with 1..max_mm substitutions at positions mis_start..mis_end counted from
/// the 3' end (1 = last base), in MFEprimer's numbering order: by number of mismatches, then by the 3'-nearest
/// position, then the next position, with the later position's substitution varying fastest.
pub fn variants(seq: &str, max_mm: usize, mis_start: usize, mis_end: usize) -> Vec<Variant> {
    let s = seq.as_bytes();
    let n = s.len();
    let mis_end = mis_end.min(n);
    // string indices, 3'-nearest first; an inosine already pairs with every base: no variants there
    let pos: Vec<usize> = (mis_start..=mis_end).map(|d| n - d).filter(|&p| s[p] != b'I').collect();
    let mut out = vec![Variant { index: 0, seq: seq.to_string() }];
    let mut idx = 1;
    let mut cur = s.to_vec();
    for k in 1..=max_mm {
        // position combinations first (3'-nearest first), then the substitutions of that combination
        for_combos(pos.len(), k, &mut |combo| {
            subs(s, &pos, combo, 0, &mut cur, &mut out, &mut idx);
        });
    }
    out
}

fn subs(s: &[u8], pos: &[usize], combo: &[usize], i: usize, cur: &mut Vec<u8>, out: &mut Vec<Variant>, idx: &mut usize) {
    if i == combo.len() {
        out.push(Variant { index: *idx, seq: String::from_utf8(cur.clone()).unwrap() });
        *idx += 1;
        return;
    }
    let p = pos[combo[i]];
    for &a in sub_order(s[p]) {
        cur[p] = a;
        subs(s, pos, combo, i + 1, cur, out, idx);
    }
    cur[p] = s[p];
}

/// Calls f with every k-subset of 0..n-1 in lexicographic order.
fn for_combos(n: usize, k: usize, f: &mut dyn FnMut(&[usize])) {
    fn rec(c: &mut Vec<usize>, i: usize, start: usize, n: usize, k: usize, f: &mut dyn FnMut(&[usize])) {
        if i == k {
            f(c);
            return;
        }
        if n < k - i {
            return;
        }
        for v in start..=n - (k - i) {
            c[i] = v;
            rec(c, i + 1, v + 1, n, k, f);
        }
    }
    let mut c = vec![0; k];
    rec(&mut c, 0, 0, n, k, f);
}

/// One concrete (non-degenerate) sequence of an oligo with the 3'-window variants searched for it.
#[derive(Clone, Debug)]
#[allow(dead_code)]
pub struct Concrete {
    /// NAME.k, or NAME when the oligo is not degenerate
    pub name: String,
    pub seq: String,
    /// [0] is seq itself; the rest are named name.i
    pub variants: Vec<Variant>,
}

impl Concrete {
    /// The name MFEprimer gives variant v.
    pub fn variant_name(&self, v: &Variant) -> String {
        if v.index == 0 {
            self.name.clone()
        } else {
            format!("{}.{}", self.name, v.index)
        }
    }
}

/// Expands an oligo into the concrete sequences and variants MFEprimer searches. The expansions of one oligo share a
/// set of sequences already generated: an expansion that is itself already generated (a variant of an earlier
/// expansion) gets no variants, and the others skip variants already generated; the kept variants are numbered 1, 2,
/// ... in generation order.
pub fn build(o: &Oligo, max_mm: usize, mis_start: usize, mis_end: usize) -> Result<Vec<Concrete>, String> {
    let seqs = expand(&o.seq).map_err(|e| format!("{}: {}", o.name, e))?;
    let mut seen: HashSet<String> = HashSet::new();
    let mut out = Vec::with_capacity(seqs.len());
    let multi = seqs.len() > 1;
    for (i, s) in seqs.into_iter().enumerate() {
        let name = if multi { format!("{}.{}", o.name, i + 1) } else { o.name.clone() };
        let mut c = Concrete { name, seq: s.clone(), variants: vec![Variant { index: 0, seq: s.clone() }] };
        if seen.insert(s.clone()) {
            let mut idx = 1;
            for v in variants(&s, max_mm, mis_start, mis_end).into_iter().skip(1) {
                if seen.insert(v.seq.clone()) {
                    c.variants.push(Variant { index: idx, seq: v.seq });
                    idx += 1;
                }
            }
        }
        out.push(c);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn expand_order() {
        assert_eq!(expand("ARN").unwrap(), vec!["AGG", "AGA", "AGT", "AGC", "AAG", "AAA", "AAT", "AAC"]);
        assert_eq!(expand("AIC").unwrap(), vec!["AIC"]);
    }

    #[test]
    fn variant_numbering() {
        let v = variants("ACGTACGTACGT", 2, 1, 9);
        assert_eq!(v.len(), 1 + 9 * 3 + 36 * 9);
        assert_eq!(v[1].seq, "ACGTACGTACGA"); // last base T -> A first
        assert_eq!(v[4].seq, "ACGTACGTACTT"); // then the next position, G -> T
    }

    /// The primer table at the top of a real mfeprimer report (every concrete sequence and variant, by name) must be
    /// reproduced exactly from the FR.fa it was run on. PE_REPORT points at a <set>_mm2.txt next to <set>_FR.fa, or
    /// EDITIONS at the pipeline editions; nothing confidential is copied into the repository.
    #[test]
    fn variants_match_report() {
        let path = std::env::var("PE_REPORT").unwrap_or_else(|_| {
            format!("{}/ebov/output/primers/P001/primer_qc/P001_mm2.txt", std::env::var("EDITIONS").unwrap_or_default())
        });
        let p = std::path::Path::new(&path);
        let base = p.file_name().and_then(|b| b.to_str()).unwrap_or("");
        let Some(cut) = base.find("_mm") else { return };
        let fr = p.with_file_name(format!("{}_FR.fa", &base[..cut]));
        let (Ok(oligos), Ok(report)) = (read_fasta(fr.to_str().unwrap()), std::fs::read_to_string(p)) else {
            eprintln!("skip: no report / FR.fa (set EDITIONS or PE_REPORT)");
            return;
        };
        let mut want = Vec::new();
        for o in &oligos {
            for c in build(o, 2, 1, 9).unwrap() {
                for v in &c.variants {
                    want.push((c.variant_name(v), v.seq.clone()));
                }
            }
        }
        let mut got = Vec::new();
        for l in report.lines() {
            if l.starts_with("Descriptions") {
                break;
            }
            let f: Vec<&str> = l.split_whitespace().collect();
            if f.len() < 2 {
                continue;
            }
            if oligos.iter().any(|o| f[0] == o.name || f[0].starts_with(&format!("{}.", o.name))) {
                got.push((f[0].to_string(), f[1].to_string()));
            }
        }
        assert_eq!(got, want);
    }
}
