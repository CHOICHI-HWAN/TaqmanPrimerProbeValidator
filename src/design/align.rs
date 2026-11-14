//! The lineage alignment: majority consensus, per-genome mismatch prefix sums, greedy IUPAC codes and pass bits.

use super::oligo::{bit, code_of, iupac, is_acgt};
use std::io::{BufRead, BufReader};

/// A set of genomes (bit g = genome g of the alignment).
pub type Bits = Vec<u64>;

pub fn new_bits(n: usize) -> Bits {
    vec![0; n.div_ceil(64)]
}

#[inline]
pub fn set_bit(b: &mut Bits, i: usize) {
    b[i >> 6] |= 1 << (i & 63);
}

pub fn count(b: &[u64]) -> usize {
    b.iter().map(|w| w.count_ones() as usize).sum()
}

/// |a & b & c &^ not|
pub fn and3_count(a: &[u64], b: &[u64], c: &[u64], not: Option<&[u64]>) -> usize {
    let mut n = 0;
    for i in 0..a.len() {
        let mut w = a[i] & b[i] & c[i];
        if let Some(x) = not {
            w &= !x[i];
        }
        n += w.count_ones() as usize;
    }
    n
}

pub fn and3(a: &[u64], b: &[u64], c: &[u64]) -> Bits {
    (0..a.len()).map(|i| a[i] & b[i] & c[i]).collect()
}

/// Oligo kind: b'F' (3' end = right end of the window), b'R' (3' end = left end), b'P' (perfect match).
pub type Kind = u8;

pub struct Aligned {
    #[allow(dead_code)]
    pub ids: Vec<String>,
    pub g: usize,
    pub l: usize,
    pub a: Vec<Vec<u8>>,
    /// majority ACGT base where >= half the genomes have it, '-' otherwise
    pub cons: Vec<u8>,
    cs: Vec<Vec<i32>>,
    gapcs: Vec<i32>,
}

/// Reads an aligned FASTA (upper-cased).
pub fn read_aligned(path: &str) -> Result<(Vec<String>, Vec<Vec<u8>>), String> {
    let f = std::fs::File::open(path).map_err(|e| format!("{path}: {e}"))?;
    let (mut ids, mut rows): (Vec<String>, Vec<Vec<u8>>) = (Vec::new(), Vec::new());
    for ln in BufReader::new(f).lines() {
        let ln = ln.map_err(|e| e.to_string())?;
        let ln = ln.trim();
        if let Some(h) = ln.strip_prefix('>') {
            ids.push(h.split_whitespace().next().unwrap_or("x").to_string());
            rows.push(Vec::new());
        } else if !ln.is_empty() {
            if let Some(r) = rows.last_mut() {
                r.extend(ln.bytes().map(|b| b.to_ascii_uppercase()));
            }
        }
    }
    Ok((ids, rows))
}

impl Aligned {
    pub fn new(ids: Vec<String>, rows: Vec<Vec<u8>>) -> Result<Aligned, String> {
        if rows.is_empty() {
            return Err("empty alignment".into());
        }
        let (g, l) = (rows.len(), rows[0].len());
        for (i, r) in rows.iter().enumerate() {
            if r.len() != l {
                return Err(format!("{}: row length {}, expected {}", ids[i], r.len(), l));
            }
        }
        let mut cons = vec![b'-'; l];
        for c in 0..l {
            let mut n = [0usize; 4];
            for r in &rows {
                match r[c] {
                    b'A' => n[0] += 1,
                    b'C' => n[1] += 1,
                    b'G' => n[2] += 1,
                    b'T' => n[3] += 1,
                    _ => {}
                }
            }
            let mut top = 0;
            for i in 1..4 {
                if n[i] > n[top] {
                    top = i;
                }
            }
            if n[top] as f64 >= 0.5 * g as f64 {
                cons[c] = b"ACGT"[top];
            }
        }
        let cs = rows
            .iter()
            .map(|r| {
                let mut p = vec![0i32; l + 1];
                for c in 0..l {
                    p[c + 1] = p[c] + (r[c] != cons[c]) as i32;
                }
                p
            })
            .collect();
        let mut gapcs = vec![0i32; l + 1];
        for c in 0..l {
            gapcs[c + 1] = gapcs[c] + (cons[c] == b'-') as i32;
        }
        Ok(Aligned { ids, g, l, a: rows, cons, cs, gapcs })
    }

    /// The k-window at s has no '-' consensus column.
    pub fn conserved(&self, s: isize, k: usize) -> bool {
        s >= 0 && s as usize + k <= self.l && self.gapcs[s as usize + k] == self.gapcs[s as usize]
    }

    pub fn windows(&self, k: usize) -> Vec<usize> {
        if k > self.l {
            return Vec::new();
        }
        (0..=self.l - k).filter(|&s| self.gapcs[s + k] == self.gapcs[s]).collect()
    }

    /// Genomes with <= thr mismatches against the consensus window.
    pub fn upper_bound(&self, s: usize, k: usize, thr: usize) -> usize {
        self.cs.iter().filter(|p| (p[s + k] - p[s]) as usize <= thr).count()
    }

    /// Conserved columns in [s, e).
    pub fn amplicon_len(&self, s: usize, e: usize) -> usize {
        self.cons[s..e].iter().filter(|&&c| c != b'-').count()
    }

    /// Greedily turns up to amb consensus positions of the window [s, s+k) into two-base IUPAC codes, each where it
    /// adds the most genomes PASSING the rule (then the most brought within reach, then the most fixed); the codes in
    /// the order chosen. A prefix of the result is the greedy answer for a smaller amb.
    pub fn eval_codes(&self, s: usize, k: usize, kind: Kind, amb: usize, pass_mm: usize, tp_n: usize) -> Vec<(usize, u8)> {
        let g = self.g;
        let mut bad = vec![false; g * k];
        let mut nbad = vec![0usize; g];
        let mut endbad = vec![0usize; g];
        let mut end = vec![false; k];
        if kind == b'F' {
            for j in k.saturating_sub(tp_n)..k {
                end[j] = true;
            }
        } else if kind == b'R' {
            for j in 0..tp_n.min(k) {
                end[j] = true;
            }
        }
        for (gi, r) in self.a.iter().enumerate() {
            for j in 0..k {
                if r[s + j] != self.cons[s + j] {
                    bad[gi * k + j] = true;
                }
            }
        }
        let limit = if kind == b'P' { 0 } else { pass_mm };
        let mut coded = vec![false; k];
        let mut codes = Vec::new();
        let mut vals: Vec<u8> = Vec::with_capacity(g);
        for step in 0..amb {
            for gi in 0..g {
                let (mut nb, mut eb) = (0, 0);
                for j in 0..k {
                    if bad[gi * k + j] {
                        nb += 1;
                        if end[j] {
                            eb += 1;
                        }
                    }
                }
                nbad[gi] = nb;
                endbad[gi] = eb;
            }
            let left = amb - step - 1;
            let mut best: Option<(usize, u8, [i64; 4])> = None;
            for j in 0..k {
                vals.clear();
                for gi in 0..g {
                    if bad[gi * k + j] {
                        vals.push(self.a[gi][s + j]);
                    }
                }
                if vals.is_empty() || coded[j] || !is_acgt(self.cons[s + j]) {
                    continue;
                }
                for b in py_set_order(&vals) {
                    if !is_acgt(b) {
                        continue;
                    }
                    let (mut okc, mut nearc, mut fixc) = (0i64, 0i64, 0i64);
                    for gi in 0..g {
                        let fix = bad[gi * k + j] && self.a[gi][s + j] == b;
                        let mut nb = nbad[gi];
                        let mut fi = 0;
                        if fix {
                            nb -= 1;
                            fixc += 1;
                            if end[j] {
                                fi = 1;
                            }
                        }
                        if nb <= limit && (kind == b'P' || endbad[gi] - fi == 0) {
                            okc += 1;
                        }
                        if fix && nb <= limit + left {
                            nearc += 1;
                        }
                    }
                    let key = [okc, nearc, fixc, -(j as i64)];
                    if best.as_ref().map_or(true, |x| key > x.2) {
                        best = Some((j, b, key));
                    }
                }
            }
            let Some((bj, bb, key)) = best else { break };
            if key[0] + key[1] == 0 {
                break;
            }
            coded[bj] = true;
            codes.push((bj, bb));
            for gi in 0..g {
                if self.a[gi][s + bj] == bb {
                    bad[gi * k + bj] = false;
                }
            }
        }
        codes
    }

    /// The consensus window with the given codes applied (genome sense).
    pub fn oligo_with(&self, s: usize, k: usize, codes: &[(usize, u8)]) -> Vec<u8> {
        let mut w = self.cons[s..s + k].to_vec();
        for &(j, b) in codes {
            w[j] = code_of(bit(w[j]) | bit(b));
        }
        w
    }

    fn tab(oligo: &[u8]) -> Vec<[bool; 256]> {
        oligo
            .iter()
            .map(|&c| {
                let mut t = [false; 256];
                let b = iupac(c);
                for &x in if b.is_empty() { b"ACGT" as &[u8] } else { b } {
                    t[x as usize] = true;
                }
                t
            })
            .collect()
    }

    /// Genomes the genome-sense oligo passes at window s (P: perfect; F/R: <= pass_mm mismatches, none in the tp_n
    /// bases at the 3' end).
    pub fn pass_bits(&self, oligo: &[u8], s: usize, kind: Kind, pass_mm: usize, tp_n: usize) -> Bits {
        let k = oligo.len();
        let tab = Self::tab(oligo);
        let mut out = new_bits(self.g);
        for (gi, r) in self.a.iter().enumerate() {
            let row = &r[s..s + k];
            let (mut nb, mut end_bad) = (0, false);
            for j in 0..k {
                if !tab[j][row[j] as usize] {
                    nb += 1;
                    if (kind == b'F' && j + tp_n >= k) || (kind == b'R' && j < tp_n) {
                        end_bad = true;
                    }
                }
            }
            let ok = if kind == b'P' { nb == 0 } else { nb <= pass_mm && !end_bad };
            if ok {
                set_bit(&mut out, gi);
            }
        }
        out
    }

    /// The window where most genomes match the genome-sense oligo with <= 3 mismatches (IUPAC-aware).
    pub fn locate(&self, oligo: &[u8]) -> Option<(usize, usize)> {
        let k = oligo.len();
        if k > self.l {
            return None;
        }
        let tab = Self::tab(oligo);
        let (mut best, mut best_n) = (0usize, -1i64);
        for s in 0..=self.l - k {
            let mut n = 0i64;
            for r in &self.a {
                let row = &r[s..s + k];
                let mut mm = 0;
                for j in 0..k {
                    if !tab[j][row[j] as usize] {
                        mm += 1;
                        if mm > 3 {
                            break;
                        }
                    }
                }
                if mm <= 3 {
                    n += 1;
                }
            }
            if n > best_n {
                best = s;
                best_n = n;
            }
        }
        (best_n > 0).then_some((best, best_n as usize))
    }
}

/// The distinct values of vals in the iteration order of CPython 3.12's set(vals) (small ints hash to themselves), so
/// ties between bases break exactly like the Python engine.
pub fn py_set_order(vals: &[u8]) -> Vec<u8> {
    let mut table: Vec<Option<usize>> = vec![None; 8];
    let mut mask = 7usize;
    let mut fill = 0usize;
    fn insert_clean(t: &mut [Option<usize>], m: usize, h: usize) {
        let (mut i, mut perturb) = (h & m, h);
        loop {
            if t[i].is_none() {
                t[i] = Some(h);
                return;
            }
            if i + 9 <= m {
                for j in 1..=9 {
                    if t[i + j].is_none() {
                        t[i + j] = Some(h);
                        return;
                    }
                }
            }
            perturb >>= 5;
            i = (i * 5 + 1 + perturb) & m;
        }
    }
    for &b in vals {
        let h = b as usize;
        let (mut i, mut perturb) = (h & mask, h);
        let mut found = false;
        'probe: loop {
            let probes = if i + 9 <= mask { 9 } else { 0 };
            for j in 0..=probes {
                match table[i + j] {
                    None => {
                        table[i + j] = Some(h);
                        fill += 1;
                        break 'probe;
                    }
                    Some(v) if v == h => {
                        found = true;
                        break 'probe;
                    }
                    _ => {}
                }
            }
            perturb >>= 5;
            i = (i * 5 + 1 + perturb) & mask;
        }
        if found || fill * 5 < mask * 3 {
            continue;
        }
        let mut size = 8;
        while size <= fill * 4 {
            size <<= 1;
        }
        let mut nt = vec![None; size];
        for e in table.iter().flatten() {
            insert_clean(&mut nt, size - 1, *e);
        }
        table = nt;
        mask = size - 1;
    }
    table.iter().flatten().map(|&v| v as u8).collect()
}
