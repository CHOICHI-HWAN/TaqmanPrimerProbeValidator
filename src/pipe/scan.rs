//! `chichipes scan`: the whole-genome brute-force scans of the validation pipeline (mfe_common.py _oligo_sites, scan_fr,
//! scan_probe), result-for-result the same as the pipeline's go/scan, which it replaces. One JSON request on stdin, one JSON
//! answer on stdout (byte-identical to go/scan's); genomes are scanned in parallel and answered in the order given.
//!
//!   {"mode":"fr",    "f":"..","r":"..","three_prime_n":0,"max_mm":2,"max_size":2000,"genomes":[{"id":"x","seq":".."}]}
//!   {"mode":"probe", "probe":"..","max_mm":8,"genomes":[{"id":"x","seq":"..","amp_start":null|int,"amp_end":null|int}]}
//!   {"mode":"sites", "oligo":"..","max_mm":2,"genomes":[...]}
//!   {"mode":"multi", "oligos":[{"seq":"..","max_mm":2},...],"genomes":[...]}

use super::{out_int, out_str};
use crate::design::cover::parallel;
use crate::json::{self, Value};

/// IUPAC codes of the pipeline's scan (I = any base); unknown letters accept ACGT, as in Python.
fn iupac(c: u8) -> &'static [u8] {
    match c {
        b'A' => b"A",
        b'C' => b"C",
        b'G' => b"G",
        b'T' => b"T",
        b'R' => b"AG",
        b'Y' => b"CT",
        b'S' => b"GC",
        b'W' => b"AT",
        b'K' => b"GT",
        b'M' => b"AC",
        b'B' => b"CGT",
        b'D' => b"AGT",
        b'H' => b"ACT",
        b'V' => b"ACG",
        _ => b"ACGT",
    }
}

fn complement(c: u8) -> u8 {
    match c {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'I' => b'I',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'D' => b'H',
        b'H' => b'D',
        b'B' => b'V',
        b'V' => b'B',
        _ => b'N',
    }
}

pub fn rev_comp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&c| complement(c.to_ascii_uppercase())).collect()
}

fn up(s: &[u8]) -> Vec<u8> {
    s.to_ascii_uppercase()
}

#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Site {
    pub start: usize,
    pub end: usize,
    pub plus: bool,
    pub mm: usize,
}

/// _oligo_sites(): every window with <= max_mm mismatches, "+" strand first, then "-", each by ascending start.
pub fn oligo_sites(oligo: &[u8], genome: &[u8], max_mm: usize) -> Vec<Site> {
    let mut sites = Vec::new();
    let n = genome.len();
    for (plus, q) in [(true, oligo.to_vec()), (false, rev_comp(oligo))] {
        let l = q.len();
        let mut tab = vec![[false; 256]; l];
        for (j, &c) in q.iter().enumerate() {
            for &b in iupac(c) {
                tab[j][b as usize] = true;
            }
        }
        if l > n {
            continue;
        }
        for i in 0..=n - l {
            let mut m = 0;
            let mut j = 0;
            while j < l {
                if !tab[j][genome[i + j] as usize] {
                    m += 1;
                    if m > max_mm {
                        break;
                    }
                }
                j += 1;
            }
            if j == l {
                sites.push(Site { start: i, end: i + l, plus, mm: m });
            }
        }
    }
    sites
}

fn iupac_mismatch(q: &[u8], t: &[u8]) -> usize {
    let n = q.len().min(t.len());
    let m = (0..n).filter(|&i| !iupac(q[i].to_ascii_uppercase()).contains(&t[i].to_ascii_uppercase())).count();
    m + q.len().abs_diff(t.len())
}

fn terminal_ok(query: &[u8], target: &[u8], is_fp: bool, n: usize) -> bool {
    let n = n.min(query.len()).min(target.len());
    if n == 0 {
        return true;
    }
    if is_fp {
        iupac_mismatch(&query[query.len() - n..], &target[target.len() - n..]) == 0
    } else {
        iupac_mismatch(&query[..n], &target[..n]) == 0
    }
}

pub struct Fr {
    pub f_mm: usize,
    pub r_mm: usize,
    pub f3p: bool,
    pub r3p: bool,
    pub amp_s: usize,
    pub amp_e: usize,
    pub f: Site,
    pub r: Site,
}

/// scan_fr(): best F/R pair by (3' ends intact, total mismatches, shorter amplicon); the first of equal keys wins.
pub fn scan_fr(f_seq: &[u8], r_seq: &[u8], genome: &[u8], three_prime_n: usize, max_mm: usize, max_size: usize) -> Option<Fr> {
    let g = up(genome);
    let (fu, ru) = (up(f_seq), up(r_seq));
    let (fs, rs) = (oligo_sites(&fu, &g, max_mm), oligo_sites(&ru, &g, max_mm));
    let (frc, rrc) = (rev_comp(&fu), rev_comp(&ru));
    let mut best: Option<(( usize, usize, usize), Fr)> = None;
    for f in &fs {
        for r in &rs {
            let (left, right) = if f.plus && !r.plus && f.start <= r.start && r.end > f.start && r.end - f.start <= max_size {
                (f, r)
            } else if r.plus && !f.plus && r.start <= f.start && f.end > r.start && f.end - r.start <= max_size {
                (r, f)
            } else {
                continue;
            };
            let fq = if f.plus { &fu } else { &frc };
            let rq = if r.plus { &ru } else { &rrc };
            let fok = terminal_ok(fq, &g[f.start..f.end], f.plus, three_prime_n);
            let rok = terminal_ok(rq, &g[r.start..r.end], r.plus, three_prime_n);
            let key = (if fok && rok { 0 } else { 1 }, f.mm + r.mm, right.end - left.start);
            if best.as_ref().map_or(true, |b| key < b.0) {
                best = Some((
                    key,
                    Fr { f_mm: f.mm, r_mm: r.mm, f3p: fok, r3p: rok, amp_s: left.start + 1, amp_e: right.end, f: *f, r: *r },
                ));
            }
        }
    }
    best.map(|b| b.1)
}

/// scan_probe(): the best window (fewest mismatches, first wins) inside [amp_s, amp_e] (1-based inclusive) when given.
pub fn scan_probe(probe: &[u8], genome: &[u8], amp: Option<(i64, i64)>, max_mm: usize) -> Option<(usize, usize, usize)> {
    let mut best: Option<(usize, usize, usize)> = None;
    for s in oligo_sites(&up(probe), &up(genome), max_mm) {
        if let Some((a, b)) = amp {
            if !(a <= s.start as i64 + 1 && s.end as i64 <= b) {
                continue;
            }
        }
        if best.map_or(true, |x| s.mm < x.0) {
            best = Some((s.mm, s.start + 1, s.end));
        }
    }
    best
}

fn write_site(o: &mut Vec<u8>, s: &Site) {
    o.extend_from_slice(b"{\"start\":");
    out_int(o, s.start as i64);
    o.extend_from_slice(b",\"end\":");
    out_int(o, s.end as i64);
    o.extend_from_slice(if s.plus { b",\"strand\":\"+\",\"mm\":" } else { b",\"strand\":\"-\",\"mm\":" });
    out_int(o, s.mm as i64);
    o.push(b'}');
}

fn write_sites(o: &mut Vec<u8>, v: &[Site]) {
    o.push(b'[');
    for (i, s) in v.iter().enumerate() {
        if i > 0 {
            o.push(b',');
        }
        write_site(o, s);
    }
    o.push(b']');
}

fn write_pos(o: &mut Vec<u8>, s: &Site) {
    o.extend_from_slice(b"{\"start\":");
    out_int(o, s.start as i64 + 1);
    o.extend_from_slice(b",\"end\":");
    out_int(o, s.end as i64);
    o.extend_from_slice(if s.plus { b",\"strand\":\"+\"}" } else { b",\"strand\":\"-\"}" });
}

fn int_field(v: &Value, k: &str) -> usize {
    v.get(k).as_f64() as usize
}

fn str_field(v: &Value, k: &str) -> Vec<u8> {
    match v.get(k) {
        Value::Str(s) => s.as_bytes().to_vec(),
        _ => Vec::new(),
    }
}

fn opt_int(v: &Value) -> Option<i64> {
    match v {
        Value::Num(x) => Some(*x as i64),
        _ => None,
    }
}

/// Answers one request (the JSON text go/scan would write, newline included).
pub fn run(input: &str) -> Result<Vec<u8>, String> {
    let req = json::parse(input).map_err(|e| format!("bad request: {e}"))?;
    let mode = String::from_utf8(str_field(&req, "mode")).unwrap_or_default();
    if !matches!(mode.as_str(), "fr" | "probe" | "sites" | "multi") {
        return Err(format!("error: unknown mode {mode:?}"));
    }
    let genomes = req.get("genomes").as_array();
    let (f, r, probe, oligo) = (str_field(&req, "f"), str_field(&req, "r"), str_field(&req, "probe"), str_field(&req, "oligo"));
    let (tpn, max_mm, max_size) = (int_field(&req, "three_prime_n"), int_field(&req, "max_mm"), int_field(&req, "max_size"));
    let oligos: Vec<(Vec<u8>, usize)> =
        req.get("oligos").as_array().iter().map(|o| (up(&str_field(o, "seq")), int_field(o, "max_mm"))).collect();
    let workers = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
    let parts = parallel(genomes.len(), workers, |i| {
        let g = &genomes[i];
        let seq = str_field(g, "seq");
        let mut o = Vec::new();
        o.extend_from_slice(b"{\"id\":");
        json::write_go_string(&mut o, &str_field(g, "id"));
        o.extend_from_slice(b",\"hit\":");
        match mode.as_str() {
            "fr" => match scan_fr(&f, &r, &seq, tpn, max_mm, max_size) {
                Some(h) => {
                    o.extend_from_slice(b"{\"f_mismatch\":");
                    out_int(&mut o, h.f_mm as i64);
                    o.extend_from_slice(b",\"r_mismatch\":");
                    out_int(&mut o, h.r_mm as i64);
                    o.extend_from_slice(b",\"total\":");
                    out_int(&mut o, (h.f_mm + h.r_mm) as i64);
                    o.extend_from_slice(if h.f3p { b",\"f_3p_ok\":true" } else { b",\"f_3p_ok\":false" });
                    o.extend_from_slice(if h.r3p { b",\"r_3p_ok\":true" } else { b",\"r_3p_ok\":false" });
                    o.extend_from_slice(b",\"amp_start\":");
                    out_str(&mut o, &h.amp_s.to_string());
                    o.extend_from_slice(b",\"amp_end\":");
                    out_str(&mut o, &h.amp_e.to_string());
                    o.extend_from_slice(b",\"f_tm\":\"\",\"f_dg\":\"\",\"r_tm\":\"\",\"r_dg\":\"\",\"f_pos\":");
                    write_pos(&mut o, &h.f);
                    o.extend_from_slice(b",\"r_pos\":");
                    write_pos(&mut o, &h.r);
                    o.push(b'}');
                }
                None => o.extend_from_slice(b"null"),
            },
            "probe" => {
                let amp = match (opt_int(g.get("amp_start")), opt_int(g.get("amp_end"))) {
                    (Some(a), Some(b)) => Some((a, b)),
                    _ => None,
                };
                match scan_probe(&probe, &seq, amp, max_mm) {
                    Some((mm, s, e)) => {
                        o.extend_from_slice(b"{\"mismatch\":");
                        out_int(&mut o, mm as i64);
                        o.extend_from_slice(b",\"start\":");
                        out_str(&mut o, &s.to_string());
                        o.extend_from_slice(b",\"end\":");
                        out_str(&mut o, &e.to_string());
                        o.extend_from_slice(b",\"tm\":\"\",\"dg\":\"\"}");
                    }
                    None => o.extend_from_slice(b"null"),
                }
            }
            "sites" => {
                o.extend_from_slice(b"null");
                let s = oligo_sites(&up(&oligo), &up(&seq), max_mm);
                if !s.is_empty() {
                    o.extend_from_slice(b",\"sites\":");
                    write_sites(&mut o, &s);
                }
            }
            _ => {
                o.extend_from_slice(b"null");
                if !oligos.is_empty() {
                    let gs = up(&seq);
                    o.extend_from_slice(b",\"multi\":[");
                    for (k, (q, mm)) in oligos.iter().enumerate() {
                        if k > 0 {
                            o.push(b',');
                        }
                        write_sites(&mut o, &oligo_sites(q, &gs, *mm));
                    }
                    o.push(b']');
                }
            }
        }
        o.push(b'}');
        o
    });
    let mut out = b"{\"results\":[".to_vec();
    for (i, p) in parts.iter().enumerate() {
        if i > 0 {
            out.push(b',');
        }
        out.extend_from_slice(p);
    }
    out.extend_from_slice(b"]}\n");
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sites_and_pairs() {
        let g = b"TTTTACGTACGGTTAGCAAAAAAAAAAAAAAAAGGCATTGCAAAA";
        let s = oligo_sites(b"ACGTACGG", g, 0);
        assert_eq!(s, vec![Site { start: 4, end: 12, plus: true, mm: 0 }]);
        // R binds the minus strand: its reverse complement TTGCC... sits at 33
        let fr = scan_fr(b"ACGTACGG", b"TTGCAATGCC", g, 2, 1, 2000).unwrap();
        assert_eq!((fr.amp_s, fr.amp_e, fr.f_mm, fr.r_mm), (5, 43, 0, 0));
        assert_eq!(scan_probe(b"TTAGCA", g, Some((5, 43)), 2), Some((0, 13, 18)));
        assert_eq!(scan_probe(b"TTAGCA", g, Some((20, 43)), 0), None);
    }
}
