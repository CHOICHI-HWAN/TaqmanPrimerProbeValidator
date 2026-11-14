//! `chichipes probedesign`: the rule-based ("own") path of the pipeline's scripts/probe_designer.py, bit-for-bit the same as
//! the pipeline's go/probedesign, which it replaces. One JSON request on stdin, the candidates as JSON on stdout.
//!
//!   {"consensus": "...", "f_col": [a,b], "r_col": [a,b], "dots": ["...", ...],
//!    "min_len": 18, "max_len": 30, "tm_range": null|[lo,hi], "gc_range": [30,80], "top": 10, "max_iupac": 2, "tm_above": null|[lo,hi]}

use super::{go_float, out_int, out_str};
use crate::design::oligo::{gc_percent, longest_run, py_sum, resolutions, tm_nn};
use crate::json::{self, Value};

const PROBE_TM_TARGET: f64 = 6.0;

/// tm_nn() of a plain ACGT oligo; None where the Python raises KeyError (any other letter).
fn tm_checked(seq: &[u8]) -> Option<f64> {
    let s = seq.to_ascii_uppercase();
    if s.is_empty() || (s.len() > 1 && s.iter().any(|c| !matches!(c, b'A' | b'C' | b'G' | b'T'))) {
        return None;
    }
    Some(tm_nn(&s))
}

fn mean_tm_gc(res: &[Vec<u8>]) -> Option<(f64, f64)> {
    let mut tms = Vec::with_capacity(res.len());
    let mut gcs = Vec::with_capacity(res.len());
    for r in res {
        tms.push(tm_checked(r)?);
        gcs.push(gc_percent(r));
    }
    let n = res.len() as f64;
    Some((py_sum(&tms) / n, py_sum(&gcs) / n))
}

/// rule_problems(): true when the oligo is NOT a usable probe.
fn rule_problems(seq: &[u8], tm: (f64, f64), gc: (f64, f64)) -> bool {
    let res = match resolutions(seq, 64) {
        Some(r) if !r.is_empty() && r.len() <= 64 => r,
        _ => return true,
    };
    let Some((t, g)) = mean_tm_gc(&res) else { return true };
    if !(tm.0 <= t && t <= tm.1) || !(gc.0 <= g && g <= gc.1) {
        return true;
    }
    res.iter().any(|r| r[0] == b'G' || longest_run(r) >= 4)
}

fn mask_bit(c: u8) -> u8 {
    match c {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' => 8,
        _ => 0,
    }
}

fn code_by_mask(m: u8) -> u8 {
    match m {
        5 => b'R',
        10 => b'Y',
        6 => b'S',
        9 => b'W',
        12 => b'K',
        3 => b'M',
        _ => 0,
    }
}

fn acgt(c: u8) -> bool {
    matches!(c, b'A' | b'C' | b'G' | b'T')
}

/// degenerate_window(): greedily add up to max_iupac two-base codes where the most genomes still need them.
fn degenerate_window(window: &[u8], row_mm: &[Vec<(usize, u8)>], max_iupac: usize) -> (Vec<u8>, Vec<usize>) {
    let mut seq = window.to_vec();
    let mut allowed = vec![0u8; window.len()];
    let mut left: Vec<Vec<(usize, u8)>> = row_mm.to_vec();
    for k in (1..=max_iupac).rev() {
        let mut order: Vec<(usize, u8)> = Vec::new();
        let mut score: Vec<usize> = Vec::new();
        for m in &left {
            if !m.is_empty() && m.len() <= k {
                for &(pos, base) in m {
                    if acgt(base) && base != window[pos] && allowed[pos] == 0 && acgt(window[pos]) {
                        match order.iter().position(|&x| x == (pos, base)) {
                            Some(i) => score[i] += 1,
                            None => {
                                order.push((pos, base));
                                score.push(1);
                            }
                        }
                    }
                }
            }
        }
        if order.is_empty() {
            break;
        }
        let mut b = 0;
        for i in 1..order.len() {
            if score[i] > score[b] || (score[i] == score[b] && order[i].0 < order[b].0) {
                b = i;
            }
        }
        let (pos, base) = order[b];
        allowed[pos] = mask_bit(window[pos]) | mask_bit(base);
        seq[pos] = code_by_mask(allowed[pos]);
        for m in left.iter_mut() {
            if m.contains(&(pos, base)) {
                m.retain(|&x| x != (pos, base));
            }
        }
    }
    (seq, left.iter().map(|m| m.len()).collect())
}

#[derive(Clone)]
struct Cand {
    seq: Vec<u8>,
    plus: bool,
    col_start: i64,
    col_end: i64,
    length: i64,
    n_iupac: usize,
    tm: f64,
    gc: f64,
    cov0: usize,
    cov1: usize,
    n: usize,
    source: &'static str,
}

fn round1(x: f64) -> f64 {
    format!("{x:.1}").parse().unwrap()
}

fn pair(v: &Value) -> Option<(f64, f64)> {
    match v {
        Value::Arr(a) if a.len() >= 2 => Some((a[0].as_f64(), a[1].as_f64())),
        _ => None,
    }
}

fn int_pair(v: &Value) -> (i64, i64) {
    let a = v.as_array();
    (a.first().map_or(0.0, |x| x.as_f64()) as i64, a.get(1).map_or(0.0, |x| x.as_f64()) as i64)
}

fn pick(c: &mut [Cand], quota: usize, primer_tm: f64) -> Vec<Cand> {
    c.sort_by(|a, b| {
        b.cov0
            .cmp(&a.cov0)
            .then(b.cov1.cmp(&a.cov1))
            .then(a.n_iupac.cmp(&b.n_iupac))
            .then_with(|| {
                let da = (a.tm - (primer_tm + PROBE_TM_TARGET)).abs();
                let db = (b.tm - (primer_tm + PROBE_TM_TARGET)).abs();
                da.partial_cmp(&db).unwrap()
            })
            .then(a.length.cmp(&b.length))
    });
    let mut picked: Vec<Cand> = Vec::new();
    for x in c.iter() {
        if picked.len() >= quota {
            break;
        }
        if !picked.iter().any(|p| p.plus == x.plus && (p.col_start - x.col_start).abs() <= 3) {
            picked.push(x.clone());
        }
    }
    picked
}

/// Answers one request (the JSON text go/probedesign would write, newline included).
pub fn run(input: &str) -> Result<Vec<u8>, (i32, String)> {
    let r = json::parse(input).map_err(|e| (2, format!("bad request: {e}")))?;
    let cons = match r.get("consensus") {
        Value::Str(s) => s.to_ascii_uppercase().into_bytes(),
        _ => Vec::new(),
    };
    let (fc, rc) = (int_pair(r.get("f_col")), int_pair(r.get("r_col")));
    let dots: Vec<Vec<u8>> = r
        .get("dots")
        .as_array()
        .iter()
        .map(|d| match d {
            Value::Str(s) => s.as_bytes().to_vec(),
            _ => Vec::new(),
        })
        .collect();
    let gi = |k: &str| r.get(k).as_f64() as i64;
    let (min_len, max_len, top, max_iupac) = (gi("min_len"), gi("max_len"), gi("top").max(0) as usize, gi("max_iupac").max(0) as usize);
    let gc_range = pair(r.get("gc_range")).unwrap_or((0.0, 0.0));
    let mut pt = f64::NEG_INFINITY;
    for c in [fc, rc] {
        let hi = ((c.1 + 1) as usize).min(cons.len());
        let lo = (c.0.max(0) as usize).min(hi);
        let w: Vec<u8> = cons[lo..hi].iter().copied().filter(|&b| b != b'-').collect();
        if w.is_empty() {
            continue;
        }
        let t = tm_checked(&w).ok_or((1, "error: primer window is not plain ACGT".to_string()))?;
        pt = pt.max(t);
    }
    if pt == f64::NEG_INFINITY {
        return Err((1, "error: the consensus has no bases under the F/R primer windows".into()));
    }
    let above = pair(r.get("tm_above")).unwrap_or((3.0, 10.0));
    let tm_range = pair(r.get("tm_range")).unwrap_or((pt + above.0, pt + above.1));
    let lo = fc.1.min(rc.1) + 1;
    let hi = fc.0.max(rc.0) - 1;
    let n_rows = dots.len();
    let lengths: Vec<i64> = (min_len..=max_len).collect();
    let parts = crate::design::cover::parallel(lengths.len(), lengths.len().max(1), |li| {
        let length = lengths[li];
        let (mut out, mut rejected) = (Vec::new(), 0usize);
        let mut s = lo;
        while s <= hi - length + 1 {
            let e = ((s + length) as usize).min(cons.len());
            let su = s as usize;
            let window = &cons[su.min(e)..e];
            let row_mm: Vec<Vec<(usize, u8)>> = dots
                .iter()
                .map(|d| (su..e.min(d.len())).filter(|&i| d[i] != b'.').map(|i| (i - su, d[i])).collect())
                .collect();
            let (window, left) = degenerate_window(window, &row_mm, max_iupac);
            let n_iupac = window.iter().filter(|&&c| !acgt(c)).count();
            let cov0 = left.iter().filter(|&&m| m == 0).count();
            let cov1 = left.iter().filter(|&&m| m <= 1).count();
            for plus in [true, false] {
                // go/probedesign's complement table has no I (it becomes N), unlike the scan's
                let oligo = if plus {
                    window.clone()
                } else {
                    super::scan::rev_comp(&window).into_iter().map(|c| if c == b'I' { b'N' } else { c }).collect()
                };
                if rule_problems(&oligo, tm_range, gc_range) {
                    rejected += 1;
                    continue;
                }
                let res = resolutions(&oligo, usize::MAX).unwrap_or_default();
                let (tm, gc) = mean_tm_gc(&res).unwrap_or((f64::NAN, f64::NAN));
                out.push(Cand {
                    seq: oligo,
                    plus,
                    col_start: s + 1,
                    col_end: e as i64,
                    length,
                    n_iupac,
                    tm: round1(tm),
                    gc: round1(gc),
                    cov0,
                    cov1,
                    n: n_rows,
                    source: "own",
                });
            }
            s += 1;
        }
        (out, rejected)
    });
    let mut all: Vec<Cand> = Vec::new();
    let mut rejected = 0;
    for (o, rj) in parts {
        all.extend(o);
        rejected += rj;
    }
    let n_total = all.len();
    let picked = pick(&mut all, top, pt);
    let mut uniq: Vec<Cand> = Vec::new();
    for c in picked {
        if let Some(u) = uniq.iter_mut().find(|u| u.seq == c.seq) {
            u.source = "own+primer3";
        } else {
            uniq.push(c);
        }
    }
    uniq.sort_by(|a, b| b.cov0.cmp(&a.cov0).then(b.cov1.cmp(&a.cov1)).then(a.n_iupac.cmp(&b.n_iupac)).then(a.length.cmp(&b.length)));
    let mut o = b"{\"primer_tm\":".to_vec();
    go_float(&mut o, pt);
    o.extend_from_slice(b",\"tm_range\":[");
    go_float(&mut o, tm_range.0);
    o.push(b',');
    go_float(&mut o, tm_range.1);
    o.extend_from_slice(b"],\"window\":[");
    out_int(&mut o, lo + 1);
    o.push(b',');
    out_int(&mut o, hi + 1);
    o.extend_from_slice(b"],\"n_total\":");
    out_int(&mut o, n_total as i64);
    o.extend_from_slice(b",\"n_rejected\":");
    out_int(&mut o, rejected as i64);
    o.extend_from_slice(b",\"n_rows\":");
    out_int(&mut o, n_rows as i64);
    o.extend_from_slice(b",\"candidates\":[");
    for (i, c) in uniq.iter().enumerate() {
        if i > 0 {
            o.push(b',');
        }
        o.extend_from_slice(b"{\"seq\":");
        json::write_go_string(&mut o, &c.seq);
        o.extend_from_slice(if c.plus { b",\"strand\":\"+\"" } else { b",\"strand\":\"-\"" });
        for (k, v) in [("col_start", c.col_start), ("col_end", c.col_end), ("length", c.length), ("n_iupac", c.n_iupac as i64)] {
            o.extend_from_slice(format!(",\"{k}\":").as_bytes());
            out_int(&mut o, v);
        }
        o.extend_from_slice(b",\"tm\":");
        go_float(&mut o, c.tm);
        o.extend_from_slice(b",\"gc\":");
        go_float(&mut o, c.gc);
        for (k, v) in [("msa_cov0", c.cov0), ("msa_cov1", c.cov1), ("msa_n", c.n)] {
            o.extend_from_slice(format!(",\"{k}\":").as_bytes());
            out_int(&mut o, v as i64);
        }
        o.extend_from_slice(b",\"source\":");
        out_str(&mut o, c.source);
        o.push(b'}');
    }
    o.extend_from_slice(b"]}\n");
    Ok(o)
}
