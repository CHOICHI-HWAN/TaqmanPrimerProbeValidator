//! Coverage-first design (port of the panels' scripts/redesign_coverage.py, same designs for the same input).

use super::align::{and3, and3_count, count, new_bits, Aligned, Bits, Kind};
use super::chem::Chem;
use super::oligo::{longest_run, n_ambig, rev_comp, stats};
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::sync::atomic::{AtomicUsize, Ordering as AO};
use std::sync::{Arc, Mutex};

pub type RejectFn = Arc<dyn Fn(&[u8], &[u8], &[u8]) -> bool + Send + Sync>;

/// One coverage-first run's parameters (design_from_alignment's arguments).
#[derive(Clone)]
pub struct Params {
    pub ambig: usize,
    pub probe_ambig: usize,
    pub pass_mm: usize,
    pub tp_n: usize,
    pub amp: (usize, usize),
    pub primer_tm: (f64, f64),
    pub gc: (f64, f64),
    pub probe_above: (f64, f64),
    pub max_dtm: f64,
    pub len_primer: (usize, usize),
    pub len_probe: (usize, usize),
    pub pool_bound: usize,
    pub pool_keep: usize,
    pub top: usize,
    pub probe: Chem,
    pub workers: usize,
    /// vetoes a triple in the greedy pick (e.g. a non-target product with the probe); the next best is taken
    pub reject: Option<RejectFn>,
    /// large genomes (> 20,000 windows of a primer length): only primer windows whose consensus GC / rough Tm can pass
    /// (within 5) are evaluated, as the MPXV / VARV editions' Python engine does
    pub large_prefilter: bool,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            ambig: 2,
            probe_ambig: 2,
            pass_mm: 2,
            tp_n: 2,
            amp: (70, 200),
            primer_tm: (52.0, 66.0),
            gc: (30.0, 75.0),
            probe_above: (2.0, 14.0),
            max_dtm: 5.0,
            len_primer: (18, 25),
            len_probe: (20, 30),
            pool_bound: 4000,
            pool_keep: 1200,
            top: 3,
            probe: Chem::default(),
            workers: crate::spec::num_cpu(),
            reject: None,
            large_prefilter: false,
        }
    }
}

/// One evaluated oligo window.
#[derive(Clone)]
pub struct Cand {
    pub cnt: usize,
    pub s: usize,
    pub e: usize,
    /// 5'->3' as ordered
    pub seq: Vec<u8>,
    pub strand: u8,
    pub tm: f64,
    pub gc: f64,
    pub bits: Arc<Bits>,
    pub amb: usize,
}

/// One chosen set (fields of the Python engine's designs).
#[derive(Clone)]
pub struct Design {
    pub name: String,
    pub f: Vec<u8>,
    pub p: Vec<u8>,
    pub r: Vec<u8>,
    pub f_tm: f64,
    pub p_tm: f64,
    pub r_tm: f64,
    pub f_gc: f64,
    pub p_gc: f64,
    pub r_gc: f64,
    pub amplicon_len: usize,
    pub est_cov: usize,
    pub est_cov_pct: f64,
    pub est_added: usize,
    pub p_strand: u8,
    pub fs: usize,
    pub re: usize,
    #[allow(dead_code)]
    pub bits: Bits,
}

/// Python's round(x, 1): the exact binary value rounded to one decimal, ties to even (56.25 -> 56.2).
pub fn round1(x: f64) -> f64 {
    format!("{x:.1}").parse().unwrap_or(x)
}

/// fn(i) for i in 0..n on `workers` threads, results in index order.
pub fn parallel<T: Send, F: Fn(usize) -> T + Sync>(n: usize, workers: usize, f: F) -> Vec<T> {
    let next = AtomicUsize::new(0);
    let out: Mutex<Vec<(usize, T)>> = Mutex::new(Vec::with_capacity(n));
    std::thread::scope(|sc| {
        for _ in 0..workers.max(1).min(n.max(1)) {
            sc.spawn(|| {
                let mut local = Vec::new();
                loop {
                    let i = next.fetch_add(1, AO::Relaxed);
                    if i >= n {
                        break;
                    }
                    local.push((i, f(i)));
                }
                out.lock().unwrap().extend(local);
            });
        }
    });
    let mut v = out.into_inner().unwrap();
    v.sort_by_key(|x| x.0);
    v.into_iter().map(|x| x.1).collect()
}

/// The strands a genome-sense oligo is ordered on.
pub fn orientations(kind: Kind, oligo: &[u8]) -> Vec<(u8, Vec<u8>)> {
    match kind {
        b'P' => vec![(b'+', oligo.to_vec()), (b'-', rev_comp(oligo))],
        b'R' => vec![(b'+', rev_comp(oligo))],
        _ => vec![(b'+', oligo.to_vec())],
    }
}

/// The per-oligo rules of the Python engine (and the probe format): (Tm, GC) or None.
pub fn screen(kind: Kind, seq: &[u8], p: &Params) -> Option<(f64, f64)> {
    let (res, tm, gc) = stats(seq)?;
    if kind == b'P' {
        if res.iter().any(|r| r[0] == b'G' || longest_run(r) >= 4) || !p.probe.length_ok(seq.len()) {
            return None;
        }
        return Some((tm + p.probe.offset, gc));
    }
    if tm < p.primer_tm.0 || tm > p.primer_tm.1 || gc < p.gc.0 || gc > p.gc.1 || res.iter().any(|r| longest_run(r) >= 5) {
        return None;
    }
    Some((tm, gc))
}

pub fn sort_cands(v: &mut [Cand]) {
    v.sort_by(|a, b| b.cnt.cmp(&a.cnt).then(a.amb.cmp(&b.amb)).then(a.s.cmp(&b.s)));
}

/// The most promising windows of one kind, evaluated (the Python engine's candidates()).
pub fn candidates(al: &Aligned, kind: Kind, p: &Params) -> Vec<Cand> {
    let (amb, lr, thr) = if kind == b'P' {
        (p.probe_ambig, p.len_probe, p.probe_ambig)
    } else {
        (p.ambig, p.len_primer, p.pass_mm + p.ambig)
    };
    let mut pool: Vec<(usize, usize, usize)> = Vec::new(); // ub, s, k
    for k in lr.0..=lr.1 {
        let mut starts = al.windows(k);
        if p.large_prefilter && kind != b'P' && starts.len() > 20000 {
            // numpy: tm_est = 64.9 + 41.0 * (n_gc - 16.4) / k; gc = 100.0 * n_gc / k
            starts.retain(|&s| {
                let n_gc = al.cons[s..s + k].iter().filter(|&&c| c == b'G' || c == b'C').count() as f64;
                let kf = k as f64;
                let tm_est = 64.9 + 41.0 * (n_gc - 16.4) / kf;
                let gcp = 100.0 * n_gc / kf;
                tm_est >= p.primer_tm.0 - 5.0 && tm_est <= p.primer_tm.1 + 5.0 && gcp >= p.gc.0 - 5.0 && gcp <= p.gc.1 + 5.0
            });
        }
        for s in starts {
            pool.push((al.upper_bound(s, k, thr), s, k));
        }
    }
    pool.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
    pool.truncate(p.pool_bound);
    let per = parallel(pool.len(), p.workers, |i| {
        let (_, s, k) = pool[i];
        let codes = al.eval_codes(s, k, kind, amb, p.pass_mm, p.tp_n);
        let oligo = al.oligo_with(s, k, &codes);
        let bits = Arc::new(al.pass_bits(&oligo, s, kind, p.pass_mm, p.tp_n));
        let cnt = count(&bits);
        let mut v = Vec::new();
        for (strand, seq) in orientations(kind, &oligo) {
            if let Some((tm, gc)) = screen(kind, &seq, p) {
                let amb = n_ambig(&seq);
                v.push(Cand { cnt, s, e: s + k, seq, strand, tm, gc, bits: bits.clone(), amb });
            }
        }
        v
    });
    let mut out: Vec<Cand> = per.into_iter().flatten().collect();
    sort_cands(&mut out);
    out.truncate(p.pool_keep);
    out
}

/// A ranked F/R/P combination (indices into the candidate lists).
#[derive(Clone)]
pub struct Triple {
    pub new: usize,
    pub total: usize,
    pub f: usize,
    pub r: usize,
    pub p: usize,
    key: [f64; 5],
}

fn cmp_key(a: &[f64; 5], b: &[f64; 5]) -> Ordering {
    for i in 0..5 {
        match a[i].partial_cmp(&b[i]).unwrap_or(Ordering::Equal) {
            Ordering::Equal => continue,
            o => return o,
        }
    }
    Ordering::Equal
}

struct MinItem(Triple);
impl PartialEq for MinItem {
    fn eq(&self, o: &Self) -> bool {
        cmp_key(&self.0.key, &o.0.key) == Ordering::Equal
    }
}
impl Eq for MinItem {}
impl PartialOrd for MinItem {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for MinItem {
    // reversed: BinaryHeap pops the smallest key first
    fn cmp(&self, o: &Self) -> Ordering {
        cmp_key(&o.0.key, &self.0.key)
    }
}

/// Combinations ranked by the genomes they pass that `covered` does not (then all they pass, fewer IUPAC codes,
/// closer primer Tm, found first); the Python engine's best_triples() with the same pruning.
pub fn best_triples(f: &[Cand], r: &[Cand], pc: &[Cand], covered: &[u64], p: &Params, want: usize) -> Vec<Triple> {
    let mut h: BinaryHeap<MinItem> = BinaryHeap::new();
    let (mut floor, mut seq) = (0usize, 0usize);
    let mut fp = vec![0u64; covered.len()];
    for (pi, pr) in pc.iter().enumerate() {
        if pr.cnt <= floor {
            break;
        }
        for (fi, fc) in f.iter().enumerate() {
            if !(fc.e <= pr.s && pr.e - fc.s <= p.amp.1) {
                continue;
            }
            if fc.cnt <= floor {
                break;
            }
            let mut nfp = 0;
            for i in 0..fp.len() {
                fp[i] = fc.bits[i] & pr.bits[i] & !covered[i];
                nfp += fp[i].count_ones() as usize;
            }
            if nfp <= floor {
                continue;
            }
            for (ri, rc) in r.iter().enumerate() {
                if !(rc.s >= pr.e && rc.e - pr.s <= p.amp.1) {
                    continue;
                }
                if rc.cnt <= floor {
                    break;
                }
                let amp = rc.e - fc.s;
                if !(p.amp.0 <= amp && amp <= p.amp.1) || (fc.tm - rc.tm).abs() > p.max_dtm {
                    continue;
                }
                let hi = fc.tm.max(rc.tm);
                if !(hi + p.probe_above.0 <= pr.tm && pr.tm <= hi + p.probe_above.1) {
                    continue;
                }
                let mut nw = 0;
                for i in 0..fp.len() {
                    nw += (fp[i] & rc.bits[i]).count_ones() as usize;
                }
                if nw <= floor {
                    continue;
                }
                let tot = and3_count(&fc.bits, &pr.bits, &rc.bits, None);
                seq += 1;
                let key = [nw as f64, tot as f64, -((fc.amb + rc.amb + pr.amb) as f64), -(fc.tm - rc.tm).abs(), -(seq as f64)];
                let t = Triple { new: nw, total: tot, f: fi, r: ri, p: pi, key };
                if h.len() < want {
                    h.push(MinItem(t));
                } else if cmp_key(&h.peek().unwrap().0.key, &t.key) == Ordering::Less {
                    h.pop();
                    h.push(MinItem(t));
                }
                if h.len() >= want {
                    let m = h.peek().unwrap().0.key[0] as usize;
                    floor = if m > 0 { m - 1 } else { 0 };
                }
            }
        }
    }
    let mut out: Vec<Triple> = h.into_iter().map(|x| x.0).collect();
    out.sort_by(|a, b| cmp_key(&b.key, &a.key));
    out
}

pub struct Run {
    pub n: usize,
    pub status: String,
    pub designs: Vec<Design>,
    pub message: String,
    pub nf: usize,
    pub nr: usize,
    pub np: usize,
}

/// The coverage-first greedy panel: each round takes the best triple for the genomes the earlier rounds miss (a new
/// triple must not repeat an earlier amplicon within 6 bases at both ends).
pub fn design_from_alignment(al: &Aligned, p: &Params) -> Run {
    let f = candidates(al, b'F', p);
    let r = candidates(al, b'R', p);
    let pc = candidates(al, b'P', p);
    let mut run = Run { n: al.g, status: "ok".into(), designs: Vec::new(), message: String::new(), nf: f.len(), nr: r.len(), np: pc.len() };
    if f.is_empty() || r.is_empty() || pc.is_empty() {
        run.status = "none".into();
        run.message = "not enough F/R/probe candidates for the Tm, GC and IUPAC limits".into();
        return run;
    }
    let mut covered = new_bits(al.g);
    for i in 0..p.top {
        let cand = best_triples(&f, &r, &pc, &covered, p, 12);
        let pick = cand.iter().find(|c| {
            let dup = run.designs.iter().any(|d| f[c.f].s.abs_diff(d.fs) < 6 && r[c.r].e.abs_diff(d.re) < 6);
            !dup && p.reject.as_ref().map_or(true, |rej| !rej(&f[c.f].seq, &pc[c.p].seq, &r[c.r].seq))
        });
        let Some(t) = pick else { break };
        let (fc, rc, pr) = (&f[t.f], &r[t.r], &pc[t.p]);
        let b = and3(&fc.bits, &rc.bits, &pr.bits);
        for w in 0..covered.len() {
            covered[w] |= b[w];
        }
        run.designs.push(Design {
            name: format!("COV_{i}"),
            f: fc.seq.clone(),
            p: pr.seq.clone(),
            r: rc.seq.clone(),
            f_tm: round1(fc.tm),
            p_tm: round1(pr.tm),
            r_tm: round1(rc.tm),
            f_gc: round1(fc.gc),
            p_gc: round1(pr.gc),
            r_gc: round1(rc.gc),
            amplicon_len: al.amplicon_len(fc.s, rc.e),
            est_cov: t.total,
            est_cov_pct: round1(100.0 * t.total as f64 / al.g as f64),
            est_added: t.new,
            p_strand: pr.strand,
            fs: fc.s,
            re: rc.e,
            bits: b,
        });
    }
    if run.designs.is_empty() {
        run.status = "none".into();
        run.message = "no F/R/probe combination fits (check amplicon length and Tm limits, or every candidate amplifies a non-target)".into();
    }
    run
}
