//! Exhaustive (brute-force) searches: every shifted / resized / degenerate variant of a set, and the best set in
//! every amplicon-sized window along the genome.

use super::align::{and3, and3_count, count, new_bits, Aligned, Kind};
use super::cover::{best_triples, orientations, parallel, screen, sort_cands, Cand, Params};
use super::oligo::{n_ambig, rev_comp, stats};
use super::qc::Set;
use std::collections::{HashMap, HashSet};
use std::sync::Arc;

/// For every conserved window start in starts and every length in lens: the consensus with 0..amb greedy IUPAC codes
/// (each level one variant) that pass the per-oligo rules. tm_ok (shift) replaces the primer Tm window and GC limits.
fn variants(al: &Aligned, kind: Kind, starts: &[isize], lens: (usize, usize), amb: usize, p: &Params, tm_ok: Option<&(dyn Fn(f64) -> bool + Sync)>) -> Vec<Cand> {
    let mut jobs = Vec::new();
    for &s in starts {
        for k in lens.0..=lens.1 {
            if al.conserved(s, k) {
                jobs.push((s as usize, k));
            }
        }
    }
    let mut q = p.clone();
    q.primer_tm = (-1e9, 1e9);
    q.gc = (0.0, 100.0);
    let pp = if tm_ok.is_some() { &q } else { p };
    let per = parallel(jobs.len(), p.workers, |i| {
        let (s, k) = jobs[i];
        let codes = al.eval_codes(s, k, kind, amb, p.pass_mm, p.tp_n);
        let mut seen: HashSet<Vec<u8>> = HashSet::new();
        let mut v = Vec::new();
        for lv in 0..=codes.len() {
            let oligo = al.oligo_with(s, k, &codes[..lv]);
            if !seen.insert(oligo.clone()) {
                continue;
            }
            let bits = Arc::new(al.pass_bits(&oligo, s, kind, p.pass_mm, p.tp_n));
            let cnt = count(&bits);
            for (strand, seq) in orientations(kind, &oligo) {
                let Some((tm, gc)) = screen(kind, &seq, pp) else { continue };
                if kind != b'P' && tm_ok.is_some_and(|f| !f(tm)) {
                    continue;
                }
                let amb = n_ambig(&seq);
                v.push(Cand { cnt, s, e: s + k, seq, strand, tm, gc, bits: bits.clone(), amb });
            }
        }
        v
    });
    per.into_iter().flatten().collect()
}

fn range_step(lo: isize, hi: isize, step: usize) -> Vec<isize> {
    let mut out = Vec::new();
    let mut s = lo;
    while s <= hi {
        out.push(s);
        s += step.max(1) as isize;
    }
    out
}

/// A set to optimise, as ordered (F and P in any orientation; R as ordered).
pub struct Seed {
    pub name: String,
    pub f: Vec<u8>,
    pub p: Vec<u8>,
    pub r: Vec<u8>,
}

pub struct ShiftOptions {
    pub shift: usize,
    pub step: usize,
    pub keep: usize,
    pub primer_amb: usize,
    pub probe_amb: usize,
    /// a primer variant's Tm stays within this of the original
    pub dtm: f64,
}

#[derive(Clone, Default)]
struct Located {
    s: usize,
    seq: Vec<u8>,
    n: usize,
}

fn locate_either(al: &Aligned, o: &[u8]) -> (Located, Located) {
    let rc = rev_comp(o);
    let plus = al.locate(o).map(|(s, n)| Located { s, seq: o.to_vec(), n }).unwrap_or_default();
    let minus = al.locate(&rc).map(|(s, n)| Located { s, seq: rc.clone(), n }).unwrap_or_default();
    (plus, minus)
}

/// Every variant of the set's F, R and probe: each start moved up to opt.shift bases either way in opt.step steps,
/// every length of the length ranges, 0..amb greedy IUPAC codes; the opt.keep best of each are combined in full.
/// Returns the best set for each change class (P, F, R, FR, FP, RP, FRP) that beats the original, and the original.
pub fn shift(al: &Aligned, seed: &Seed, p: &Params, opt: &ShiftOptions) -> Result<(Vec<Set>, Set), String> {
    let (fp, fm) = locate_either(al, &seed.f);
    let (rp, rm) = locate_either(al, &seed.r);
    // one primer sits in genome sense (forward), the other on the minus strand; manifests sometimes swap the
    // columns, so both assignments are tried and the one more genomes match wins
    let s1 = if fp.n > 0 && rm.n > 0 { fp.n + rm.n } else { 0 };
    let s2 = if rp.n > 0 && fm.n > 0 { rp.n + fm.n } else { 0 };
    let (fl, rl, sf, sr) = if s1 == 0 && s2 == 0 {
        return Err(format!("{}: primers not located in the alignment", seed.name));
    } else if s1 >= s2 {
        (fp, rm, seed.f.clone(), seed.r.clone())
    } else {
        (rp, fm, seed.r.clone(), seed.f.clone())
    };
    if rl.s < fl.s {
        return Err(format!("{}: primers located but not facing each other", seed.name));
    }
    let (pp, pm) = locate_either(al, &seed.p);
    let pl = if pm.n > pp.n { pm } else { pp };
    if pl.n == 0 {
        return Err(format!("{}: probe not located in the alignment", seed.name));
    }
    let orig = |kind: Kind, l: &Located, ordered: &[u8]| -> Cand {
        let b = al.pass_bits(&l.seq, l.s, kind, p.pass_mm, p.tp_n);
        let (tm, gc) = stats(ordered).map(|x| (x.1, x.2)).unwrap_or((0.0, 0.0));
        let tm = if kind == b'P' { tm + p.probe.offset } else { tm };
        Cand { cnt: count(&b), s: l.s, e: l.s + l.seq.len(), seq: ordered.to_vec(), strand: b'+', tm, gc, bits: Arc::new(b), amb: n_ambig(ordered) }
    };
    let (f0, r0, p0) = (orig(b'F', &fl, &sf), orig(b'R', &rl, &sr), orig(b'P', &pl, &seed.p));
    let base = Set {
        name: seed.name.clone(),
        f: sf.clone(),
        p: seed.p.clone(),
        r: sr.clone(),
        cov: and3_count(&f0.bits, &r0.bits, &p0.bits, None),
        fs: f0.s,
        re: r0.e,
        amplicon: al.amplicon_len(f0.s, r0.e),
        ..Default::default()
    };
    let (sh, st) = (opt.shift as isize, opt.step);
    let (ft, rt, dtm) = (f0.tm, r0.tm, opt.dtm);
    let near_f = move |t: f64| (t - ft).abs() <= dtm;
    let near_r = move |t: f64| (t - rt).abs() <= dtm;
    let mut fs = variants(al, b'F', &range_step(f0.s as isize - sh, f0.s as isize + sh, st), p.len_primer, opt.primer_amb, p, Some(&near_f));
    let mut rs = variants(al, b'R', &range_step(r0.s as isize - sh, r0.s as isize + sh, st), p.len_primer, opt.primer_amb, p, Some(&near_r));
    let mut ps = variants(al, b'P', &range_step(f0.s as isize - sh, r0.e as isize + sh, st), p.len_probe, opt.probe_amb, p, None);
    let prune = |cs: &mut Vec<Cand>, o: &Cand| {
        cs.sort_by(|a, b| b.cnt.cmp(&a.cnt).then(a.amb.cmp(&b.amb)).then(a.s.abs_diff(o.s).cmp(&b.s.abs_diff(o.s))));
        cs.truncate(opt.keep);
        cs.push(o.clone());
    };
    prune(&mut fs, &f0);
    prune(&mut rs, &r0);
    prune(&mut ps, &p0);
    let same = |a: &Cand, b: &Cand| a.s == b.s && a.e == b.e && a.seq == b.seq;
    let mut classes: HashMap<String, ([f64; 3], usize, usize, usize, usize)> = HashMap::new();
    let mut fr = new_bits(al.g);
    for (fi, f) in fs.iter().enumerate() {
        for (ri, r) in rs.iter().enumerate() {
            let amp = r.e as isize - f.s as isize;
            if amp < p.amp.0 as isize || amp > p.amp.1 as isize || r.s < f.e {
                continue;
            }
            let (fc, rc) = (!same(f, &f0), !same(r, &r0));
            let mut nfr = 0;
            for w in 0..fr.len() {
                fr[w] = f.bits[w] & r.bits[w];
                nfr += fr[w].count_ones() as usize;
            }
            if nfr <= base.cov {
                continue;
            }
            let hi = f.tm.max(r.tm);
            for (qi, q) in ps.iter().enumerate() {
                if q.cnt <= base.cov || q.s < f.e || q.e > r.s {
                    continue;
                }
                let qc = !same(q, &p0);
                if !fc && !rc && !qc {
                    continue;
                }
                // as optimize_user_sets.py: each primer within dtm of its original, the probe at least
                // probe_above.0 above the higher primer (no upper limit)
                if q.tm < hi + p.probe_above.0 {
                    continue;
                }
                let mut cnt = 0;
                for w in 0..fr.len() {
                    cnt += (fr[w] & q.bits[w]).count_ones() as usize;
                }
                if cnt <= base.cov {
                    continue;
                }
                let cls = format!("{}{}{}", if fc { "F" } else { "" }, if rc { "R" } else { "" }, if qc { "P" } else { "" });
                let dist = f.s.abs_diff(f0.s) + r.e.abs_diff(r0.e) + q.s.abs_diff(p0.s);
                let key = [cnt as f64, -((f.amb + r.amb + q.amb) as f64), -(dist as f64)];
                let better = match classes.get(&cls) {
                    None => true,
                    Some(b) => key.partial_cmp(&b.0) == Some(std::cmp::Ordering::Greater),
                };
                if better {
                    classes.insert(cls, (key, fi, ri, qi, cnt));
                }
            }
        }
    }
    let mut out = Vec::new();
    for cls in ["P", "F", "R", "FR", "FP", "RP", "FRP"] {
        let Some(&(_, fi, ri, qi, cnt)) = classes.get(cls) else { continue };
        let (f, r, q) = (&fs[fi], &rs[ri], &ps[qi]);
        out.push(Set {
            name: seed.name.clone(),
            change: cls.into(),
            base_cov: base.cov,
            cov: cnt,
            added: cnt - base.cov,
            f: f.seq.clone(),
            p: q.seq.clone(),
            r: r.seq.clone(),
            fs: f.s,
            re: r.e,
            amplicon: al.amplicon_len(f.s, r.e),
            bits: and3(&f.bits, &r.bits, &q.bits),
            ..Default::default()
        });
    }
    Ok((out, base))
}

/// The best set in one amplicon-sized window.
pub struct TileRow {
    pub start: usize,
    pub set: Option<Set>,
}

/// Evaluates EVERY conserved window of every length (0..amb IUPAC codes each), then for each region [x, x+amp_max)
/// with x stepping by `step` finds the best triple inside it: a coverage map along the genome.
pub fn tile(al: &Aligned, p: &Params, step: usize, keep: usize) -> Vec<TileRow> {
    let all = |kind: Kind, lens: (usize, usize), amb: usize| {
        let hi = al.l as isize - lens.0 as isize;
        let mut cs = variants(al, kind, &range_step(0, hi, 1), lens, amb, p, None);
        sort_cands(&mut cs);
        cs
    };
    let f = all(b'F', p.len_primer, p.ambig);
    let r = all(b'R', p.len_primer, p.ambig);
    let pc = all(b'P', p.len_probe, p.probe_ambig);
    let starts: Vec<usize> = (0..al.l).step_by(step.max(1)).collect();
    let empty = new_bits(al.g);
    parallel(starts.len(), p.workers, |i| {
        let x = starts[i];
        let lim = x + p.amp.1;
        let pick = |cs: &[Cand]| -> Vec<Cand> { cs.iter().filter(|c| c.s >= x && c.e <= lim).take(keep).cloned().collect() };
        let (fl, rl, ql) = (pick(&f), pick(&r), pick(&pc));
        if fl.is_empty() || rl.is_empty() || ql.is_empty() {
            return TileRow { start: x, set: None };
        }
        let t = best_triples(&fl, &rl, &ql, &empty, p, 1);
        let Some(b) = t.first() else { return TileRow { start: x, set: None } };
        let (fc, rc, qc) = (&fl[b.f], &rl[b.r], &ql[b.p]);
        TileRow {
            start: x,
            set: Some(Set {
                name: format!("TILE_{x}"),
                f: fc.seq.clone(),
                p: qc.seq.clone(),
                r: rc.seq.clone(),
                cov: b.total,
                fs: fc.s,
                re: rc.e,
                amplicon: al.amplicon_len(fc.s, rc.e),
                bits: and3(&fc.bits, &rc.bits, &qc.bits),
                ..Default::default()
            }),
        }
    })
}

/// From tile results, sets that together cover the most genomes (each next one adds the most).
pub fn greedy_panel(rows: &[TileRow], top: usize, n: usize) -> Vec<Set> {
    let mut covered = new_bits(n);
    let mut out: Vec<Set> = Vec::new();
    let mut used: HashSet<(usize, usize)> = HashSet::new();
    while out.len() < top {
        let (mut bi, mut bn): (Option<usize>, usize) = (None, 0);
        for (i, r) in rows.iter().enumerate() {
            let Some(s) = &r.set else { continue };
            if used.contains(&(s.fs, s.re)) {
                continue;
            }
            let add: usize = (0..covered.len()).map(|w| (s.bits[w] & !covered[w]).count_ones() as usize).sum();
            let better = add > bn || (add == bn && add > 0 && bi.is_some_and(|b| s.cov > rows[b].set.as_ref().unwrap().cov));
            if better {
                bi = Some(i);
                bn = add;
            }
        }
        let Some(b) = bi else { break };
        if bn == 0 {
            break;
        }
        let mut s = rows[b].set.clone().unwrap();
        s.added = bn;
        s.name = format!("TILE_{}", out.len());
        used.insert((s.fs, s.re));
        for w in 0..covered.len() {
            covered[w] |= s.bits[w];
        }
        out.push(s);
    }
    out
}
