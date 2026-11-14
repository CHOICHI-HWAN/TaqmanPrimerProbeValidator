//! Non-target sequences (other viruses, the host): a set is cross-reactive on a record when two of its primers bind
//! it (<= mm mismatches, 3' bases intact) facing each other within the amplicon limit; the probe binding inside that
//! product makes it a false signal. Searched in overlapping chunks, keeping the best hit of each oligo per chunk.

use super::cover::Params;
use super::oligo::iupac;
use super::qc::Set;
use crate::{reads, seqio};
use std::collections::{HashMap, HashSet};

pub struct Exclusion {
    recs: Vec<seqio::Record>,
    mm: usize,
    tp: usize,
    amp_max: usize,
}

struct XHit {
    role: u8,
    /// extends towards higher coordinates
    fwd: bool,
    start: usize,
    end: usize,
}

impl Exclusion {
    pub fn load(path: &str, p: &Params) -> Result<Exclusion, String> {
        let mut r = seqio::open(path)?;
        let mut recs = Vec::new();
        while let Some(mut rec) = r.next()? {
            rec.seq.make_ascii_uppercase();
            recs.push(rec);
        }
        if recs.is_empty() {
            return Err(format!("{path}: no sequences"));
        }
        Ok(Exclusion { recs, mm: p.pass_mm, tp: p.tp_n, amp_max: p.amp.1.max(400) })
    }

    /// Sets s.off_target, e.g. "2 products (1 with probe) on 2 records; oligo hits F3 R1 P0", and flags the set.
    pub fn check(&self, s: &mut Set) {
        let olig = vec![
            reads::Oligo { name: "F".into(), set: String::new(), role: "F".into(), seq: s.f.clone() },
            reads::Oligo { name: "R".into(), set: String::new(), role: "R".into(), seq: s.r.clone() },
            reads::Oligo { name: "P".into(), set: String::new(), role: "P".into(), seq: s.p.clone() },
        ];
        let seqs: Vec<(u8, Vec<u8>)> = olig.iter().map(|o| (o.role.as_bytes()[0], o.seq.clone())).collect();
        let srch = match reads::new_searcher(olig, self.mm, 0, false) {
            Ok(x) => x,
            Err(e) => {
                s.off_target = format!("error: {e}");
                return;
            }
        };
        let mut hits_by: HashMap<u8, usize> = HashMap::new();
        let (mut products, mut with_probe) = (0usize, 0usize);
        let mut recs_hit: HashSet<&str> = HashSet::new();
        let chunk = 4000usize;
        let step = chunk - self.amp_max - 60;
        for rec in &self.recs {
            let mut seen: HashSet<(u8, usize, bool)> = HashSet::new();
            let mut hs: Vec<XHit> = Vec::new();
            let mut off = 0;
            while off < rec.seq.len() {
                let end = rec.seq.len().min(off + chunk);
                for h in srch.scan(&rec.seq[off..end]) {
                    let (role, oseq) = (&seqs[h.oligo].0, &seqs[h.oligo].1);
                    if *role != b'P' {
                        let n = oseq.len();
                        if (n - self.tp.min(n)..n).any(|j| !iupac(oseq[j]).contains(&h.bases[j])) {
                            continue;
                        }
                    }
                    if seen.insert((*role, off + h.start, h.minus)) {
                        hs.push(XHit { role: *role, fwd: !h.minus, start: off + h.start, end: off + h.end });
                    }
                }
                if end == rec.seq.len() {
                    break;
                }
                off += step;
            }
            for h in &hs {
                *hits_by.entry(h.role).or_default() += 1;
            }
            hs.sort_by_key(|h| h.start);
            for a in hs.iter().filter(|a| a.role != b'P' && a.fwd) {
                for b in hs.iter().filter(|b| b.role != b'P' && !b.fwd) {
                    if b.end <= a.start || b.end - a.start > self.amp_max {
                        continue;
                    }
                    products += 1;
                    recs_hit.insert(&rec.id);
                    if hs.iter().any(|q| q.role == b'P' && q.start >= a.start && q.end <= b.end) {
                        with_probe += 1;
                    }
                }
            }
        }
        let hb = |r: u8| hits_by.get(&r).copied().unwrap_or(0);
        s.off_target = format!(
            "{products} products ({with_probe} with probe) on {} records; oligo hits F{} R{} P{}",
            recs_hit.len(),
            hb(b'F'),
            hb(b'R'),
            hb(b'P')
        );
        if with_probe > 0 {
            s.flags.push("CROSS-REACTIVE (non-target product with probe)".into());
        } else if products > 0 {
            s.flags.push("non-target product (no probe)".into());
        }
        s.penalty += 10.0 * with_probe as f64 + 2.0 * products as f64 + 0.2 * (hb(b'F') + hb(b'R') + hb(b'P')) as f64;
    }

    /// A non-target product with the probe inside (for Params::reject).
    pub fn rejects(&self, f: &[u8], p: &[u8], r: &[u8]) -> bool {
        let mut s = Set { f: f.to_vec(), p: p.to_vec(), r: r.to_vec(), ..Default::default() };
        self.check(&mut s);
        s.flags.iter().any(|x| x.starts_with("CROSS-REACTIVE"))
    }
}
