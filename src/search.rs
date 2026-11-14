//! `mfeprimer search`: every oligo (degenerate ones expanded, named NAME.k) is anchored by its 5' and 3' k-mers on
//! either strand, and the span between them must be within the allowed edits (Levenshtein) of the query.

use crate::primers::{self, Oligo};
use crate::seqio::Record;
use crate::spec::{encode_all, kmer_mask, rev_comp, upper};
use crate::thermo::{self, Cond};
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::sync::atomic::{AtomicUsize, Ordering};

/// One hit of `mfeprimer search`.
pub struct SearchHit {
    pub rec: usize,
    /// 0-based, plus strand
    pub start: usize,
    /// exclusive
    pub end: usize,
    pub q: usize,
    /// '+': the query reads on the plus strand
    pub strand: u8,
    pub mismatch: usize,
    /// 0 for a gapped hit (span length differs from the query)
    pub tm: f64,
    pub dg: f64,
    /// the span in the query's orientation
    pub hit: Vec<u8>,
}

pub struct Query {
    pub name: String,
    pub seq: Vec<u8>,
}

/// The anchors of one query: 5' and 3' k-mers, and those of its reverse complement (an inosine gives one key per
/// base).
struct Anchors {
    a: Vec<u32>,
    b: Vec<u32>,
    ra: Vec<u32>,
    rb: Vec<u32>,
}

pub fn search(db: &[Record], oligos: &[Oligo], k: usize, max_mm: usize, cpu: usize, c: Cond) -> Result<(Vec<Query>, Vec<SearchHit>), String> {
    let mut qs = Vec::new();
    for o in oligos {
        let seqs = primers::expand(&o.seq).map_err(|e| format!("{}: {}", o.name, e))?;
        let multi = seqs.len() > 1;
        for (i, s) in seqs.into_iter().enumerate() {
            let name = if multi { format!("{}.{}", o.name, i + 1) } else { o.name.clone() };
            if s.len() >= k {
                qs.push(Query { name, seq: s.into_bytes() });
            }
        }
    }
    let an: Vec<Anchors> = qs
        .iter()
        .map(|q| {
            let n = q.seq.len();
            let rcq = rev_comp(&q.seq);
            // on the plus strand the reverse complement starts with rc(3' k-mer)
            Anchors {
                a: encode_all(&q.seq[..k]),
                b: encode_all(&q.seq[n - k..]),
                ra: encode_all(&rcq[..k]),
                rb: encode_all(&rcq[n - k..]),
            }
        })
        .collect();
    let mut want: HashMap<u32, usize> = HashMap::new(); // wanted k-mer -> slot in the per-genome position lists
    for x in &an {
        for ks in [&x.a, &x.b, &x.ra, &x.rb] {
            for &w in ks {
                let l = want.len();
                want.entry(w).or_insert(l);
            }
        }
    }
    let dense: Option<Vec<i32>> = (k <= 11).then(|| {
        let mut d = vec![-1i32; 1 << (2 * k)];
        for (&w, &slot) in &want {
            d[w as usize] = slot as i32;
        }
        d
    });
    let ctx = Ctx { qs: &qs, an: &an, want: &want, dense: dense.as_deref(), k, max_mm, c };
    let cpu = cpu.max(1);
    let next = AtomicUsize::new(0);
    let mut per: Vec<(usize, Vec<SearchHit>)> = std::thread::scope(|sc| {
        let hs: Vec<_> = (0..cpu)
            .map(|_| {
                sc.spawn(|| {
                    let mut mine = Vec::new();
                    loop {
                        let i = next.fetch_add(1, Ordering::Relaxed);
                        if i >= db.len() {
                            break;
                        }
                        mine.push((i, ctx.search_one(i, &db[i].seq)));
                    }
                    mine
                })
            })
            .collect();
        hs.into_iter().flat_map(|h| h.join().unwrap()).collect()
    });
    per.sort_by_key(|x| x.0);
    Ok((qs, per.into_iter().flat_map(|x| x.1).collect()))
}

struct Ctx<'a> {
    qs: &'a [Query],
    an: &'a [Anchors],
    want: &'a HashMap<u32, usize>,
    dense: Option<&'a [i32]>,
    k: usize,
    max_mm: usize,
    c: Cond,
}

impl Ctx<'_> {
    #[inline]
    fn slot(&self, w: u32) -> Option<usize> {
        match self.dense {
            Some(d) => {
                let s = d[w as usize];
                (s >= 0).then_some(s as usize)
            }
            None => self.want.get(&w).copied(),
        }
    }

    /// The start of every exact occurrence of each wanted k-mer (A/C/G/T only), by slot.
    fn kmer_positions(&self, g: &[u8]) -> Vec<Vec<usize>> {
        let k = self.k;
        let mut pos = vec![Vec::new(); self.want.len()];
        let mask = kmer_mask(k);
        let mut w = 0u32;
        let mut valid = 0;
        for (i, &b) in g.iter().enumerate() {
            let c = thermo::code(b);
            if c < 0 {
                valid = 0;
                w = 0;
                continue;
            }
            w = (w << 2 | c as u32) & mask;
            valid += 1;
            if valid >= k {
                if let Some(s) = self.slot(w) {
                    pos[s].push(i + 1 - k);
                }
            }
        }
        pos
    }

    /// The occurrences of several anchor keys merged (one key unless the anchor holds inosine).
    fn positions<'p>(&self, pos: &'p [Vec<usize>], keys: &[u32]) -> std::borrow::Cow<'p, [usize]> {
        if keys.len() == 1 {
            return std::borrow::Cow::Borrowed(&pos[self.want[&keys[0]]]);
        }
        let mut out: Vec<usize> = keys.iter().flat_map(|w| pos[self.want[w]].iter().copied()).collect();
        out.sort_unstable();
        std::borrow::Cow::Owned(out)
    }

    fn search_one(&self, rec: usize, g: &[u8]) -> Vec<SearchHit> {
        let k = self.k;
        let pos = self.kmer_positions(g);
        let up: Vec<u8> = g.iter().map(|&b| upper(b)).collect();
        let mut hits = Vec::new();
        for (qi, q) in self.qs.iter().enumerate() {
            let n = q.seq.len();
            let a = &self.an[qi];
            for strand in [b'+', b'-'] {
                let (lead, tail) = if strand == b'+' {
                    (self.positions(&pos, &a.a), self.positions(&pos, &a.b))
                } else {
                    (self.positions(&pos, &a.ra), self.positions(&pos, &a.rb))
                };
                for &s in lead.iter() {
                    for &t in tail.iter() {
                        let end = t + k;
                        if end < s {
                            continue;
                        }
                        let span = end - s;
                        if span + self.max_mm < n || span > n + self.max_mm || span < k {
                            continue;
                        }
                        let hit = if strand == b'-' { rev_comp(&up[s..end]) } else { up[s..end].to_vec() };
                        let mm = levenshtein(&q.seq, &hit);
                        if mm > self.max_mm {
                            continue;
                        }
                        let (mut tm, mut dg) = (0.0, 0.0);
                        if span == n {
                            if let Some((t, d)) = thermo::search(&q.seq, &hit, self.c) {
                                tm = t;
                                dg = d;
                            }
                        }
                        hits.push(SearchHit { rec, start: s, end, q: qi, strand, mismatch: mm, tm, dg, hit });
                    }
                }
            }
        }
        hits.sort_by_key(|h| h.start); // stable
        hits
    }
}

fn levenshtein(a: &[u8], b: &[u8]) -> usize {
    let mut prev: Vec<usize> = (0..=b.len()).collect();
    let mut cur = vec![0; b.len() + 1];
    for i in 1..=a.len() {
        cur[0] = i;
        for j in 1..=b.len() {
            let mut d = prev[j - 1];
            if a[i - 1] != b[j - 1] && !(a[i - 1] == b'I' && thermo::code(b[j - 1]) >= 0) {
                // inosine faces any base
                d += 1;
            }
            d = d.min(prev[j] + 1).min(cur[j - 1] + 1);
            cur[j] = d;
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    prev[b.len()]
}

/// Writes hits in `mfeprimer search` format.
pub fn write_search(path: &str, db: &[Record], qs: &[Query], hits: &[SearchHit], max_mm: usize) -> io::Result<()> {
    let mut w = BufWriter::with_capacity(1 << 20, File::create(path)?);
    write!(w, "# mfeprimer search\n# max allowed mismatch: {max_mm}\n\n")?;
    writeln!(w, "#Chrom\tStart(0-based)\tEnd\tQueryName\tSize\tStrand\tMisMatch\tDeltaG\tTm\tQuerySeq\tHitSeq")?;
    for h in hits {
        let q = &qs[h.q];
        write!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}\t",
            db[h.rec].id,
            h.start,
            h.end,
            q.name,
            h.end - h.start,
            h.strand as char,
            h.mismatch,
            h.dg,
            h.tm
        )?;
        w.write_all(&q.seq)?;
        w.write_all(b"\t")?;
        w.write_all(&h.hit)?;
        w.write_all(b"\n")?;
    }
    w.flush()
}
