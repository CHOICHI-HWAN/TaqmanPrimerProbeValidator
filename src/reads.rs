//! Checks sequencing reads (FASTQ or FASTA, single-end, paired-end or long reads) for the oligos of primer sets: which
//! reads carry each primer / probe, with how many mismatches, and which bases the reads have at every oligo position.
//!
//! Search: every oligo (and its reverse complement) is cut into disjoint seeds of length s; with at most m mismatches
//! (or m edits) at least one of m+1 disjoint pieces is untouched, so it appears exactly in the read. Each seed hit is
//! verified against the whole oligo: Hamming distance by default, or a semi-global edit distance (indels allowed, for
//! error-prone long reads) with indel. IUPAC codes and inosine (I) in the oligo match every base they stand for; a
//! read N is a mismatch. An oligo must lie wholly inside the read.
//!
//! Unlike the Go engine, hits are taken in oligo order, so when both mates of a pair carry an oligo with the same
//! number of mismatches the first mate's hit is always the one counted, and .hits.tsv lines come in input order.

use crate::primers;
use crate::seqio;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fmt::Write as _;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::mpsc;

/// The read bases each oligo letter accepts (A=1 C=2 G=4 T=8).
static MASK: [u8; 256] = {
    let mut t = [0u8; 256];
    let m: [(u8, u8); 17] = [
        (b'A', 1), (b'C', 2), (b'G', 4), (b'T', 8), (b'U', 8), (b'R', 1 | 4), (b'Y', 2 | 8), (b'S', 2 | 4),
        (b'W', 1 | 8), (b'K', 4 | 8), (b'M', 1 | 2), (b'B', 2 | 4 | 8), (b'D', 1 | 4 | 8), (b'H', 1 | 2 | 8),
        (b'V', 1 | 2 | 4), (b'N', 15), (b'I', 15),
    ];
    let mut i = 0;
    while i < m.len() {
        t[m[i].0 as usize] = m[i].1;
        t[(m[i].0 + 32) as usize] = m[i].1;
        i += 1;
    }
    t
};

static READ_BIT: [u8; 256] = {
    let mut t = [0u8; 256];
    let b = *b"ACGT";
    let mut i = 0;
    while i < 4 {
        t[b[i] as usize] = 1 << i;
        t[(b[i] + 32) as usize] = 1 << i;
        i += 1;
    }
    t
};

static CODE2: [i8; 256] = {
    let mut t = [-1i8; 256];
    let b = *b"ACGT";
    let mut i = 0;
    while i < 4 {
        t[b[i] as usize] = i as i8;
        t[(b[i] + 32) as usize] = i as i8;
        i += 1;
    }
    t
};

static COMPL: [u8; 256] = {
    let mut t = [b'N'; 256];
    let p: [&[u8; 2]; 17] = [
        b"AT", b"CG", b"GC", b"TA", b"RY", b"YR", b"KM", b"MK", b"SS", b"WW", b"BV", b"VB", b"DH", b"HD", b"NN", b"II",
        b"UA",
    ];
    let mut i = 0;
    while i < p.len() {
        t[p[i][0] as usize] = p[i][1];
        t[(p[i][0] + 32) as usize] = p[i][1];
        i += 1;
    }
    t
};

fn rev_comp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| COMPL[b as usize]).collect()
}

/// One primer or probe; set and role ("F", "R", "P" or "") come from a NAME_F / NAME_R / NAME_P name.
pub struct Oligo {
    pub name: String,
    pub set: String,
    pub role: String,
    pub seq: Vec<u8>,
}

/// Turns input oligos into Oligos, splitting SET_F / SET_R / SET_P names.
pub fn from_primers(input: &[primers::Oligo]) -> Result<Vec<Oligo>, String> {
    let mut out = Vec::new();
    for o in input {
        let seq = o.seq.to_uppercase().into_bytes();
        if let Some(&b) = seq.iter().find(|&&b| MASK[b as usize] == 0) {
            return Err(format!("{}: unsupported base '{}'", o.name, b as char));
        }
        let mut x = Oligo { name: o.name.clone(), set: String::new(), role: String::new(), seq };
        if let Some(i) = o.name.rfind('_').filter(|&i| i > 0) {
            let r = o.name[i + 1..].to_uppercase();
            if matches!(r.as_str(), "F" | "R" | "P") {
                x.set = o.name[..i].to_string();
                x.role = r;
            }
        }
        out.push(x);
    }
    Ok(out)
}

struct Pattern {
    oligo: usize,
    /// the oligo's reverse complement: the oligo reads on the read's minus strand
    minus: bool,
    /// as searched on the read
    seq: Vec<u8>,
}

/// The seed index of a set of oligos.
pub struct Searcher {
    pub oligos: Vec<Oligo>,
    pub max_mm: usize,
    pub indel: bool,
    seed: usize,
    pats: Vec<Pattern>,
    index: HashMap<u32, Vec<(u32, u32)>>,
}

const MAX_SEED_EXPANSIONS: usize = 1 << 12;

/// Indexes the oligos. seed 0 picks the longest seed that still gives every oligo m+1 disjoint pieces.
pub fn new_searcher(oligos: Vec<Oligo>, max_mm: usize, seed: usize, indel: bool) -> Result<Searcher, String> {
    if oligos.is_empty() {
        return Err("no oligos".into());
    }
    let mut seed = seed;
    if seed == 0 {
        seed = 12;
        for o in &oligos {
            seed = seed.min(o.seq.len() / (max_mm + 1));
        }
    }
    if seed < 4 {
        return Err(format!("an oligo is too short for {max_mm} mismatches (seed would be {seed} bases; need >= 4)"));
    }
    seed = seed.min(15);
    let mut s = Searcher { oligos, max_mm, indel, seed, pats: Vec::new(), index: HashMap::new() };
    for oi in 0..s.oligos.len() {
        for minus in [false, true] {
            let seq = if minus { rev_comp(&s.oligos[oi].seq) } else { s.oligos[oi].seq.clone() };
            let pi = s.pats.len() as u32;
            let mut indexed = 0;
            let mut off = 0;
            while off + seed <= seq.len() {
                let keys = expand_seed(&seq[off..off + seed]);
                for &k in &keys {
                    s.index.entry(k).or_default().push((pi, off as u32));
                }
                if !keys.is_empty() {
                    indexed += 1;
                }
                off += seed;
            }
            if indexed < max_mm + 1 {
                return Err(format!(
                    "{}: too degenerate to search with {} mismatches (only {} usable seeds)",
                    s.oligos[oi].name, max_mm, indexed
                ));
            }
            s.pats.push(Pattern { oligo: oi, minus, seq });
        }
    }
    Ok(s)
}

/// The 2-bit codes of every plain sequence a seed stands for (none if there are too many).
fn expand_seed(seed: &[u8]) -> Vec<u32> {
    let mut keys = vec![0u32];
    for &b in seed {
        let m = MASK[b as usize];
        let mut next = Vec::new();
        for &k in &keys {
            for c in 0..4u32 {
                if m & (1 << c) != 0 {
                    next.push(k << 2 | c);
                }
            }
        }
        if next.len() > MAX_SEED_EXPANSIONS {
            return Vec::new();
        }
        keys = next;
    }
    keys
}

/// The best occurrence of one oligo in one read.
#[derive(Clone)]
pub struct Hit {
    pub oligo: usize,
    pub minus: bool,
    /// 0-based on the read
    pub start: usize,
    /// exclusive
    pub end: usize,
    /// substitutions (Hamming) or edits
    pub mismatch: usize,
    /// the read base facing each oligo position, 5'->3' in the oligo's orientation ('-' = deleted)
    pub bases: Vec<u8>,
}

impl Searcher {
    /// The best hit of each oligo found in seq (at most one per oligo; fewer mismatches win), in oligo order.
    pub fn scan(&self, seq: &[u8]) -> Vec<Hit> {
        let mut best: BTreeMap<usize, Hit> = BTreeMap::new();
        let mut tried: HashSet<(u32, i64)> = HashSet::new();
        let k = self.seed;
        let kmask = ((1u64 << (2 * k)) - 1) as u32;
        let mut w = 0u32;
        let mut valid = 0;
        for i in 0..seq.len() {
            let c = CODE2[seq[i] as usize];
            if c < 0 {
                valid = 0;
                w = 0;
                continue;
            }
            w = (w << 2 | c as u32) & kmask;
            valid += 1;
            if valid < k {
                continue;
            }
            let kstart = (i + 1 - k) as i64;
            let Some(refs) = self.index.get(&w) else { continue };
            for &(pat, off) in refs {
                let p = &self.pats[pat as usize];
                let cand = kstart - off as i64;
                if !tried.insert((pat, cand)) {
                    continue;
                }
                let h = if self.indel { self.verify_edit(p, seq, cand) } else { self.verify_hamming(p, seq, cand) };
                let Some(h) = h else { continue };
                match best.get(&p.oligo) {
                    Some(old) if h.mismatch >= old.mismatch => {}
                    _ => {
                        best.insert(p.oligo, h);
                    }
                }
            }
        }
        best.into_values().collect()
    }

    fn verify_hamming(&self, p: &Pattern, seq: &[u8], start: i64) -> Option<Hit> {
        let len = p.seq.len();
        if start < 0 || start as usize + len > seq.len() {
            return None;
        }
        let start = start as usize;
        let mut mm = 0;
        for j in 0..len {
            if MASK[p.seq[j] as usize] & READ_BIT[seq[start + j] as usize] == 0 {
                mm += 1;
                if mm > self.max_mm {
                    return None;
                }
            }
        }
        let b = &seq[start..start + len];
        let bases = if p.minus { rev_comp(b) } else { b.to_vec() };
        Some(Hit { oligo: p.oligo, minus: p.minus, start, end: start + len, mismatch: mm, bases })
    }

    /// Aligns the whole oligo against the read around start (free read ends, edits <= max_mm).
    fn verify_edit(&self, p: &Pattern, seq: &[u8], start: i64) -> Option<Hit> {
        let m = self.max_mm as i64;
        let len = p.seq.len() as i64;
        let lo = (start - m).max(0);
        let hi = (start + len + m).min(seq.len() as i64);
        if hi - lo < len - m {
            return None;
        }
        let (lo, hi) = (lo as usize, hi as usize);
        let t = &seq[lo..hi];
        let (ll, n) = (p.seq.len(), t.len());
        // d[i][j]: edits aligning oligo[:i] with a read stretch ending at t[j-1]; the read start is free (d[0][j] = 0)
        let w = n + 1;
        let mut d = vec![0usize; (ll + 1) * w];
        for i in 0..=ll {
            d[i * w] = i;
        }
        for i in 1..=ll {
            let ob = MASK[p.seq[i - 1] as usize];
            for j in 1..=n {
                let mut sub = d[(i - 1) * w + j - 1];
                if ob & READ_BIT[t[j - 1] as usize] == 0 {
                    sub += 1;
                }
                let mut v = sub;
                v = v.min(d[(i - 1) * w + j] + 1); // oligo base with no read base (deletion in the read)
                v = v.min(d[i * w + j - 1] + 1); // extra read base (insertion in the read)
                d[i * w + j] = v;
            }
        }
        let mut bj: Option<usize> = None;
        for j in 0..=n {
            if d[ll * w + j] <= self.max_mm && bj.is_none_or(|b| d[ll * w + j] < d[ll * w + b]) {
                bj = Some(j);
            }
        }
        let bj = bj?;
        let mut bases = vec![0u8; ll];
        let (mut i, mut j) = (ll, bj);
        while i > 0 {
            let ob = MASK[p.seq[i - 1] as usize];
            if j > 0 && d[i * w + j] == d[(i - 1) * w + j - 1] + usize::from(ob & READ_BIT[t[j - 1] as usize] == 0) {
                bases[i - 1] = t[j - 1];
                i -= 1;
                j -= 1;
            } else if d[i * w + j] == d[(i - 1) * w + j] + 1 {
                bases[i - 1] = b'-';
                i -= 1;
            } else {
                j -= 1;
            }
        }
        if p.minus {
            // bases were collected along the reverse complement: put them in the oligo's orientation
            bases.reverse();
            for b in bases.iter_mut() {
                if *b != b'-' {
                    *b = COMPL[*b as usize];
                }
            }
        }
        Some(Hit { oligo: p.oligo, minus: p.minus, start: lo + j, end: lo + bj, mismatch: d[ll * w + bj], bases })
    }
}

/// Counts over all reads (single-end) or read pairs (paired-end: a pair counts once, as a fragment, and carries an
/// oligo when either mate does).
pub struct Stats {
    pub units: i64,
    pub bases: i64,
    /// [oligo][mismatches] units whose best hit has that many mismatches
    pub mm: Vec<Vec<i64>>,
    /// units where the oligo reads on the read's plus strand
    pub plus: Vec<i64>,
    /// [oligo][position] A, C, G, T, other, deleted
    pub pos: Vec<Vec<[i64; 6]>>,
    pub sets: HashMap<String, SetCount>,
    set_order: Vec<String>,
}

/// Units carrying the oligos of one primer set.
#[derive(Default, Clone, Copy)]
pub struct SetCount {
    pub f: i64,
    pub r: i64,
    pub p: i64,
    pub fr: i64,
    pub frp: i64,
    pub any: i64,
}

fn new_stats(s: &Searcher) -> Stats {
    let mut st = Stats {
        units: 0,
        bases: 0,
        mm: vec![vec![0; s.max_mm + 1]; s.oligos.len()],
        plus: vec![0; s.oligos.len()],
        pos: s.oligos.iter().map(|o| vec![[0; 6]; o.seq.len()]).collect(),
        sets: HashMap::new(),
        set_order: Vec::new(),
    };
    for o in &s.oligos {
        if !o.set.is_empty() && !st.sets.contains_key(&o.set) {
            st.sets.insert(o.set.clone(), SetCount::default());
            st.set_order.push(o.set.clone());
        }
    }
    st
}

impl Stats {
    fn add(&mut self, o: &Stats) {
        self.units += o.units;
        self.bases += o.bases;
        for i in 0..self.mm.len() {
            for k in 0..self.mm[i].len() {
                self.mm[i][k] += o.mm[i][k];
            }
            self.plus[i] += o.plus[i];
            for p in 0..self.pos[i].len() {
                for b in 0..6 {
                    self.pos[i][p][b] += o.pos[i][p][b];
                }
            }
        }
        for (k, v) in &o.sets {
            let x = self.sets.get_mut(k).unwrap();
            x.f += v.f;
            x.r += v.r;
            x.p += v.p;
            x.fr += v.fr;
            x.frp += v.frp;
            x.any += v.any;
        }
    }

    /// The number of units holding oligo i.
    pub fn found(&self, i: usize) -> i64 {
        self.mm[i].iter().sum()
    }
}

struct Unit {
    ids: Vec<String>,
    seqs: Vec<Vec<u8>>,
}

#[inline]
fn base_idx(b: u8) -> usize {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        b'-' => 5,
        _ => 4,
    }
}

impl Searcher {
    /// Adds one unit's hits (all mates together) to st, and writes hit lines when hits is Some.
    fn count(&self, st: &mut Stats, u: &Unit, hits: Option<&mut String>) {
        st.units += 1;
        let mut best: BTreeMap<usize, Hit> = BTreeMap::new();
        let mut hits = hits;
        for (m, seq) in u.seqs.iter().enumerate() {
            st.bases += seq.len() as i64;
            for h in self.scan(seq) {
                if let Some(out) = hits.as_deref_mut() {
                    let _ = writeln!(
                        out,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        u.ids[m],
                        m + 1,
                        self.oligos[h.oligo].name,
                        if h.minus { '-' } else { '+' },
                        h.start + 1,
                        h.end,
                        h.mismatch,
                        String::from_utf8_lossy(&h.bases)
                    );
                }
                match best.get(&h.oligo) {
                    Some(old) if h.mismatch >= old.mismatch => {}
                    _ => {
                        best.insert(h.oligo, h);
                    }
                }
            }
        }
        let mut sets: HashMap<&str, (bool, bool, bool)> = HashMap::new();
        for (&oi, h) in &best {
            st.mm[oi][h.mismatch] += 1;
            if !h.minus {
                st.plus[oi] += 1;
            }
            for (p, &b) in h.bases.iter().enumerate() {
                st.pos[oi][p][base_idx(b)] += 1;
            }
            let o = &self.oligos[oi];
            if o.set.is_empty() {
                continue;
            }
            let x = sets.entry(o.set.as_str()).or_default();
            match o.role.as_str() {
                "F" => x.0 = true,
                "R" => x.1 = true,
                "P" => x.2 = true,
                _ => {}
            }
        }
        for (name, (f, r, p)) in sets {
            let c = st.sets.get_mut(name).unwrap();
            c.any += 1;
            if f {
                c.f += 1;
            }
            if r {
                c.r += 1;
            }
            if p {
                c.p += 1;
            }
            if f && r {
                c.fr += 1;
                if p {
                    c.frp += 1;
                }
            }
        }
    }

    /// Scans single-end files (each read is a unit) or, with r2 set, paired files read in lockstep (each pair is a
    /// unit). Hit lines go to hits_out when it is Some, in input order.
    pub fn run(&self, r1: &[String], r2: &[String], cpu: usize, hits_out: Option<File>) -> Result<Stats, String> {
        if !r2.is_empty() && r2.len() != r1.len() {
            return Err(format!("{} R1 files but {} R2 files", r1.len(), r2.len()));
        }
        let cpu = cpu.max(1);
        let want_hits = hits_out.is_some();
        let (btx, brx) = mpsc::sync_channel::<(usize, Vec<Unit>)>(cpu * 2);
        let brx = std::sync::Mutex::new(brx);
        let (htx, hrx) = mpsc::sync_channel::<(usize, String)>(cpu * 2);
        std::thread::scope(|sc| {
            let workers: Vec<_> = (0..cpu)
                .map(|_| {
                    let htx = htx.clone();
                    let brx = &brx;
                    sc.spawn(move || {
                        let mut st = new_stats(self);
                        loop {
                            let job = brx.lock().unwrap().recv();
                            let Ok((seq, units)) = job else { break };
                            let mut hb = String::new();
                            for u in &units {
                                self.count(&mut st, u, want_hits.then_some(&mut hb));
                            }
                            if want_hits {
                                let _ = htx.send((seq, hb));
                            }
                        }
                        st
                    })
                })
                .collect();
            drop(htx);
            let writer = hits_out.map(|f| {
                sc.spawn(move || -> std::io::Result<()> {
                    let mut bw = BufWriter::with_capacity(1 << 20, f);
                    writeln!(bw, "#read_id\tmate\toligo\tstrand\tstart\tend\tmismatch\tread_bases")?;
                    let mut pending: BTreeMap<usize, String> = BTreeMap::new();
                    let mut next = 0;
                    for (seq, s) in hrx {
                        pending.insert(seq, s);
                        while let Some(s) = pending.remove(&next) {
                            bw.write_all(s.as_bytes())?;
                            next += 1;
                        }
                    }
                    bw.flush()
                })
            });
            let read_err = feed(r1, r2, &btx);
            drop(btx);
            let mut total = new_stats(self);
            for w in workers {
                total.add(&w.join().unwrap());
            }
            if let Some(w) = writer {
                w.join().unwrap().map_err(|e| e.to_string())?;
            }
            read_err.map(|_| total)
        })
    }
}

fn feed(r1: &[String], r2: &[String], out: &mpsc::SyncSender<(usize, Vec<Unit>)>) -> Result<(), String> {
    // a batch closes at 4096 units or 256 kb, so a few thousand long reads still spread over every thread
    const BATCH: usize = 4096;
    const BATCH_BASES: usize = 1 << 18;
    let mut seq = 0;
    let mut send = |buf: Vec<Unit>| {
        let _ = out.send((seq, buf));
        seq += 1;
    };
    for fi in 0..r1.len() {
        let mut a = seqio::open(&r1[fi])?;
        let mut b = if r2.is_empty() { None } else { Some(seqio::open(&r2[fi])?) };
        let mut buf = Vec::with_capacity(BATCH);
        let mut bases = 0;
        loop {
            let Some(x) = a.next()? else {
                if let Some(b) = b.as_mut() {
                    if !matches!(b.next(), Ok(None)) {
                        return Err(format!("{} has more reads than {}", r2[fi], r1[fi]));
                    }
                }
                break;
            };
            let mut u = Unit { ids: vec![x.id], seqs: vec![x.seq] };
            if let Some(b) = b.as_mut() {
                let y = match b.next() {
                    Ok(Some(y)) => y,
                    Ok(None) => return Err(format!("{} ends before {} (EOF)", r2[fi], r1[fi])),
                    Err(e) => return Err(format!("{} ends before {} ({})", r2[fi], r1[fi], e)),
                };
                if mate_id(&u.ids[0]) != mate_id(&y.id) {
                    return Err(format!("read pair out of order: {} / {}", u.ids[0], y.id));
                }
                u.ids.push(y.id);
                u.seqs.push(y.seq);
            }
            bases += u.seqs.iter().map(Vec::len).sum::<usize>();
            buf.push(u);
            if buf.len() == BATCH || bases >= BATCH_BASES {
                send(std::mem::replace(&mut buf, Vec::with_capacity(BATCH)));
                bases = 0;
            }
        }
        if !buf.is_empty() {
            send(buf);
        }
    }
    Ok(())
}

/// Strips a /1 or /2 suffix so the two mates of a pair compare equal.
fn mate_id(id: &str) -> &str {
    let b = id.as_bytes();
    let n = b.len();
    if n > 2 && b[n - 2] == b'/' && (b[n - 1] == b'1' || b[n - 1] == b'2') {
        &id[..n - 2]
    } else {
        id
    }
}

/// Writes <prefix>.summary.tsv, .sets.tsv (when the oligos are named SET_F/_R/_P) and .positions.tsv, and returns a
/// short human-readable report.
pub fn write_reports(prefix: &str, s: &Searcher, st: &Stats, paired: bool) -> Result<String, String> {
    let unit_name = if paired { "read_pairs" } else { "reads" };
    let pct = |n: i64| if st.units == 0 { 0.0 } else { 100.0 * n as f64 / st.units as f64 };
    let wr = |p: String, s: &str| std::fs::write(&p, s).map_err(|e| format!("{p}: {e}"));
    let mut sum = String::new();
    let _ = write!(sum, "#oligo\tset\trole\tlength\t{unit_name}\twith_oligo\tpct");
    for k in 0..=s.max_mm {
        let _ = write!(sum, "\tmm{k}");
    }
    sum.push_str("\tplus_strand\tminus_strand\n");
    let mut rep = String::new();
    let mode = if s.indel { "edits (indels allowed)" } else { "substitutions" };
    let _ = write!(
        rep,
        "chichipes reads: {} {}, {} bases, up to {} {} per oligo\n\n",
        st.units,
        unit_name.replace('_', " "),
        st.bases,
        s.max_mm,
        mode
    );
    let _ = writeln!(rep, "{:<32} {:>12} {:>9}", "oligo", "with oligo", "%");
    for (i, o) in s.oligos.iter().enumerate() {
        let n = st.found(i);
        let _ = write!(sum, "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}", o.name, o.set, o.role, o.seq.len(), st.units, n, pct(n));
        for v in &st.mm[i] {
            let _ = write!(sum, "\t{v}");
        }
        let _ = writeln!(sum, "\t{}\t{}", st.plus[i], n - st.plus[i]);
        let _ = writeln!(rep, "{:<32} {:>12} {:>9.4}", o.name, n, pct(n));
    }
    wr(format!("{prefix}.summary.tsv"), &sum)?;
    if !st.set_order.is_empty() {
        let mut b = String::new();
        let _ = writeln!(b, "#set\t{unit_name}\tF\tR\tP\tF_and_R\tF_R_and_P\tany\tpct_F_and_R\tpct_F_R_and_P");
        let _ = write!(rep, "\n{:<24} {:>10} {:>10} {:>10} {:>10} {:>10}\n", "set", "F", "R", "P", "F+R", "F+R+P");
        let mut names = st.set_order.clone();
        names.sort();
        for n in &names {
            let c = st.sets[n];
            let _ = writeln!(
                b,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}",
                n,
                st.units,
                c.f,
                c.r,
                c.p,
                c.fr,
                c.frp,
                c.any,
                pct(c.fr),
                pct(c.frp)
            );
            let _ = writeln!(rep, "{:<24} {:>10} {:>10} {:>10} {:>10} {:>10}", n, c.f, c.r, c.p, c.fr, c.frp);
        }
        wr(format!("{prefix}.sets.tsv"), &b)?;
    }
    let mut pb = String::from("#oligo\tpos\toligo_base\tdepth\tA\tC\tG\tT\tother\tdeleted\tmatch_pct\n");
    for (i, o) in s.oligos.iter().enumerate() {
        for (p, c) in st.pos[i].iter().enumerate() {
            let depth: i64 = c.iter().sum();
            let mut matched = 0;
            for (b, &v) in c.iter().take(4).enumerate() {
                if MASK[o.seq[p] as usize] & (1 << b) != 0 {
                    matched += v;
                }
            }
            let mp = if depth > 0 { 100.0 * matched as f64 / depth as f64 } else { 0.0 };
            let _ = writeln!(
                pb,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}",
                o.name,
                p + 1,
                o.seq[p] as char,
                depth,
                c[0],
                c[1],
                c[2],
                c[3],
                c[4],
                c[5],
                mp
            );
        }
    }
    wr(format!("{prefix}.positions.tsv"), &pb)?;
    wr(format!("{prefix}.txt"), &rep)?;
    Ok(rep)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn searcher(seq: &str, mm: usize, indel: bool) -> Searcher {
        let o = Oligo { name: "S_F".into(), set: "S".into(), role: "F".into(), seq: seq.as_bytes().to_vec() };
        new_searcher(vec![o], mm, 0, indel).unwrap()
    }

    #[test]
    fn scan_hamming_both_strands() {
        let s = searcher("ACGTTGCAGTGCATCGAGTC", 2, false);
        let read = "TTTTTACGTTGCAGTGAATCGAGTCGGGGG"; // one mismatch
        let h = s.scan(read.as_bytes());
        assert!(h.len() == 1 && h[0].mismatch == 1 && !h[0].minus && h[0].start == 5);
        let rc = rev_comp(read.as_bytes());
        let h = s.scan(&rc);
        assert!(h.len() == 1 && h[0].mismatch == 1 && h[0].minus && h[0].bases == b"ACGTTGCAGTGAATCGAGTC");
        assert!(s.scan(b"ACGTTGCAGTGAATCGAG").is_empty()); // oligo not wholly inside the read
    }

    #[test]
    fn scan_inosine_and_iupac() {
        let s = searcher("ACGTTGCAITGCATCGRGTC", 0, false);
        for x in ["A", "C", "G", "T"] {
            for y in ["A", "G"] {
                let h = s.scan(format!("ACGTTGCA{x}TGCATCG{y}GTC").as_bytes());
                assert!(h.len() == 1 && h[0].mismatch == 0, "I={x} R={y}");
            }
        }
        assert!(s.scan(b"ACGTTGCAATGCATCGCGTC").is_empty()); // R does not take C
    }

    #[test]
    fn scan_indel() {
        let s = searcher("ACGTTGCAGTGCATCGAGTC", 2, true);
        let h = s.scan(b"CCCCACGTTGCAGTGCAATCGAGTCCCCC"); // one inserted base
        assert!(h.len() == 1 && h[0].mismatch == 1);
        let h = s.scan(b"CCCCACGTTGCAGTCATCGAGTCCCCC"); // one deleted base
        assert!(h.len() == 1 && h[0].mismatch == 1);
    }
}
