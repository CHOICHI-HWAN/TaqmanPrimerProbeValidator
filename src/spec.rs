//! Finds primer binding sites and the amplicons they form across a genome database — the clean-room counterpart of
//! `mfeprimer --virus`.

use crate::primers::{self, Oligo};
use crate::seqio::Record;
use crate::thermo::{self, Cond};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::sync::atomic::{AtomicUsize, Ordering};

/// The mfeprimer flags the pipeline uses.
#[derive(Clone, Debug)]
pub struct Options {
    pub k: usize,
    pub mis_match: usize,
    pub mis_start: usize,
    pub mis_end: usize,
    pub min_size: i64,
    pub max_size: i64,
    pub min_tm: f64,
    pub max_tm: f64,
    pub cpu: usize,
    pub cond: Cond,
}

pub fn num_cpu() -> usize {
    std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1)
}

impl Default for Options {
    /// mfeprimer's defaults.
    fn default() -> Self {
        Options {
            k: 9,
            mis_match: 0,
            mis_start: 1,
            mis_end: 9,
            min_size: 0,
            max_size: 2000,
            min_tm: 30.0,
            max_tm: 100.0,
            cpu: num_cpu(),
            cond: thermo::DEFAULT,
        }
    }
}

/// One searchable sequence: a variant of a concrete primer.
pub struct Probe {
    pub name: String,
    /// 5'->3'
    pub seq: Vec<u8>,
    /// reverse complement
    pub rcs: Vec<u8>,
}

/// One primer binding site on a genome, in plus-strand coordinates (0-based start, length n).
#[derive(Clone, Copy)]
pub struct Site {
    pub pos: usize,
    pub n: usize,
    /// the primer binds the plus strand: its reverse complement is read at pos
    pub rev: bool,
    pub p: u32,
    pub tm: f64,
    pub dg: f64,
}

/// How close (°C) a binding Tm may lie to the -t cut-off to be listed as a borderline case.
pub const BORDER_BAND: f64 = 0.02;

/// One binding whose Tm lies within BORDER_BAND of the cut-off.
pub struct Border {
    pub site: Site,
    pub kept: bool,
}

/// One product: forward site f (binds the minus strand) and reverse site r (binds the plus strand).
pub struct Amplicon {
    pub rec: usize,
    pub f: Site,
    pub r: Site,
}

impl Amplicon {
    pub fn start(&self) -> usize {
        self.f.pos
    }
    pub fn end(&self) -> usize {
        self.r.pos + self.r.n
    }
    pub fn size(&self) -> i64 {
        self.end() as i64 - self.start() as i64
    }
}

static COMP: [u8; 256] = {
    let mut t = [b'N'; 256];
    let pairs: [&[u8; 2]; 16] =
        [b"AT", b"CG", b"GC", b"TA", b"RY", b"YR", b"KM", b"MK", b"SS", b"WW", b"BV", b"VB", b"DH", b"HD", b"NN", b"II"];
    let mut i = 0;
    while i < pairs.len() {
        let p = pairs[i];
        t[p[0] as usize] = p[1];
        t[(p[0] + 32) as usize] = p[1] + 32;
        i += 1;
    }
    t
};

#[inline]
pub fn comp(b: u8) -> u8 {
    COMP[b as usize]
}

pub fn rev_comp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| COMP[b as usize]).collect()
}

#[inline]
pub fn upper(b: u8) -> u8 {
    b.to_ascii_uppercase()
}

/// Packs an oligo window into 2 bits per base; an inosine gives one key per base it can face (4^m keys for m
/// inosines), so the window is found whatever the genome has there. Empty when it holds another letter.
pub fn encode_all(s: &[u8]) -> Vec<u32> {
    let mut keys = vec![0u32];
    for &b in s {
        if b == b'I' || b == b'i' {
            let mut next = Vec::with_capacity(4 * keys.len());
            for &v in &keys {
                for c in 0..4u32 {
                    next.push(v << 2 | c);
                }
            }
            keys = next;
            continue;
        }
        let c = thermo::code(b);
        if c < 0 {
            return Vec::new();
        }
        for k in keys.iter_mut() {
            *k = *k << 2 | c as u32;
        }
    }
    keys
}

pub fn kmer_mask(k: usize) -> u32 {
    ((1u64 << (2 * k)) - 1) as u32
}

/// k-mer -> list of u32 values: a dense table for short k, a hash map otherwise.
pub enum KmerIndex {
    Dense { start: Vec<u32>, items: Vec<u32> },
    Sparse(HashMap<u32, Vec<u32>>),
}

impl KmerIndex {
    pub fn build(k: usize, pairs: &[(u32, u32)]) -> Self {
        if k <= 11 {
            let size = 1usize << (2 * k);
            let mut count = vec![0u32; size + 1];
            for &(w, _) in pairs {
                count[w as usize + 1] += 1;
            }
            for i in 0..size {
                count[i + 1] += count[i];
            }
            let mut fill = count.clone();
            let mut items = vec![0u32; pairs.len()];
            for &(w, v) in pairs {
                items[fill[w as usize] as usize] = v;
                fill[w as usize] += 1;
            }
            KmerIndex::Dense { start: count, items }
        } else {
            let mut m: HashMap<u32, Vec<u32>> = HashMap::new();
            for &(w, v) in pairs {
                m.entry(w).or_default().push(v);
            }
            KmerIndex::Sparse(m)
        }
    }
    #[inline]
    pub fn get(&self, w: u32) -> &[u32] {
        match self {
            KmerIndex::Dense { start, items } => {
                &items[start[w as usize] as usize..start[w as usize + 1] as usize]
            }
            KmerIndex::Sparse(m) => m.get(&w).map(|v| &v[..]).unwrap_or(&[]),
        }
    }
}

/// The probes indexed by their 3' window.
pub struct Engine {
    opt: Options,
    pub probes: Vec<Probe>,
    /// 3' window (as written 5'->3') -> probes
    fwd: KmerIndex,
    /// reverse complement of the 3' window -> probes
    rev: KmerIndex,
}

/// Bindings of one genome.
#[derive(Default)]
struct SiteSet {
    fs: Vec<Site>,
    rs: Vec<Site>,
    border: Vec<Border>,
}

impl Engine {
    /// Builds the search set from the input oligos.
    pub fn new(oligos: &[Oligo], opt: Options) -> Result<Engine, String> {
        let mut probes = Vec::new();
        let (mut fwd, mut rev) = (Vec::new(), Vec::new());
        for o in oligos {
            let cs = primers::build(o, opt.mis_match, opt.mis_start, opt.mis_end)?;
            // an expansion equal to an earlier variant is searched once, under the first name
            let mut seen: HashSet<&str> = HashSet::new();
            for c in &cs {
                for v in &c.variants {
                    if !seen.insert(v.seq.as_str()) {
                        continue;
                    }
                    let seq = v.seq.as_bytes().to_vec();
                    let p = Probe { name: c.variant_name(v), rcs: rev_comp(&seq), seq };
                    let pi = probes.len() as u32;
                    if p.seq.len() >= opt.k {
                        for w in encode_all(&p.seq[p.seq.len() - opt.k..]) {
                            fwd.push((w, pi));
                        }
                        for w in encode_all(&p.rcs[..opt.k]) {
                            rev.push((w, pi));
                        }
                    }
                    probes.push(p);
                }
            }
        }
        let (fwd, rev) = (KmerIndex::build(opt.k, &fwd), KmerIndex::build(opt.k, &rev));
        Ok(Engine { opt, probes, fwd, rev })
    }

    /// Every binding site on one genome that passes the Tm window.
    fn sites(&self, g: &[u8]) -> SiteSet {
        let k = self.opt.k;
        let n = g.len();
        let mask = kmer_mask(k);
        let mut out = SiteSet::default();
        let mut s = Vec::with_capacity(64);
        let mut w = 0u32;
        let mut valid = 0usize; // consecutive A/C/G/T bases ending at the current position
        let c = self.opt.cond;
        for i in 0..n {
            let x = thermo::code(g[i]);
            if x < 0 {
                valid = 0;
                w = 0;
                continue;
            }
            w = (w << 2 | x as u32) & mask;
            valid += 1;
            if valid < k {
                continue;
            }
            let end = i + 1; // the window is g[end-k..end]
            for &pi in self.fwd.get(w) {
                let p = &self.probes[pi as usize];
                let pl = p.seq.len();
                if end < pl {
                    continue;
                }
                let pos = end - pl;
                let l5 = if pos > 0 { upper(g[pos - 1]) } else { b'-' }; // the site starts the genome
                s.clear();
                s.extend(g[pos..end].iter().map(|&b| upper(b)));
                let r3 = if end < n { upper(g[end]) } else { b'A' };
                let bind = thermo::binding(&p.seq, &s, l5, r3, false, c);
                self.keep(&mut out.fs, &mut out.border, bind, Site { pos, n: pl, rev: false, p: pi, tm: 0.0, dg: 0.0 });
            }
            let start = end - k; // the reverse complement of the primer's 3' window starts here
            for &pi in self.rev.get(w) {
                let p = &self.probes[pi as usize];
                let pl = p.seq.len();
                if start + pl > n || start < 1 {
                    // (a site at the very start cannot be the right end of an amplicon)
                    continue;
                }
                s.clear(); // the site read on the minus strand, 5'->3' like the primer
                s.extend(g[start..start + pl].iter().rev().map(|&b| comp(upper(b))));
                let l5 = if start + pl < n { comp(upper(g[start + pl])) } else { b'-' }; // the site ends the genome
                let r3 = comp(upper(g[start - 1]));
                let bind = thermo::binding(&p.seq, &s, l5, r3, true, c);
                self.keep(&mut out.rs, &mut out.border, bind, Site { pos: start, n: pl, rev: true, p: pi, tm: 0.0, dg: 0.0 });
            }
        }
        out
    }

    #[inline]
    fn keep(&self, list: &mut Vec<Site>, border: &mut Vec<Border>, bind: Option<(f64, f64)>, mut site: Site) {
        let ok = bind.is_some();
        let (tm, dg) = bind.unwrap_or((0.0, 0.0));
        site.tm = tm;
        site.dg = dg;
        let keep = ok && tm >= self.opt.min_tm && tm <= self.opt.max_tm;
        if keep {
            list.push(site);
        }
        if ok && (tm - self.opt.min_tm).abs() <= BORDER_BAND {
            border.push(Border { site, kept: keep });
        }
    }

    /// Pairs the sites of one genome into products within the size limits.
    fn amplicons(&self, rec: usize, fs: &[Site], rs: &mut [Site]) -> Vec<Amplicon> {
        // the sites come in position order already; a stable sort leaves them as the Go engine's sort does
        rs.sort_by_key(|s| s.pos);
        let mut out = Vec::new();
        for f in fs {
            for r in rs.iter() {
                if r.pos < f.pos {
                    continue;
                }
                let a = Amplicon { rec, f: *f, r: *r };
                let sz = a.size();
                if sz > self.opt.max_size {
                    break;
                }
                // the two primers may overlap by fewer than k bases (larger overlaps are never reported)
                if sz < self.opt.min_size || f.pos as i64 + f.n as i64 - r.pos as i64 >= self.opt.k as i64 {
                    continue;
                }
                out.push(a);
            }
        }
        out
    }

    /// Searches every record in parallel and returns the amplicons in database order.
    pub fn run(&self, db: &[Record]) -> Result_ {
        let cpu = self.opt.cpu.max(1);
        let next = AtomicUsize::new(0);
        let mut per: Vec<(usize, Vec<Amplicon>, Vec<Border>)> = std::thread::scope(|sc| {
            let hs: Vec<_> = (0..cpu)
                .map(|_| {
                    sc.spawn(|| {
                        let mut mine = Vec::new();
                        loop {
                            let i = next.fetch_add(1, Ordering::Relaxed);
                            if i >= db.len() {
                                break;
                            }
                            let mut ss = self.sites(&db[i].seq);
                            let amps = self.amplicons(i, &ss.fs, &mut ss.rs);
                            mine.push((i, amps, ss.border));
                        }
                        mine
                    })
                })
                .collect();
            hs.into_iter().flat_map(|h| h.join().unwrap()).collect()
        });
        per.sort_by_key(|x| x.0);
        let mut res = Result_ { amplicons: Vec::new(), border: Vec::with_capacity(db.len()), failed: Vec::new(), total: db.len() };
        for (i, a, b) in per {
            if a.is_empty() {
                res.failed.push(i);
            }
            res.amplicons.extend(a);
            res.border.push(b);
        }
        res
    }

    pub fn probe(&self, s: &Site) -> &Probe {
        &self.probes[s.p as usize]
    }
}

/// The outcome over a whole database.
pub struct Result_ {
    pub amplicons: Vec<Amplicon>,
    /// per record: bindings near the Tm cut-off
    pub border: Vec<Vec<Border>>,
    /// records without any amplicon, in database order
    pub failed: Vec<usize>,
    pub total: usize,
}

fn create(path: &str) -> io::Result<BufWriter<File>> {
    Ok(BufWriter::with_capacity(1 << 20, File::create(path)?))
}

/// Writes the files the pipeline reads, in mfeprimer's formats: <out> (report with the coverage line),
/// <out>.spec.tsv, <out>.amp.fasta (first max_amp products) and <out>.virus.Failed.txt.
pub fn write_outputs(e: &Engine, out: &str, db_path: &str, db: &[Record], res: &Result_, max_amp: usize) -> io::Result<()> {
    let spec_path = format!("{out}.spec.tsv");
    let _ = std::fs::remove_file(&spec_path);
    // mfeprimer writes no spec.tsv at all when nothing amplifies
    let mut sw: Box<dyn Write> =
        if res.amplicons.is_empty() { Box::new(io::sink()) } else { Box::new(create(&spec_path)?) };
    writeln!(sw, "#1-based coordinate, should be idential to blat output if align primer sequences to genome")?;
    writeln!(
        sw,
        "#name\tchrom\tampStart\tampEnd\tampGC\tampSize\tfpName\tfpStart\tfpEnd\tfpSeq\tfpTm\tfpGC\tfpDg\t\
         rpName\trpEnd\trpStart\trpSeq\trpTm\trpGC\trpDg\tnote"
    )?;
    let mut aw = create(&format!("{out}.amp.fasta"))?;
    let mut prod = Vec::new();
    for (i, a) in res.amplicons.iter().enumerate() {
        let g = &db[a.rec].seq;
        let (fp, rp) = (e.probe(&a.f), e.probe(&a.r));
        let fseq = &fp.seq[..];
        let rplus = &rp.rcs[..]; // the reverse primer as it reads on the plus strand
        let overlap = a.r.pos < a.f.pos + a.f.n;
        let mid: &[u8] = if overlap { &[] } else { &g[a.f.pos + a.f.n..a.r.pos] };
        prod.clear();
        if overlap {
            // the two primers overlap: no middle; GC from the genome region (docs/FINDINGS.md, known difference)
            prod.extend_from_slice(&g[a.start()..a.end()]);
        } else {
            // the product as amplified: primer sequences at both ends
            prod.extend_from_slice(fseq);
            prod.extend_from_slice(mid);
            prod.extend_from_slice(rplus);
        }
        let id = &db[a.rec].id;
        write!(
            sw,
            "Amp_{}\t{}\t{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t",
            i + 1,
            id,
            a.start() + 1,
            a.end(),
            thermo::gc(&prod),
            a.size(),
            fp.name,
            a.f.pos + 1,
            a.f.pos + a.f.n
        )?;
        sw.write_all(fseq)?;
        write!(sw, "\t{:.2}\t{:.2}\t{:.2}\t{}\t{}\t{}\t", a.f.tm, thermo::gc(fseq), a.f.dg, rp.name, a.r.pos + 1, a.r.pos + a.r.n)?;
        sw.write_all(rplus)?;
        write!(sw, "\t{:.2}\t{:.2}\t{:.2}\t\n", a.r.tm, thermo::gc(&rp.seq), a.r.dg)?;
        if i < max_amp {
            write!(aw, ">Amp_{} {}+{}==>{}:{}-{} {}\n", i + 1, fp.name, rp.name, id, a.start(), a.end(), db_path)?;
            if overlap {
                aw.write_all(&prod.to_ascii_uppercase())?;
            } else {
                aw.write_all(fseq)?;
                aw.write_all(&mid.to_ascii_lowercase())?;
                aw.write_all(rplus)?;
            }
            aw.write_all(b"\n")?;
        }
    }
    sw.flush()?;
    drop(sw);
    aw.flush()?;
    drop(aw);

    let mut fw = create(&format!("{out}.virus.Failed.txt"))?;
    for &i in &res.failed {
        write!(fw, "{} \n", db[i].id)?;
    }
    fw.flush()?;
    drop(fw);

    write_border(e, &format!("{out}.borderline.tsv"), db, res)?;

    let mut rep = File::create(out)?;
    let cov = if res.total > 0 { 100.0 * (res.total - res.failed.len()) as f64 / res.total as f64 } else { 0.0 };
    let mut s = String::from("chichipes report (mfeprimer-compatible outputs)\n\n");
    s.push_str("Primer Coverate Rate (virus mode, e.g., to check whether these primers can amplify all the variants of SARS-CoV-2 genomes)\n");
    s.push_str(&format!(
        "Coverage {:.2}%, total records in database: {}, amplicon number: {}, missed number: {}.\n",
        cov,
        res.total,
        res.amplicons.len(),
        res.failed.len()
    ));
    rep.write_all(s.as_bytes())
}

/// Lists the bindings whose Tm lies within BORDER_BAND of the -t cut-off (engine-only side file).
fn write_border(e: &Engine, path: &str, db: &[Record], res: &Result_) -> io::Result<()> {
    let mut w = create(path)?;
    writeln!(w, "# bindings with Tm within {BORDER_BAND:.2} C of the cut-off; decision may differ from mfeprimer")?;
    writeln!(w, "#chrom\tstart\tend\tstrand\tname\tseq\ttm\tdg\tdecision")?;
    for (i, bs) in res.border.iter().enumerate() {
        for b in bs {
            let s = &b.site;
            let p = e.probe(s);
            write!(w, "{}\t{}\t{}\t{}\t{}\t", db[i].id, s.pos + 1, s.pos + s.n, if s.rev { "-" } else { "+" }, p.name)?;
            w.write_all(&p.seq)?;
            write!(w, "\t{:.3}\t{:.3}\t{}\n", s.tm, s.dg, if b.kept { "kept" } else { "dropped" })?;
        }
    }
    w.flush()
}
