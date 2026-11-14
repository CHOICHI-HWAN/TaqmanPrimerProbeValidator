//! Oligo, primer-binding, probe-search, dimer and hairpin thermodynamics (docs/FINDINGS.md §3, §6–8).
//!
//! Perfect duplexes use the SantaLucia & Hicks (2004) nearest-neighbour model with the von Ahsen (2001) Mg2+/dNTP
//! salt equivalent; mismatch stacks, dangling ends and tri/tetraloop bonuses are SantaLucia & Hicks values as
//! distributed with primer3, with the few entries where MFEprimer's table differs replaced by the measured ones.
//! Deoxyinosine stacks are Watkins & SantaLucia (2005), Table 2.
//!
//! Every sum is taken in the same order as the Go engine so the numbers agree bit for bit.

use crate::golog::go_log;
use std::collections::HashMap;
use std::sync::OnceLock;

const GAS_R: f64 = 1.987;
/// The hairpin Tm and the report oligo tables use this (the search and `thermo` use 1.987).
const HP_GAS_R: f64 = 1.9872;
const T37: f64 = 310.15;

/// The reaction buffer (mfeprimer's --mono/--diva/--dntp in mM, --oligo in nM).
#[derive(Clone, Copy, Debug)]
pub struct Cond {
    pub mono: f64,
    pub diva: f64,
    pub dntp: f64,
    pub oligo: f64,
}

/// mfeprimer's defaults.
pub const DEFAULT: Cond = Cond { mono: 50.0, diva: 1.5, dntp: 0.25, oligo: 50.0 };

impl Cond {
    /// The monovalent-equivalent salt concentration in M.
    pub fn na_eq(&self) -> f64 {
        (self.mono + 120.0 * (self.diva - self.dntp).max(0.0).sqrt()) / 1000.0
    }
}

/// Code maps A,C,G,T (either case) to 0..3 and anything else to -1.
#[inline]
pub fn code(b: u8) -> i32 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => -1,
    }
}

/// The code of deoxyinosine (A, C, G, T are 0..3, N/ambiguity 4).
pub const INO: i32 = 5;

/// Code for an oligo base: A/C/G/T 0..3, I INO, anything else -1.
#[inline]
pub fn ocode(b: u8) -> i32 {
    if b == b'I' || b == b'i' {
        INO
    } else {
        code(b)
    }
}

/// Whether the oligo holds an I.
pub fn has_inosine(s: &[u8]) -> bool {
    s.iter().any(|&b| b == b'I' || b == b'i')
}

/// Code for a target base, every IUPAC ambiguity code (R, Y, ..., N) as 4: mfeprimer scores them all like N.
#[inline]
fn site_code(b: u8) -> i32 {
    let c = code(b);
    if c >= 0 {
        return c;
    }
    match b | 0x20 {
        b'r' | b'y' | b'k' | b'm' | b's' | b'w' | b'b' | b'd' | b'h' | b'v' | b'n' => 4,
        _ => -1,
    }
}

/// site_code for a flanking base, with '-' (the site touches the genome end) as 5.
#[inline]
fn flank_code(b: u8) -> i32 {
    if b == b'-' {
        5
    } else {
        site_code(b)
    }
}

/// The code of the complementary base (N stays N).
#[inline]
fn comp_code(c: i32) -> i32 {
    if c == 4 {
        4
    } else {
        3 - c
    }
}

struct Tables {
    nn_h: [[f64; 4]; 4],
    nn_s: [[f64; 4]; 4],
    /// 5'ab3'/3'cd5' -> dH, dS
    mm: [[[[[f64; 2]; 4]; 4]; 4]; 4],
    /// [half][pair 3' base][pair 5' base][dangling base]
    dangle: [[[[[f64; 2]; 4]; 4]; 4]; 2],
    bonus: HashMap<Vec<u8>, f64>,
    loop_l: Vec<i64>,
    loop_g: Vec<f64>,
    ino: HashMap<[i32; 4], [f64; 2]>,
}

fn tables() -> &'static Tables {
    static T: OnceLock<Tables> = OnceLock::new();
    T.get_or_init(build_tables)
}

/// Watkins & SantaLucia (2005) Table 2: "top/bottom" with top 5'->3' over bottom 3'->5', dH (kcal/mol), dS (cal/K/mol).
const WATKINS_2005: [(&str, [f64; 2]); 36] = [
    // I·C
    ("AI/TC", [-8.9, -25.5]), ("TI/AC", [-5.9, -17.4]), ("AC/TI", [-8.8, -25.4]), ("TC/AI", [-4.9, -13.9]),
    ("CI/GC", [-5.4, -13.7]), ("GI/CC", [-6.8, -19.1]), ("CC/GI", [-8.3, -23.8]), ("GC/CI", [-5.0, -12.6]),
    // I·A
    ("AI/TA", [-8.3, -25.0]), ("TI/AA", [-3.4, -11.2]), ("AA/TI", [-0.7, -2.6]), ("TA/AI", [-1.3, -4.6]),
    ("CI/GA", [2.6, 8.9]), ("GI/CA", [-7.8, -21.1]), ("CA/GI", [-7.0, -20.0]), ("GA/CI", [-7.6, -20.2]),
    // I·T
    ("AI/TT", [0.49, -0.7]), ("TI/AT", [-6.5, -22.0]), ("AT/TI", [-5.6, -18.7]), ("TT/AI", [-0.8, -4.3]),
    ("CI/GT", [-1.0, -2.4]), ("GI/CT", [-3.5, -10.6]), ("CT/GI", [0.1, -1.0]), ("GT/CI", [-4.3, -12.1]),
    // I·G
    ("AI/TG", [-4.9, -15.8]), ("TI/AG", [-1.9, -8.5]), ("AG/TI", [0.1, -1.8]), ("TG/AI", [1.0, 1.0]),
    ("CI/GG", [7.1, 21.3]), ("GI/CG", [-1.1, -3.2]), ("CG/GI", [5.8, 16.9]), ("GG/CI", [-7.6, -22.0]),
    // I·I (oligo against oligo)
    ("AI/TI", [-3.3, -11.9]), ("TI/AI", [0.1, -2.3]), ("CI/GI", [1.3, 3.0]), ("GI/CI", [-0.5, -1.3]),
];

/// The printed dG37 of each Watkins row, kept to check the transcription of dH/dS.
#[cfg(test)]
const WATKINS_DG37: [(&str, f64); 36] = [
    ("AI/TC", -0.96), ("TI/AC", -0.46), ("AC/TI", -0.89), ("TC/AI", -0.59), ("CI/GC", -1.14), ("GI/CC", -0.86),
    ("CC/GI", -0.88), ("GC/CI", -1.07), ("AI/TA", -0.51), ("TI/AA", 0.09), ("AA/TI", 0.12), ("TA/AI", 0.12),
    ("CI/GA", -0.18), ("GI/CA", -1.24), ("CA/GI", -0.77), ("GA/CI", -1.33), ("AI/TT", 0.71), ("TI/AT", 0.36),
    ("AT/TI", 0.22), ("TT/AI", 0.54), ("CI/GT", -0.26), ("GI/CT", -0.19), ("CT/GI", 0.41), ("GT/CI", -0.54),
    ("AI/TG", 0.02), ("TI/AG", 0.76), ("AG/TI", 0.65), ("TG/AI", 0.70), ("CI/GG", 0.47), ("GI/CG", -0.10),
    ("CG/GI", 0.54), ("GG/CI", -0.74), ("AI/TI", 0.40), ("TI/AI", 0.81), ("CI/GI", 0.36), ("GI/CI", -0.09),
];

fn build_tables() -> Tables {
    let mut t = Tables {
        nn_h: [[0.0; 4]; 4],
        nn_s: [[0.0; 4]; 4],
        mm: [[[[[0.0; 2]; 4]; 4]; 4]; 4],
        dangle: [[[[[0.0; 2]; 4]; 4]; 4]; 2],
        bonus: HashMap::new(),
        loop_l: Vec::new(),
        loop_g: Vec::new(),
        ino: HashMap::new(),
    };
    let nn: [(&[u8; 2], [f64; 2]); 10] = [
        (b"AA", [-7.6, -21.3]), (b"AT", [-7.2, -20.4]), (b"TA", [-7.2, -21.3]), (b"CA", [-8.5, -22.7]),
        (b"GT", [-8.4, -22.4]), (b"CT", [-7.8, -21.0]), (b"GA", [-8.2, -22.2]), (b"CG", [-10.6, -27.2]),
        (b"GC", [-9.8, -24.4]), (b"GG", [-8.0, -19.9]),
    ];
    for (k, v) in nn {
        let (a, b) = (code(k[0]) as usize, code(k[1]) as usize);
        t.nn_h[a][b] = v[0];
        t.nn_s[a][b] = v[1];
        let (ra, rb) = (3 - b, 3 - a); // the same stack read on the other strand
        t.nn_h[ra][rb] = v[0];
        t.nn_s[ra][rb] = v[1];
    }
    let j = crate::json::parse(include_str!("../data/hairpin_tables.json")).expect("hairpin_tables.json");
    for (k, v) in j.get("mm").as_object() {
        let k = k.as_bytes(); // "ab/cd"
        let v = v.as_array();
        t.mm[code(k[0]) as usize][code(k[1]) as usize][code(k[3]) as usize][code(k[4]) as usize] =
            [v[0].as_f64(), v[1].as_f64()];
    }
    for (k, v) in j.get("dangle").as_object() {
        let k = k.as_bytes(); // "half:abc"
        let v = v.as_array();
        t.dangle[(k[0] - b'0') as usize][code(k[2]) as usize][code(k[3]) as usize][code(k[4]) as usize] =
            [v[0].as_f64(), v[1].as_f64()];
    }
    for (k, v) in j.get("bonus").as_object() {
        t.bonus.insert(k.as_bytes().to_vec(), v.as_f64());
    }
    t.loop_l = j.get("loop_sizes").as_array().iter().map(|v| v.as_f64() as i64).collect();
    t.loop_g = j.get("loop_g").as_array().iter().map(|v| v.as_f64()).collect();
    for (k, v) in WATKINS_2005 {
        let k = k.as_bytes();
        let (t0, t1, b0, b1) = (ocode(k[0]), ocode(k[1]), ocode(k[3]), ocode(k[4]));
        t.ino.insert([t0, t1, b0, b1], v);
        t.ino.insert([b1, b0, t1, t0], v); // the same stack read on the other strand
    }
    t
}

fn is_self_comp(f: &[u8]) -> bool {
    let n = f.len();
    (0..n).all(|i| code(f[i]) == 3 - code(f[n - 1 - i]))
}

/// dH and salt-corrected dS of f against its exact complement. f must be A/C/G/T (or I).
fn perfect_hs(f: &[u8], na: f64) -> (f64, f64) {
    if has_inosine(f) {
        return perfect_hs_ino(f, na);
    }
    let t = tables();
    let n = f.len();
    let (mut h, mut s) = (0.0, 0.0);
    for i in 0..n.saturating_sub(1) {
        let (a, b) = (code(f[i]) as usize, code(f[i + 1]) as usize);
        h += t.nn_h[a][b];
        s += t.nn_s[a][b];
    }
    h += 0.2;
    s += -5.7;
    for e in [f[0], f[n - 1]] {
        let c = code(e);
        if c == 0 || c == 3 {
            h += 2.2;
            s += 6.9;
        }
    }
    if is_self_comp(f) {
        s += -1.4;
    }
    s += 0.368 * (n as f64 - 1.0) * go_log(na);
    (h, s)
}

/// Tm (C) and dG37 (kcal/mol) of an oligo against its exact complement (== `mfeprimer thermo`).
#[allow(dead_code)]
pub fn oligo(f: &[u8], c: Cond) -> (f64, f64) {
    let (h, s) = perfect_hs(f, c.na_eq());
    let x = if is_self_comp(f) { 1.0 } else { 4.0 };
    let tm = 1000.0 * h / (s + GAS_R * go_log(c.oligo * 1e-9 / x)) - 273.15;
    (tm, h - T37 * s / 1000.0)
}

/// oligo() as printed in the oligo table of the dimer and hairpin reports, which use R = 1.9872.
pub fn oligo_report(f: &[u8], c: Cond) -> (f64, f64) {
    let (h, s) = perfect_hs(f, c.na_eq());
    let x = if is_self_comp(f) { 1.0 } else { 4.0 };
    let tm = 1000.0 * h / (s + HP_GAS_R * go_log(c.oligo * 1e-9 / x)) - 273.15;
    (tm, h - T37 * s / 1000.0)
}

/// The GC percentage of seq.
pub fn gc(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let n = seq.iter().filter(|&&b| matches!(b, b'G' | b'C' | b'g' | b'c')).count();
    100.0 * n as f64 / seq.len() as f64
}

// ---- stacks ----

/// A stack 5'-t0 t1-3' / 3'-b0 b1-5' by codes (4 = N): N anywhere gives nothing (the stack still counts for salt),
/// a perfect stack the NN values, otherwise the mismatch-stack table.
#[inline]
fn t_stack(t0: i32, t1: i32, b0: i32, b1: i32) -> (f64, f64) {
    if t0 > 3 || t1 > 3 || b0 > 3 || b1 > 3 {
        return (0.0, 0.0);
    }
    let t = tables();
    if b0 == 3 - t0 && b1 == 3 - t1 {
        return (t.nn_h[t0 as usize][t1 as usize], t.nn_s[t0 as usize][t1 as usize]);
    }
    let v = t.mm[t0 as usize][t1 as usize][b0 as usize][b1 as usize];
    (v[0], v[1])
}

/// The stack 5'-t0 t1-3' / 3'-b0 b1-5' by codes, inosine included. A stack the paper has no value for (I next to a
/// mismatch or an N, two inosines in a row on one strand) adds nothing but still counts for salt.
#[inline]
fn stack(t0: i32, t1: i32, b0: i32, b1: i32) -> (f64, f64) {
    if t0 == INO || t1 == INO || b0 == INO || b1 == INO {
        if t0 == 4 || t1 == 4 || b0 == 4 || b1 == 4 {
            return (0.0, 0.0);
        }
        let v = tables().ino.get(&[t0, t1, b0, b1]).copied().unwrap_or([0.0, 0.0]);
        return (v[0], v[1]);
    }
    t_stack(t0, t1, b0, b1)
}

/// perfect_hs for an oligo holding inosine, against its most stable complement (I·C).
fn perfect_hs_ino(f: &[u8], na: f64) -> (f64, f64) {
    let n = f.len();
    let (mut h, mut s) = (0.2, -5.7);
    for i in 0..n.saturating_sub(1) {
        let (a, b) = (ocode(f[i]), ocode(f[i + 1]));
        let (dh, ds) = stack(a, b, ino_partner(a), ino_partner(b));
        h += dh;
        s += ds;
    }
    for e in [f[0], f[n - 1]] {
        let c = code(e);
        if c == 0 || c == 3 {
            h += 2.2;
            s += 6.9;
        }
    }
    s += 0.368 * (n as f64 - 1.0) * go_log(na);
    (h, s)
}

fn ino_partner(c: i32) -> i32 {
    if c == INO {
        1 // C
    } else {
        3 - c
    }
}

// ---- primer binding and probe search ----

/// The oligo base pairs with the site base (inosine with any of A/C/G/T, never with N).
#[inline]
fn col_match(f: u8, s: i32) -> bool {
    let o = ocode(f);
    o == s || (o == INO && (0..=3).contains(&s))
}

/// The number of 5' columns mfeprimer drops: a mismatch run (N counts) touching the 5' end, or the 3'-most run of
/// three or more. The duplex then starts at the last column of that run.
fn bind_cut(f: &[u8], s: &[i32]) -> usize {
    let n = f.len();
    let mut d = 0;
    let mut i = 0;
    while i < n {
        if col_match(f[i], s[i]) {
            i += 1;
            continue;
        }
        let mut j = i;
        while j < n && !col_match(f[j], s[j]) {
            j += 1;
        }
        if j - i >= 3 || i == 0 {
            d = j;
        }
        i = j;
    }
    d
}

/// Tm and dG of primer f bound to site s (same length, both 5'->3' in the primer's orientation, i.e. a perfect site
/// has s == f). l5 is the template-side base beyond the primer's 5' end and r3 the one beyond its 3' end (both in the
/// primer's orientation, '-' at the genome end). rev is true when the primer binds the plus strand.
///
/// Measured rules: the dropped 5' columns (bind_cut) leave a blunt end with no terminal term; otherwise each end gets
/// the dangling template base (half a salt unit; an N flank gives only that), or a terminal-AT term at the genome end.
/// A forward primer reads its 5' dangling base uncomplemented and takes its own 3' base as the 3' dangling base; a
/// reverse primer is read on the plus strand, which turns every mismatch into the complementary one.
pub fn binding(f: &[u8], s: &[u8], l5: u8, r3: u8, rev: bool, c: Cond) -> Option<(f64, f64)> {
    let n0 = f.len();
    if n0 != s.len() || n0 < 2 {
        return None;
    }
    let mut buf = [0i32; 128];
    let mut heap;
    let sc: &mut [i32] = if n0 <= buf.len() {
        &mut buf[..n0]
    } else {
        heap = vec![0i32; n0];
        &mut heap[..]
    };
    for i in 0..n0 {
        let x = site_code(s[i]);
        if ocode(f[i]) < 0 || x < 0 {
            return None;
        }
        sc[i] = x;
    }
    let (lc, rc) = (flank_code(l5), flank_code(r3));
    if lc < 0 || rc < 0 {
        return None;
    }
    let d = bind_cut(f, sc);
    let (f, sc) = if d > 0 { (&f[d - 1..], &sc[d - 1..]) } else { (f, &sc[..]) };
    let n = f.len();
    let t = tables();
    let (mut h, mut sd, mut units) = (0.2, -5.7, 0.0);
    for x in 0..n - 1 {
        let y = n - 2 - x;
        let (dh, ds) = if rev && ocode(f[y]) != INO && ocode(f[y + 1]) != INO {
            // rc(primer) over the plus-strand site
            t_stack(3 - code(f[y + 1]), 3 - code(f[y]), sc[y + 1], sc[y])
        } else if rev {
            // a stack holding inosine is read as the physical duplex
            stack(ocode(f[y]), ocode(f[y + 1]), comp_code(sc[y]), comp_code(sc[y + 1]))
        } else {
            stack(ocode(f[x]), ocode(f[x + 1]), comp_code(sc[x]), comp_code(sc[x + 1]))
        };
        h += dh;
        sd += ds;
        units += 1.0;
    }
    let (f0, fl) = (ocode(f[0]), ocode(f[n - 1]));
    if d > 0 {
        // blunt cut end
    } else if f0 == INO {
        // no dangle values for an inosine end; terminal-AT does not apply
        if lc <= 4 {
            units += 0.5;
        }
    } else if lc <= 3 {
        let v = if rev {
            t.dangle[1][(3 - f0) as usize][f0 as usize][(3 - lc) as usize]
        } else {
            t.dangle[0][(3 - f0) as usize][f0 as usize][lc as usize]
        };
        h += v[0];
        sd += v[1];
        units += 0.5;
    } else if lc == 4 {
        units += 0.5;
    } else if f0 == 0 || f0 == 3 {
        h += 2.2;
        sd += 6.9;
    }
    if fl == INO {
        if !rev || rc <= 4 {
            units += 0.5;
        }
    } else if !rev {
        let v = t.dangle[1][fl as usize][(3 - fl) as usize][fl as usize];
        h += v[0];
        sd += v[1];
        units += 0.5;
    } else if rc <= 3 {
        let v = t.dangle[0][fl as usize][(3 - fl) as usize][(3 - rc) as usize];
        h += v[0];
        sd += v[1];
        units += 0.5;
    } else if rc == 4 {
        units += 0.5;
    } else if fl == 0 || fl == 3 {
        // not seen in the data; by symmetry with the 5' end
        h += 2.2;
        sd += 6.9;
    }
    sd += 0.368 * units * go_log(c.na_eq());
    Some((1000.0 * h / (sd + HP_GAS_R * go_log(c.oligo * 1e-9 / 4.0)) - 273.15, h - T37 * sd / 1000.0))
}

/// Tm and dG the way `mfeprimer search` scores a hit of the same length as the query: the perfect duplex of q with
/// every stack that holds a mismatch replaced by its table value (read as written), no end or flank terms and no
/// truncation of mismatch runs. hit may hold N / ambiguity codes (scored as N).
pub fn search(q: &[u8], hit: &[u8], c: Cond) -> Option<(f64, f64)> {
    let n = q.len();
    if n != hit.len() || n < 2 {
        return None;
    }
    for i in 0..n {
        if ocode(q[i]) < 0 || site_code(hit[i]) < 0 {
            return None;
        }
    }
    if has_inosine(q) {
        return Some(search_ino(q, hit, c));
    }
    let t = tables();
    let (mut h, mut sd) = perfect_hs(q, c.na_eq());
    for i in 0..n - 1 {
        let (f1, f2, s1, s2) = (code(q[i]), code(q[i + 1]), site_code(hit[i]), site_code(hit[i + 1]));
        if f1 == s1 && f2 == s2 {
            continue;
        }
        let (mh, ms) = t_stack(f1, f2, comp_code(s1), comp_code(s2));
        h += mh - t.nn_h[f1 as usize][f2 as usize];
        sd += ms - t.nn_s[f1 as usize][f2 as usize];
    }
    let x = if is_self_comp(q) { 1.0 } else { 4.0 };
    Some((1000.0 * h / (sd + GAS_R * go_log(c.oligo * 1e-9 / x)) - 273.15, h - T37 * sd / 1000.0))
}

/// search() for a query holding inosine: every stack read from the actual pairs, with the same initiation,
/// terminal-AT and salt terms as the perfect duplex.
fn search_ino(q: &[u8], hit: &[u8], c: Cond) -> (f64, f64) {
    let n = q.len();
    let (mut h, mut sd) = (0.2, -5.7);
    for i in 0..n - 1 {
        let (dh, ds) =
            stack(ocode(q[i]), ocode(q[i + 1]), comp_code(site_code(hit[i])), comp_code(site_code(hit[i + 1])));
        h += dh;
        sd += ds;
    }
    for e in [q[0], q[n - 1]] {
        let x = code(e);
        if x == 0 || x == 3 {
            h += 2.2;
            sd += 6.9;
        }
    }
    sd += 0.368 * (n as f64 - 1.0) * go_log(c.na_eq());
    (1000.0 * h / (sd + GAS_R * go_log(c.oligo * 1e-9 / 4.0)) - 273.15, h - T37 * sd / 1000.0)
}

// ---- dimer and hairpin ----

/// The stack 5'-t0 t1-3' / 3'-b0 b1-5' by bases: a nearest-neighbour stack when both columns pair, else a mismatch
/// stack (zero for two mismatches in a row, bar the few mfeprimer has values for).
#[inline]
fn hp_stack(t0: u8, t1: u8, b0: u8, b1: u8) -> (f64, f64) {
    let (a, b, c, d) = (ocode(t0), ocode(t1), ocode(b0), ocode(b1));
    if a == INO || b == INO || c == INO || d == INO {
        return stack(a, b, c, d);
    }
    let t = tables();
    if c == 3 - a && d == 3 - b {
        return (t.nn_h[a as usize][b as usize], t.nn_s[a as usize][b as usize]);
    }
    let v = t.mm[a as usize][b as usize][c as usize][d as usize];
    (v[0], v[1])
}

fn hp_loop_term(l: i64) -> f64 {
    let t = tables();
    if l <= 10 {
        for (i, &x) in t.loop_l.iter().enumerate() {
            if x == l {
                return t.loop_g[i];
            }
        }
        return 0.0;
    }
    for (i, &x) in t.loop_l.iter().enumerate() {
        // longer loops take the next listed size up; none past the table
        if x > l {
            return t.loop_g[i];
        }
    }
    0.0
}

/// dG37 (kcal/mol) and Tm (C) of the hairpin of s with the given pairs (outermost first, (i, j) with i < j, all on
/// one diagonal i+j apart from the gaps between them).
pub fn hairpin_energy(s: &[u8], pairs: &[[usize; 2]], c: Cond) -> (f64, f64) {
    let n = s.len();
    let (mut h, mut sd, mut units) = (0.2, -5.7, 0.0);
    for k in 0..pairs.len() - 1 {
        let (i, j, i2) = (pairs[k][0], pairs[k][1], pairs[k + 1][0]);
        for x in 0..i2 - i {
            // mismatched columns between two pairs are read as mismatch stacks
            let (dh, ds) = hp_stack(s[i + x], s[i + x + 1], s[j - x], s[j - x - 1]);
            h += dh;
            sd += ds;
            units += 1.0;
        }
    }
    let (i, j) = (pairs[0][0], pairs[0][1]);
    if i > 0 && j + 1 < n {
        // both flanks: a terminal mismatch stack
        let (dh, ds) = hp_stack(s[i - 1], s[i], s[j + 1], s[j]);
        h += dh;
        sd += ds;
        units += 1.0;
    } else if j + 1 < n {
        let v = dangle(0, s[j], s[i], s[j + 1]);
        h += v[0];
        sd += v[1];
        units += 0.5;
    } else if i > 0 {
        let v = dangle(1, s[j], s[i], s[i - 1]);
        h += v[0];
        sd += v[1];
        units += 0.5;
    }
    let last = pairs[pairs.len() - 1];
    let (i, j) = (last[0], last[1]);
    // the closing stack is read with the closing 5' base twice (the loop base after it is never looked at)
    let (dh, ds) = hp_stack(s[i], s[i], s[j], s[j - 1]);
    h += dh;
    sd += ds;
    units += 1.0;
    sd += 0.368 * units * go_log(c.na_eq());
    let bonus = tables().bonus.get(&s[i..=j]).copied().unwrap_or(0.0);
    let dg = h - T37 * sd / 1000.0 + hp_loop_term(j as i64 - i as i64 - 1) + bonus;
    let tm = 1000.0 * h / (sd + HP_GAS_R * go_log(c.oligo * 1e-9 / 4.0)) - 273.15;
    (dg, tm)
}

/// dG37 (kcal/mol) and Tm (C) of oligo a paired with oligo b over a's columns i..j, b written reversed (br = b
/// 3'->5') and shifted so that br column c-shift faces a column c; both span ends are pairs. A duplex whose two
/// strands read the same over the span and whose overhanging neighbours mirror (always so for an oligo with itself)
/// gets the symmetry term and its Tm uses Ct instead of Ct/4. Same tables as the hairpin.
pub fn dimer_dg(a: &[u8], br: &[u8], shift: isize, i: usize, j: usize, c: Cond) -> Option<(f64, f64)> {
    let at = |seq: &[u8], k: isize| -> u8 {
        if k < 0 || k >= seq.len() as isize {
            0
        } else {
            seq[k as usize]
        }
    };
    let (ii, jj) = (i as isize, j as isize);
    for x in ii - 1..=jj + 1 {
        if x >= 0 && x < a.len() as isize && ocode(a[x as usize]) < 0 {
            return None;
        }
        let y = x - shift;
        if y >= 0 && y < br.len() as isize && ocode(br[y as usize]) < 0 {
            return None;
        }
    }
    let b = |k: isize| br[(k - shift) as usize];
    let (mut h, mut sd, mut units) = (0.2, -5.7, 0.0);
    for x in ii..jj {
        let (dh, ds) = hp_stack(a[x as usize], a[x as usize + 1], b(x), b(x + 1));
        h += dh;
        sd += ds;
        units += 1.0;
    }
    // left end: 5'-aL a[i]-3' / 3'-bL br[i]-5'
    let (al, bl) = (at(a, ii - 1), at(br, ii - 1 - shift));
    if al != 0 && bl != 0 {
        let (dh, ds) = hp_stack(al, a[i], bl, b(ii));
        h += dh;
        sd += ds;
        units += 1.0;
    } else if al != 0 {
        let v = dangle(1, b(ii), a[i], al);
        h += v[0];
        sd += v[1];
        units += 0.5;
    } else if bl != 0 {
        let v = dangle(0, b(ii), a[i], bl);
        h += v[0];
        sd += v[1];
        units += 0.5;
    }
    // right end: 5'-a[j] aR-3' / 3'-br[j] bR-5'
    let (ar, brr) = (at(a, jj + 1), at(br, jj + 1 - shift));
    if ar != 0 && brr != 0 {
        let (dh, ds) = hp_stack(a[j], ar, b(jj), brr);
        h += dh;
        sd += ds;
        units += 1.0;
    } else if brr != 0 {
        let v = dangle(1, a[j], b(jj), brr);
        h += v[0];
        sd += v[1];
        units += 0.5;
    } else if ar != 0 {
        let v = dangle(0, a[j], b(jj), ar);
        h += v[0];
        sd += v[1];
        units += 0.5;
    }
    let mut x = 4.0;
    // the overhanging neighbours mirror each other too; top a[i..j] against the bottom strand read 5'->3'
    let sym = al == brr && bl == ar && (ii..=jj).all(|k| a[k as usize] == b(ii + jj - k));
    if sym {
        sd += -1.4;
        x = 1.0;
    }
    sd += 0.368 * units * go_log(c.na_eq());
    Some((h - T37 * sd / 1000.0, 1000.0 * h / (sd + HP_GAS_R * go_log(c.oligo * 1e-9 / x)) - 273.15))
}

/// The dangle table by bases; an inosine in it has no published value and adds nothing (its half salt unit is still
/// counted by the caller).
fn dangle(half: usize, p3: u8, p5: u8, d: u8) -> [f64; 2] {
    let (a, b, c) = (ocode(p3), ocode(p5), ocode(d));
    if a < 0 || b < 0 || c < 0 || a == INO || b == INO || c == INO {
        return [0.0, 0.0];
    }
    tables().dangle[half][a as usize][b as usize][c as usize]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn watkins_transcription() {
        let rows: HashMap<&str, [f64; 2]> = WATKINS_2005.iter().copied().collect();
        for (k, dg) in WATKINS_DG37 {
            let v = rows[k];
            let calc = v[0] - T37 * v[1] / 1000.0;
            assert!((calc - dg).abs() <= 0.051, "{k}: {calc} vs {dg}");
        }
    }

    fn bind_tm(f: &str, s: &str, rev: bool) -> f64 {
        binding(f.as_bytes(), s.as_bytes(), b'A', b'C', rev, DEFAULT).expect("binding").0
    }

    /// An internal inosine pairs best with C, then A, then T and G (roughly equal), and always more weakly than a
    /// Watson-Crick G·C. The site is written in the primer's orientation, so the template base facing the I is the
    /// complement of the site letter: I·C is site G, I·A site T, I·T site A, I·G site C.
    #[test]
    fn inosine_binding_order() {
        let at = |x: &str| format!("ACGTTGCA{x}TGCATCGAGTC");
        for rev in [false, true] {
            let wc = bind_tm(&at("G"), &at("G"), rev);
            let ic = bind_tm(&at("I"), &at("G"), rev);
            let ia = bind_tm(&at("I"), &at("T"), rev);
            let it = bind_tm(&at("I"), &at("A"), rev);
            let ig = bind_tm(&at("I"), &at("C"), rev);
            assert!(wc > ic && ic > ia && ia > it && ia > ig, "rev={rev}: {wc} {ic} {ia} {it} {ig}");
        }
    }

    /// Probe-search scoring of a query with inosine against a site with C equals the oligo's own (I·C) duplex.
    #[test]
    fn search_inosine_self() {
        let q = b"ACGTTGCAITGCATCGAGTC";
        let (tm, _) = search(q, b"ACGTTGCAGTGCATCGAGTC", DEFAULT).unwrap();
        let (own, _) = oligo(q, DEFAULT);
        assert!((tm - own).abs() <= 0.3, "search Tm {tm} vs oligo Tm {own}");
    }

    // ---- against MFEprimer reference runs in work/ (git-ignored; skipped when absent) ----

    use crate::json::{self, glob, work_dir};

    fn read_tsv(path: &std::path::Path) -> Option<Vec<Vec<String>>> {
        let text = std::fs::read_to_string(path).ok()?;
        Some(
            text.lines()
                .filter(|l| !l.is_empty() && !l.starts_with('#'))
                .map(|l| l.split('\t').map(String::from).collect())
                .collect(),
        )
    }

    /// `mfeprimer thermo` on random oligos under several buffers (printed to 2 decimals).
    #[test]
    fn oligo_matches_mfeprimer_thermo() {
        let conds = [
            ("rand.thermo.tsv", DEFAULT),
            ("r_mono50diva0dntp0.thermo.tsv", Cond { mono: 50.0, diva: 0.0, dntp: 0.0, oligo: 50.0 }),
            ("r_mono100diva0dntp0.thermo.tsv", Cond { mono: 100.0, diva: 0.0, dntp: 0.0, oligo: 50.0 }),
            ("r_mono50diva3dntp0.thermo.tsv", Cond { mono: 50.0, diva: 3.0, dntp: 0.0, oligo: 50.0 }),
            ("r_mono50diva1.5dntp0.25oligo200.thermo.tsv", Cond { mono: 50.0, diva: 1.5, dntp: 0.25, oligo: 200.0 }),
        ];
        for (file, c) in conds {
            let Some(rows) = read_tsv(&work_dir().join(file)) else {
                eprintln!("skip: {file} not available");
                return;
            };
            let (mut etm, mut edg) = (0f64, 0f64);
            for r in &rows {
                let (tm, dg) = oligo(r[1].as_bytes(), c);
                etm = etm.max((tm - r[4].parse::<f64>().unwrap()).abs());
                edg = edg.max((dg - r[5].parse::<f64>().unwrap()).abs());
            }
            assert!(etm <= 0.011 && edg <= 0.011, "{file}: max |dTm| {etm:.4} |ddG| {edg:.4}");
        }
    }

    /// Primer bindings from MFEprimer's spec.tsv (harness/synth.py experiments; printed to 2 decimals).
    #[test]
    fn binding_matches_spec() {
        for exp in ["exp_holdout", "exp_R_rand", "exp_tandem", "exp_R_tandem", "exp_t5"] {
            let Some(rows) = read_tsv(&work_dir().join(exp).join("joined.tsv")) else {
                eprintln!("skip: {exp} not available");
                return;
            };
            let col = |h: &str| rows[0].iter().position(|x| x == h).unwrap();
            // a forward-only experiment has no "side" column
            let side = rows[0].iter().position(|x| x == "side");
            let (found, fp, site, l5, r3, ftm, fdg) =
                (col("found"), col("fpSeq"), col("site"), col("l5"), col("r3"), col("fpTm"), col("fpDg"));
            let (mut etm, mut edg) = (0f64, 0f64);
            for r in &rows[1..] {
                if r[found] != "1" {
                    continue;
                }
                let (tm, dg) = binding(r[fp].as_bytes(), r[site].as_bytes(), r[l5].as_bytes()[0], r[r3].as_bytes()[0], side.is_some_and(|c| r[c] == "R"), DEFAULT)
                    .unwrap_or_else(|| panic!("{exp}: binding not computable for {r:?}"));
                etm = etm.max((tm - r[ftm].parse::<f64>().unwrap()).abs());
                edg = edg.max((dg - r[fdg].parse::<f64>().unwrap()).abs());
            }
            assert!(etm <= 0.03 && edg <= 0.015, "{exp}: max |dTm| {etm:.4} |ddG| {edg:.4}");
        }
    }

    /// work/bind_all.json: full-precision bindings from `mfeprimer -j` reruns (harness/bind_json.py).
    #[test]
    fn binding_matches_mfeprimer_json() {
        let Ok(text) = std::fs::read_to_string(work_dir().join("bind_all.json")) else {
            eprintln!("skip: bind_all.json not available");
            return;
        };
        let rows = json::parse(&text).unwrap();
        let mut bad = 0;
        for r in rows.as_array() {
            let (f, s) = (r.get("F").as_str(), r.get("S").as_str());
            let (l5, r3) = (r.get("L5").as_str().as_bytes()[0], r.get("R3").as_str().as_bytes()[0]);
            let got = binding(f.as_bytes(), s.as_bytes(), l5, r3, r.get("Side").as_str() == "R", DEFAULT);
            match got {
                Some((tm, dg)) if (dg - r.get("Dg").as_f64()).abs() <= 1e-9 && (tm - r.get("Tm").as_f64()).abs() <= 1e-6 => {}
                _ => bad += 1,
            }
        }
        // one record sits in a poly-G run where the harness cannot place the site unambiguously
        assert!(bad <= 1, "{bad} of {} bindings differ", rows.as_array().len());
    }

    /// `mfeprimer dimer -j` on harness dimer sets (work/dh/dj*.json): full-precision dG and Tm.
    #[test]
    fn dimer_dg_matches_mfeprimer() {
        let files = glob(&work_dir().join("dh"), |n| n.starts_with("dj") && n.ends_with(".json"));
        if files.is_empty() {
            eprintln!("skip: work/dh/dj*.json not available");
            return;
        }
        let mut n = 0;
        for f in files {
            let d = json::parse(&std::fs::read_to_string(&f).unwrap()).unwrap();
            for x in d.get("DimerList").as_array() {
                let a = x.get("S1").get("Seq").as_str().as_bytes();
                let br = x.get("S2").get("Seq").as_str().as_bytes(); // S2 is already written 3'->5'
                let (s1d, ad, s2d) =
                    (x.get("S1Dangling").as_f64() as isize, x.get("AseqDangling").as_f64() as isize, x.get("S2Dangling").as_f64() as isize);
                let (shift, off) = (s2d - s1d, ad - s1d);
                let cols: Vec<isize> =
                    x.get("Aseq").as_str().chars().enumerate().filter(|(_, c)| *c == ':' || *c == '.').map(|(k, _)| off + k as isize).collect();
                let (i, j) = (cols[0] as usize, *cols.last().unwrap() as usize);
                let (dg, tm) = dimer_dg(a, br, shift, i, j, DEFAULT).expect("dimer");
                assert!(
                    (dg - x.get("Dg").as_f64()).abs() <= 1e-9 && (tm - x.get("Tm").as_f64()).abs() <= 1e-6,
                    "{f:?}: shift {shift} [{i},{j}]: dG {dg} Tm {tm}"
                );
                n += 1;
            }
        }
        assert!(n > 0);
    }

    /// `mfeprimer hairpin -j` on synth_hairpin.py sets (work/hpx/*/rows.json): full-precision dG and Tm.
    #[test]
    fn hairpin_energy_matches_mfeprimer() {
        let dirs = glob(&work_dir().join("hpx"), |_| true);
        let mut n = 0;
        for d in dirs {
            let Ok(text) = std::fs::read_to_string(d.join("rows.json")) else { continue };
            for r in json::parse(&text).unwrap().as_array() {
                if r.get("Score").is_null() {
                    continue;
                }
                let st = r.get("Struct").as_str();
                let op: Vec<usize> = st.char_indices().filter(|c| c.1 == '/').map(|c| c.0).collect();
                let cl: Vec<usize> = st.char_indices().filter(|c| c.1 == '\\').map(|c| c.0).collect();
                let pairs: Vec<[usize; 2]> = (0..op.len()).map(|k| [op[k], cl[cl.len() - 1 - k]]).collect();
                let (dg, tm) = hairpin_energy(r.get("Seq").as_str().as_bytes(), &pairs, DEFAULT);
                assert!(
                    (dg - r.get("Dg").as_f64()).abs() <= 1e-9 && (tm - r.get("Tm").as_f64()).abs() <= 1e-6,
                    "{d:?}: {} {st}: dG {dg} Tm {tm}",
                    r.get("Seq").as_str()
                );
                n += 1;
            }
        }
        if n == 0 {
            eprintln!("skip: work/hpx/*/rows.json not available");
        }
    }
}
