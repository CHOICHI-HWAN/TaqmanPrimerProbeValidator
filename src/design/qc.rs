//! Reported sets and their check with chichipes' own thermodynamics.

use super::align::Bits;
use super::chem::more_c_than_g;
use super::cover::Params;
use super::oligo::{gc_of, longest_run, n_ambig, resolutions, rev_comp};
use crate::{dimer, hairpin, primers, thermo};

/// A reported primer/probe set with its coverage and check.
#[derive(Clone, Default)]
pub struct Set {
    pub name: String,
    pub f: Vec<u8>,
    pub p: Vec<u8>,
    pub r: Vec<u8>,
    /// shift: which oligos differ from the starting set ("-" for the starting set itself)
    pub change: String,
    pub base_cov: usize,
    pub cov: usize,
    pub added: usize,
    pub amplicon: usize,
    pub fs: usize,
    pub re: usize,
    pub f_tm: f64,
    pub p_tm: f64,
    pub r_tm: f64,
    pub p_tm_label: f64,
    pub f_hp: f64,
    pub p_hp: f64,
    pub r_hp: f64,
    pub dimer3: f64,
    pub dimer_probe: f64,
    pub penalty: f64,
    pub flags: Vec<String>,
    pub off_target: String,
    pub bits: Bits,
}

pub const SET_COLUMNS: &[&str] = &[
    "name", "change", "coverage", "coverage_pct", "base_pct", "added", "amplicon", "F", "P", "R", "F_tm", "P_tm", "P_tm_label", "R_tm",
    "hairpin_dg_F/P/R", "dimer3_dg", "probe_dimer_dg", "penalty", "off_target", "flags",
];

/// Go's fmt "%.1f" / "%.2f" (round half to even on the exact binary value, like Rust's formatter).
fn f1(x: f64) -> String {
    format!("{x:.1}")
}

impl Set {
    pub fn row(&self, n: usize) -> Vec<String> {
        let pct = |c: usize| format!("{:.2}", 100.0 * c as f64 / n.max(1) as f64);
        let s = |v: &[u8]| String::from_utf8_lossy(v).into_owned();
        vec![
            self.name.clone(),
            self.change.clone(),
            self.cov.to_string(),
            pct(self.cov),
            if self.change.is_empty() { String::new() } else { pct(self.base_cov) },
            self.added.to_string(),
            self.amplicon.to_string(),
            s(&self.f),
            s(&self.p),
            s(&self.r),
            f1(self.f_tm),
            f1(self.p_tm),
            f1(self.p_tm_label),
            f1(self.r_tm),
            format!("{}/{}/{}", f1(self.f_hp), f1(self.p_hp), f1(self.r_hp)),
            f1(self.dimer3),
            f1(self.dimer_probe),
            format!("{:.2}", self.penalty),
            self.off_target.clone(),
            self.flags.join("; "),
        ]
    }
}

fn expand(seq: &[u8]) -> Vec<String> {
    primers::expand(&String::from_utf8_lossy(seq)).unwrap_or_default()
}

/// Mean chichipes Tm (SantaLucia & Hicks 2004, spec-search buffer) over the IUPAC resolutions.
fn tm04(seq: &[u8], c: thermo::Cond) -> f64 {
    let res = expand(seq);
    if res.is_empty() {
        return 0.0;
    }
    res.iter().fold(0.0, |t, r| t + thermo::oligo(r.as_bytes(), c).0) / res.len() as f64
}

fn worst_hairpin(seq: &[u8], c: thermo::Cond) -> f64 {
    let mut w = 0.0f64;
    for r in expand(seq) {
        if let Some((_, _, dg, _)) = hairpin::find_hairpin(r.as_bytes(), &hairpin::DEFAULT, c) {
            if dg < w {
                w = dg;
            }
        }
    }
    w
}

/// chichipes Tm, hairpins, dimers and the probe-format rules, plus a varVAMP-like penalty (lower is better). The
/// probe carries a 3' quencher and cannot be extended: only dimers at a primer's 3' end count as extendable.
pub fn qc(s: &mut Set, p: &Params) {
    let c = thermo::DEFAULT;
    s.flags.clear();
    s.f_tm = tm04(&s.f, c);
    s.p_tm = tm04(&s.p, c);
    s.r_tm = tm04(&s.r, c);
    s.p_tm_label = s.p_tm + p.probe.offset;
    s.f_hp = worst_hairpin(&s.f, c);
    s.p_hp = worst_hairpin(&s.p, c);
    s.r_hp = worst_hairpin(&s.r, c);
    let hi = s.f_tm.max(s.r_tm);
    if s.p_tm_label < hi + 2.0 {
        s.flags.push(format!("probe Tm {:.1} < primers+2 ({:.1})", s.p_tm_label, hi + 2.0));
    }
    if (s.f_tm - s.r_tm).abs() > 5.0 {
        s.flags.push(format!("primer Tm differ by {:.1}", (s.f_tm - s.r_tm).abs()));
    }
    if s.p.first() == Some(&b'G') {
        s.flags.push("probe 5' G (quenches the reporter)".into());
    }
    if !p.probe.length_ok(s.p.len()) {
        let (lo, hi) = p.probe.probe_len_range();
        s.flags.push(format!("probe {} nt outside the {} range [{} {}]", s.p.len(), p.probe.name, lo, hi));
    }
    if !more_c_than_g(&s.p) {
        s.flags.push("probe has more G than C (the other strand is preferred)".into());
    }
    s.dimer3 = 0.0;
    s.dimer_probe = 0.0;
    let st = |v: &[u8]| String::from_utf8_lossy(v).into_owned();
    let input = vec![
        primers::Oligo { name: "F".into(), seq: st(&s.f) },
        primers::Oligo { name: "R".into(), seq: st(&s.r) },
        primers::Oligo { name: "P".into(), seq: st(&s.p) },
    ];
    if let Ok(ex) = dimer::expand_all(&input) {
        if ex.len() <= 64 {
            for d in dimer::dimers(&ex, 3, 3, 0.0, c) {
                let (a, b) = (&ex[d.a], &ex[d.b]);
                let (ra, rb) = (&a.name[..1], &b.name[..1]);
                let ext = (ra != "P" && d.j + 1 >= a.seq.len()) || (rb != "P" && d.i as isize <= d.shift && d.shift <= d.j as isize);
                if (ra == "P" || rb == "P") && d.dg < s.dimer_probe {
                    s.dimer_probe = d.dg;
                }
                if ext && d.dg < s.dimer3 {
                    s.dimer3 = d.dg;
                }
            }
        }
    }
    if s.dimer3 < -6.0 {
        s.flags.push(format!("3' dimer {:.1} kcal/mol", s.dimer3));
    }
    let mut pen = (s.f_tm - s.r_tm).abs();
    for x in [&s.f, &s.r, &s.p] {
        let g = gc_of(x);
        if g < 40.0 {
            pen += (40.0 - g) / 5.0;
        } else if g > 60.0 {
            pen += (g - 60.0) / 5.0;
        }
        pen += 0.5 * n_ambig(x) as f64;
    }
    for h in [s.f_hp, s.r_hp, s.p_hp] {
        pen += (-h - 2.0).max(0.0);
    }
    pen += 2.0 * (-s.dimer3 - 5.0).max(0.0) + (-s.dimer_probe - 8.0).max(0.0);
    s.penalty = pen;
}

/// Puts the probe on the strand with more C than G when the other strand also has no 5' G and no 4-base run (a
/// perfect-match probe passes the same genomes on either strand; the screening Tm is strand-symmetric).
pub fn prefer_c_strand(s: &mut Set) {
    if s.p.is_empty() || more_c_than_g(&s.p) {
        return;
    }
    let rc = rev_comp(&s.p);
    if !more_c_than_g(&rc) {
        return;
    }
    let Some(res) = resolutions(&rc, 64) else { return };
    if res.is_empty() || res.iter().any(|r| r[0] == b'G' || longest_run(r) >= 4) {
        return;
    }
    s.p = rc;
}
