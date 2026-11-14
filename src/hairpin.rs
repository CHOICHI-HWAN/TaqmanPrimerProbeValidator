//! `mfeprimer hairpin`: the one hairpin reported for each oligo, and the report.

use crate::dimer::{now, pair_ok, write_cond, write_oligo_table, Expanded};
use crate::thermo::{self, Cond};
use std::fs::File;
use std::io::{self, BufWriter, Write};

/// `mfeprimer hairpin`'s settings.
#[derive(Clone, Copy)]
pub struct Options {
    /// -s: minimum score
    pub score: i64,
    /// -d: maximum dG
    pub dg: f64,
    /// -t: minimum Tm
    pub tm: f64,
    /// --mindh (accepted; no effect found)
    pub min_dh: i64,
    /// --minl, --maxl (only the default maxl is verified)
    pub min_l: i64,
    pub max_l: i64,
    /// --maxgap
    pub max_gap: i64,
}

pub const DEFAULT: Options = Options { score: 5, dg: 0.0, tm: 0.0, min_dh: 2, min_l: 3, max_l: 12, max_gap: 2 };

/// The one hairpin reported for an oligo.
pub struct Hairpin {
    pub oligo: usize,
    /// outermost first
    pub pairs: Vec<[usize; 2]>,
    pub score: i64,
    pub dg: f64,
    pub tm: f64,
}

struct Cand {
    pairs: Vec<[usize; 2]>,
    score: i64,
    /// column at which it was found, diagonal
    at: i64,
    k: i64,
}

/// Walks diagonal kk (column c pairs s[c] with s[kk-c]) from its outer end towards the loop while at least minl bases
/// would stay unpaired. Mismatched columns are gaps; when a segment collects more than maxgap of them it is reported
/// at that column with its end gaps dropped and scored by its pairs alone, and a new segment starts after it. The
/// segment that reaches the loop keeps its leading gaps, drops the loop-side pieces (gaps plus the helix after them)
/// that do not raise its score, and scores pairs - gaps.
fn hp_walk(s: &[u8], kk: i64, o: &Options, out: &mut Vec<Cand>) {
    let n = s.len() as i64;
    let mut c = 0i64;
    if kk - n + 1 > 0 {
        c = kk - n + 1;
    }
    let emit = |out: &mut Vec<Cand>, mut start: i64, mut cols: &[bool], overflow: bool, at: i64| {
        while let [rest @ .., false] = cols {
            cols = rest;
        }
        if overflow {
            while let [false, rest @ ..] = cols {
                cols = rest;
                start += 1;
            }
        } else {
            loop {
                let mut k = cols.len();
                while k > 0 && cols[k - 1] {
                    k -= 1;
                }
                if k == 0 {
                    break;
                }
                let mut g = k;
                while g > 0 && !cols[g - 1] {
                    g -= 1;
                }
                if (cols.len() - k) as i64 - (k - g) as i64 > 0 {
                    break;
                }
                cols = &cols[..g];
            }
        }
        let mut pairs = Vec::new();
        let mut gaps = 0;
        for (x, &ok) in cols.iter().enumerate() {
            if ok {
                let p = start + x as i64;
                pairs.push([p as usize, (kk - p) as usize]);
            } else {
                gaps += 1;
            }
        }
        if pairs.is_empty() {
            return;
        }
        let mut score = pairs.len() as i64;
        if !overflow {
            score -= gaps;
        }
        out.push(Cand { pairs, score, at, k: kk });
    };
    let mut start = c;
    let mut cols: Vec<bool> = Vec::new();
    let mut gaps = 0;
    while (kk - c) - c - 1 >= o.min_l {
        let ok = pair_ok(s[c as usize], s[(kk - c) as usize]);
        cols.push(ok);
        if !ok {
            gaps += 1;
        }
        if gaps > o.max_gap {
            emit(out, start, &cols[..cols.len() - 1], true, c);
            start = c + 1;
            cols.clear();
            gaps = 0;
        }
        c += 1;
    }
    emit(out, start, &cols, false, c - 1);
}

/// The hairpin mfeprimer reports for s (None when none passes the cutoffs): the first candidate with the highest
/// score, candidates taken in the order they are found (by column, then diagonal); the cutoffs are applied to that
/// one only.
pub fn find_hairpin(s: &[u8], o: &Options, c: Cond) -> Option<(Vec<[usize; 2]>, i64, f64, f64)> {
    let mut cands = Vec::new();
    let n = s.len() as i64;
    for kk in 0..=2 * n - 2 {
        hp_walk(s, kk, o, &mut cands);
    }
    cands.sort_by_key(|x| (x.at, x.k)); // stable
    let mut best: Option<usize> = None;
    for x in 0..cands.len() {
        if best.is_none_or(|b| cands[x].score > cands[b].score) {
            best = Some(x);
        }
    }
    let b = best?;
    if cands[b].score < o.score {
        return None;
    }
    let b = cands.swap_remove(b);
    let (dg, tm) = thermo::hairpin_energy(s, &b.pairs, c);
    if dg > o.dg || tm < o.tm {
        return None;
    }
    Some((b.pairs, b.score, dg, tm))
}

/// The hairpin of every expanded oligo, sorted by dG (most stable first).
pub fn hairpins(ex: &[Expanded], o: &Options, c: Cond) -> Vec<Hairpin> {
    let mut out: Vec<Hairpin> = ex
        .iter()
        .enumerate()
        .filter_map(|(i, e)| {
            find_hairpin(&e.seq, o, c).map(|(pairs, score, dg, tm)| Hairpin { oligo: i, pairs, score, dg, tm })
        })
        .collect();
    out.sort_by(|a, b| a.dg.partial_cmp(&b.dg).unwrap_or(std::cmp::Ordering::Equal));
    out
}

/// The hairpin drawn as mfeprimer does: '/' for the 5' arm, '\' for the 3' arm, '-' unpaired.
fn structure(h: &Hairpin, n: usize) -> String {
    let mut b = vec![b'-'; n];
    for p in &h.pairs {
        b[p[0]] = b'/';
        b[p[1]] = b'\\';
    }
    String::from_utf8(b).unwrap()
}

/// Writes a `mfeprimer hairpin`-style report.
pub fn write_hairpins(path: &str, input: &str, ex: &[Expanded], hs: &[Hairpin], o: &Options, c: Cond) -> io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    write!(w, "MFEprimer-3.0 Hairpin Reports ({})\n\n", now())?;
    write_oligo_table(&mut w, ex, c)?;
    write!(w, "\n\n\nHairpin List ({})\n\n", hs.len())?;
    if hs.is_empty() {
        write!(w, "No hairpins found.\n\n\n")?;
    } else {
        write!(w, "\n")?;
        for (k, h) in hs.iter().enumerate() {
            if k > 0 {
                write!(w, "\n\n\n")?;
            }
            let e = &ex[h.oligo];
            write!(
                w,
                "Hairpin {}: {}\n\n    Score: {}, Tm = {:.2} °C, Delta G = {:.2} kcal/mol\n    \n    {}\n    {}",
                k + 1,
                e.name,
                h.score,
                h.tm,
                h.dg,
                structure(h, e.seq.len()),
                String::from_utf8_lossy(&e.seq)
            )?;
        }
        write!(w, "\n\n\n\n\n")?;
    }
    write!(
        w,
        "Parameters\n\n{:>30}: {}\n{:>30}: {}\n{:>30}: {:.2}\n{:>30}: {:.2}\n{:>30}: {}\n{:>30}: {}\n{:>30}: {}\n{:>30}: {}\n",
        "Primer file",
        input,
        "Score cutoff",
        o.score,
        "Delta G cutoff",
        o.dg,
        "Tm cutoff (°C)",
        o.tm,
        "Min double helix size",
        o.min_dh,
        "Min loop size",
        o.min_l,
        "Max loop size",
        o.max_l,
        "Max gap count",
        o.max_gap
    )?;
    write_cond(&mut w, c, "hairpin")?;
    w.flush()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::json::{self, glob, work_dir};

    /// Default-option `mfeprimer hairpin -j` runs on synth_hairpin.py inputs (work/hpx/def_*/, git-ignored): every
    /// oligo's reported hairpin (or none) must be reproduced, bar the known equal-score tie-breaks (FINDINGS §8).
    #[test]
    fn hairpins_match_mfeprimer() {
        let (mut total, mut bad) = (0, 0);
        for d in glob(&work_dir().join("hpx"), |n| n.starts_with("def_")) {
            let (Ok(fa), Ok(js)) = (std::fs::read_to_string(d.join("in.fa")), std::fs::read_to_string(d.join("out.json"))) else {
                continue;
            };
            let js = json::parse(&js).unwrap();
            let want: std::collections::HashMap<&str, &str> = js
                .get("HairpinList")
                .as_array()
                .iter()
                .map(|h| (h.get("Seq").get("ID").as_str(), h.get("Aseq").as_str()))
                .collect();
            let lines: Vec<&str> = fa.trim().lines().collect();
            for k in (0..lines.len().saturating_sub(1)).step_by(2) {
                let (id, seq) = (lines[k].trim_start_matches('>'), lines[k + 1].as_bytes());
                let got = find_hairpin(seq, &DEFAULT, thermo::DEFAULT)
                    .map(|(pairs, score, dg, tm)| structure(&Hairpin { oligo: 0, pairs, score, dg, tm }, seq.len()))
                    .unwrap_or_default();
                total += 1;
                if got != want.get(id).copied().unwrap_or("") {
                    bad += 1;
                }
            }
        }
        if total == 0 {
            eprintln!("skip: work/hpx/def_* not available");
        }
        assert!(bad <= total / 2000, "{bad} of {total} oligos differ (more than the known tie-breaks)");
    }
}
