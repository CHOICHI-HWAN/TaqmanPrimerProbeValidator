//! `mfeprimer dimer`: oligo-oligo duplexes and the report.

use crate::primers::{self, Oligo};
use crate::thermo::{self, Cond};
use std::fs::File;
use std::io::{self, BufWriter, Write};

/// One concrete oligo of the input (degenerate ones expanded, named like the primer search).
#[derive(Clone)]
pub struct Expanded {
    pub name: String,
    pub seq: Vec<u8>,
}

/// Expands the input oligos in order.
pub fn expand_all(oligos: &[Oligo]) -> Result<Vec<Expanded>, String> {
    let mut out = Vec::new();
    for o in oligos {
        let seqs = primers::expand(&o.seq).map_err(|e| format!("{}: {}", o.name, e))?;
        let multi = seqs.len() > 1;
        for (i, s) in seqs.into_iter().enumerate() {
            let name = if multi { format!("{}.{}", o.name, i + 1) } else { o.name.clone() };
            out.push(Expanded { name, seq: s.into_bytes() });
        }
    }
    Ok(out)
}

/// Inosine counts as a pair with C and A (the stabilising pairs); I·T, I·G and I·I are scored as mismatches.
#[inline]
pub fn pair_ok(a: u8, b: u8) -> bool {
    matches!(
        (a, b),
        (b'A', b'T') | (b'T', b'A') | (b'C', b'G') | (b'G', b'C') | (b'I', b'C') | (b'C', b'I') | (b'I', b'A') | (b'A', b'I')
    )
}

/// One oligo-oligo duplex as `mfeprimer dimer` reports it. a is the later oligo of the (expanded) input, b the
/// earlier one (or the same one); b is written reversed under a, shifted by shift columns.
pub struct Dimer {
    pub a: usize,
    pub b: usize,
    pub shift: isize,
    /// paired span on a's columns
    pub i: usize,
    pub j: usize,
    /// pairs - mismatches inside the span
    pub score: i64,
    pub dg: f64,
}

/// For every pair of expanded oligos (and every oligo with itself) and every relative shift, the span from the first
/// to the last complementary column; it is a dimer when it has at most max_mm mismatches, a score (pairs -
/// mismatches) of at least min_score and dG <= max_dg. Sorted by dG.
pub fn dimers(ex: &[Expanded], min_score: i64, max_mm: i64, max_dg: f64, c: Cond) -> Vec<Dimer> {
    let mut out = Vec::new();
    for ai in 0..ex.len() {
        for bi in 0..=ai {
            let a = &ex[ai].seq;
            let br: Vec<u8> = ex[bi].seq.iter().rev().copied().collect();
            let (la, lb) = (a.len() as isize, br.len() as isize);
            for sh in -lb + 1..la {
                let (lo, hi) = (sh.max(0), la.min(sh + lb));
                let (mut first, mut last, mut pairs) = (-1isize, -1isize, 0i64);
                for col in lo..hi {
                    if pair_ok(a[col as usize], br[(col - sh) as usize]) {
                        if first < 0 {
                            first = col;
                        }
                        last = col;
                        pairs += 1;
                    }
                }
                if first < 0 {
                    continue;
                }
                let mm = (last - first + 1) as i64 - pairs;
                if mm > max_mm || pairs - mm < min_score {
                    continue;
                }
                let Some((dg, _)) = thermo::dimer_dg(a, &br, sh, first as usize, last as usize, c) else { continue };
                if dg > max_dg {
                    continue;
                }
                out.push(Dimer { a: ai, b: bi, shift: sh, i: first as usize, j: last as usize, score: pairs - mm, dg });
            }
        }
    }
    out.sort_by(|x, y| x.dg.partial_cmp(&y.dg).unwrap_or(std::cmp::Ordering::Equal)); // stable
    out
}

/// The three-line picture mfeprimer prints for a dimer.
fn align(ex: &[Expanded], d: &Dimer) -> String {
    let a = &ex[d.a].seq;
    let br: Vec<u8> = ex[d.b].seq.iter().rev().copied().collect();
    let (alead, blead) = if d.shift < 0 { ((-d.shift) as usize, 0) } else { (0, d.shift as usize) };
    let mut marks = vec![b' '; alead + 3 + d.j + 1];
    for col in d.i..=d.j {
        marks[alead + 3 + col] = if pair_ok(a[col], br[(col as isize - d.shift) as usize]) { b':' } else { b'.' };
    }
    format!(
        "{}5' {} 3'\n{}\n{}3' {} 5'",
        " ".repeat(alead),
        String::from_utf8_lossy(a),
        String::from_utf8_lossy(&marks),
        " ".repeat(blead),
        String::from_utf8_lossy(&br)
    )
}

pub fn now() -> String {
    chrono::Local::now().format("%Y-%m-%d %H:%M:%S").to_string()
}

pub fn write_oligo_table(w: &mut impl Write, ex: &[Expanded], c: Cond) -> io::Result<()> {
    write!(w, "Primer ID                      Sequence (5'-->3')                    Length     GC      Tm      Dg     \n")?;
    write!(w, "                                                                      (bp)      (%)    (°C)  (kcal/mol)\n\n")?;
    for e in ex {
        let (tm, dg) = thermo::oligo_report(&e.seq, c);
        writeln!(
            w,
            "{:<31}{:<36}{:>8}{:>8.2}{:>8.2}{:>10.2}",
            e.name,
            String::from_utf8_lossy(&e.seq),
            e.seq.len(),
            thermo::gc(&e.seq),
            tm,
            dg
        )?;
    }
    Ok(())
}

pub fn write_cond(w: &mut impl Write, c: Cond, what: &str) -> io::Result<()> {
    write!(
        w,
        "{:>30}: {:.2}\n{:>30}: {:.2}\n{:>30}: {:.2}\n{:>30}: {:.2}\n\nchichipes (mfeprimer-compatible {what} report)\n",
        "Monovalent cations [mM]",
        c.mono,
        "Divalent cations [mM]",
        c.diva,
        "dNTPs [mM]",
        c.dntp,
        "Annealing oligos [nM]",
        c.oligo
    )
}

/// Writes a `mfeprimer dimer`-style report.
#[allow(clippy::too_many_arguments)]
pub fn write_dimers(path: &str, input: &str, ex: &[Expanded], ds: &[Dimer], min_score: i64, max_mm: i64, max_dg: f64, c: Cond) -> io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    write!(w, "MFEprimer-3.0 Dimer Reports ({})\n\n", now())?;
    write_oligo_table(&mut w, ex, c)?;
    write!(w, "\n\n\nDimer List ({})\n\n", ds.len())?;
    if ds.is_empty() {
        write!(w, "No dimer found.\n\n\n\n")?;
    } else {
        write!(w, "\n")?;
    }
    for (k, d) in ds.iter().enumerate() {
        if k > 0 {
            write!(w, "\n\n\n")?;
        }
        write!(
            w,
            "Dimer {}: {} x {}\n\n    Score: {}, Delta G = {:.2} kcal/mol\n\n{}",
            k + 1,
            ex[d.a].name,
            ex[d.b].name,
            d.score,
            d.dg,
            align(ex, d)
        )?;
    }
    if !ds.is_empty() {
        write!(w, "\n\n\n\n\n\n")?;
    }
    write!(
        w,
        "Parameters\n\n{:>30}: {}\n{:>30}: {}\n{:>30}: {:.2}\n{:>30}: {}\n{:>30}: false\n",
        "Primer file", input, "Score cutoff", min_score, "Delta G cutoff", max_dg, "MisMatch", max_mm, "Primer mode"
    )?;
    write_cond(&mut w, c, "dimer")?;
    w.flush()
}
