//! Bases, IUPAC codes and the screening Tm of the design engine (bit-identical to the panels' Python engine).

use crate::golog::go_log;

/// IUPAC bases in the expansion order of the pipeline (mfe_common.IUPAC); empty for other letters.
pub fn iupac(b: u8) -> &'static [u8] {
    match b {
        b'A' => b"A",
        b'C' => b"C",
        b'G' => b"G",
        b'T' => b"T",
        b'R' => b"AG",
        b'Y' => b"CT",
        b'S' => b"GC",
        b'W' => b"AT",
        b'K' => b"GT",
        b'M' => b"AC",
        b'B' => b"CGT",
        b'D' => b"AGT",
        b'H' => b"ACT",
        b'V' => b"ACG",
        b'N' => b"ACGT",
        _ => b"",
    }
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'D' => b'H',
        b'H' => b'D',
        b'B' => b'V',
        b'V' => b'B',
        b'I' => b'I',
        _ => b'N',
    }
}

pub fn bit(b: u8) -> u8 {
    match b {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' => 8,
        _ => 0,
    }
}

/// The two-base IUPAC code of a mask of two bases.
pub fn code_of(mask: u8) -> u8 {
    match mask {
        5 => b'R',
        10 => b'Y',
        6 => b'S',
        9 => b'W',
        12 => b'K',
        3 => b'M',
        _ => b'N',
    }
}

/// Reverse complement (IUPAC-aware; unknown letters become N).
pub fn rev_comp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| complement(b)).collect()
}

pub fn is_acgt(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T')
}

/// Non-ACGT letters.
pub fn n_ambig(s: &[u8]) -> usize {
    s.iter().filter(|&&b| !is_acgt(b)).count()
}

fn i2(b: u8) -> usize {
    match b {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        _ => 3,
    }
}

/// SantaLucia (1998) unified NN parameters [first][second] = (dH, dS), as probe_designer.tm_nn.
const NN98: [[(f64, f64); 4]; 4] = [
    [(-7.9, -22.2), (-8.4, -22.4), (-7.8, -21.0), (-7.2, -20.4)], // AA AC AG AT
    [(-8.5, -22.7), (-8.0, -19.9), (-10.6, -27.2), (-7.8, -21.0)], // CA CC CG CT
    [(-8.2, -22.2), (-9.8, -24.4), (-8.0, -19.9), (-8.4, -22.4)], // GA GC GG GT
    [(-7.2, -21.3), (-8.2, -22.2), (-8.5, -22.7), (-7.9, -22.2)], // TA TC TG TT
];

/// Python 3.12's sum() of floats: Neumaier-compensated after the first term.
pub fn py_sum(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        return 0.0;
    }
    let (mut f, mut c) = (xs[0], 0.0f64);
    for &x in &xs[1..] {
        let t = f + x;
        if f.abs() >= x.abs() {
            c += (f - t) + x;
        } else {
            c += (x - t) + f;
        }
        f = t;
    }
    if c != 0.0 && c.is_finite() {
        f += c;
    }
    f
}

/// Screening Tm of a plain ACGT oligo (250 nM oligo, 50 mM Na+).
pub fn tm_nn(s: &[u8]) -> f64 {
    let ln_na = go_log(50.0 / 1000.0);
    let ln_conc = go_log(250.0 * 1e-9 / 4.0);
    let mut h = [0.0f64; 64];
    let mut sx = [0.0f64; 64];
    let n = s.len().saturating_sub(1).min(64);
    let mut hv: Vec<f64> = Vec::new();
    let mut sv: Vec<f64> = Vec::new();
    let (hs, ss): (&[f64], &[f64]) = if s.len() <= 65 {
        for i in 0..n {
            let v = NN98[i2(s[i])][i2(s[i + 1])];
            h[i] = v.0;
            sx[i] = v.1;
        }
        (&h[..n], &sx[..n])
    } else {
        for i in 0..s.len() - 1 {
            let v = NN98[i2(s[i])][i2(s[i + 1])];
            hv.push(v.0);
            sv.push(v.1);
        }
        (&hv, &sv)
    };
    let (mut dh, mut ds) = (py_sum(hs), py_sum(ss));
    for e in [s[0], s[s.len() - 1]] {
        if e == b'G' || e == b'C' {
            dh += 2.3;
            ds += 4.1;
        } else {
            dh += 0.1;
            ds -= 2.8;
        }
    }
    ds += 0.368 * (s.len() - 1) as f64 * ln_na;
    dh * 1000.0 / (ds + 1.987 * ln_conc) - 273.15
}

pub fn gc_percent(s: &[u8]) -> f64 {
    100.0 * s.iter().filter(|&&b| b == b'G' || b == b'C').count() as f64 / s.len() as f64
}

pub fn longest_run(s: &[u8]) -> usize {
    let (mut best, mut cur) = (1, 1);
    for i in 1..s.len() {
        if s[i] == s[i - 1] {
            cur += 1;
        } else {
            cur = 1;
        }
        best = best.max(cur);
    }
    best
}

/// Every plain oligo an IUPAC oligo stands for, or None when there would be more than limit (or a letter is unknown).
pub fn resolutions(s: &[u8], limit: usize) -> Option<Vec<Vec<u8>>> {
    let mut n = 1usize;
    for &b in s {
        n *= iupac(b).len().max(1);
        if n > limit {
            return None;
        }
    }
    let mut outs: Vec<Vec<u8>> = vec![Vec::new()];
    for &c in s {
        let b = iupac(c);
        if b.is_empty() {
            return None;
        }
        let mut next = Vec::with_capacity(outs.len() * b.len());
        for o in &outs {
            for &x in b {
                let mut v = o.clone();
                v.push(x);
                next.push(v);
            }
        }
        outs = next;
    }
    Some(outs)
}

/// (resolutions <= 64, mean screening Tm, mean GC), as the Python engine's _stats.
pub fn stats(s: &[u8]) -> Option<(Vec<Vec<u8>>, f64, f64)> {
    let res = resolutions(s, 64)?;
    if res.is_empty() {
        return None;
    }
    let ts: Vec<f64> = res.iter().map(|r| tm_nn(r)).collect();
    let gs: Vec<f64> = res.iter().map(|r| gc_percent(r)).collect();
    let n = res.len() as f64;
    let (tm, gc) = (py_sum(&ts) / n, py_sum(&gs) / n);
    Some((res, tm, gc))
}

/// Mean GC over the resolutions (50 when there are too many).
pub fn gc_of(s: &[u8]) -> f64 {
    match resolutions(s, 1 << 12) {
        Some(r) if !r.is_empty() => r.iter().fold(0.0, |t, x| t + gc_percent(x)) / r.len() as f64,
        _ => 50.0,
    }
}
