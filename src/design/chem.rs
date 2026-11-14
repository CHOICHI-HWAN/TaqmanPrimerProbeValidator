//! Hydrolysis-probe formats (length limits, Tm offset) and reporter / quencher checks.
//!
//! - linear: 5' reporter, 3' quencher (FAM/BHQ1). Quenching depends on the dye-quencher distance: probe <= 30 nt. A
//!   terminal dye or quencher changes duplex stability by about a degree (Moreira et al. 2005): no Tm offset.
//! - zen: double-quenched (internal quencher ~9 nt from the reporter + 3' quencher): probes up to 40 nt, no offset.
//! - mgb: 3' minor-groove binder + non-fluorescent quencher. The MGB raises the duplex Tm strongly (Kutyavin et al.
//!   2000), so MGB probes are 13-22 nt; the default +10 C is a rough figure to calibrate (--probe-tm-offset).
//!
//! Every format: no G at the 5' end (it quenches the reporter), no run of 4+ identical bases.

#[derive(Clone, Debug)]
pub struct Chem {
    pub name: String,
    pub offset: f64,
    pub min_len: usize,
    pub max_len: usize,
}

impl Default for Chem {
    fn default() -> Self {
        Chem { name: "linear".into(), offset: 0.0, min_len: 0, max_len: 30 }
    }
}

/// offset NaN and min/max 0 take the format's defaults.
pub fn parse_chem(name: &str, offset: f64, min_len: usize, max_len: usize) -> Result<Chem, String> {
    let mut c = match name.to_lowercase().as_str() {
        "" | "linear" => Chem { name: "linear".into(), offset: 0.0, min_len: 0, max_len: 30 },
        "zen" | "double" | "double-quenched" => Chem { name: "zen".into(), offset: 0.0, min_len: 0, max_len: 40 },
        "mgb" => Chem { name: "mgb".into(), offset: 10.0, min_len: 13, max_len: 22 },
        _ => return Err(format!("unknown probe format {name:?} (linear, zen, mgb)")),
    };
    if !offset.is_nan() {
        c.offset = offset;
    }
    if min_len > 0 {
        c.min_len = min_len;
    }
    if max_len > 0 {
        c.max_len = max_len;
    }
    Ok(c)
}

impl Chem {
    pub fn length_ok(&self, n: usize) -> bool {
        (self.min_len == 0 || n >= self.min_len) && (self.max_len == 0 || n <= self.max_len)
    }

    pub fn probe_len_range(&self) -> (usize, usize) {
        (if self.min_len > 0 { self.min_len } else { 20 }, if self.max_len > 0 { self.max_len } else { 30 })
    }
}

/// ABI strand rule for TaqMan probes: prefer the strand with more C than G.
pub fn more_c_than_g(s: &[u8]) -> bool {
    s.iter().filter(|&&b| b == b'C').count() > s.iter().filter(|&&b| b == b'G').count()
}

fn quenchers(dye: &str) -> Option<&'static [&'static str]> {
    Some(match dye {
        "FAM" => &["BHQ1", "BHQ-1", "IBFQ", "ZEN", "NFQ", "MGB", "TAMRA", "QSY7"],
        "HEX" => &["BHQ1", "BHQ-1", "IBFQ", "ZEN", "NFQ", "MGB", "QSY7"],
        "TET" => &["BHQ1", "BHQ-1", "NFQ", "MGB"],
        "ROX" => &["BHQ2", "BHQ-2", "IBRQ", "TAO", "QSY7", "QSY21"],
        "TEXASRED" => &["BHQ2", "BHQ-2", "IBRQ", "TAO"],
        "CY3" => &["BHQ2", "BHQ-2", "IBRQ", "TAO"],
        "CY5" => &["BHQ2", "BHQ-2", "BHQ3", "IBRQ", "TAO", "QSY21"],
        "CY5.5" => &["BHQ3", "BHQ2", "IBRQ"],
        _ => return None,
    })
}

/// Reporter and quenchers with vendor spellings unified ("6FAM / MGBNFQ" -> FAM, NFQ).
fn norm_label(label: &str) -> Vec<String> {
    label
        .to_uppercase()
        .split(|c| c == '/' || c == '-' || c == ' ' || c == ',')
        .filter(|x| !x.is_empty())
        .enumerate()
        .map(|(i, x)| {
            match x {
                "6FAM" | "FAM6" => "FAM",
                "SUN" | "VIC" | "HEX" | "JOE" if i == 0 => "HEX",
                "MGBNFQ" | "NFQMGB" => "NFQ",
                "IABKFQ" | "IOWABLACKFQ" => "IBFQ",
                "IABRQSP" | "IOWABLACKRQ" => "IBRQ",
                o => o,
            }
            .to_string()
        })
        .collect()
}

/// Detection channel of a label's reporter ("" when unknown).
pub fn channel(label: &str) -> String {
    norm_label(label).into_iter().next().unwrap_or_default()
}

/// A quencher that does not cover the reporter; empty when fine or unknown.
pub fn label_problems(label: &str) -> Vec<String> {
    let p = norm_label(label);
    if p.len() < 2 {
        return Vec::new();
    }
    let Some(ok) = quenchers(&p[0]) else { return Vec::new() };
    p[1..]
        .iter()
        .filter(|q| q.as_str() != "BHQ" && !ok.contains(&q.as_str()))
        .map(|q| format!("{q} does not quench {} well", p[0]))
        .collect()
}

/// Sets meant to run in one tube must use different reporters (one per detection channel).
pub fn multiplex_problems(labels: &[(String, String)]) -> Vec<String> {
    let mut l = labels.to_vec();
    l.sort();
    let mut seen: Vec<(String, String)> = Vec::new();
    let mut out = Vec::new();
    for (set, lab) in &l {
        let d = channel(lab);
        if d.is_empty() {
            continue;
        }
        if let Some((_, o)) = seen.iter().find(|(c, _)| *c == d) {
            out.push(format!("{o} and {set} share the {d} channel"));
        }
        seen.retain(|(c, _)| *c != d);
        seen.push((d, set.clone()));
        out.extend(label_problems(lab));
    }
    out
}
