//! `chichipes design cover|sweep|shift|tile|labels`

use super::align::{read_aligned, Aligned};
use super::chem::{label_problems, multiplex_problems, parse_chem};
use super::cover::{design_from_alignment, Design, Params};
use super::exclude::Exclusion;
use super::qc::{prefer_c_strand, qc, Set, SET_COLUMNS};
use super::shift::{greedy_panel, shift, tile, Seed, ShiftOptions};
use crate::flags::FlagSet;
use crate::json;
use std::collections::HashSet;
use std::io::Write;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Instant;

pub const USAGE: &str = "chichipes design: primer / probe design from an alignment

  chichipes design cover -a aligned.fa [options]            coverage-first panel (greedy set cover)
  chichipes design sweep -a aligned.fa --grid mega -o results.jsonl [--jobs N]
                                                            many parameter combinations (design_sweep.py compatible)
  chichipes design shift -a aligned.fa --sets sets.tsv [--shift 50 --step 1]
                                                            exhaustive shifted / resized / degenerate variants of sets
  chichipes design tile  -a aligned.fa [--step 10]          best set in every amplicon-sized window along the genome
  chichipes design labels --sets primers.tsv [--groups set_groups.tsv]
                                                            reporter/quencher pairs; channel clashes inside a group

Common options (cover, shift, tile):
  --ambig N --probe-ambig N   IUPAC codes allowed in a primer / the probe (2 / 2)
  --mm N --tp N               primer mismatches allowed / 3' bases that must match (2 / 2)
  --amp 70,200                amplicon length range
  --primer-tm 52,66           primer Tm window (screening Tm, SantaLucia 1998, 250 nM / 50 mM Na+)
  --probe-above 2,14          probe Tm above the higher primer Tm
  --primer-len 18,25 --probe-len 20,30
  --probe-chem linear|zen|mgb probe format (length limits, Tm offset; see README)
  --probe-tm-offset X         override the format's probe Tm offset
  --exclude nontarget.fa      reject sets that would amplify a non-target sequence; penalise oligos that bind one
  --top N                     designs to report (3)
  --cpu N                     worker threads
";

fn flags(name: &str) -> FlagSet {
    let cpu = crate::spec::num_cpu() as i64;
    let mut f = FlagSet::new(&format!("design {name}"));
    f.usage(USAGE);
    f.string("a", "aligned FASTA");
    f.string("o", "output file (default: stdout / tile.tsv)");
    f.string("exclude", "non-target FASTA");
    f.string("sets", "sets: TSV with name, F, R, P columns (or primers.tsv)");
    f.string("segment", "shift: only sets of this segment (primers.tsv)");
    f.string("groups", "labels: set groups (group_id, label, member_ids) run in one tube");
    f.string("grid", "sweep grid: wide, hightm, midtm, longoligo, mega, mega2, mega3");
    f.string("probe-chem", "probe format: linear, zen, mgb");
    f.float("probe-tm-offset", f64::NAN, "probe Tm offset (C)");
    f.string("amp", "amplicon length range (70,200)");
    f.string("primer-tm", "primer Tm range (52,66)");
    f.string("probe-above", "probe Tm above the higher primer (2,14)");
    f.string("primer-len", "primer length range (18,25)");
    f.string("probe-len", "probe length range (default 20,30 or the format's)");
    f.int("ambig", 2, "IUPAC codes per primer");
    f.int("probe-ambig", 2, "IUPAC codes in the probe");
    f.int("mm", 2, "primer mismatches allowed");
    f.int("tp", 2, "3' bases that must match");
    f.int("top", 3, "designs to report");
    f.int("cpu", cpu, "worker threads");
    f.int("shift", 50, "shift: bases each oligo may move either way");
    f.int("step", 1, "shift / tile: step in bases");
    f.int("keep", 300, "shift / tile: variants kept per oligo for the combination step");
    f.int("max-amb", -1, "shift: IUPAC codes a variant may have (default --ambig / --probe-ambig)");
    f.int("jobs", (cpu / 4).max(1), "sweep: runs in parallel");
    f.boolean_d("qc", true, "check the chosen sets with chichipes thermodynamics (Tm 2004, hairpin, dimer)");
    f
}

fn pair(s: &str, d: (f64, f64)) -> Result<(f64, f64), String> {
    if s.is_empty() {
        return Ok(d);
    }
    let p: Vec<&str> = s.split(',').collect();
    match (p.len(), p.first().and_then(|x| x.trim().parse().ok()), p.get(1).and_then(|x| x.trim().parse().ok())) {
        (2, Some(a), Some(b)) => Ok((a, b)),
        _ => Err(format!("{s:?}: expected lo,hi")),
    }
}

fn ipair(s: &str, d: (usize, usize)) -> Result<(usize, usize), String> {
    let p = pair(s, (d.0 as f64, d.1 as f64))?;
    Ok((p.0 as usize, p.1 as usize))
}

fn params(f: &FlagSet) -> Result<Params, String> {
    let mut p = Params {
        ambig: f.get_int("ambig").max(0) as usize,
        probe_ambig: f.get_int("probe-ambig").max(0) as usize,
        pass_mm: f.get_int("mm").max(0) as usize,
        tp_n: f.get_int("tp").max(0) as usize,
        top: f.get_int("top").max(0) as usize,
        workers: f.get_int("cpu").max(1) as usize,
        ..Default::default()
    };
    p.amp = ipair(&f.get_str("amp"), p.amp)?;
    p.primer_tm = pair(&f.get_str("primer-tm"), p.primer_tm)?;
    p.probe_above = pair(&f.get_str("probe-above"), p.probe_above)?;
    p.len_primer = ipair(&f.get_str("primer-len"), p.len_primer)?;
    p.probe = parse_chem(&f.get_str("probe-chem"), f.get_float("probe-tm-offset"), 0, 0)?;
    let ql = f.get_str("probe-len");
    if !ql.is_empty() {
        p.len_probe = ipair(&ql, p.len_probe)?;
    } else if p.probe.name != "linear" {
        p.len_probe = p.probe.probe_len_range();
    }
    Ok(p)
}

fn load_aligned(path: &str) -> Result<Aligned, String> {
    if path.is_empty() {
        return Err("-a aligned.fa is required".into());
    }
    let (ids, rows) = read_aligned(path)?;
    Aligned::new(ids, rows)
}

fn out_writer(path: &str) -> Result<Box<dyn Write>, String> {
    if path.is_empty() || path == "-" {
        Ok(Box::new(std::io::BufWriter::new(std::io::stdout())))
    } else {
        let f = std::fs::File::create(path).map_err(|e| format!("{path}: {e}"))?;
        Ok(Box::new(std::io::BufWriter::new(f)))
    }
}

pub fn run(args: &[String]) -> Result<(), String> {
    let Some(cmd) = args.first() else {
        print!("{USAGE}");
        return Ok(());
    };
    if cmd == "-h" || cmd == "--help" {
        print!("{USAGE}");
        return Ok(());
    }
    let mut f = flags(cmd);
    f.parse(&args[1..]);
    match cmd.as_str() {
        "cover" => cover(&f),
        "sweep" => sweep(&f),
        "shift" => shift_cmd(&f),
        "tile" => tile_cmd(&f),
        "labels" => labels(&f),
        _ => Err(format!("unknown design command {cmd:?} (cover, sweep, shift, tile, labels)")),
    }
}

fn write_sets(w: &mut dyn Write, al: &Aligned, sets: &mut [Set], p: &Params, xs: Option<&Exclusion>, do_qc: bool) -> Result<(), String> {
    let io = |e: std::io::Error| e.to_string();
    writeln!(w, "{}", SET_COLUMNS.join("\t")).map_err(io)?;
    for s in sets.iter_mut() {
        if s.change != "-" {
            prefer_c_strand(s);
        }
        if do_qc {
            qc(s, p);
        }
        if let Some(x) = xs {
            x.check(s);
        }
        writeln!(w, "{}", s.row(al.g).join("\t")).map_err(io)?;
    }
    Ok(())
}

fn cover(f: &FlagSet) -> Result<(), String> {
    let al = load_aligned(&f.get_str("a"))?;
    let mut p = params(f)?;
    let xs = match f.get_str("exclude").as_str() {
        "" => None,
        path => Some(Arc::new(Exclusion::load(path, &p)?)),
    };
    if let Some(x) = &xs {
        let x = x.clone();
        p.reject = Some(Arc::new(move |a: &[u8], b: &[u8], c: &[u8]| x.rejects(a, b, c)));
    }
    let t0 = Instant::now();
    let run = design_from_alignment(&al, &p);
    eprintln!(
        "{} genomes; candidates F {} R {} P {}; {} designs in {:.2}s",
        al.g,
        run.nf,
        run.nr,
        run.np,
        run.designs.len(),
        t0.elapsed().as_secs_f64()
    );
    if run.status != "ok" {
        eprintln!("{}", run.message);
    }
    let mut sets: Vec<Set> = run
        .designs
        .iter()
        .map(|d| Set {
            name: d.name.clone(),
            f: d.f.clone(),
            p: d.p.clone(),
            r: d.r.clone(),
            cov: d.est_cov,
            added: d.est_added,
            amplicon: d.amplicon_len,
            fs: d.fs,
            re: d.re,
            ..Default::default()
        })
        .collect();
    let mut w = out_writer(&f.get_str("o"))?;
    write_sets(&mut *w, &al, &mut sets, &p, xs.as_deref(), f.get_bool("qc"))?;
    w.flush().map_err(|e| e.to_string())
}

// ---- sweep ----

#[derive(Clone)]
enum GVal {
    I(i64),
    A(usize, usize),
    F(f64, f64),
}

fn grid(name: &str) -> Result<Vec<(&'static str, Vec<GVal>)>, String> {
    use GVal::*;
    let ii = |xs: &[i64]| xs.iter().map(|&x| I(x)).collect::<Vec<_>>();
    let ap = |xs: &[(usize, usize)]| xs.iter().map(|&(a, b)| A(a, b)).collect::<Vec<_>>();
    let pr = |xs: &[(f64, f64)]| xs.iter().map(|&(a, b)| F(a, b)).collect::<Vec<_>>();
    Ok(match name {
        "hightm" => vec![
            ("ambig", ii(&[2, 3, 4, 5, 6])),
            ("probe_ambig", ii(&[2, 3, 4])),
            ("pass_mm", ii(&[1])),
            ("tp_n", ii(&[2])),
            ("amp", ap(&[(70, 200), (70, 250)])),
            ("primer_tm", pr(&[(56.0, 66.0), (58.0, 68.0)])),
            ("probe_above", pr(&[(6.0, 14.0), (8.0, 16.0)])),
        ],
        "longoligo" => vec![
            ("ambig", ii(&[3, 4, 5, 6])),
            ("probe_ambig", ii(&[2, 3, 4])),
            ("pass_mm", ii(&[1])),
            ("tp_n", ii(&[2])),
            ("amp", ap(&[(70, 200), (70, 250)])),
            ("len_primer_max", ii(&[30])),
            ("len_probe_max", ii(&[36])),
            ("primer_tm", pr(&[(54.0, 68.0), (56.0, 68.0)])),
            ("probe_above", pr(&[(3.0, 16.0), (6.0, 16.0)])),
        ],
        "midtm" => vec![
            ("ambig", ii(&[3, 4, 5, 6])),
            ("probe_ambig", ii(&[2, 3, 4])),
            ("pass_mm", ii(&[1])),
            ("tp_n", ii(&[2])),
            ("amp", ap(&[(70, 200), (70, 250)])),
            ("primer_tm", pr(&[(54.0, 66.0), (55.0, 67.0)])),
            ("probe_above", pr(&[(3.0, 14.0), (5.0, 14.0)])),
        ],
        "mega" => vec![
            ("ambig", ii(&[1, 2, 3, 4, 5, 6])),
            ("probe_ambig", ii(&[1, 2, 3, 4, 5])),
            ("pass_mm", ii(&[1, 2, 3])),
            ("tp_n", ii(&[1, 2])),
            ("amp", ap(&[(60, 150), (70, 200), (70, 250), (70, 300), (100, 250), (150, 300)])),
            ("primer_tm", pr(&[(50.0, 64.0), (52.0, 66.0), (54.0, 68.0), (56.0, 70.0)])),
            ("probe_above", pr(&[(2.0, 12.0), (4.0, 14.0), (6.0, 16.0)])),
            ("len_primer_max", ii(&[24, 30])),
            ("len_probe_max", ii(&[28, 36])),
        ],
        "mega2" => vec![
            ("ambig", ii(&[1, 2, 3, 4, 5, 6])),
            ("probe_ambig", ii(&[1, 2, 3, 4, 5])),
            ("pass_mm", ii(&[1, 2, 3])),
            ("tp_n", ii(&[1, 2])),
            ("amp", ap(&[(70, 180), (90, 220), (110, 260), (130, 280), (150, 320), (180, 350)])),
            ("primer_tm", pr(&[(51.0, 63.0), (53.0, 65.0), (55.0, 67.0), (57.0, 69.0), (59.0, 71.0)])),
            ("probe_above", pr(&[(3.0, 13.0), (5.0, 15.0), (7.0, 17.0)])),
            ("len_primer_max", ii(&[28, 32])),
            ("len_probe_max", ii(&[32, 40])),
        ],
        "mega3" => vec![
            ("ambig", ii(&[1, 2, 3, 4, 5, 6])),
            ("probe_ambig", ii(&[1, 2, 3, 4, 5])),
            ("pass_mm", ii(&[1, 2, 3])),
            ("tp_n", ii(&[1, 2, 3])),
            ("amp", ap(&[(60, 160), (80, 210), (100, 240), (120, 270), (140, 300), (160, 330), (200, 360)])),
            ("primer_tm", pr(&[(48.0, 62.0), (60.0, 72.0), (53.0, 69.0)])),
            ("probe_above", pr(&[(2.0, 10.0), (4.0, 16.0), (8.0, 18.0)])),
            ("len_primer_max", ii(&[26, 32])),
            ("len_probe_max", ii(&[30, 42])),
        ],
        "" | "wide" => vec![
            ("ambig", ii(&[1, 2, 3, 4, 5])),
            ("probe_ambig", ii(&[1, 2, 3, 4])),
            ("pass_mm", ii(&[1, 2])),
            ("tp_n", ii(&[2])),
            ("amp", ap(&[(70, 200), (70, 300)])),
            ("primer_tm", pr(&[(52.0, 66.0), (50.0, 68.0)])),
        ],
        _ => return Err(format!("unknown grid {name:?}")),
    })
}

/// Python's repr of a float that is a whole or short decimal (50.0, 52.5).
fn py_float(x: f64) -> String {
    let s = format!("{x}");
    if s.contains('.') || s.contains('e') || s.contains("inf") || s.contains("NaN") {
        s
    } else {
        s + ".0"
    }
}

type Combo = Vec<(&'static str, GVal)>;

fn params_json(c: &Combo) -> String {
    let parts: Vec<String> = c
        .iter()
        .map(|(k, v)| {
            let v = match v {
                GVal::I(x) => x.to_string(),
                GVal::A(a, b) => format!("[{a}, {b}]"),
                GVal::F(a, b) => format!("[{}, {}]", py_float(*a), py_float(*b)),
            };
            format!("\"{k}\": {v}")
        })
        .collect();
    format!("{{{}}}", parts.join(", "))
}

/// Order-independent key of a params object (design_sweep.py's resume check: sorted keys, numbers by value).
fn canon(v: &json::Value) -> String {
    match v {
        json::Value::Obj(kv) => {
            let mut kv: Vec<&(String, json::Value)> = kv.iter().collect();
            kv.sort_by(|a, b| a.0.cmp(&b.0));
            format!("{{{}}}", kv.iter().map(|(k, v)| format!("{k}:{}", canon(v))).collect::<Vec<_>>().join(","))
        }
        json::Value::Arr(a) => format!("[{}]", a.iter().map(canon).collect::<Vec<_>>().join(",")),
        json::Value::Num(x) => format!("{x}"),
        json::Value::Str(s) => format!("{s:?}"),
        json::Value::Bool(b) => b.to_string(),
        json::Value::Null => "null".into(),
    }
}

fn apply(c: &Combo, p: &mut Params) {
    for (k, v) in c {
        match (*k, v) {
            ("ambig", GVal::I(x)) => p.ambig = *x as usize,
            ("probe_ambig", GVal::I(x)) => p.probe_ambig = *x as usize,
            ("pass_mm", GVal::I(x)) => p.pass_mm = *x as usize,
            ("tp_n", GVal::I(x)) => p.tp_n = *x as usize,
            ("len_primer_max", GVal::I(x)) => p.len_primer.1 = *x as usize,
            ("len_probe_max", GVal::I(x)) => p.len_probe.1 = *x as usize,
            ("amp", GVal::A(a, b)) => p.amp = (*a, *b),
            ("primer_tm", GVal::F(a, b)) => p.primer_tm = (*a, *b),
            ("probe_above", GVal::F(a, b)) => p.probe_above = (*a, *b),
            _ => {}
        }
    }
}

/// A float the way Go's encoding/json writes it (shortest repr; whole numbers without ".0").
fn go_num(x: f64) -> String {
    format!("{x}")
}

fn design_json(d: &Design) -> String {
    let mut out = Vec::new();
    let s = |out: &mut Vec<u8>, v: &[u8]| json::write_go_string(out, v);
    out.extend_from_slice(b"{\"name\":");
    s(&mut out, d.name.as_bytes());
    for (k, v) in [("F", &d.f), ("P", &d.p), ("R", &d.r)] {
        out.extend_from_slice(format!(",\"{k}\":").as_bytes());
        s(&mut out, v);
    }
    let nums = [("F_tm", d.f_tm), ("P_tm", d.p_tm), ("R_tm", d.r_tm), ("F_gc", d.f_gc), ("P_gc", d.p_gc), ("R_gc", d.r_gc)];
    for (k, v) in nums {
        out.extend_from_slice(format!(",\"{k}\":{}", go_num(v)).as_bytes());
    }
    out.extend_from_slice(
        format!(
            ",\"amplicon_len\":{},\"std_est\":null,\"est_cov\":{},\"est_cov_pct\":{},\"est_added\":{},\"P_strand\":\"{}\"}}",
            d.amplicon_len,
            d.est_cov,
            go_num(d.est_cov_pct),
            d.est_added,
            d.p_strand as char
        )
        .as_bytes(),
    );
    String::from_utf8(out).unwrap()
}

fn sweep(f: &FlagSet) -> Result<(), String> {
    let aln = f.get_str("a");
    let al = load_aligned(&aln)?;
    let g = grid(&f.get_str("grid"))?;
    let mut base = params(f)?;
    if f.get_int("top") == 3 {
        base.top = 10;
    }
    let out = f.get_str("o");
    if out.is_empty() {
        return Err("-o results.jsonl is required".into());
    }
    let mut todo: Vec<Combo> = vec![Vec::new()];
    for (k, vals) in &g {
        todo = todo.into_iter().flat_map(|c| vals.iter().map(move |v| {
            let mut c2 = c.clone();
            c2.push((*k, v.clone()));
            c2
        })).collect();
    }
    let mut done: HashSet<String> = HashSet::new();
    if let Ok(b) = std::fs::read_to_string(&out) {
        for ln in b.lines() {
            if let Ok(v) = json::parse(ln) {
                if let json::Value::Obj(_) = v.get("params") {
                    done.insert(canon(v.get("params")));
                }
            }
        }
    }
    let left: Vec<(Combo, String)> = todo
        .into_iter()
        .map(|c| {
            let js = params_json(&c);
            (c, js)
        })
        .filter(|(_, js)| json::parse(js).map(|v| !done.contains(&canon(&v))).unwrap_or(true))
        .collect();
    let segment = std::path::Path::new(&aln).parent().and_then(|p| p.file_name()).map(|s| s.to_string_lossy().into_owned()).unwrap_or_default();
    let jobs = f.get_int("jobs").max(1) as usize;
    let total = left.len();
    eprintln!("{} genomes; {} runs to do ({} done), {} in parallel", al.g, total, done.len(), jobs);
    let file = std::fs::OpenOptions::new().create(true).append(true).open(&out).map_err(|e| format!("{out}: {e}"))?;
    let file = Mutex::new((file, 0usize));
    let next = AtomicUsize::new(0);
    let per = (base.workers / jobs).max(1);
    let t0 = Instant::now();
    std::thread::scope(|sc| {
        for _ in 0..jobs.min(total.max(1)) {
            sc.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::Relaxed);
                if i >= total {
                    break;
                }
                let (c, js) = &left[i];
                let mut p = base.clone();
                p.workers = per;
                apply(c, &mut p);
                let t = Instant::now();
                let run = design_from_alignment(&al, &p);
                let designs: Vec<String> = run.designs.iter().map(design_json).collect();
                let mut seg = Vec::new();
                json::write_go_string(&mut seg, segment.as_bytes());
                let line = format!(
                    "{{\"segment\":{},\"params\":{},\"n\":{},\"status\":\"{}\",\"designs\":[{}],\"seconds\":{},\"engine\":\"chichipes\"}}\n",
                    String::from_utf8_lossy(&seg),
                    js,
                    run.n,
                    run.status,
                    designs.join(","),
                    go_num((t.elapsed().as_secs_f64() * 100.0).round() / 100.0)
                );
                let mut g = file.lock().unwrap();
                let _ = g.0.write_all(line.as_bytes());
                g.1 += 1;
                if g.1 % 50 == 0 || g.1 == total {
                    eprintln!("[{}/{}] {:.0}s", g.1, total, t0.elapsed().as_secs_f64());
                }
            });
        }
    });
    Ok(())
}

// ---- shift / tile / labels ----

fn read_tsv(path: &str) -> Result<Vec<Vec<(String, String)>>, String> {
    let b = std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
    let mut lines = b.trim_end_matches('\n').split('\n');
    let hdr: Vec<String> = lines.next().unwrap_or("").split('\t').map(|h| h.trim().to_string()).collect();
    Ok(lines
        .map(|ln| {
            let f: Vec<&str> = ln.split('\t').collect();
            hdr.iter().enumerate().filter(|(i, _)| *i < f.len()).map(|(i, h)| (h.clone(), f[i].trim().to_string())).collect()
        })
        .collect())
}

fn col<'a>(r: &'a [(String, String)], names: &[&str]) -> Option<&'a str> {
    for n in names {
        if let Some((_, v)) = r.iter().find(|(h, _)| h.eq_ignore_ascii_case(n)) {
            return Some(v);
        }
    }
    None
}

fn read_seeds(path: &str, segment: &str) -> Result<Vec<Seed>, String> {
    let rows = read_tsv(path)?;
    let mut out = Vec::new();
    for (i, r) in rows.iter().enumerate() {
        let up = |v: Option<&str>| v.unwrap_or("").to_uppercase().into_bytes();
        if !segment.is_empty() {
            if let Some(s) = col(r, &["segment"]) {
                if !s.eq_ignore_ascii_case(segment) {
                    continue;
                }
            }
        }
        let (fw, rv, pr) = (up(col(r, &["F", "forward_primer"])), up(col(r, &["R", "reverse_primer"])), up(col(r, &["P", "probe"])));
        if fw.is_empty() || rv.is_empty() || pr.is_empty() {
            continue;
        }
        let name = col(r, &["name", "primer_name", "primer_id"]).map(String::from).unwrap_or(format!("set{}", i + 1));
        out.push(Seed { name, f: fw, p: pr, r: rv });
    }
    if out.is_empty() && rows.is_empty() {
        return Err(format!("{path}: no sets"));
    }
    Ok(out)
}

fn shift_cmd(f: &FlagSet) -> Result<(), String> {
    let al = load_aligned(&f.get_str("a"))?;
    let p = params(f)?;
    let sets = f.get_str("sets");
    if sets.is_empty() {
        return Err("--sets sets.tsv is required".into());
    }
    let seeds = read_seeds(&sets, &f.get_str("segment"))?;
    let xs = match f.get_str("exclude").as_str() {
        "" => None,
        path => Some(Exclusion::load(path, &p)?),
    };
    let ma = f.get_int("max-amb");
    let opt = ShiftOptions {
        shift: f.get_int("shift").max(0) as usize,
        step: f.get_int("step").max(1) as usize,
        keep: f.get_int("keep").max(1) as usize,
        primer_amb: if ma >= 0 { ma as usize } else { p.ambig },
        probe_amb: if ma >= 0 { ma as usize } else { p.probe_ambig.max(3) },
        dtm: 3.0,
    };
    let mut w = out_writer(&f.get_str("o"))?;
    let io = |e: std::io::Error| e.to_string();
    writeln!(w, "{}", SET_COLUMNS.join("\t")).map_err(io)?;
    for s in &seeds {
        let t0 = Instant::now();
        let (res, mut base) = match shift(&al, s, &p, &opt) {
            Ok(x) => x,
            Err(e) => {
                eprintln!("{e}");
                continue;
            }
        };
        eprintln!(
            "{}: now {}/{} ({:.1}%); {} improved classes ({:.1}s)",
            s.name,
            base.cov,
            al.g,
            100.0 * base.cov as f64 / al.g as f64,
            res.len(),
            t0.elapsed().as_secs_f64()
        );
        base.change = "-".into();
        base.base_cov = base.cov;
        let mut all = vec![base];
        all.extend(res);
        for x in all.iter_mut() {
            if x.change != "-" {
                prefer_c_strand(x);
            }
            if f.get_bool("qc") {
                qc(x, &p);
            }
            if let Some(e) = &xs {
                e.check(x);
            }
            writeln!(w, "{}", x.row(al.g).join("\t")).map_err(io)?;
        }
    }
    w.flush().map_err(io)
}

fn tile_cmd(f: &FlagSet) -> Result<(), String> {
    let al = load_aligned(&f.get_str("a"))?;
    let p = params(f)?;
    let step = f.get_int("step").max(1) as usize;
    let t0 = Instant::now();
    let rows = tile(&al, &p, step, f.get_int("keep").max(1) as usize);
    eprintln!("{} genomes, {} windows (step {}) in {:.1}s", al.g, rows.len(), step, t0.elapsed().as_secs_f64());
    let mut out = f.get_str("o");
    if out.is_empty() {
        out = "tile.tsv".into();
    }
    let io = |e: std::io::Error| e.to_string();
    let mut w = out_writer(&out)?;
    writeln!(w, "window_start\twindow_end\tcoverage\tcoverage_pct\tf_start\tr_end\tamplicon\tF\tP\tR").map_err(io)?;
    for r in &rows {
        match &r.set {
            None => writeln!(w, "{}\t{}\t0\t0.00\t\t\t\t\t\t", r.start, r.start + p.amp.1).map_err(io)?,
            Some(s) => writeln!(
                w,
                "{}\t{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{}\t{}",
                r.start,
                r.start + p.amp.1,
                s.cov,
                100.0 * s.cov as f64 / al.g as f64,
                s.fs,
                s.re,
                s.amplicon,
                String::from_utf8_lossy(&s.f),
                String::from_utf8_lossy(&s.p),
                String::from_utf8_lossy(&s.r)
            )
            .map_err(io)?,
        }
    }
    w.flush().map_err(io)?;
    let xs = match f.get_str("exclude").as_str() {
        "" => None,
        path => Some(Exclusion::load(path, &p)?),
    };
    let mut panel = greedy_panel(&rows, f.get_int("top").max(0) as usize, al.g);
    let pp = format!("{}.panel.tsv", out.strip_suffix(".tsv").unwrap_or(&out));
    let mut pw = out_writer(&pp)?;
    write_sets(&mut *pw, &al, &mut panel, &p, xs.as_ref(), f.get_bool("qc"))?;
    pw.flush().map_err(io)
}

fn labels(f: &FlagSet) -> Result<(), String> {
    let sets = f.get_str("sets");
    if sets.is_empty() {
        return Err("--sets primers.tsv is required".into());
    }
    let rows = read_tsv(&sets)?;
    let mut label: Vec<(String, String)> = Vec::new();
    let mut bad = 0;
    for r in &rows {
        let id = col(r, &["primer_id"]).filter(|x| !x.is_empty()).or(col(r, &["name"])).unwrap_or("").to_string();
        let l = col(r, &["fluorophore"]).unwrap_or("").to_string();
        if l.is_empty() {
            continue;
        }
        for pr in label_problems(&l) {
            println!("{} {}: {}", id, col(r, &["primer_name"]).unwrap_or(""), pr);
            bad += 1;
        }
        label.push((id, l));
    }
    let groups = f.get_str("groups");
    if !groups.is_empty() {
        for g in read_tsv(&groups)? {
            let members: Vec<String> = col(&g, &["member_ids"]).unwrap_or("").split(';').map(|x| x.trim().to_string()).collect();
            let m: Vec<(String, String)> = label.iter().filter(|(id, _)| members.contains(id)).cloned().collect();
            for pr in multiplex_problems(&m) {
                if pr.contains("share") {
                    println!("group {} {} (if run in one tube): {}", col(&g, &["group_id"]).unwrap_or(""), col(&g, &["label"]).unwrap_or(""), pr);
                    bad += 1;
                }
            }
        }
    }
    eprintln!("{bad} label problems");
    Ok(())
}

