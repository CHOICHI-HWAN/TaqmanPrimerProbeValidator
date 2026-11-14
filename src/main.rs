//! chichipes: clean-room primer engine with an mfeprimer-compatible command line (the subset the pipeline uses).
//!
//!     chichipes -d db.fa -i primers.fa -o out.txt [--fasta] [--virus] [--misMatch 2] [--cpu 8] [-t 30] [-T 100] [-s 0] [-S 2000]
//!     chichipes index -i db.fa        (writes the .fai/.json the pipeline expects; the search needs no index)
//!     chichipes dimer -i oligos.fa -o dimer.txt / chichipes hairpin -i oligos.fa -o hairpin.txt
//!     chichipes reads -i oligos.fa -r reads.fq.gz -o prefix   (or -1 R1 -2 R2; see `chichipes reads -h`)
//!     chichipes version

mod dimer;
mod flags;
mod golog;
mod hairpin;
mod json;
mod primers;
mod reads;
mod search;
mod seqio;
mod spec;
mod thermo;

use flags::FlagSet;
use std::process::exit;
use std::time::Instant;

/// Set at release build time: CHICHIPES_VERSION=v0.4.0 cargo build --release
fn version() -> String {
    format!("chichipes {} (mfeprimer-4.2.4 compatible subset)", option_env!("CHICHIPES_VERSION").unwrap_or("dev"))
}

fn check<T, E: std::fmt::Display>(r: Result<T, E>) -> T {
    r.unwrap_or_else(|e| {
        eprintln!("error: {e}");
        exit(1)
    })
}

fn cond_flags(fs: &mut FlagSet) {
    fs.float("mono", 50.0, "");
    fs.float("diva", 1.5, "");
    fs.float("dntp", 0.25, "");
    fs.float("oligo", 50.0, "");
}

fn cond(fs: &FlagSet) -> thermo::Cond {
    thermo::Cond { mono: fs.get_float("mono"), diva: fs.get_float("diva"), dntp: fs.get_float("dntp"), oligo: fs.get_float("oligo") }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() > 1 {
        let rest = &args[2..];
        match args[1].as_str() {
            "version" | "-v" | "--version" | "-version" => {
                println!("{}", version());
                return;
            }
            "index" => return run_index(rest),
            "search" => return run_search(rest),
            "dimer" => return run_dimer(rest),
            "hairpin" => return run_hairpin(rest),
            "reads" => return run_reads(rest),
            _ => {}
        }
    }
    let d = spec::Options::default();
    let mut fs = FlagSet::new("chichipes");
    fs.string("d", "database FASTA");
    fs.string("i", "primer FASTA");
    fs.string("o", "output report");
    fs.boolean("fasta", "write amplicons in FASTA (always written)");
    fs.boolean("virus", "virus coverage mode (always on)");
    fs.int("misMatch", 0, "max mismatches in the 3' window");
    fs.int("misStart", 1, "");
    fs.int("misEnd", 9, "");
    fs.int("k", 9, "");
    fs.int("cpu", d.cpu as i64, "");
    fs.alias("c", "cpu");
    fs.float("t", 30.0, "");
    fs.float("T", 100.0, "");
    fs.int("s", 0, "");
    fs.int("S", 2000, "");
    cond_flags(&mut fs);
    fs.int("M", 10000, "max amplicons written to the FASTA");
    let mut a = &args[1..];
    if a.first().map(String::as_str) == Some("spec") {
        a = &a[1..];
    }
    fs.parse(a);
    let (db, input, out) = (fs.get_str("d"), fs.get_str("i"), fs.get_str("o"));
    if db.is_empty() || input.is_empty() || out.is_empty() {
        eprintln!("need -d, -i and -o");
        exit(2);
    }
    let opt = spec::Options {
        k: fs.get_int("k") as usize,
        mis_match: fs.get_int("misMatch") as usize,
        mis_start: fs.get_int("misStart") as usize,
        mis_end: fs.get_int("misEnd") as usize,
        min_size: fs.get_int("s"),
        max_size: fs.get_int("S"),
        min_tm: fs.get_float("t"),
        max_tm: fs.get_float("T"),
        cpu: fs.get_int("cpu").max(1) as usize,
        cond: cond(&fs),
    };
    let t0 = Instant::now();
    let oligos = check(primers::read_fasta(&input));
    let e = check(spec::Engine::new(&oligos, opt));
    let recs = check(seqio::read_all(&db));
    let res = e.run(&recs);
    check(spec::write_outputs(&e, &out, &db, &recs, &res, fs.get_int("M").max(0) as usize));
    eprintln!(
        "{} records, {} amplicons, {} missed, {:.1}s",
        res.total,
        res.amplicons.len(),
        res.failed.len(),
        t0.elapsed().as_secs_f64()
    );
}

/// chichipes search db.fa query.fa out.txt [-m 2] [-k 9] [--cpu N]
fn run_search(args: &[String]) {
    let mut fs = FlagSet::new("search");
    fs.int("m", 2, "max mismatches (edit distance)");
    fs.alias("mismatch", "m");
    fs.int("k", 9, "");
    fs.int("cpu", spec::num_cpu() as i64, "");
    cond_flags(&mut fs);
    let mut pos = Vec::new();
    let mut rest = args.to_vec();
    while !rest.is_empty() {
        // positional arguments may come before or after the flags
        rest = fs.parse(&rest);
        if !rest.is_empty() {
            pos.push(rest.remove(0));
        }
    }
    if pos.len() != 3 {
        eprintln!("usage: chichipes search db.fa query.fa out.txt [-m N]");
        exit(2);
    }
    let t0 = Instant::now();
    let recs = check(seqio::read_all(&pos[0]));
    let oligos = check(primers::read_fasta(&pos[1]));
    let (mm, k) = (fs.get_int("m").max(0) as usize, fs.get_int("k") as usize);
    let (qs, hits) = check(search::search(&recs, &oligos, k, mm, fs.get_int("cpu").max(1) as usize, cond(&fs)));
    check(search::write_search(&pos[2], &recs, &qs, &hits, mm));
    eprintln!("{} records, {} hits, {:.1}s", recs.len(), hits.len(), t0.elapsed().as_secs_f64());
}

/// chichipes dimer -i oligos.fa -o dimer.txt [-s 5] [-m 2] [-d 0] [--cpu N]
fn run_dimer(args: &[String]) {
    let mut fs = FlagSet::new("dimer");
    fs.string("i", "oligo FASTA");
    fs.string("o", "report");
    fs.int("s", 5, "score cutoff");
    fs.int("m", 2, "max mismatches");
    fs.float("d", 0.0, "dG cutoff");
    fs.int("cpu", 1, "");
    fs.int("c", 1, "");
    fs.boolean("p", "only dimers with 3' end binding (not implemented)");
    cond_flags(&mut fs);
    fs.parse(args);
    let c = cond(&fs);
    let input = fs.get_str("i");
    let oligos = check(primers::read_fasta(&input));
    let ex = check(dimer::expand_all(&oligos));
    let (s, m, d) = (fs.get_int("s"), fs.get_int("m"), fs.get_float("d"));
    let ds = dimer::dimers(&ex, s, m, d, c);
    check(dimer::write_dimers(&fs.get_str("o"), &input, &ex, &ds, s, m, d, c));
}

/// chichipes hairpin -i oligos.fa -o hairpin.txt [-s 5] [-d 0] [-t 0] [--mindh 2] [--minl 3] [--maxl 12] [--maxgap 2]
fn run_hairpin(args: &[String]) {
    let mut fs = FlagSet::new("hairpin");
    let d = hairpin::DEFAULT;
    fs.string("i", "oligo FASTA");
    fs.string("o", "report");
    fs.int("s", d.score, "score cutoff");
    fs.float("d", d.dg, "dG cutoff");
    fs.float("t", d.tm, "Tm cutoff");
    fs.int("mindh", d.min_dh, "");
    fs.int("minl", d.min_l, "");
    fs.int("maxl", d.max_l, "");
    fs.int("maxgap", d.max_gap, "");
    fs.int("cpu", 1, "");
    fs.int("c", 1, "");
    cond_flags(&mut fs);
    fs.parse(args);
    let o = hairpin::Options {
        score: fs.get_int("s"),
        dg: fs.get_float("d"),
        tm: fs.get_float("t"),
        min_dh: fs.get_int("mindh"),
        min_l: fs.get_int("minl"),
        max_l: fs.get_int("maxl"),
        max_gap: fs.get_int("maxgap"),
    };
    let c = cond(&fs);
    let input = fs.get_str("i");
    let oligos = check(primers::read_fasta(&input));
    let ex = check(dimer::expand_all(&oligos));
    let hs = hairpin::hairpins(&ex, &o, c);
    check(hairpin::write_hairpins(&fs.get_str("o"), &input, &ex, &hs, &o, c));
}

struct FaiEntry {
    id: Vec<u8>,
    desc: Vec<u8>,
    size: i64,
    offset: i64,
    line_bases: i64,
    line_width: i64,
}

/// Writes the files `mfeprimer index` leaves next to the database, because the pipeline's Snakemake rules declare them
/// as outputs and the web app reads the .fai. The search itself never needs them: .fai (samtools layout) and .json are
/// real, .primerqc/.primerqc.fai are empty placeholders for mfeprimer's own k-mer index.
fn run_index(args: &[String]) {
    let mut fs = FlagSet::new("index");
    fs.string("i", "database FASTA");
    fs.int("cpu", 1, "");
    fs.int("c", 1, "");
    fs.int("k", 9, "");
    fs.parse(args);
    let input = fs.get_str("i");
    if input.is_empty() {
        eprintln!("usage: chichipes index -i db.fa");
        exit(2);
    }
    let entries = check(scan_fai(&input));
    let mut fai = Vec::new();
    for e in &entries {
        fai.extend_from_slice(&e.id);
        fai.extend_from_slice(format!("\t{}\t{}\t{}\t{}\n", e.size, e.offset, e.line_bases, e.line_width).as_bytes());
    }
    // Go's encoding/json writes map keys sorted as strings: "0", "1", "10", ...
    let mut keys: Vec<(String, usize)> = (0..entries.len()).map(|i| (i.to_string(), i)).collect();
    keys.sort();
    let mut js = vec![b'{'];
    for (n, (k, i)) in keys.iter().enumerate() {
        if n > 0 {
            js.push(b',');
        }
        let e = &entries[*i];
        json::write_go_string(&mut js, k.as_bytes());
        js.extend_from_slice(b":{\"desc\":");
        json::write_go_string(&mut js, &e.desc);
        js.extend_from_slice(b",\"id\":");
        json::write_go_string(&mut js, &e.id);
        js.extend_from_slice(format!(",\"size\":{}}}", e.size).as_bytes());
    }
    js.push(b'}');
    check(std::fs::write(format!("{input}.fai"), fai));
    check(std::fs::write(format!("{input}.json"), js));
    check(std::fs::write(format!("{input}.log"), format!("{}\nindexed {} records of {}\n", version(), entries.len(), input)));
    check(std::fs::write(format!("{input}.primerqc"), b""));
    check(std::fs::write(format!("{input}.primerqc.fai"), b""));
}

fn scan_fai(path: &str) -> Result<Vec<FaiEntry>, String> {
    let data = std::fs::read(path).map_err(|e| format!("open {path}: {e}"))?;
    let mut out: Vec<FaiEntry> = Vec::new();
    let mut pos = 0i64;
    for line in data.split_inclusive(|&b| b == b'\n') {
        let n = line.len() as i64;
        let mut body = line;
        while let [rest @ .., b'\n' | b'\r'] = body {
            body = rest;
        }
        if let Some(h) = body.strip_prefix(b">") {
            let (id, desc) = match h.iter().position(|&b| b == b' ') {
                Some(i) => (&h[..i], &h[i + 1..]),
                None => (h, &b""[..]),
            };
            out.push(FaiEntry { id: id.to_vec(), desc: desc.to_vec(), size: 0, offset: pos + n, line_bases: 0, line_width: 0 });
        } else if let Some(cur) = out.last_mut() {
            if !body.is_empty() {
                if cur.line_bases == 0 {
                    cur.line_bases = body.len() as i64;
                    cur.line_width = n;
                }
                cur.size += body.len() as i64;
            }
        }
        pos += n;
    }
    Ok(out)
}

const READS_USAGE: &str = "usage:
  chichipes reads -i oligos.fa -r reads.fastq.gz [-r more.fq.gz] -o prefix [-m 2] [--indel] [--cpu N] [--hits]
  chichipes reads -i oligos.fa -1 R1.fastq.gz -2 R2.fastq.gz -o prefix   (paired-end: a pair counts once)

Checks which reads carry each oligo (FASTQ or FASTA, gzip or not). Name oligos SET_F / SET_R / SET_P to also get
per-set counts (reads with F and R, with F, R and P). -m is mismatches per oligo (substitutions; with --indel, edits,
for error-prone long reads). Writes <prefix>.summary.tsv, .sets.tsv, .positions.tsv (bases the reads have at every
oligo position), .txt and, with --hits, .hits.tsv (one line per read and oligo found).
";

/// chichipes reads ...
fn run_reads(args: &[String]) {
    let mut fs = FlagSet::new("reads");
    fs.usage(READS_USAGE);
    fs.string("i", "oligo FASTA (primers / probes)");
    fs.list("r", "single-end or long-read file(s), repeatable or comma-separated");
    fs.list("1", "paired-end R1 file(s)");
    fs.list("2", "paired-end R2 file(s), same order as -1");
    fs.string("o", "output prefix");
    fs.int("m", 2, "max mismatches (or edits with --indel) per oligo");
    fs.boolean("indel", "allow insertions/deletions (ONT, PacBio CLR)");
    fs.int("seed", 0, "seed length (0 = automatic)");
    fs.int("cpu", spec::num_cpu() as i64, "threads");
    fs.boolean("hits", "also write <prefix>.hits.tsv");
    fs.parse(args);
    let (se, r1, r2) = (fs.get_list("r"), fs.get_list("1"), fs.get_list("2"));
    let (input, out) = (fs.get_str("i"), fs.get_str("o"));
    let paired = !r1.is_empty();
    if input.is_empty() || out.is_empty() || se.is_empty() == !paired || (paired && r2.is_empty()) {
        eprint!("{READS_USAGE}");
        exit(2);
    }
    let t0 = Instant::now();
    let ol = check(primers::read_fasta(&input));
    let oligos = check(reads::from_primers(&ol));
    let s = check(reads::new_searcher(oligos, fs.get_int("m").max(0) as usize, fs.get_int("seed").max(0) as usize, fs.get_bool("indel")));
    let hf = fs.get_bool("hits").then(|| check(std::fs::File::create(format!("{out}.hits.tsv"))));
    let cpu = fs.get_int("cpu").max(1) as usize;
    let st = if paired { s.run(&r1, &r2, cpu, hf) } else { s.run(&se, &[], cpu, hf) };
    let st = check(st);
    let rep = check(reads::write_reports(&out, &s, &st, paired));
    print!("{rep}");
    eprintln!("{:.1}s", t0.elapsed().as_secs_f64());
}
