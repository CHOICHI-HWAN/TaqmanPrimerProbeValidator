//! Engines of the primer-validation pipeline that used to be separate Go programs (go/scan, go/probedesign). Each reads one
//! JSON request on stdin and writes the same JSON answer the Go program wrote, so the pipeline's Python only changes the
//! command it starts.

pub mod probedesign;
pub mod scan;

use std::io::{Read, Write};

pub fn out_int(o: &mut Vec<u8>, v: i64) {
    o.extend_from_slice(v.to_string().as_bytes());
}

pub fn out_str(o: &mut Vec<u8>, s: &str) {
    crate::json::write_go_string(o, s.as_bytes());
}

/// A float the way Go's encoding/json writes it (shortest round-trip; exponent form below 1e-6 or from 1e21).
pub fn go_float(o: &mut Vec<u8>, x: f64) {
    let a = x.abs();
    if a != 0.0 && (a < 1e-6 || a >= 1e21) {
        // Go: strconv 'e' with -1 precision, then "e-07" -> "e-07" (two-digit exponent kept)
        let s = format!("{x:e}");
        let (m, e) = s.split_once('e').unwrap();
        let (sign, digits) = if let Some(d) = e.strip_prefix('-') { ("-", d) } else { ("+", e) };
        let digits = if digits.len() < 2 { format!("0{digits}") } else { digits.to_string() };
        o.extend_from_slice(format!("{m}e{sign}{digits}").as_bytes());
    } else {
        o.extend_from_slice(format!("{x}").as_bytes());
    }
}

fn stdin_text() -> Result<String, String> {
    let mut s = String::new();
    std::io::stdin().read_to_string(&mut s).map_err(|e| format!("bad request: {e}"))?;
    Ok(s)
}

fn emit(out: &[u8]) -> i32 {
    let mut w = std::io::stdout().lock();
    if w.write_all(out).and_then(|_| w.flush()).is_err() {
        return 1;
    }
    0
}

/// `chichipes scan` (stdin JSON -> stdout JSON, as go/scan).
pub fn run_scan(_args: &[String]) -> i32 {
    let input = match stdin_text() {
        Ok(s) => s,
        Err(e) => {
            eprintln!("{e}");
            return 2;
        }
    };
    match scan::run(&input) {
        Ok(o) => emit(&o),
        Err(e) => {
            eprintln!("{e}");
            if e.starts_with("bad request") {
                2
            } else {
                1
            }
        }
    }
}

/// `chichipes probedesign` (stdin JSON -> stdout JSON, as go/probedesign).
pub fn run_probedesign(_args: &[String]) -> i32 {
    let input = match stdin_text() {
        Ok(s) => s,
        Err(e) => {
            eprintln!("{e}");
            return 2;
        }
    };
    match probedesign::run(&input) {
        Ok(o) => emit(&o),
        Err((code, e)) => {
            eprintln!("{e}");
            code
        }
    }
}
