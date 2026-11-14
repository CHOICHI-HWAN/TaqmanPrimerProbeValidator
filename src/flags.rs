//! Command-line flags parsed like Go's `flag` package (the Go engine's and MFEprimer's syntax): `-name value`,
//! `--name value`, `-name=value`, boolean `-name`; parsing stops at the first non-flag argument or at `--`.

use std::collections::HashMap;
use std::process::exit;

#[derive(Clone)]
enum Val {
    Int(i64),
    Float(f64),
    Bool(bool),
    Str(String),
    List(Vec<String>),
}

pub struct FlagSet {
    name: String,
    usage: String,
    order: Vec<(String, String)>, // name, help
    vals: HashMap<String, Val>,
    alias: HashMap<String, String>,
}

impl FlagSet {
    pub fn new(name: &str) -> Self {
        FlagSet { name: name.into(), usage: String::new(), order: Vec::new(), vals: HashMap::new(), alias: HashMap::new() }
    }
    pub fn usage(&mut self, u: &str) {
        self.usage = u.into();
    }
    fn def(&mut self, name: &str, v: Val, help: &str) {
        self.order.push((name.into(), help.into()));
        self.vals.insert(name.into(), v);
    }
    pub fn int(&mut self, name: &str, d: i64, help: &str) {
        self.def(name, Val::Int(d), help)
    }
    pub fn float(&mut self, name: &str, d: f64, help: &str) {
        self.def(name, Val::Float(d), help)
    }
    pub fn boolean(&mut self, name: &str, help: &str) {
        self.def(name, Val::Bool(false), help)
    }
    pub fn string(&mut self, name: &str, help: &str) {
        self.def(name, Val::Str(String::new()), help)
    }
    /// A repeatable flag whose values may also be comma-separated.
    pub fn list(&mut self, name: &str, help: &str) {
        self.def(name, Val::List(Vec::new()), help)
    }
    /// A second name that sets the same value (Go's two flag.IntVar calls on one variable).
    pub fn alias(&mut self, alias: &str, of: &str) {
        self.order.push((alias.into(), String::new()));
        self.alias.insert(alias.into(), of.into());
    }

    pub fn get_int(&self, n: &str) -> i64 {
        match &self.vals[n] {
            Val::Int(v) => *v,
            _ => unreachable!(),
        }
    }
    pub fn get_float(&self, n: &str) -> f64 {
        match &self.vals[n] {
            Val::Float(v) => *v,
            _ => unreachable!(),
        }
    }
    pub fn get_bool(&self, n: &str) -> bool {
        match &self.vals[n] {
            Val::Bool(v) => *v,
            _ => unreachable!(),
        }
    }
    pub fn get_str(&self, n: &str) -> String {
        match &self.vals[n] {
            Val::Str(v) => v.clone(),
            _ => unreachable!(),
        }
    }
    pub fn get_list(&self, n: &str) -> Vec<String> {
        match &self.vals[n] {
            Val::List(v) => v.clone(),
            _ => unreachable!(),
        }
    }

    fn print_usage(&self) {
        if self.usage.is_empty() {
            eprintln!("Usage of {}:", self.name);
        } else {
            eprint!("{}", self.usage);
        }
        let mut names: Vec<&(String, String)> = self.order.iter().collect();
        names.sort_by(|a, b| a.0.cmp(&b.0));
        for (n, h) in names {
            eprintln!("  -{n}\n    \t{h}");
        }
    }

    fn fail(&self, msg: &str) -> ! {
        eprintln!("{msg}");
        self.print_usage();
        exit(2)
    }

    /// Parses flags and returns the remaining arguments.
    pub fn parse(&mut self, args: &[String]) -> Vec<String> {
        let mut i = 0;
        while i < args.len() {
            let s = &args[i];
            if s.len() < 2 || !s.starts_with('-') {
                break;
            }
            let mut name = &s[1..];
            if let Some(rest) = name.strip_prefix('-') {
                if rest.is_empty() {
                    i += 1;
                    break;
                }
                name = rest;
            }
            if name.is_empty() || name.starts_with('-') || name.starts_with('=') {
                self.fail(&format!("bad flag syntax: {s}"));
            }
            i += 1;
            let (name, value) = match name.split_once('=') {
                Some((n, v)) => (n.to_string(), Some(v.to_string())),
                None => (name.to_string(), None),
            };
            let key = self.alias.get(&name).cloned().unwrap_or(name.clone());
            let Some(cur) = self.vals.get(&key).cloned() else {
                if name == "help" || name == "h" {
                    self.print_usage();
                    exit(0);
                }
                self.fail(&format!("flag provided but not defined: -{name}"));
            };
            let new = match cur {
                Val::Bool(_) => match value.as_deref() {
                    None => Val::Bool(true),
                    Some("1" | "t" | "T" | "true" | "TRUE" | "True") => Val::Bool(true),
                    Some("0" | "f" | "F" | "false" | "FALSE" | "False") => Val::Bool(false),
                    Some(v) => self.fail(&format!("invalid boolean value {v:?} for -{name}: parse error")),
                },
                other => {
                    let v = match value {
                        Some(v) => v,
                        None => {
                            if i >= args.len() {
                                self.fail(&format!("flag needs an argument: -{name}"));
                            }
                            i += 1;
                            args[i - 1].clone()
                        }
                    };
                    match other {
                        Val::Int(_) => match parse_int(&v) {
                            Some(x) => Val::Int(x),
                            None => self.fail(&format!("invalid value {v:?} for flag -{name}: parse error")),
                        },
                        Val::Float(_) => match v.parse::<f64>() {
                            Ok(x) => Val::Float(x),
                            Err(_) => self.fail(&format!("invalid value {v:?} for flag -{name}: parse error")),
                        },
                        Val::Str(_) => Val::Str(v),
                        Val::List(mut l) => {
                            l.extend(v.split(',').map(str::trim).filter(|p| !p.is_empty()).map(String::from));
                            Val::List(l)
                        }
                        Val::Bool(_) => unreachable!(),
                    }
                }
            };
            self.vals.insert(key, new);
        }
        args[i..].to_vec()
    }
}

/// strconv.ParseInt(s, 0, 64): decimal, or 0x / 0o / 0b / leading-0 octal, with an optional sign and underscores.
fn parse_int(s: &str) -> Option<i64> {
    let (neg, body) = match s.as_bytes().first() {
        Some(b'-') => (true, &s[1..]),
        Some(b'+') => (false, &s[1..]),
        _ => (false, s),
    };
    let lower = body.to_ascii_lowercase();
    let (radix, digits) = if let Some(d) = lower.strip_prefix("0x") {
        (16, d.to_string())
    } else if let Some(d) = lower.strip_prefix("0o") {
        (8, d.to_string())
    } else if let Some(d) = lower.strip_prefix("0b") {
        (2, d.to_string())
    } else if lower.len() > 1 && lower.starts_with('0') {
        (8, lower[1..].to_string())
    } else {
        (10, lower.clone())
    };
    let digits = if radix != 10 || digits.contains('_') && body != digits { digits.replace('_', "") } else { digits };
    let v = i64::from_str_radix(&digits, radix).ok()?;
    Some(if neg { -v } else { v })
}
