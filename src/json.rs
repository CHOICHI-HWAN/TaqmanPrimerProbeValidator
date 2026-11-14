//! A small JSON reader for the embedded parameter table, and Go-compatible string escaping for `index`.
//! Numbers go through `str::parse::<f64>`, which rounds correctly like Go's encoding/json.

#[allow(dead_code)]
pub enum Value {
    Null,
    Bool(bool),
    Num(f64),
    Str(String),
    Arr(Vec<Value>),
    Obj(Vec<(String, Value)>),
}

impl Value {
    pub fn get(&self, key: &str) -> &Value {
        match self {
            // exact name first, then any case (Go's encoding/json matches field names case-insensitively)
            Value::Obj(kv) => kv
                .iter()
                .find(|(k, _)| k == key)
                .or_else(|| kv.iter().find(|(k, _)| k.eq_ignore_ascii_case(key)))
                .map(|(_, v)| v)
                .unwrap_or(&Value::Null),
            _ => &Value::Null,
        }
    }
    pub fn as_object(&self) -> &[(String, Value)] {
        match self {
            Value::Obj(kv) => kv,
            _ => &[],
        }
    }
    pub fn as_array(&self) -> &[Value] {
        match self {
            Value::Arr(v) => v,
            _ => &[],
        }
    }
    pub fn as_f64(&self) -> f64 {
        match self {
            Value::Num(x) => *x,
            _ => 0.0,
        }
    }
    #[cfg(test)]
    pub fn as_str(&self) -> &str {
        match self {
            Value::Str(s) => s,
            _ => "",
        }
    }
    #[cfg(test)]
    pub fn is_null(&self) -> bool {
        matches!(self, Value::Null)
    }
}

/// The git-ignored folder with MFEprimer reference runs (CHICHIPES_WORK, default ./work); tests skip without it.
#[cfg(test)]
pub fn work_dir() -> std::path::PathBuf {
    std::env::var_os("CHICHIPES_WORK")
        .map(Into::into)
        .unwrap_or_else(|| std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("work"))
}

#[cfg(test)]
pub fn glob(dir: &std::path::Path, pred: impl Fn(&str) -> bool) -> Vec<std::path::PathBuf> {
    let mut v: Vec<_> = std::fs::read_dir(dir)
        .into_iter()
        .flatten()
        .flatten()
        .map(|e| e.path())
        .filter(|p| p.file_name().and_then(|n| n.to_str()).is_some_and(&pred))
        .collect();
    v.sort();
    v
}

pub fn parse(s: &str) -> Result<Value, String> {
    let mut p = Parser { b: s.as_bytes(), i: 0 };
    let v = p.value()?;
    p.ws();
    if p.i != p.b.len() {
        return Err(format!("trailing data at {}", p.i));
    }
    Ok(v)
}

struct Parser<'a> {
    b: &'a [u8],
    i: usize,
}

impl Parser<'_> {
    fn ws(&mut self) {
        while self.i < self.b.len() && matches!(self.b[self.i], b' ' | b'\t' | b'\n' | b'\r') {
            self.i += 1;
        }
    }
    fn eat(&mut self, c: u8) -> Result<(), String> {
        self.ws();
        if self.b.get(self.i) == Some(&c) {
            self.i += 1;
            Ok(())
        } else {
            Err(format!("expected '{}' at {}", c as char, self.i))
        }
    }
    fn value(&mut self) -> Result<Value, String> {
        self.ws();
        match self.b.get(self.i) {
            Some(b'{') => {
                self.i += 1;
                let mut kv = Vec::new();
                self.ws();
                if self.b.get(self.i) == Some(&b'}') {
                    self.i += 1;
                    return Ok(Value::Obj(kv));
                }
                loop {
                    self.ws();
                    let k = self.string()?;
                    self.eat(b':')?;
                    kv.push((k, self.value()?));
                    self.ws();
                    match self.b.get(self.i) {
                        Some(b',') => self.i += 1,
                        Some(b'}') => {
                            self.i += 1;
                            return Ok(Value::Obj(kv));
                        }
                        _ => return Err(format!("bad object at {}", self.i)),
                    }
                }
            }
            Some(b'[') => {
                self.i += 1;
                let mut v = Vec::new();
                self.ws();
                if self.b.get(self.i) == Some(&b']') {
                    self.i += 1;
                    return Ok(Value::Arr(v));
                }
                loop {
                    v.push(self.value()?);
                    self.ws();
                    match self.b.get(self.i) {
                        Some(b',') => self.i += 1,
                        Some(b']') => {
                            self.i += 1;
                            return Ok(Value::Arr(v));
                        }
                        _ => return Err(format!("bad array at {}", self.i)),
                    }
                }
            }
            Some(b'"') => Ok(Value::Str(self.string()?)),
            Some(b't') if self.b[self.i..].starts_with(b"true") => {
                self.i += 4;
                Ok(Value::Bool(true))
            }
            Some(b'f') if self.b[self.i..].starts_with(b"false") => {
                self.i += 5;
                Ok(Value::Bool(false))
            }
            Some(b'n') if self.b[self.i..].starts_with(b"null") => {
                self.i += 4;
                Ok(Value::Null)
            }
            _ => {
                let st = self.i;
                while self.i < self.b.len() && matches!(self.b[self.i], b'-' | b'+' | b'.' | b'e' | b'E' | b'0'..=b'9') {
                    self.i += 1;
                }
                let t = std::str::from_utf8(&self.b[st..self.i]).unwrap();
                t.parse::<f64>().map(Value::Num).map_err(|_| format!("bad number at {st}"))
            }
        }
    }
    fn string(&mut self) -> Result<String, String> {
        if self.b.get(self.i) != Some(&b'"') {
            return Err(format!("expected string at {}", self.i));
        }
        self.i += 1;
        let mut out = String::new();
        loop {
            let st = self.i;
            while self.i < self.b.len() && self.b[self.i] != b'"' && self.b[self.i] != b'\\' {
                self.i += 1;
            }
            out.push_str(std::str::from_utf8(&self.b[st..self.i]).map_err(|e| e.to_string())?);
            match self.b.get(self.i) {
                Some(b'"') => {
                    self.i += 1;
                    return Ok(out);
                }
                Some(b'\\') => {
                    let c = *self.b.get(self.i + 1).ok_or("bad escape")?;
                    self.i += 2;
                    match c {
                        b'n' => out.push('\n'),
                        b't' => out.push('\t'),
                        b'r' => out.push('\r'),
                        b'b' => out.push('\u{8}'),
                        b'f' => out.push('\u{c}'),
                        b'u' => {
                            let h = std::str::from_utf8(&self.b[self.i..self.i + 4]).map_err(|e| e.to_string())?;
                            let v = u32::from_str_radix(h, 16).map_err(|e| e.to_string())?;
                            self.i += 4;
                            out.push(char::from_u32(v).unwrap_or('\u{fffd}'));
                        }
                        c => out.push(c as char),
                    }
                }
                _ => return Err("unterminated string".into()),
            }
        }
    }
}

/// Writes s as a JSON string the way Go's encoding/json does (HTML-safe: <, >, & escaped; invalid UTF-8 as U+FFFD).
pub fn write_go_string(out: &mut Vec<u8>, s: &[u8]) {
    const HEX: &[u8; 16] = b"0123456789abcdef";
    out.push(b'"');
    let text = String::from_utf8_lossy(s);
    for ch in text.chars() {
        match ch {
            '\\' => out.extend_from_slice(b"\\\\"),
            '"' => out.extend_from_slice(b"\\\""),
            '\u{8}' => out.extend_from_slice(b"\\b"),
            '\u{c}' => out.extend_from_slice(b"\\f"),
            '\n' => out.extend_from_slice(b"\\n"),
            '\r' => out.extend_from_slice(b"\\r"),
            '\t' => out.extend_from_slice(b"\\t"),
            '<' | '>' | '&' | '\u{0}'..='\u{1f}' => {
                let v = ch as u32;
                out.extend_from_slice(b"\\u00");
                out.push(HEX[(v >> 4) as usize]);
                out.push(HEX[(v & 15) as usize]);
            }
            '\u{2028}' => out.extend_from_slice(b"\\u2028"),
            '\u{2029}' => out.extend_from_slice(b"\\u2029"),
            c => {
                let mut b = [0u8; 4];
                out.extend_from_slice(c.encode_utf8(&mut b).as_bytes());
            }
        }
    }
    out.push(b'"');
}
