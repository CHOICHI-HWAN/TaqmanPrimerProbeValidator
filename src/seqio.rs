//! FASTA databases (read whole) and FASTA/FASTQ read files (streamed, gzip or not).

use flate2::bufread::MultiGzDecoder;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};

/// One FASTA entry. seq keeps the original letters (case included).
pub struct Record {
    pub id: String,
    pub seq: Vec<u8>,
}

#[inline]
fn is_space(b: u8) -> bool {
    matches!(b, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c)
}

fn trim(mut s: &[u8]) -> &[u8] {
    while let [f, rest @ ..] = s {
        if is_space(*f) {
            s = rest;
        } else {
            break;
        }
    }
    while let [rest @ .., l] = s {
        if is_space(*l) {
            s = rest;
        } else {
            break;
        }
    }
    s
}

fn trim_eol(mut s: &[u8]) -> &[u8] {
    while let [rest @ .., l] = s {
        if *l == b'\n' || *l == b'\r' {
            s = rest;
        } else {
            break;
        }
    }
    s
}

/// The header up to the first whitespace.
fn header_id(h: &[u8]) -> Option<String> {
    h.split(|&b| is_space(b)).find(|w| !w.is_empty()).map(|w| String::from_utf8_lossy(w).into_owned())
}

/// Reads every record of a FASTA file (the ID is the header up to the first whitespace).
pub fn read_all(path: &str) -> Result<Vec<Record>, String> {
    let data = std::fs::read(path).map_err(|e| format!("open {path}: {e}"))?;
    let mut out: Vec<Record> = Vec::new();
    for line in data.split(|&b| b == b'\n') {
        let line = trim_eol(line);
        if line.is_empty() {
            continue;
        }
        if line[0] == b'>' {
            let id = header_id(&line[1..]).ok_or("empty FASTA header")?;
            out.push(Record { id, seq: Vec::new() });
        } else if let Some(cur) = out.last_mut() {
            cur.seq.extend_from_slice(trim(line));
        }
    }
    Ok(out)
}

/// Streams records from a FASTA or FASTQ file, gzip-compressed or not (both detected from the content).
pub struct Reader {
    br: Box<dyn BufRead + Send>,
    fastq: bool,
    path: String,
    line: usize,
    buf: Vec<u8>,
}

/// Opens path ("-" is standard input) for streaming.
pub fn open(path: &str) -> Result<Reader, String> {
    let f: Box<dyn Read + Send> = if path == "-" {
        Box::new(io::stdin())
    } else {
        Box::new(File::open(path).map_err(|e| format!("open {path}: {e}"))?)
    };
    let mut br = BufReader::with_capacity(1 << 20, f);
    let magic = br.fill_buf().map_err(|e| format!("{path}: {e}"))?;
    let mut br: Box<dyn BufRead + Send> = if magic.len() >= 2 && magic[0] == 0x1f && magic[1] == 0x8b {
        Box::new(BufReader::with_capacity(1 << 20, MultiGzDecoder::new(br)))
    } else {
        Box::new(br)
    };
    let mut fastq = false;
    loop {
        let b = match br.fill_buf() {
            Ok(b) if !b.is_empty() => b[0],
            Ok(_) => break, // empty file: next() returns None
            Err(e) => return Err(format!("{path}: {e}")),
        };
        match b {
            b'@' => {
                fastq = true;
                break;
            }
            b'>' => break,
            b'\n' | b'\r' | b' ' | b'\t' => br.consume(1),
            _ => return Err(format!("{path}: neither FASTA ('>') nor FASTQ ('@')")),
        }
    }
    Ok(Reader { br, fastq, path: path.to_string(), line: 0, buf: Vec::new() })
}

impl Reader {
    /// The next line without its line end, or None at the end of the file.
    fn read_line(&mut self) -> Result<Option<Vec<u8>>, String> {
        self.buf.clear();
        let n = self.br.read_until(b'\n', &mut self.buf).map_err(|e| format!("{}: {e}", self.path))?;
        self.line += 1;
        if n == 0 {
            return Ok(None);
        }
        Ok(Some(trim_eol(&self.buf).to_vec()))
    }

    fn peek(&mut self) -> Option<u8> {
        match self.br.fill_buf() {
            Ok(b) if !b.is_empty() => Some(b[0]),
            _ => None,
        }
    }

    /// The next record (ID = header up to the first whitespace), or None at the end.
    pub fn next(&mut self) -> Result<Option<Record>, String> {
        if self.fastq {
            self.next_fastq()
        } else {
            self.next_fasta()
        }
    }

    fn next_fastq(&mut self) -> Result<Option<Record>, String> {
        let h = loop {
            match self.read_line()? {
                None => return Ok(None),
                Some(l) if !l.is_empty() => break l,
                Some(_) => {}
            }
        };
        if h[0] != b'@' {
            return Err(format!("{}:{}: FASTQ record does not start with '@'", self.path, self.line));
        }
        let seq = self.read_line()?.ok_or_else(|| format!("{}:{}: truncated FASTQ record", self.path, self.line))?;
        match self.read_line()? {
            Some(p) if !p.is_empty() && p[0] == b'+' => {}
            _ => {
                return Err(format!(
                    "{}:{}: expected '+' line (multi-line FASTQ is not supported)",
                    self.path, self.line
                ))
            }
        }
        if self.read_line()?.is_none() {
            return Err(format!("{}:{}: truncated FASTQ record", self.path, self.line));
        }
        Ok(Some(Record { id: header_id(&h[1..]).unwrap_or_default(), seq }))
    }

    fn next_fasta(&mut self) -> Result<Option<Record>, String> {
        let h = loop {
            match self.read_line()? {
                None => return Ok(None),
                Some(l) if !l.is_empty() && l[0] == b'>' => break l,
                Some(_) => {}
            }
        };
        let mut rec = Record { id: header_id(&h[1..]).unwrap_or_default(), seq: Vec::new() };
        loop {
            match self.peek() {
                None | Some(b'>') => return Ok(Some(rec)),
                _ => {}
            }
            match self.read_line()? {
                Some(l) => rec.seq.extend_from_slice(trim(&l)),
                None => return Ok(Some(rec)),
            }
        }
    }
}
