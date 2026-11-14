package seqio

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"fmt"
	"io"
	"os"
)

// Reader streams records from a FASTA or FASTQ file, gzip-compressed or not (both detected from the content), so
// read files of any size can be scanned without loading them.
type Reader struct {
	br    *bufio.Reader
	close func() error
	fastq bool
	next  []byte // a FASTA header line already read
	path  string
	line  int
}

// Open opens path ("-" is standard input) for streaming.
func Open(path string) (*Reader, error) {
	var f io.ReadCloser = os.Stdin
	if path != "-" {
		fh, err := os.Open(path)
		if err != nil {
			return nil, err
		}
		f = fh
	}
	br := bufio.NewReaderSize(f, 1<<20)
	closer := f.Close
	if magic, _ := br.Peek(2); len(magic) == 2 && magic[0] == 0x1f && magic[1] == 0x8b {
		gz, err := gzip.NewReader(br)
		if err != nil {
			f.Close()
			return nil, fmt.Errorf("%s: %w", path, err)
		}
		gz.Multistream(true)
		br = bufio.NewReaderSize(gz, 1<<20)
		closer = func() error { gz.Close(); return f.Close() }
	}
	r := &Reader{br: br, close: closer, path: path}
	for {
		b, err := br.Peek(1)
		if err != nil {
			return r, nil // empty file: Next returns io.EOF
		}
		switch b[0] {
		case '@':
			r.fastq = true
			return r, nil
		case '>':
			return r, nil
		case '\n', '\r', ' ', '\t':
			br.ReadByte()
		default:
			f.Close()
			return nil, fmt.Errorf("%s: neither FASTA ('>') nor FASTQ ('@')", path)
		}
	}
}

// Close releases the file.
func (r *Reader) Close() error { return r.close() }

func (r *Reader) readLine() ([]byte, error) {
	l, err := r.br.ReadBytes('\n')
	r.line++
	if err == io.EOF && len(l) > 0 {
		err = nil
	}
	return bytes.TrimRight(l, "\r\n"), err
}

// Next returns the next record (ID = header up to the first whitespace), or io.EOF.
func (r *Reader) Next() (Record, error) {
	if r.fastq {
		return r.nextFastq()
	}
	return r.nextFasta()
}

func headerID(h []byte) string {
	if f := bytes.Fields(h); len(f) > 0 {
		return string(f[0])
	}
	return ""
}

func (r *Reader) nextFastq() (Record, error) {
	var h []byte
	for {
		l, err := r.readLine()
		if err != nil {
			return Record{}, err
		}
		if len(l) > 0 {
			h = l
			break
		}
	}
	if h[0] != '@' {
		return Record{}, fmt.Errorf("%s:%d: FASTQ record does not start with '@'", r.path, r.line)
	}
	seq, err := r.readLine()
	if err != nil {
		return Record{}, fmt.Errorf("%s:%d: truncated FASTQ record", r.path, r.line)
	}
	plus, err := r.readLine()
	if err != nil || len(plus) == 0 || plus[0] != '+' {
		return Record{}, fmt.Errorf("%s:%d: expected '+' line (multi-line FASTQ is not supported)", r.path, r.line)
	}
	if _, err := r.readLine(); err != nil {
		return Record{}, fmt.Errorf("%s:%d: truncated FASTQ record", r.path, r.line)
	}
	return Record{ID: headerID(h[1:]), Seq: append([]byte(nil), seq...)}, nil
}

func (r *Reader) nextFasta() (Record, error) {
	h := r.next
	r.next = nil
	for h == nil {
		l, err := r.readLine()
		if err != nil {
			return Record{}, err
		}
		if len(l) > 0 && l[0] == '>' {
			h = append([]byte(nil), l...)
		}
	}
	rec := Record{ID: headerID(h[1:])}
	for {
		b, err := r.br.Peek(1)
		if err != nil || b[0] == '>' {
			return rec, nil
		}
		l, err := r.readLine()
		rec.Seq = append(rec.Seq, bytes.TrimSpace(l)...)
		if err != nil {
			return rec, nil
		}
	}
}
