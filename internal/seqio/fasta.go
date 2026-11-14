// Package seqio reads FASTA databases.
package seqio

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"os"
)

// Record is one FASTA entry. Seq keeps the original letters (case included).
type Record struct {
	ID  string
	Seq []byte
}

// ReadAll reads every record of a FASTA file (the ID is the header up to the first whitespace).
func ReadAll(path string) ([]Record, error) {
	fh, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer fh.Close()
	return Read(fh)
}

// Read parses FASTA from r.
func Read(r io.Reader) ([]Record, error) {
	br := bufio.NewReaderSize(r, 1<<20)
	var out []Record
	var cur *Record
	for {
		line, err := br.ReadBytes('\n')
		line = bytes.TrimRight(line, "\r\n")
		if len(line) > 0 {
			if line[0] == '>' {
				f := bytes.Fields(line[1:])
				if len(f) == 0 {
					return nil, fmt.Errorf("empty FASTA header")
				}
				out = append(out, Record{ID: string(f[0])})
				cur = &out[len(out)-1]
			} else if cur != nil {
				cur.Seq = append(cur.Seq, bytes.TrimSpace(line)...)
			}
		}
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
	}
	return out, nil
}
