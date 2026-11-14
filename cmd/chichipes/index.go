package main

import (
	"bufio"
	"encoding/json"
	"flag"
	"fmt"
	"os"
	"strconv"
	"strings"
)

type faiEntry struct {
	id, desc                string
	size, offset, lineBases int64
	lineWidth               int64
}

// runIndex writes the files `mfeprimer index` leaves next to the database, because the pipeline's Snakemake rules
// declare them as outputs and the web app reads the .fai. The search itself never needs them: .fai (samtools
// layout) and .json are real, .primerqc/.primerqc.fai are empty placeholders for mfeprimer's own k-mer index.
func runIndex(args []string) {
	fs := flag.NewFlagSet("index", flag.ExitOnError)
	in := fs.String("i", "", "database FASTA")
	fs.Int("cpu", 1, "")
	fs.Int("c", 1, "")
	fs.Int("k", 9, "")
	check(fs.Parse(args))
	if *in == "" {
		fmt.Fprintln(os.Stderr, "usage: chichipes index -i db.fa")
		os.Exit(2)
	}
	entries, err := scanFai(*in)
	check(err)
	var fai strings.Builder
	js := make(map[string]map[string]any, len(entries))
	for i, e := range entries {
		fmt.Fprintf(&fai, "%s\t%d\t%d\t%d\t%d\n", e.id, e.size, e.offset, e.lineBases, e.lineWidth)
		js[strconv.Itoa(i)] = map[string]any{"id": e.id, "desc": e.desc, "size": e.size}
	}
	b, err := json.Marshal(js)
	check(err)
	check(os.WriteFile(*in+".fai", []byte(fai.String()), 0o644))
	check(os.WriteFile(*in+".json", b, 0o644))
	check(os.WriteFile(*in+".log", []byte(fmt.Sprintf("%s\nindexed %d records of %s\n", version, len(entries), *in)), 0o644))
	check(os.WriteFile(*in+".primerqc", nil, 0o644))
	check(os.WriteFile(*in+".primerqc.fai", nil, 0o644))
}

func scanFai(path string) ([]faiEntry, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	r := bufio.NewReaderSize(f, 1<<20)
	var out []faiEntry
	var cur *faiEntry
	var pos int64
	for {
		line, err := r.ReadString('\n')
		if len(line) > 0 {
			n := int64(len(line))
			body := strings.TrimRight(line, "\r\n")
			if strings.HasPrefix(body, ">") {
				id, desc, _ := strings.Cut(body[1:], " ")
				out = append(out, faiEntry{id: id, desc: desc, offset: pos + n})
				cur = &out[len(out)-1]
			} else if cur != nil && len(body) > 0 {
				if cur.lineBases == 0 {
					cur.lineBases, cur.lineWidth = int64(len(body)), n
				}
				cur.size += int64(len(body))
			}
			pos += n
		}
		if err != nil {
			break
		}
	}
	return out, nil
}
