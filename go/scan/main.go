// scan: Go port of the whole-genome brute-force scans in scripts/mfe_common.py (_oligo_sites, scan_fr, scan_probe), kept
// result-for-result comparable with them. One JSON request on stdin, one JSON answer on stdout; genomes are scanned in parallel
// and answered in the order they were given.
//
//	{"mode":"fr",    "f":"..","r":"..","three_prime_n":0,"max_mm":2,"max_size":2000,"genomes":[{"id":"x","seq":".."}]}
//	{"mode":"probe", "probe":"..","max_mm":8,"genomes":[{"id":"x","seq":"..","amp_start":null|int,"amp_end":null|int}]}
//	{"mode":"sites", "oligo":"..","max_mm":2,"genomes":[...]}      (raw sites; used by the checks)
//
// `scan -selftest` runs the built-in unit checks.
package main

import (
	"encoding/json"
	"fmt"
	"os"
	"runtime"
	"strconv"
	"strings"
	"sync"
)

var iupac = map[byte]string{
	'A': "A", 'C': "C", 'G': "G", 'T': "T", 'R': "AG", 'Y': "CT", 'S': "GC", 'W': "AT", 'K': "GT", 'M': "AC",
	'B': "CGT", 'D': "AGT", 'H': "ACT", 'V': "ACG", 'N': "ACGT",
}

var complement = map[byte]byte{
	'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K',
	'D': 'H', 'H': 'D', 'B': 'V', 'V': 'B',
}

func upper(c byte) byte {
	if c >= 'a' && c <= 'z' {
		return c - 32
	}
	return c
}

func up(s string) string { return strings.ToUpper(s) }

func revComp(s string) string {
	out := make([]byte, len(s))
	for i := 0; i < len(s); i++ {
		c, ok := complement[upper(s[len(s)-1-i])]
		if !ok {
			c = 'N'
		}
		out[i] = c
	}
	return string(out)
}

// allowedTable[j][b] is true when genome byte b satisfies the oligo code at position j (unknown codes accept ACGT, as in Python).
func allowedTable(q string) [][256]bool {
	t := make([][256]bool, len(q))
	for j := 0; j < len(q); j++ {
		set, ok := iupac[q[j]]
		if !ok {
			set = "ACGT"
		}
		for k := 0; k < len(set); k++ {
			t[j][set[k]] = true
		}
	}
	return t
}

type Site struct {
	Start  int    `json:"start"` // 0-based
	End    int    `json:"end"`   // exclusive
	Strand string `json:"strand"`
	MM     int    `json:"mm"`
}

// oligoSites is _oligo_sites(): every window with <= maxMM mismatches, "+" strand first, then "-", each by ascending start.
func oligoSites(oligo, genome string, maxMM int) []Site {
	var sites []Site
	n := len(genome)
	for _, st := range []struct {
		strand string
		q      string
	}{{"+", oligo}, {"-", revComp(oligo)}} {
		tab := allowedTable(st.q)
		L := len(st.q)
		for i := 0; i+L <= n; i++ {
			m := 0
			j := 0
			for ; j < L; j++ {
				if !tab[j][genome[i+j]] {
					m++
					if m > maxMM {
						break
					}
				}
			}
			if j == L {
				sites = append(sites, Site{i, i + L, st.strand, m})
			}
		}
	}
	return sites
}

func iupacMismatch(q, t string) int {
	n := len(q)
	if len(t) < n {
		n = len(t)
	}
	m := 0
	for i := 0; i < n; i++ {
		set, ok := iupac[upper(q[i])]
		if !ok {
			set = "ACGT"
		}
		if !strings.Contains(set, string(upper(t[i]))) {
			m++
		}
	}
	d := len(q) - len(t)
	if d < 0 {
		d = -d
	}
	return m + d
}

func terminalOK(query, target string, isFP bool, n int) bool {
	if len(query) < n {
		n = len(query)
	}
	if len(target) < n {
		n = len(target)
	}
	if n <= 0 {
		return true
	}
	var q, t string
	if isFP {
		q, t = query[len(query)-n:], target[len(target)-n:]
	} else {
		q, t = query[:n], target[:n]
	}
	return iupacMismatch(q, t) == 0
}

type Pos struct {
	Start  int    `json:"start"`
	End    int    `json:"end"`
	Strand string `json:"strand"`
}

type FR struct {
	FMM   int    `json:"f_mismatch"`
	RMM   int    `json:"r_mismatch"`
	Total int    `json:"total"`
	F3P   bool   `json:"f_3p_ok"`
	R3P   bool   `json:"r_3p_ok"`
	AmpS  string `json:"amp_start"`
	AmpE  string `json:"amp_end"`
	FTm   string `json:"f_tm"` // no mfeprimer thermodynamics for a scan hit: always "" (as in Python)
	FDg   string `json:"f_dg"`
	RTm   string `json:"r_tm"`
	RDg   string `json:"r_dg"`
	FPos  Pos    `json:"f_pos"`
	RPos  Pos    `json:"r_pos"`
}

// scanFR is scan_fr(): best F/R pair by (3' ends intact, total mismatches, shorter amplicon); the first of equal keys wins.
func scanFR(fSeq, rSeq, genome string, threePrimeN, maxMM, maxSize int) *FR {
	g := up(genome)
	fU, rU := up(fSeq), up(rSeq)
	fs, rs := oligoSites(fU, g, maxMM), oligoSites(rU, g, maxMM)
	var best *FR
	var bestKey [3]int
	for _, f := range fs {
		for _, r := range rs {
			var left, right Site
			switch {
			case f.Strand == "+" && r.Strand == "-" && f.Start <= r.Start && r.End > f.Start && r.End-f.Start <= maxSize:
				left, right = f, r
			case r.Strand == "+" && f.Strand == "-" && r.Start <= f.Start && f.End > r.Start && f.End-r.Start <= maxSize:
				left, right = r, f
			default:
				continue
			}
			fq, rq := fU, rU
			if f.Strand == "-" {
				fq = revComp(fU)
			}
			if r.Strand == "-" {
				rq = revComp(rU)
			}
			fok := terminalOK(fq, g[f.Start:f.End], f.Strand == "+", threePrimeN)
			rok := terminalOK(rq, g[r.Start:r.End], r.Strand == "+", threePrimeN)
			k0 := 1
			if fok && rok {
				k0 = 0
			}
			key := [3]int{k0, f.MM + r.MM, right.End - left.Start}
			if best == nil || key[0] < bestKey[0] || (key[0] == bestKey[0] && (key[1] < bestKey[1] || (key[1] == bestKey[1] && key[2] < bestKey[2]))) {
				bestKey = key
				best = &FR{f.MM, r.MM, f.MM + r.MM, fok, rok, strconv.Itoa(left.Start + 1), strconv.Itoa(right.End), "", "", "", "",
					Pos{f.Start + 1, f.End, f.Strand}, Pos{r.Start + 1, r.End, r.Strand}}
			}
		}
	}
	return best
}

type Probe struct {
	MM    int    `json:"mismatch"`
	Start string `json:"start"`
	End   string `json:"end"`
	Tm    string `json:"tm"`
	Dg    string `json:"dg"`
}

// scanProbe is scan_probe(): the best window (fewest mismatches, first wins) inside [ampStart, ampEnd] (1-based inclusive) if given.
func scanProbe(probe, genome string, ampS, ampE *int, maxMM int) *Probe {
	var best *Probe
	for _, s := range oligoSites(up(probe), up(genome), maxMM) {
		if ampS != nil && !(*ampS <= s.Start+1 && s.End <= *ampE) {
			continue
		}
		if best == nil || s.MM < best.MM {
			best = &Probe{s.MM, strconv.Itoa(s.Start + 1), strconv.Itoa(s.End), "", ""}
		}
	}
	return best
}

type Genome struct {
	ID       string `json:"id"`
	Seq      string `json:"seq"`
	AmpStart *int   `json:"amp_start"`
	AmpEnd   *int   `json:"amp_end"`
}

type Request struct {
	Mode        string   `json:"mode"`
	F           string   `json:"f"`
	R           string   `json:"r"`
	Probe       string   `json:"probe"`
	Oligo       string   `json:"oligo"`
	ThreePrimeN int      `json:"three_prime_n"`
	MaxMM       int      `json:"max_mm"`
	MaxSize     int      `json:"max_size"`
	Genomes     []Genome `json:"genomes"`
}

type Result struct {
	ID    string `json:"id"`
	Hit   any    `json:"hit"`
	Sites []Site `json:"sites,omitempty"`
}

func run(req Request) ([]Result, error) {
	switch req.Mode {
	case "fr", "probe", "sites":
	default:
		return nil, fmt.Errorf("unknown mode %q", req.Mode)
	}
	out := make([]Result, len(req.Genomes))
	workers := runtime.NumCPU()
	jobs := make(chan int)
	var wg sync.WaitGroup
	for w := 0; w < workers; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for i := range jobs {
				g := req.Genomes[i]
				res := Result{ID: g.ID}
				switch req.Mode {
				case "fr":
					if h := scanFR(req.F, req.R, g.Seq, req.ThreePrimeN, req.MaxMM, req.MaxSize); h != nil {
						res.Hit = h
					}
				case "probe":
					if h := scanProbe(req.Probe, g.Seq, g.AmpStart, g.AmpEnd, req.MaxMM); h != nil {
						res.Hit = h
					}
				case "sites":
					res.Sites = oligoSites(up(req.Oligo), up(g.Seq), req.MaxMM)
				}
				out[i] = res
			}
		}()
	}
	for i := range req.Genomes {
		jobs <- i
	}
	close(jobs)
	wg.Wait()
	return out, nil
}

func main() {
	if len(os.Args) > 1 && os.Args[1] == "-selftest" {
		os.Exit(selftest())
	}
	var req Request
	if err := json.NewDecoder(os.Stdin).Decode(&req); err != nil {
		fmt.Fprintln(os.Stderr, "bad request:", err)
		os.Exit(2)
	}
	res, err := run(req)
	if err != nil {
		fmt.Fprintln(os.Stderr, "error:", err)
		os.Exit(1)
	}
	_ = json.NewEncoder(os.Stdout).Encode(map[string]any{"results": res})
}
