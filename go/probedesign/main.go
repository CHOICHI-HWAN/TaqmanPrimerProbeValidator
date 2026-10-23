// probedesign: Go port of the "own" (rule-based) path of scripts/probe_designer.py, kept bit-for-bit comparable with it.
// It reads one JSON request on stdin and writes the candidates as JSON on stdout.
//
//	{"consensus": "...", "f_col": [a,b], "r_col": [a,b], "dots": ["...", ...],
//	 "min_len": 18, "max_len": 30, "tm_range": null|[lo,hi], "gc_range": [30,80], "top": 10, "max_iupac": 2, "tm_above": null|[lo,hi]}
//
// Run `probedesign -selftest` for the built-in unit checks.
package main

import (
	"encoding/json"
	"fmt"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"
	"sync"
)

var nn = map[string][2]float64{
	"AA": {-7.9, -22.2}, "TT": {-7.9, -22.2}, "AT": {-7.2, -20.4}, "TA": {-7.2, -21.3}, "CA": {-8.5, -22.7},
	"TG": {-8.5, -22.7}, "GT": {-8.4, -22.4}, "AC": {-8.4, -22.4}, "CT": {-7.8, -21.0}, "AG": {-7.8, -21.0},
	"GA": {-8.2, -22.2}, "TC": {-8.2, -22.2}, "CG": {-10.6, -27.2}, "GC": {-9.8, -24.4}, "GG": {-8.0, -19.9},
	"CC": {-8.0, -19.9},
}

var iupac = map[byte]string{
	'A': "A", 'C': "C", 'G': "G", 'T': "T", 'R': "AG", 'Y': "CT", 'S': "GC", 'W': "AT", 'K': "GT", 'M': "AC",
	'B': "CGT", 'D': "AGT", 'H': "ACT", 'V': "ACG", 'N': "ACGT",
}

var complement = map[byte]byte{
	'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K',
	'D': 'H', 'H': 'D', 'B': 'V', 'V': 'B',
}

// tmNN is tm_nn(): nearest-neighbour Tm of a plain ACGT oligo (250 nM oligo, 50 mM Na+). ok=false where Python raises KeyError.
func tmNN(seq string) (float64, bool) {
	seq = strings.ToUpper(seq)
	if len(seq) < 1 {
		return 0, false
	}
	dhs, dss := make([]float64, 0, len(seq)), make([]float64, 0, len(seq))
	for i := 0; i+1 < len(seq); i++ {
		p, found := nn[seq[i:i+2]]
		if !found {
			return 0, false
		}
		dhs, dss = append(dhs, p[0]), append(dss, p[1])
	}
	dh, ds := pySum(dhs), pySum(dss)
	for _, e := range []byte{seq[0], seq[len(seq)-1]} {
		if e == 'G' || e == 'C' {
			dh, ds = dh+2.3, ds+4.1
		} else {
			dh, ds = dh+0.1, ds-2.8
		}
	}
	ds += 0.368 * float64(len(seq)-1) * math.Log(50.0/1000.0)
	return dh*1000.0/(ds+1.987*math.Log(250*1e-9/4.0)) - 273.15, true
}

// pySum reproduces CPython's sum() over floats (3.12+: Neumaier-compensated), so means and sums agree to the last bit.
func pySum(xs []float64) float64 {
	var sum, c float64
	for _, x := range xs {
		t := sum + x
		if math.Abs(sum) >= math.Abs(x) {
			c += (sum - t) + x
		} else {
			c += (x - t) + sum
		}
		sum = t
	}
	if c != 0 && !math.IsInf(c, 0) && !math.IsNaN(c) {
		sum += c
	}
	return sum
}

func gcPercent(s string) float64 {
	n := 0
	for i := 0; i < len(s); i++ {
		if s[i] == 'G' || s[i] == 'C' {
			n++
		}
	}
	return 100.0 * float64(n) / float64(len(s))
}

func longestRun(s string) int {
	best, cur := 1, 1
	for i := 1; i < len(s); i++ {
		if s[i] == s[i-1] {
			cur++
		} else {
			cur = 1
		}
		if cur > best {
			best = cur
		}
	}
	return best
}

// resolutions is resolutions(): every ACGT oligo an IUPAC oligo stands for; empty if any character is not an IUPAC code.
func resolutions(seq string, limit int) []string {
	if limit > 0 {
		n := 1
		for i := 0; i < len(seq); i++ {
			k := len(iupac[upper(seq[i])])
			if k < 1 {
				k = 1
			}
			if n *= k; n > limit {
				return nil
			}
		}
	}
	outs := []string{""}
	for i := 0; i < len(seq); i++ {
		b, ok := iupac[upper(seq[i])]
		if !ok {
			return nil
		}
		next := make([]string, 0, len(outs)*len(b))
		for _, o := range outs {
			for j := 0; j < len(b); j++ {
				next = append(next, o+string(b[j]))
			}
		}
		outs = next
	}
	return outs
}

func upper(c byte) byte {
	if c >= 'a' && c <= 'z' {
		return c - 32
	}
	return c
}

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

// meanTmGC averages Tm and GC over the resolutions (as the Python does for IUPAC positions).
func meanTmGC(res []string) (float64, float64, bool) {
	tms, gcs := make([]float64, 0, len(res)), make([]float64, 0, len(res))
	for _, r := range res {
		t, ok := tmNN(r)
		if !ok {
			return 0, 0, false
		}
		tms, gcs = append(tms, t), append(gcs, gcPercent(r))
	}
	n := float64(len(res))
	return pySum(tms) / n, pySum(gcs) / n, true
}

// ruleProblems is rule_problems(): true when the oligo is NOT a usable probe.
func ruleProblems(seq string, tmLo, tmHi, gcLo, gcHi float64) bool {
	res := resolutions(seq, 64)
	if len(res) == 0 || len(res) > 64 {
		return true
	}
	tm, gc, ok := meanTmGC(res)
	if !ok {
		return true
	}
	if !(tmLo <= tm && tm <= tmHi) || !(gcLo <= gc && gc <= gcHi) {
		return true
	}
	for _, r := range res {
		if r[0] == 'G' || longestRun(r) >= 4 {
			return true
		}
	}
	return false
}

var maskBit = map[byte]uint8{'A': 1, 'C': 2, 'G': 4, 'T': 8}

// codeByMask[mask] is the IUPAC letter for a set of exactly two bases (mask = bits of A,C,G,T).
var codeByMask = func() map[uint8]byte {
	m := map[uint8]byte{}
	for c, v := range iupac {
		if len(v) == 2 {
			m[maskBit[v[0]]|maskBit[v[1]]] = c
		}
	}
	return m
}()

type mm struct {
	pos  int
	base byte
}

// degenerateWindow is degenerate_window(): greedily add up to maxIUPAC two-base codes where the most genomes still need them.
func degenerateWindow(window string, rowMM [][]mm, maxIUPAC int) (string, []int) {
	seq := []byte(window)
	allowed := make([]uint8, len(window)) // 0 = untouched (only the consensus base), else the two-base mask
	left := make([][]mm, len(rowMM))
	copy(left, rowMM)
	acgt := func(c byte) bool { return c == 'A' || c == 'C' || c == 'G' || c == 'T' }
	type key struct {
		pos  int
		base byte
	}
	for k := maxIUPAC; k > 0; k-- {
		var order []key
		score := map[key]int{}
		for _, m := range left {
			if len(m) > 0 && len(m) <= k {
				for _, x := range m {
					if acgt(x.base) && x.base != window[x.pos] && allowed[x.pos] == 0 && acgt(window[x.pos]) {
						kk := key{x.pos, x.base}
						if _, seen := score[kk]; !seen {
							order = append(order, kk)
						}
						score[kk]++
					}
				}
			}
		}
		if len(order) == 0 {
			break
		}
		best := order[0] // Python max(): first of the equal maxima in insertion order
		for _, kk := range order[1:] {
			if score[kk] > score[best] || (score[kk] == score[best] && kk.pos < best.pos) {
				best = kk
			}
		}
		allowed[best.pos] = maskBit[window[best.pos]] | maskBit[best.base]
		seq[best.pos] = codeByMask[allowed[best.pos]]
		for i, m := range left {
			hit := false
			for _, x := range m {
				if x.pos == best.pos && x.base == best.base {
					hit = true
					break
				}
			}
			if !hit {
				continue
			}
			out := make([]mm, 0, len(m))
			for _, x := range m {
				if !(x.pos == best.pos && x.base == best.base) {
					out = append(out, x)
				}
			}
			left[i] = out
		}
	}
	cnt := make([]int, len(left))
	for i, m := range left {
		cnt[i] = len(m)
	}
	return string(seq), cnt
}

type Request struct {
	Consensus string      `json:"consensus"`
	FCol      [2]int      `json:"f_col"`
	RCol      [2]int      `json:"r_col"`
	Dots      []string    `json:"dots"`
	MinLen    int         `json:"min_len"`
	MaxLen    int         `json:"max_len"`
	TmRange   *[2]float64 `json:"tm_range"`
	GcRange   [2]float64  `json:"gc_range"`
	Top       int         `json:"top"`
	MaxIUPAC  int         `json:"max_iupac"`
	TmAbove   *[2]float64 `json:"tm_above"`
}

type Cand struct {
	Seq      string  `json:"seq"`
	Strand   string  `json:"strand"`
	ColStart int     `json:"col_start"`
	ColEnd   int     `json:"col_end"`
	Length   int     `json:"length"`
	NIUPAC   int     `json:"n_iupac"`
	Tm       float64 `json:"tm"`
	GC       float64 `json:"gc"`
	Cov0     int     `json:"msa_cov0"`
	Cov1     int     `json:"msa_cov1"`
	N        int     `json:"msa_n"`
	Source   string  `json:"source"`
}

type Result struct {
	PrimerTm   float64    `json:"primer_tm"`
	TmRange    [2]float64 `json:"tm_range"`
	Window     [2]int     `json:"window"`
	NTotal     int        `json:"n_total"`
	NRejected  int        `json:"n_rejected"`
	NRows      int        `json:"n_rows"`
	Candidates []Cand     `json:"candidates"`
}

const probeTmTarget = 6.0

func round1(x float64) float64 {
	v, _ := strconv.ParseFloat(strconv.FormatFloat(x, 'f', 1, 64), 64)
	return v
}

func design(r Request) (*Result, error) {
	cons := strings.ToUpper(r.Consensus)
	pt := math.Inf(-1)
	for _, c := range [][2]int{r.FCol, r.RCol} {
		hi := c[1] + 1
		if hi > len(cons) {
			hi = len(cons)
		}
		w := strings.ReplaceAll(cons[c[0]:hi], "-", "")
		if w == "" {
			continue
		}
		t, ok := tmNN(w)
		if !ok {
			return nil, fmt.Errorf("primer window is not plain ACGT")
		}
		pt = math.Max(pt, t)
	}
	if math.IsInf(pt, -1) {
		return nil, fmt.Errorf("the consensus has no bases under the F/R primer windows")
	}
	above := [2]float64{3.0, 10.0}
	if r.TmAbove != nil {
		above = *r.TmAbove
	}
	tmRange := [2]float64{pt + above[0], pt + above[1]}
	if r.TmRange != nil {
		tmRange = *r.TmRange
	}
	lo := minInt(r.FCol[1], r.RCol[1]) + 1
	hi := maxInt(r.FCol[0], r.RCol[0]) - 1
	nRows := len(r.Dots)
	type part struct {
		out      []Cand
		rejected int
	}
	nLen := r.MaxLen - r.MinLen + 1
	if nLen < 0 {
		nLen = 0
	}
	parts := make([]part, nLen)
	var wg sync.WaitGroup
	for li := 0; li < nLen; li++ {
		wg.Add(1)
		go func(li int) { // every window length is independent; results are joined in length order, so the output is unchanged
			defer wg.Done()
			length := r.MinLen + li
			p := &parts[li]
			for s := lo; s <= hi-length+1; s++ {
				e := s + length
				if e > len(cons) {
					e = len(cons)
				}
				window := cons[s:e]
				rowMM := make([][]mm, nRows)
				for ri, d := range r.Dots {
					end := e
					if end > len(d) {
						end = len(d)
					}
					for i := s; i < end; i++ {
						if d[i] != '.' {
							rowMM[ri] = append(rowMM[ri], mm{i - s, d[i]})
						}
					}
				}
				window, left := degenerateWindow(window, rowMM, r.MaxIUPAC)
				nIUPAC := 0
				for i := 0; i < len(window); i++ {
					if c := window[i]; !(c == 'A' || c == 'C' || c == 'G' || c == 'T') {
						nIUPAC++
					}
				}
				cov0, cov1 := 0, 0
				for _, m := range left {
					if m == 0 {
						cov0++
					}
					if m <= 1 {
						cov1++
					}
				}
				for _, strand := range []string{"+", "-"} {
					oligo := window
					if strand == "-" {
						oligo = revComp(window)
					}
					if ruleProblems(oligo, tmRange[0], tmRange[1], r.GcRange[0], r.GcRange[1]) {
						p.rejected++
						continue
					}
					res := resolutions(oligo, 0)
					tm, gc, _ := meanTmGC(res)
					p.out = append(p.out, Cand{Seq: oligo, Strand: strand, ColStart: s + 1, ColEnd: e, Length: length, NIUPAC: nIUPAC,
						Tm: round1(tm), GC: round1(gc), Cov0: cov0, Cov1: cov1, N: nRows, Source: "own"})
				}
			}
		}(li)
	}
	wg.Wait()
	var out []Cand
	rejected := 0
	for _, p := range parts {
		out = append(out, p.out...)
		rejected += p.rejected
	}
	picked := pick(out, r.Top, pt)
	merged := map[string]int{}
	var uniq []Cand
	for _, c := range picked {
		if i, seen := merged[c.Seq]; seen {
			uniq[i].Source = "own+primer3"
		} else {
			merged[c.Seq] = len(uniq)
			uniq = append(uniq, c)
		}
	}
	sort.SliceStable(uniq, func(i, j int) bool {
		a, b := uniq[i], uniq[j]
		if a.Cov0 != b.Cov0 {
			return a.Cov0 > b.Cov0
		}
		if a.Cov1 != b.Cov1 {
			return a.Cov1 > b.Cov1
		}
		if a.NIUPAC != b.NIUPAC {
			return a.NIUPAC < b.NIUPAC
		}
		return a.Length < b.Length
	})
	if uniq == nil {
		uniq = []Cand{}
	}
	return &Result{PrimerTm: pt, TmRange: tmRange, Window: [2]int{lo + 1, hi + 1}, NTotal: len(out), NRejected: rejected, NRows: nRows, Candidates: uniq}, nil
}

// pick is _pick(): best first, no near-duplicate (same strand, start within 3 columns) of a better one.
func pick(c []Cand, quota int, primerTm float64) []Cand {
	sort.SliceStable(c, func(i, j int) bool {
		a, b := c[i], c[j]
		if a.Cov0 != b.Cov0 {
			return a.Cov0 > b.Cov0
		}
		if a.Cov1 != b.Cov1 {
			return a.Cov1 > b.Cov1
		}
		if a.NIUPAC != b.NIUPAC {
			return a.NIUPAC < b.NIUPAC
		}
		da, db := math.Abs(a.Tm-(primerTm+probeTmTarget)), math.Abs(b.Tm-(primerTm+probeTmTarget))
		if da != db {
			return da < db
		}
		return a.Length < b.Length
	})
	var picked []Cand
	for _, x := range c {
		if len(picked) >= quota {
			break
		}
		ok := true
		for _, p := range picked {
			if p.Strand == x.Strand && absInt(p.ColStart-x.ColStart) <= 3 {
				ok = false
				break
			}
		}
		if ok {
			picked = append(picked, x)
		}
	}
	return picked
}

func minInt(a, b int) int {
	if a < b {
		return a
	}
	return b
}
func maxInt(a, b int) int {
	if a > b {
		return a
	}
	return b
}
func absInt(a int) int {
	if a < 0 {
		return -a
	}
	return a
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
	res, err := design(req)
	if err != nil {
		fmt.Fprintln(os.Stderr, "error:", err)
		os.Exit(1)
	}
	_ = json.NewEncoder(os.Stdout).Encode(res)
}
