package design

import (
	"bufio"
	"fmt"
	"math/bits"
	"os"
	"strings"
)

// Bits is a set of genomes (bit g = genome g of the alignment).
type Bits []uint64

func newBits(n int) Bits { return make(Bits, (n+63)/64) }

func (b Bits) set(i int)      { b[i>>6] |= 1 << (i & 63) }
func (b Bits) Has(i int) bool { return b[i>>6]>>(i&63)&1 == 1 }
func (b Bits) Count() int {
	n := 0
	for _, w := range b {
		n += bits.OnesCount64(w)
	}
	return n
}
func (b Bits) Clone() Bits { return append(Bits(nil), b...) }
func (b Bits) Or(o Bits) {
	for i := range b {
		b[i] |= o[i]
	}
}
func (b Bits) Any() bool {
	for _, w := range b {
		if w != 0 {
			return true
		}
	}
	return false
}

func popcount(w uint64) int { return bits.OnesCount64(w) }

// and3Count is |a & b & c &^ not| (not may be nil).
func and3Count(a, b, c, not Bits) int {
	n := 0
	for i := range a {
		w := a[i] & b[i] & c[i]
		if not != nil {
			w &^= not[i]
		}
		n += bits.OnesCount64(w)
	}
	return n
}

// Aligned is a lineage alignment with its majority consensus and per-genome mismatch prefix sums.
type Aligned struct {
	IDs   []string
	G, L  int
	A     [][]byte
	Cons  []byte    // majority ACGT base where >= half the genomes have it, '-' otherwise
	cs    [][]int32 // cs[g][c] = mismatches of genome g against Cons in columns [0, c)
	gapcs []int32
}

// ReadAligned reads an aligned FASTA (upper-cased).
func ReadAligned(path string) ([]string, []string, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, nil, err
	}
	defer f.Close()
	var ids []string
	var rows []*strings.Builder
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 1<<20), 1<<30)
	for sc.Scan() {
		ln := strings.TrimSpace(sc.Text())
		if strings.HasPrefix(ln, ">") {
			id := strings.Fields(ln[1:] + " x")[0]
			ids = append(ids, id)
			rows = append(rows, &strings.Builder{})
		} else if ln != "" && len(rows) > 0 {
			rows[len(rows)-1].WriteString(strings.ToUpper(ln))
		}
	}
	out := make([]string, len(rows))
	for i, r := range rows {
		out[i] = r.String()
	}
	return ids, out, sc.Err()
}

// NewAligned builds the consensus and prefix sums.
func NewAligned(ids, rows []string) (*Aligned, error) {
	if len(rows) == 0 {
		return nil, fmt.Errorf("empty alignment")
	}
	al := &Aligned{IDs: ids, G: len(rows), L: len(rows[0])}
	for i, r := range rows {
		if len(r) != al.L {
			return nil, fmt.Errorf("%s: row length %d, expected %d", ids[i], len(r), al.L)
		}
		al.A = append(al.A, []byte(r))
	}
	al.Cons = make([]byte, al.L)
	for c := 0; c < al.L; c++ {
		var n [4]int
		for g := 0; g < al.G; g++ {
			switch al.A[g][c] {
			case 'A':
				n[0]++
			case 'C':
				n[1]++
			case 'G':
				n[2]++
			case 'T':
				n[3]++
			}
		}
		top := 0
		for i := 1; i < 4; i++ {
			if n[i] > n[top] {
				top = i
			}
		}
		if float64(n[top]) >= 0.5*float64(al.G) {
			al.Cons[c] = "ACGT"[top]
		} else {
			al.Cons[c] = '-'
		}
	}
	al.cs = make([][]int32, al.G)
	for g := 0; g < al.G; g++ {
		p := make([]int32, al.L+1)
		for c := 0; c < al.L; c++ {
			p[c+1] = p[c]
			if al.A[g][c] != al.Cons[c] {
				p[c+1]++
			}
		}
		al.cs[g] = p
	}
	al.gapcs = make([]int32, al.L+1)
	for c := 0; c < al.L; c++ {
		al.gapcs[c+1] = al.gapcs[c]
		if al.Cons[c] == '-' {
			al.gapcs[c+1]++
		}
	}
	return al, nil
}

// Conserved reports whether the k-window at s has no '-' consensus column.
func (al *Aligned) Conserved(s, k int) bool {
	return s >= 0 && s+k <= al.L && al.gapcs[s+k] == al.gapcs[s]
}

// Windows lists the starts whose k-window is entirely conserved columns.
func (al *Aligned) Windows(k int) []int {
	var out []int
	for s := 0; s+k <= al.L; s++ {
		if al.gapcs[s+k] == al.gapcs[s] {
			out = append(out, s)
		}
	}
	return out
}

// upperBound: genomes with <= thr mismatches against the consensus window.
func (al *Aligned) upperBound(s, k, thr int) int {
	n := 0
	for g := 0; g < al.G; g++ {
		if int(al.cs[g][s+k]-al.cs[g][s]) <= thr {
			n++
		}
	}
	return n
}

// AmpliconLen is the number of conserved columns in [s, e).
func (al *Aligned) AmpliconLen(s, e int) int {
	n := 0
	for c := s; c < e; c++ {
		if al.Cons[c] != '-' {
			n++
		}
	}
	return n
}

// pySetOrder returns the distinct values of vals in the iteration order of Python's set(vals) (CPython 3.12, small
// ints hash to themselves), so ties between bases are broken exactly like the Python engine.
func pySetOrder(vals []byte) []byte {
	type slot struct {
		used bool
		v    int
	}
	table := make([]slot, 8)
	mask, fill := 7, 0
	insertClean := func(t []slot, m, h int) {
		i, perturb := h&m, h
		for {
			if !t[i].used {
				t[i] = slot{true, h}
				return
			}
			if i+9 <= m {
				for j := 1; j <= 9; j++ {
					if !t[i+j].used {
						t[i+j] = slot{true, h}
						return
					}
				}
			}
			perturb >>= 5
			i = (i*5 + 1 + perturb) & m
		}
	}
	for _, b := range vals {
		h := int(b)
		i, perturb := h&mask, h
		found := false
	probe:
		for {
			probes := 0
			if i+9 <= mask {
				probes = 9
			}
			for j := 0; j <= probes; j++ {
				e := table[i+j]
				if !e.used {
					table[i+j] = slot{true, h}
					fill++
					break probe
				}
				if e.v == h {
					found = true
					break probe
				}
			}
			perturb >>= 5
			i = (i*5 + 1 + perturb) & mask
		}
		if found || fill*5 < mask*3 {
			continue
		}
		size := 8
		for size <= fill*4 {
			size <<= 1
		}
		nt := make([]slot, size)
		for _, e := range table {
			if e.used {
				insertClean(nt, size-1, e.v)
			}
		}
		table, mask = nt, size-1
	}
	var out []byte
	for _, e := range table {
		if e.used {
			out = append(out, byte(e.v))
		}
	}
	return out
}

// Kind of an oligo: 'F' (3' end = right end of the window), 'R' (3' end = left end), 'P' (perfect match).
type Kind = byte

// evalCodes greedily turns up to amb consensus positions of the window [s, s+k) into two-base IUPAC codes, each where
// it adds the most genomes PASSING the rule (then the most brought within reach, then the most fixed), and returns the
// codes in the order chosen. A prefix of the result is the greedy answer for a smaller amb.
func (al *Aligned) evalCodes(s, k int, kind Kind, amb, passMM, tpN int) [][2]int {
	G := al.G
	bad := make([]bool, G*k)
	nbad := make([]int, G)
	endbad := make([]int, G)
	end := make([]bool, k)
	if kind == 'F' {
		for j := k - tpN; j < k; j++ {
			end[j] = true
		}
	} else if kind == 'R' {
		for j := 0; j < tpN && j < k; j++ {
			end[j] = true
		}
	}
	for g := 0; g < G; g++ {
		row := al.A[g][s : s+k]
		for j := 0; j < k; j++ {
			if row[j] != al.Cons[s+j] {
				bad[g*k+j] = true
			}
		}
	}
	limit := passMM
	if kind == 'P' {
		limit = 0
	}
	coded := make([]bool, k)
	var codes [][2]int
	vals := make([]byte, 0, G)
	for step := 0; step < amb; step++ {
		for g := 0; g < G; g++ {
			nb, eb := 0, 0
			for j := 0; j < k; j++ {
				if bad[g*k+j] {
					nb++
					if end[j] {
						eb++
					}
				}
			}
			nbad[g], endbad[g] = nb, eb
		}
		left := amb - step - 1
		bestJ, bestB := -1, byte(0)
		var bestKey [4]int
		for j := 0; j < k; j++ {
			vals = vals[:0]
			for g := 0; g < G; g++ {
				if bad[g*k+j] {
					vals = append(vals, al.A[g][s+j])
				}
			}
			if len(vals) == 0 || coded[j] || !isACGT(al.Cons[s+j]) {
				continue
			}
			for _, b := range pySetOrder(vals) {
				if !isACGT(b) {
					continue
				}
				var okc, nearc, fixc int
				for g := 0; g < G; g++ {
					fix := bad[g*k+j] && al.A[g][s+j] == b
					nb := nbad[g]
					fi := 0
					if fix {
						nb--
						fixc++
						if end[j] {
							fi = 1
						}
					}
					if nb <= limit && (kind == 'P' || endbad[g]-fi == 0) {
						okc++
					}
					if fix && nb <= limit+left {
						nearc++
					}
				}
				key := [4]int{okc, nearc, fixc, -j}
				if bestJ < 0 || greater(key, bestKey) {
					bestJ, bestB, bestKey = j, b, key
				}
			}
		}
		if bestJ < 0 || bestKey[0]+bestKey[1] == 0 {
			break
		}
		coded[bestJ] = true
		codes = append(codes, [2]int{bestJ, int(bestB)})
		for g := 0; g < G; g++ {
			if al.A[g][s+bestJ] == bestB {
				bad[g*k+bestJ] = false
			}
		}
	}
	return codes
}

func greater(a, b [4]int) bool {
	for i := range a {
		if a[i] != b[i] {
			return a[i] > b[i]
		}
	}
	return false
}

// oligoWith is the consensus window with the given codes applied (in genome sense).
func (al *Aligned) oligoWith(s, k int, codes [][2]int) string {
	w := []byte(string(al.Cons[s : s+k]))
	for _, c := range codes {
		w[c[0]] = codeOf[bit(w[c[0]])|bit(byte(c[1]))]
	}
	return string(w)
}

// PassBits: genomes the genome-sense oligo passes at window s under the rule (P: perfect; F/R: <= passMM mismatches,
// none in the tpN bases at the 3' end).
func (al *Aligned) PassBits(oligo string, s int, kind Kind, passMM, tpN int) Bits {
	k := len(oligo)
	var tab [][256]bool = make([][256]bool, k)
	for j := 0; j < k; j++ {
		b := iupac[oligo[j]]
		if b == "" {
			b = "ACGT"
		}
		for i := 0; i < len(b); i++ {
			tab[j][b[i]] = true
		}
	}
	out := newBits(al.G)
	for g := 0; g < al.G; g++ {
		row := al.A[g][s : s+k]
		nb, endBad := 0, false
		for j := 0; j < k; j++ {
			if !tab[j][row[j]] {
				nb++
				if (kind == 'F' && j >= k-tpN) || (kind == 'R' && j < tpN) {
					endBad = true
				}
			}
		}
		if kind == 'P' {
			if nb == 0 {
				out.set(g)
			}
		} else if nb <= passMM && !endBad {
			out.set(g)
		}
	}
	return out
}

// Locate finds the window of the alignment where most genomes match the genome-sense oligo with <= 3 mismatches
// (IUPAC-aware); ok is false when no genome does.
func (al *Aligned) Locate(oligo string) (start, genomes int, ok bool) {
	k := len(oligo)
	tab := make([][256]bool, k)
	for j := 0; j < k; j++ {
		b := iupac[oligo[j]]
		if b == "" {
			b = "ACGT"
		}
		for i := 0; i < len(b); i++ {
			tab[j][b[i]] = true
		}
	}
	best, bestN := 0, -1
	for s := 0; s+k <= al.L; s++ {
		n := 0
		for g := 0; g < al.G; g++ {
			row, mm := al.A[g][s:s+k], 0
			for j := 0; j < k && mm <= 3; j++ {
				if !tab[j][row[j]] {
					mm++
				}
			}
			if mm <= 3 {
				n++
			}
		}
		if n > bestN {
			best, bestN = s, n
		}
	}
	return best, bestN, bestN > 0
}
