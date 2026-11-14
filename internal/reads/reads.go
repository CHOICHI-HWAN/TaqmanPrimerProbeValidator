// Package reads checks sequencing reads (FASTQ or FASTA, single-end, paired-end or long reads) for the oligos of
// primer sets: which reads carry each primer / probe, with how many mismatches, and which bases the reads have at
// every oligo position.
//
// Search: every oligo (and its reverse complement) is cut into disjoint seeds of length s; with at most m
// mismatches (or m edits) at least one of m+1 disjoint pieces is untouched, so it appears exactly in the read. Each
// seed hit is verified against the whole oligo: Hamming distance by default, or a semi-global edit distance
// (indels allowed, for error-prone long reads) with Indel. IUPAC codes and inosine (I) in the oligo match every base
// they stand for; a read N is a mismatch. An oligo must lie wholly inside the read.
package reads

import (
	"fmt"
	"strings"

	"chichipes/internal/primers"
)

// mask of the read bases each oligo letter accepts (A=1 C=2 G=4 T=8)
var mask = func() (t [256]uint8) {
	m := map[byte]string{'A': "A", 'C': "C", 'G': "G", 'T': "T", 'U': "T", 'R': "AG", 'Y': "CT", 'S': "CG", 'W': "AT",
		'K': "GT", 'M': "AC", 'B': "CGT", 'D': "AGT", 'H': "ACT", 'V': "ACG", 'N': "ACGT", 'I': "ACGT"}
	bit := map[byte]uint8{'A': 1, 'C': 2, 'G': 4, 'T': 8}
	for k, v := range m {
		for i := 0; i < len(v); i++ {
			t[k] |= bit[v[i]]
			t[k+32] |= bit[v[i]]
		}
	}
	return
}()

var readBit = func() (t [256]uint8) {
	for b, v := range map[byte]uint8{'A': 1, 'C': 2, 'G': 4, 'T': 8} {
		t[b], t[b+32] = v, v
	}
	return
}()

var code2 = func() (t [256]int8) {
	for i := range t {
		t[i] = -1
	}
	for i, b := range []byte("ACGT") {
		t[b], t[b+32] = int8(i), int8(i)
	}
	return
}()

var compl = func() (t [256]byte) {
	for i := range t {
		t[i] = 'N'
	}
	for _, p := range []string{"AT", "CG", "GC", "TA", "RY", "YR", "KM", "MK", "SS", "WW", "BV", "VB", "DH", "HD", "NN", "II", "UA"} {
		t[p[0]], t[p[0]+32] = p[1], p[1]
	}
	return
}()

func revComp(s []byte) []byte {
	out := make([]byte, len(s))
	for i, b := range s {
		out[len(s)-1-i] = compl[b]
	}
	return out
}

// Oligo is one primer or probe; Set and Role ("F", "R", "P" or "") come from a NAME_F / NAME_R / NAME_P name.
type Oligo struct {
	Name, Set, Role string
	Seq             []byte
}

// FromPrimers turns input oligos into Oligos, splitting SET_F / SET_R / SET_P names.
func FromPrimers(in []primers.Oligo) ([]Oligo, error) {
	var out []Oligo
	for _, o := range in {
		seq := []byte(strings.ToUpper(o.Seq))
		for _, b := range seq {
			if mask[b] == 0 {
				return nil, fmt.Errorf("%s: unsupported base %q", o.Name, b)
			}
		}
		x := Oligo{Name: o.Name, Seq: seq}
		if i := strings.LastIndexByte(o.Name, '_'); i > 0 {
			switch r := strings.ToUpper(o.Name[i+1:]); r {
			case "F", "R", "P":
				x.Set, x.Role = o.Name[:i], r
			}
		}
		out = append(out, x)
	}
	return out, nil
}

type pattern struct {
	oligo  int
	minus  bool   // the oligo's reverse complement: the oligo reads on the read's minus strand
	seq    []byte // as searched on the read
	length int
}

type seedRef struct {
	pat int32
	off int32
}

// Searcher holds the seed index of a set of oligos.
type Searcher struct {
	Oligos []Oligo
	MaxMM  int
	Indel  bool
	seed   int
	pats   []pattern
	index  map[uint32][]seedRef
}

const maxSeedExpansions = 1 << 12

// NewSearcher indexes the oligos. seed 0 picks the longest seed that still gives every oligo m+1 disjoint pieces.
func NewSearcher(oligos []Oligo, maxMM, seed int, indel bool) (*Searcher, error) {
	if len(oligos) == 0 {
		return nil, fmt.Errorf("no oligos")
	}
	if seed == 0 {
		seed = 12
		for _, o := range oligos {
			if s := len(o.Seq) / (maxMM + 1); s < seed {
				seed = s
			}
		}
	}
	if seed < 4 {
		return nil, fmt.Errorf("an oligo is too short for %d mismatches (seed would be %d bases; need >= 4)", maxMM, seed)
	}
	if seed > 15 {
		seed = 15
	}
	s := &Searcher{Oligos: oligos, MaxMM: maxMM, Indel: indel, seed: seed, index: map[uint32][]seedRef{}}
	for oi, o := range oligos {
		for _, minus := range []bool{false, true} {
			seq := o.Seq
			if minus {
				seq = revComp(o.Seq)
			}
			pi := int32(len(s.pats))
			s.pats = append(s.pats, pattern{oligo: oi, minus: minus, seq: seq, length: len(seq)})
			indexed := 0
			for off := 0; off+seed <= len(seq); off += seed {
				keys := expandSeed(seq[off : off+seed])
				for _, k := range keys {
					s.index[k] = append(s.index[k], seedRef{pi, int32(off)})
				}
				if len(keys) > 0 {
					indexed++
				}
			}
			if indexed < maxMM+1 {
				return nil, fmt.Errorf("%s: too degenerate to search with %d mismatches (only %d usable seeds)", o.Name, maxMM, indexed)
			}
		}
	}
	return s, nil
}

// expandSeed lists the 2-bit codes of every plain sequence a seed stands for (none if there are too many).
func expandSeed(seed []byte) []uint32 {
	keys := []uint32{0}
	for _, b := range seed {
		m := mask[b]
		var next []uint32
		for _, k := range keys {
			for c := uint32(0); c < 4; c++ {
				if m&(1<<c) != 0 {
					next = append(next, k<<2|c)
				}
			}
		}
		if len(next) > maxSeedExpansions {
			return nil
		}
		keys = next
	}
	return keys
}

// Hit is the best occurrence of one oligo in one read.
type Hit struct {
	Oligo    int
	Minus    bool
	Start    int    // 0-based on the read
	End      int    // exclusive
	Mismatch int    // substitutions (Hamming) or edits
	Bases    []byte // the read base facing each oligo position, 5'->3' in the oligo's orientation ('-' = deleted)
}

// Scan returns the best hit of each oligo found in seq (at most one per oligo; fewer mismatches win).
func (s *Searcher) Scan(seq []byte) []Hit {
	best := map[int]Hit{}
	tried := map[int64]bool{}
	k := s.seed
	kmask := uint32(1)<<(2*uint(k)) - 1
	var w uint32
	valid := 0
	for i := 0; i < len(seq); i++ {
		c := code2[seq[i]]
		if c < 0 {
			valid, w = 0, 0
			continue
		}
		w = (w<<2 | uint32(c)) & kmask
		if valid++; valid < k {
			continue
		}
		kstart := i + 1 - k
		for _, ref := range s.index[w] {
			p := &s.pats[ref.pat]
			cand := kstart - int(ref.off)
			key := int64(ref.pat)<<32 | int64(uint32(cand))
			if tried[key] {
				continue
			}
			tried[key] = true
			var h Hit
			var ok bool
			if s.Indel {
				h, ok = s.verifyEdit(p, seq, cand)
			} else {
				h, ok = s.verifyHamming(p, seq, cand)
			}
			if !ok {
				continue
			}
			if old, seen := best[p.oligo]; !seen || h.Mismatch < old.Mismatch {
				best[p.oligo] = h
			}
		}
	}
	out := make([]Hit, 0, len(best))
	for _, h := range best {
		out = append(out, h)
	}
	return out
}

func (s *Searcher) verifyHamming(p *pattern, seq []byte, start int) (Hit, bool) {
	if start < 0 || start+p.length > len(seq) {
		return Hit{}, false
	}
	mm := 0
	for j := 0; j < p.length; j++ {
		if mask[p.seq[j]]&readBit[seq[start+j]] == 0 {
			if mm++; mm > s.MaxMM {
				return Hit{}, false
			}
		}
	}
	h := Hit{Oligo: p.oligo, Minus: p.minus, Start: start, End: start + p.length, Mismatch: mm}
	h.Bases = orient(seq[start:start+p.length], p.minus)
	return h, true
}

// orient returns read bases in the oligo's 5'->3' orientation.
func orient(b []byte, minus bool) []byte {
	if !minus {
		return append([]byte(nil), b...)
	}
	return revComp(b)
}

// verifyEdit aligns the whole oligo against the read around start (free read ends, edits <= MaxMM).
func (s *Searcher) verifyEdit(p *pattern, seq []byte, start int) (Hit, bool) {
	m := s.MaxMM
	lo, hi := start-m, start+p.length+m
	if lo < 0 {
		lo = 0
	}
	if hi > len(seq) {
		hi = len(seq)
	}
	if hi-lo < p.length-m {
		return Hit{}, false
	}
	t := seq[lo:hi]
	L, n := p.length, len(t)
	// d[i][j]: edits aligning oligo[:i] with a read stretch ending at t[j-1]; the read start is free (d[0][j] = 0)
	d := make([][]int, L+1)
	for i := range d {
		d[i] = make([]int, n+1)
		d[i][0] = i
	}
	for i := 1; i <= L; i++ {
		ob := mask[p.seq[i-1]]
		for j := 1; j <= n; j++ {
			sub := d[i-1][j-1]
			if ob&readBit[t[j-1]] == 0 {
				sub++
			}
			v := sub
			if x := d[i-1][j] + 1; x < v { // oligo base with no read base (deletion in the read)
				v = x
			}
			if x := d[i][j-1] + 1; x < v { // extra read base (insertion in the read)
				v = x
			}
			d[i][j] = v
		}
	}
	bj := -1
	for j := 0; j <= n; j++ {
		if d[L][j] <= m && (bj < 0 || d[L][j] < d[L][bj]) {
			bj = j
		}
	}
	if bj < 0 {
		return Hit{}, false
	}
	bases := make([]byte, L)
	i, j := L, bj
	for i > 0 {
		ob := mask[p.seq[i-1]]
		switch {
		case j > 0 && d[i][j] == d[i-1][j-1]+boolInt(ob&readBit[t[j-1]] == 0):
			bases[i-1] = t[j-1]
			i, j = i-1, j-1
		case d[i][j] == d[i-1][j]+1:
			bases[i-1] = '-'
			i--
		default:
			j--
		}
	}
	h := Hit{Oligo: p.oligo, Minus: p.minus, Start: lo + j, End: lo + bj, Mismatch: d[L][bj]}
	if p.minus { // bases were collected along the reverse complement: put them in the oligo's orientation
		for a, b := 0, L-1; a < b; a, b = a+1, b-1 {
			bases[a], bases[b] = bases[b], bases[a]
		}
		for x := range bases {
			if bases[x] != '-' {
				bases[x] = compl[bases[x]]
			}
		}
	}
	h.Bases = bases
	return h, true
}

func boolInt(b bool) int {
	if b {
		return 1
	}
	return 0
}
