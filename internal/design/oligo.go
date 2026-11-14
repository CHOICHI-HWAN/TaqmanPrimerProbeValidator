// Package design picks primer / probe sets from an alignment so that the genomes passing the pipeline's own rule are
// as many as possible (coverage-first), and searches shifted / lengthened / degenerate variants of sets exhaustively.
//
// The coverage-first part is a port of the pipeline's scripts/redesign_coverage.py and gives the same designs for the
// same alignment and parameters (tests compare them); the screening Tm there is SantaLucia (1998) at 250 nM / 50 mM
// Na+, exactly as in the Python engine. The final check of a chosen set uses chichipes' own thermodynamics
// (SantaLucia & Hicks 2004, the spec-search buffer).
package design

import "math"

// IUPAC bases in the expansion order of the pipeline (mfe_common.IUPAC).
var iupac = func() (t [256]string) {
	for k, v := range map[byte]string{'A': "A", 'C': "C", 'G': "G", 'T': "T", 'R': "AG", 'Y': "CT", 'S': "GC", 'W': "AT",
		'K': "GT", 'M': "AC", 'B': "CGT", 'D': "AGT", 'H': "ACT", 'V': "ACG", 'N': "ACGT"} {
		t[k] = v
	}
	return
}()

var complement = func() (t [256]byte) {
	for i := range t {
		t[i] = 'N'
	}
	for _, p := range []string{"AT", "TA", "CG", "GC", "NN", "RY", "YR", "SS", "WW", "KM", "MK", "DH", "HD", "BV", "VB", "II"} {
		t[p[0]] = p[1]
	}
	return
}()

func bit(b byte) int {
	switch b {
	case 'A':
		return 1
	case 'C':
		return 2
	case 'G':
		return 4
	case 'T':
		return 8
	}
	return 0
}

// codeOf[mask of two bases] is the two-base IUPAC code.
var codeOf = func() (t [16]byte) {
	for _, c := range []byte("RYSWKM") {
		m := 0
		for i := 0; i < len(iupac[c]); i++ {
			m |= bit(iupac[c][i])
		}
		t[m] = c
	}
	return
}()

// RevComp is the reverse complement (IUPAC-aware; unknown letters become N).
func RevComp(s string) string {
	b := make([]byte, len(s))
	for i := 0; i < len(s); i++ {
		b[len(s)-1-i] = complement[s[i]]
	}
	return string(b)
}

func isACGT(b byte) bool { return b == 'A' || b == 'C' || b == 'G' || b == 'T' }

// NAmbig counts non-ACGT letters.
func NAmbig(s string) int {
	n := 0
	for i := 0; i < len(s); i++ {
		if !isACGT(s[i]) {
			n++
		}
	}
	return n
}

// SantaLucia (1998) unified NN parameters, as in the Python engine's probe_designer.tm_nn.
var nn98 = map[string][2]float64{"AA": {-7.9, -22.2}, "TT": {-7.9, -22.2}, "AT": {-7.2, -20.4}, "TA": {-7.2, -21.3},
	"CA": {-8.5, -22.7}, "TG": {-8.5, -22.7}, "GT": {-8.4, -22.4}, "AC": {-8.4, -22.4}, "CT": {-7.8, -21.0},
	"AG": {-7.8, -21.0}, "GA": {-8.2, -22.2}, "TC": {-8.2, -22.2}, "CG": {-10.6, -27.2}, "GC": {-9.8, -24.4},
	"GG": {-8.0, -19.9}, "CC": {-8.0, -19.9}}

var nn98idx = func() (t [4][4][2]float64) {
	for k, v := range nn98 {
		t[bit2(k[0])][bit2(k[1])] = v
	}
	return
}()

func bit2(b byte) int {
	switch b {
	case 'A':
		return 0
	case 'C':
		return 1
	case 'G':
		return 2
	}
	return 3
}

var (
	lnNa   = math.Log(50.0 / 1000.0)
	lnConc = math.Log(250.0 * 1e-9 / 4.0)
)

// pySum adds floats like Python 3.12's sum() (Neumaier-compensated after the first term), so screening Tm values are
// bit-identical to the Python engine's and ties break the same way.
func pySum(xs []float64) float64 {
	if len(xs) == 0 {
		return 0
	}
	f, c := xs[0], 0.0
	for _, x := range xs[1:] {
		t := f + x
		if math.Abs(f) >= math.Abs(x) {
			c += (f - t) + x
		} else {
			c += (x - t) + f
		}
		f = t
	}
	if c != 0 && !math.IsInf(c, 0) && !math.IsNaN(c) {
		f += c
	}
	return f
}

// TmNN is the screening Tm of a plain ACGT oligo (250 nM oligo, 50 mM Na+).
func TmNN(s string) float64 {
	var hs, ss [64]float64
	h, sx := hs[:0], ss[:0]
	for i := 0; i+1 < len(s); i++ {
		v := nn98idx[bit2(s[i])][bit2(s[i+1])]
		h, sx = append(h, v[0]), append(sx, v[1])
	}
	dh, ds := pySum(h), pySum(sx)
	for _, e := range []byte{s[0], s[len(s)-1]} {
		if e == 'G' || e == 'C' {
			dh, ds = dh+2.3, ds+4.1
		} else {
			dh, ds = dh+0.1, ds-2.8
		}
	}
	ds += 0.368 * float64(len(s)-1) * lnNa
	return dh*1000.0/(ds+1.987*lnConc) - 273.15
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
		best = max(best, cur)
	}
	return best
}

// Resolutions lists every plain oligo an IUPAC oligo stands for, or nil when there would be more than limit.
func Resolutions(s string, limit int) []string {
	n := 1
	for i := 0; i < len(s); i++ {
		n *= max(1, len(iupac[s[i]]))
		if n > limit {
			return nil
		}
	}
	outs := []string{""}
	for i := 0; i < len(s); i++ {
		b := iupac[s[i]]
		if b == "" {
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

// stats: resolutions (<= 64), mean screening Tm and mean GC, as the Python engine's _stats.
func stats(s string) (res []string, tm, gc float64, ok bool) {
	res = Resolutions(s, 64)
	if len(res) == 0 {
		return nil, 0, 0, false
	}
	ts, gs := make([]float64, len(res)), make([]float64, len(res))
	for i, r := range res {
		ts[i], gs[i] = TmNN(r), gcPercent(r)
	}
	tm, gc = pySum(ts), pySum(gs)
	return res, tm / float64(len(res)), gc / float64(len(res)), true
}
