package thermo

import "math"

// Deoxyinosine (I) in an oligo: pairs with every base, unevenly (I·C > I·A > I·T ≈ I·G). Its stacks are the
// nearest-neighbour values of Watkins & SantaLucia (2005) Nucleic Acids Res. 33:6258, Table 2 (1 M NaCl, internal
// pairs). A stack the paper has no value for (I next to a mismatch or an N, two inosines in a row on one strand) adds
// nothing but still counts for salt, the convention already used for N; a dangling end on an inosine adds only its
// half salt unit, and an inosine end gets no terminal-AT term.

// Ino is the code of deoxyinosine (A, C, G, T are 0..3, N/ambiguity 4).
const Ino = 5

// OCode is Code for an oligo base: A/C/G/T 0..3, I Ino, anything else -1.
func OCode(b byte) int {
	if b == 'I' || b == 'i' {
		return Ino
	}
	return Code(b)
}

// HasInosine reports whether the oligo holds an I.
func HasInosine(s []byte) bool {
	for _, b := range s {
		if b == 'I' || b == 'i' {
			return true
		}
	}
	return false
}

// compO is the complement code on the oligo side: I reads as I on the other strand, N as N.
func compO(c int) int {
	if c >= 4 {
		return c
	}
	return 3 - c
}

// watkins2005 is Table 2: "top/bottom" with top 5'->3' over bottom 3'->5', and dH (kcal/mol), dS (cal/K/mol).
var watkins2005 = map[string][2]float64{
	// I·C
	"AI/TC": {-8.9, -25.5}, "TI/AC": {-5.9, -17.4}, "AC/TI": {-8.8, -25.4}, "TC/AI": {-4.9, -13.9},
	"CI/GC": {-5.4, -13.7}, "GI/CC": {-6.8, -19.1}, "CC/GI": {-8.3, -23.8}, "GC/CI": {-5.0, -12.6},
	// I·A
	"AI/TA": {-8.3, -25.0}, "TI/AA": {-3.4, -11.2}, "AA/TI": {-0.7, -2.6}, "TA/AI": {-1.3, -4.6},
	"CI/GA": {2.6, 8.9}, "GI/CA": {-7.8, -21.1}, "CA/GI": {-7.0, -20.0}, "GA/CI": {-7.6, -20.2},
	// I·T
	"AI/TT": {0.49, -0.7}, "TI/AT": {-6.5, -22.0}, "AT/TI": {-5.6, -18.7}, "TT/AI": {-0.8, -4.3},
	"CI/GT": {-1.0, -2.4}, "GI/CT": {-3.5, -10.6}, "CT/GI": {0.1, -1.0}, "GT/CI": {-4.3, -12.1},
	// I·G
	"AI/TG": {-4.9, -15.8}, "TI/AG": {-1.9, -8.5}, "AG/TI": {0.1, -1.8}, "TG/AI": {1.0, 1.0},
	"CI/GG": {7.1, 21.3}, "GI/CG": {-1.1, -3.2}, "CG/GI": {5.8, 16.9}, "GG/CI": {-7.6, -22.0},
	// I·I (oligo against oligo)
	"AI/TI": {-3.3, -11.9}, "TI/AI": {0.1, -2.3}, "CI/GI": {1.3, 3.0}, "GI/CI": {-0.5, -1.3},
}

// watkinsDG37 is the printed dG37 of each row, kept to check the transcription of dH/dS (inosine_test.go).
var watkinsDG37 = map[string]float64{
	"AI/TC": -0.96, "TI/AC": -0.46, "AC/TI": -0.89, "TC/AI": -0.59, "CI/GC": -1.14, "GI/CC": -0.86, "CC/GI": -0.88, "GC/CI": -1.07,
	"AI/TA": -0.51, "TI/AA": 0.09, "AA/TI": 0.12, "TA/AI": 0.12, "CI/GA": -0.18, "GI/CA": -1.24, "CA/GI": -0.77, "GA/CI": -1.33,
	"AI/TT": 0.71, "TI/AT": 0.36, "AT/TI": 0.22, "TT/AI": 0.54, "CI/GT": -0.26, "GI/CT": -0.19, "CT/GI": 0.41, "GT/CI": -0.54,
	"AI/TG": 0.02, "TI/AG": 0.76, "AG/TI": 0.65, "TG/AI": 0.70, "CI/GG": 0.47, "GI/CG": -0.10, "CG/GI": 0.54, "GG/CI": -0.74,
	"AI/TI": 0.40, "TI/AI": 0.81, "CI/GI": 0.36, "GI/CI": -0.09,
}

var inoStacks = map[[4]int][2]float64{}

func init() {
	for k, v := range watkins2005 {
		t0, t1, b0, b1 := OCode(k[0]), OCode(k[1]), OCode(k[3]), OCode(k[4])
		inoStacks[[4]int{t0, t1, b0, b1}] = v
		inoStacks[[4]int{b1, b0, t1, t0}] = v // the same stack read on the other strand
	}
}

// stack is the stack 5'-t0 t1-3' / 3'-b0 b1-5' by codes, inosine included.
func stack(t0, t1, b0, b1 int) (h, s float64) {
	if t0 == Ino || t1 == Ino || b0 == Ino || b1 == Ino {
		if t0 == 4 || t1 == 4 || b0 == 4 || b1 == 4 {
			return 0, 0
		}
		v := inoStacks[[4]int{t0, t1, b0, b1}]
		return v[0], v[1]
	}
	return tStack(t0, t1, b0, b1)
}

// perfectHSIno is perfectHS for an oligo holding inosine, against its most stable complement (I·C).
func perfectHSIno(f []byte, na float64) (h, s float64) {
	n := len(f)
	h, s = 0.2, -5.7
	for i := 0; i+1 < n; i++ {
		a, b := OCode(f[i]), OCode(f[i+1])
		dh, ds := stack(a, b, inoPartner(a), inoPartner(b))
		h, s = h+dh, s+ds
	}
	for _, e := range []byte{f[0], f[n-1]} {
		if c := Code(e); c == 0 || c == 3 {
			h, s = h+2.2, s+6.9
		}
	}
	s += 0.368 * float64(n-1) * math.Log(na)
	return h, s
}

func inoPartner(c int) int {
	if c == Ino {
		return 1 // C
	}
	return 3 - c
}
