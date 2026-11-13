// Package thermo computes oligo and primer-binding thermodynamics.
//
// Perfect duplexes use the SantaLucia & Hicks (2004) nearest-neighbour model with the von Ahsen (2001) Mg2+/dNTP salt
// equivalent. Binding to a target with mismatches adds parameters fitted against MFEprimer 4.2.4 outputs in black-box
// experiments (docs/FINDINGS.md, harness/fit_binding.py); they are embedded from params.json.
package thermo

import (
	_ "embed"
	"math"
)

const (
	gasR = 1.987
	t37  = 310.15
)

// Cond is the reaction buffer (mfeprimer's --mono/--diva/--dntp in mM, --oligo in nM).
type Cond struct{ Mono, Diva, DNTP, Oligo float64 }

// Default matches mfeprimer's defaults.
var Default = Cond{Mono: 50, Diva: 1.5, DNTP: 0.25, Oligo: 50}

// NaEq is the monovalent-equivalent salt concentration in M.
func (c Cond) NaEq() float64 {
	return (c.Mono + 120*math.Sqrt(math.Max(c.Diva-c.DNTP, 0))) / 1000
}

// Code maps A,C,G,T (either case) to 0..3 and anything else to -1.
func Code(b byte) int {
	switch b {
	case 'A', 'a':
		return 0
	case 'C', 'c':
		return 1
	case 'G', 'g':
		return 2
	case 'T', 't':
		return 3
	}
	return -1
}

const bases = "ACGT"

// siteCode is Code for a target base, with every IUPAC ambiguity code (R, Y, ..., N) as 4: mfeprimer scores them all
// like N (docs/FINDINGS.md). -1 for anything else.
func siteCode(b byte) int {
	if c := Code(b); c >= 0 {
		return c
	}
	switch b | 0x20 {
	case 'r', 'y', 'k', 'm', 's', 'w', 'b', 'd', 'h', 'v', 'n':
		return 4
	}
	return -1
}

// flankCode is siteCode for a flanking base, with '-' (the site touches the genome end) as 5.
func flankCode(b byte) int {
	if b == '-' {
		return 5
	}
	return siteCode(b)
}

// compCode is the code of the complementary base (N stays N).
func compCode(c int) int {
	if c == 4 {
		return 4
	}
	return 3 - c
}

// nearest-neighbour dH (kcal/mol), dS (cal/K/mol), indexed [5' base][3' base]
var nnH, nnS [4][4]float64

func init() {
	table := map[string][2]float64{
		"AA": {-7.6, -21.3}, "AT": {-7.2, -20.4}, "TA": {-7.2, -21.3}, "CA": {-8.5, -22.7}, "GT": {-8.4, -22.4},
		"CT": {-7.8, -21.0}, "GA": {-8.2, -22.2}, "CG": {-10.6, -27.2}, "GC": {-9.8, -24.4}, "GG": {-8.0, -19.9},
	}
	for k, v := range table {
		a, b := Code(k[0]), Code(k[1])
		nnH[a][b], nnS[a][b] = v[0], v[1]
		ra, rb := 3-b, 3-a // the same stack read on the other strand
		nnH[ra][rb], nnS[ra][rb] = v[0], v[1]
	}
}

func isSelfComp(f []byte) bool {
	n := len(f)
	for i := 0; i < n; i++ {
		if Code(f[i]) != 3-Code(f[n-1-i]) {
			return false
		}
	}
	return true
}

// perfectHS is dH and salt-corrected dS of f against its exact complement. f must be pure A/C/G/T.
func perfectHS(f []byte, na float64) (h, s float64) {
	n := len(f)
	for i := 0; i+1 < n; i++ {
		a, b := Code(f[i]), Code(f[i+1])
		h += nnH[a][b]
		s += nnS[a][b]
	}
	h += 0.2
	s += -5.7
	for _, e := range []byte{f[0], f[n-1]} {
		if c := Code(e); c == 0 || c == 3 {
			h += 2.2
			s += 6.9
		}
	}
	if isSelfComp(f) {
		s += -1.4
	}
	s += 0.368 * float64(n-1) * math.Log(na)
	return h, s
}

// Oligo returns Tm (C) and dG37 (kcal/mol) of an oligo against its exact complement (== `mfeprimer thermo`).
func Oligo(f []byte, c Cond) (tm, dg float64) {
	h, s := perfectHS(f, c.NaEq())
	x := 4.0
	if isSelfComp(f) {
		x = 1
	}
	tm = 1000*h/(s+gasR*math.Log(c.Oligo*1e-9/x)) - 273.15
	return tm, h - t37*s/1000
}

// GC is the GC percentage of seq.
func GC(seq []byte) float64 {
	if len(seq) == 0 {
		return 0
	}
	n := 0
	for _, b := range seq {
		if b == 'G' || b == 'C' || b == 'g' || b == 'c' {
			n++
		}
	}
	return 100 * float64(n) / float64(len(seq))
}

var _ = bases
