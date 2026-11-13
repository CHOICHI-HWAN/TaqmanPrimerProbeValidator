package thermo

import (
	_ "embed"
	"encoding/json"
	"math"
)

// Dimer and hairpin energies as `mfeprimer dimer` / `hairpin` compute them (docs/FINDINGS.md §7, §8). The perfect
// stacks are the binding model's; mismatch stacks, dangling ends and tri/tetraloop bonuses are SantaLucia & Hicks
// 2004 values (as distributed with primer3, with the few entries where mfeprimer's table differs replaced by the
// measured ones); the hairpin loop-size term is SantaLucia's hairpin table. Ends: a terminal mismatch stack when both
// strands overhang, a dangling end (half a salt unit) when one does, nothing at a blunt end.

//go:embed hairpin_tables.json
var hairpinJSON []byte

var (
	hpMM     [4][4][4][4][2]float64 // 5'ab3'/3'cd5' -> dH, dS
	hpDangle [2][4][4][4][2]float64 // [half][pair 3' base][pair 5' base][dangling base]
	hpBonus  map[string]float64
	hpLoopL  []int
	hpLoopG  []float64
)

const hpGasR = 1.9872 // the hairpin Tm and the report oligo tables use this (the search and `thermo` use 1.987)

func init() {
	var t struct {
		MM     map[string][2]float64 `json:"mm"`
		Dangle map[string][2]float64 `json:"dangle"`
		Bonus  map[string]float64    `json:"bonus"`
		LoopL  []int                 `json:"loop_sizes"`
		LoopG  []float64             `json:"loop_g"`
	}
	if err := json.Unmarshal(hairpinJSON, &t); err != nil {
		panic(err)
	}
	for k, v := range t.MM { // "ab/cd"
		hpMM[Code(k[0])][Code(k[1])][Code(k[3])][Code(k[4])] = v
	}
	for k, v := range t.Dangle { // "half:abc"
		hpDangle[k[0]-'0'][Code(k[2])][Code(k[3])][Code(k[4])] = v
	}
	hpBonus, hpLoopL, hpLoopG = t.Bonus, t.LoopL, t.LoopG
}

// hpStack is the stack 5'-t0 t1-3' / 3'-b0 b1-5': a nearest-neighbour stack when both columns pair, else a mismatch
// stack (zero for two mismatches in a row, bar the few mfeprimer has values for).
func hpStack(t0, t1, b0, b1 byte) (h, s float64) {
	a, b, c, d := Code(t0), Code(t1), Code(b0), Code(b1)
	if c == 3-a && d == 3-b {
		return nnH[a][b], nnS[a][b]
	}
	v := hpMM[a][b][c][d]
	return v[0], v[1]
}

func hpLoopTerm(L int) float64 {
	if L <= 10 {
		for i, l := range hpLoopL {
			if l == L {
				return hpLoopG[i]
			}
		}
		return 0
	}
	for i, l := range hpLoopL { // longer loops take the next listed size up; none past the table
		if l > L {
			return hpLoopG[i]
		}
	}
	return 0
}

// HairpinEnergy returns dG37 (kcal/mol) and Tm (C) of the hairpin of s with the given pairs (outermost first,
// (i, j) with i < j, all on one diagonal i+j apart from the gaps between them).
func HairpinEnergy(s []byte, pairs [][2]int, c Cond) (dg, tm float64) {
	n := len(s)
	h, sd, units := 0.2, -5.7, 0.0
	for k := 0; k+1 < len(pairs); k++ {
		i, j, i2 := pairs[k][0], pairs[k][1], pairs[k+1][0]
		for x := 0; x < i2-i; x++ { // mismatched columns between two pairs are read as mismatch stacks
			dh, ds := hpStack(s[i+x], s[i+x+1], s[j-x], s[j-x-1])
			h, sd, units = h+dh, sd+ds, units+1
		}
	}
	i, j := pairs[0][0], pairs[0][1]
	switch {
	case i > 0 && j+1 < n: // both flanks: a terminal mismatch stack
		dh, ds := hpStack(s[i-1], s[i], s[j+1], s[j])
		h, sd, units = h+dh, sd+ds, units+1
	case j+1 < n:
		v := hpDangle[0][Code(s[j])][Code(s[i])][Code(s[j+1])]
		h, sd, units = h+v[0], sd+v[1], units+0.5
	case i > 0:
		v := hpDangle[1][Code(s[j])][Code(s[i])][Code(s[i-1])]
		h, sd, units = h+v[0], sd+v[1], units+0.5
	}
	i, j = pairs[len(pairs)-1][0], pairs[len(pairs)-1][1]
	// the closing stack is read with the closing 5' base twice (the loop base after it is never looked at)
	dh, ds := hpStack(s[i], s[i], s[j], s[j-1])
	h, sd, units = h+dh, sd+ds, units+1
	sd += 0.368 * units * math.Log(c.NaEq())
	dg = h - t37*sd/1000 + hpLoopTerm(j-i-1) + hpBonus[string(s[i:j+1])]
	tm = 1000*h/(sd+hpGasR*math.Log(c.Oligo*1e-9/4)) - 273.15
	return dg, tm
}

// OligoReport is Oligo as printed in the oligo table of the dimer and hairpin reports, which use R = 1.9872.
func OligoReport(f []byte, c Cond) (tm, dg float64) {
	h, s := perfectHS(f, c.NaEq())
	x := 4.0
	if isSelfComp(f) {
		x = 1
	}
	tm = 1000*h/(s+hpGasR*math.Log(c.Oligo*1e-9/x)) - 273.15
	return tm, h - t37*s/1000
}

// DimerDG is dG37 (kcal/mol) and Tm (C) of oligo a paired with oligo b over a's columns i..j, b written reversed
// (br = b 3'->5') and shifted so that br column c-shift faces a column c; both span ends are pairs. A duplex whose
// two strands read the same over the span and whose overhanging neighbours mirror (always so for an oligo with
// itself) gets the symmetry term and its Tm uses Ct instead of Ct/4. Same tables as the hairpin.
func DimerDG(a, br []byte, shift, i, j int, c Cond) (dg, tm float64, ok bool) {
	for x := i - 1; x <= j+1; x++ {
		if x >= 0 && x < len(a) && Code(a[x]) < 0 {
			return 0, 0, false
		}
		if y := x - shift; y >= 0 && y < len(br) && Code(br[y]) < 0 {
			return 0, 0, false
		}
	}
	h, sd, units := 0.2, -5.7, 0.0
	for x := i; x < j; x++ {
		dh, ds := hpStack(a[x], a[x+1], br[x-shift], br[x+1-shift])
		h, sd, units = h+dh, sd+ds, units+1
	}
	at := func(seq []byte, k int) byte {
		if k < 0 || k >= len(seq) {
			return 0
		}
		return seq[k]
	}
	// left end: 5'-aL a[i]-3' / 3'-bL br[i]-5'
	aL, bL := at(a, i-1), at(br, i-1-shift)
	switch {
	case aL != 0 && bL != 0:
		dh, ds := hpStack(aL, a[i], bL, br[i-shift])
		h, sd, units = h+dh, sd+ds, units+1
	case aL != 0:
		v := hpDangle[1][Code(br[i-shift])][Code(a[i])][Code(aL)]
		h, sd, units = h+v[0], sd+v[1], units+0.5
	case bL != 0:
		v := hpDangle[0][Code(br[i-shift])][Code(a[i])][Code(bL)]
		h, sd, units = h+v[0], sd+v[1], units+0.5
	}
	// right end: 5'-a[j] aR-3' / 3'-br[j] bR-5'
	aR, bR := at(a, j+1), at(br, j+1-shift)
	switch {
	case aR != 0 && bR != 0:
		dh, ds := hpStack(a[j], aR, br[j-shift], bR)
		h, sd, units = h+dh, sd+ds, units+1
	case bR != 0:
		v := hpDangle[1][Code(a[j])][Code(br[j-shift])][Code(bR)]
		h, sd, units = h+v[0], sd+v[1], units+0.5
	case aR != 0:
		v := hpDangle[0][Code(a[j])][Code(br[j-shift])][Code(aR)]
		h, sd, units = h+v[0], sd+v[1], units+0.5
	}
	x := 4.0
	sym := aL == bR && bL == aR      // the overhanging neighbours mirror each other too
	for k := i; sym && k <= j; k++ { // top a[i..j] against the bottom strand read 5'->3'
		if a[k] != br[i+j-k-shift] {
			sym = false
			break
		}
	}
	if sym {
		sd += -1.4
		x = 1
	}
	sd += 0.368 * units * math.Log(c.NaEq())
	return h - t37*sd/1000, 1000*h/(sd+hpGasR*math.Log(c.Oligo*1e-9/x)) - 273.15, true
}
