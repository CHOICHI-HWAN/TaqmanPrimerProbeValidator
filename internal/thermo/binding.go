package thermo

import "math"

// Primer search and probe search energies (docs/FINDINGS.md §3, §6): the SantaLucia tables of dimer_hairpin.go, no
// fitted values.

// tStack is a stack 5'-t0 t1-3' / 3'-b0 b1-5' by codes (4 = N): N anywhere gives nothing (the stack still counts
// for salt), a perfect stack the NN values, otherwise the mismatch-stack table.
func tStack(t0, t1, b0, b1 int) (h, s float64) {
	if t0 > 3 || t1 > 3 || b0 > 3 || b1 > 3 {
		return 0, 0
	}
	if b0 == 3-t0 && b1 == 3-t1 {
		return nnH[t0][t1], nnS[t0][t1]
	}
	v := hpMM[t0][t1][b0][b1]
	return v[0], v[1]
}

// bindCut is the number of 5' columns mfeprimer drops: a mismatch run (N counts) touching the 5' end, or the 3'-most
// run of three or more. The duplex then starts at the last column of that run.
func bindCut(f []byte, s []int) int {
	n, d := len(f), 0
	for i := 0; i < n; {
		if Code(f[i]) == s[i] {
			i++
			continue
		}
		j := i
		for j < n && Code(f[j]) != s[j] {
			j++
		}
		if j-i >= 3 || i == 0 {
			d = j
		}
		i = j
	}
	return d
}

// Binding returns Tm and dG of primer f bound to site s (same length, both 5'->3' in the primer's orientation, i.e. a
// perfect site has s == f). l5 is the template-side base beyond the primer's 5' end and r3 the one beyond its 3' end
// (both in the primer's orientation, '-' at the genome end). rev is true when the primer binds the plus strand.
// Target bases may be N or any ambiguity code (scored as N); f must be A/C/G/T.
//
// Measured rules: the dropped 5' columns (bindCut) leave a blunt end with no terminal term; otherwise each end gets
// the dangling template base (half a salt unit; an N flank gives only that), or a terminal-AT term at the genome end.
// A forward primer reads its 5' dangling base uncomplemented and takes its own 3' base as the 3' dangling base; a
// reverse primer is read on the plus strand, which turns every mismatch into the complementary one.
func Binding(f, s []byte, l5, r3 byte, rev bool, c Cond) (tm, dg float64, ok bool) {
	n := len(f)
	if n != len(s) || n < 2 {
		return 0, 0, false
	}
	sc := make([]int, n)
	for i := 0; i < n; i++ {
		if Code(f[i]) < 0 || siteCode(s[i]) < 0 {
			return 0, 0, false
		}
		sc[i] = siteCode(s[i])
	}
	lc, rc := flankCode(l5), flankCode(r3)
	if lc < 0 || rc < 0 {
		return 0, 0, false
	}
	d := bindCut(f, sc)
	if d > 0 {
		f, sc = f[d-1:], sc[d-1:]
		n = len(f)
	}
	h, sd, units := 0.2, -5.7, 0.0
	for x := 0; x+1 < n; x++ {
		var dh, ds float64
		if rev { // rc(primer) over the plus-strand site
			y := n - 2 - x
			dh, ds = tStack(3-Code(f[y+1]), 3-Code(f[y]), sc[y+1], sc[y])
		} else {
			dh, ds = tStack(Code(f[x]), Code(f[x+1]), compCode(sc[x]), compCode(sc[x+1]))
		}
		h, sd, units = h+dh, sd+ds, units+1
	}
	f0, fl := Code(f[0]), Code(f[n-1])
	switch {
	case d > 0: // blunt cut end
	case lc <= 3:
		var v [2]float64
		if rev {
			v = hpDangle[1][3-f0][f0][3-lc]
		} else {
			v = hpDangle[0][3-f0][f0][lc]
		}
		h, sd, units = h+v[0], sd+v[1], units+0.5
	case lc == 4:
		units += 0.5
	case f0 == 0 || f0 == 3:
		h, sd = h+2.2, sd+6.9
	}
	switch {
	case !rev:
		v := hpDangle[1][fl][3-fl][fl]
		h, sd, units = h+v[0], sd+v[1], units+0.5
	case rc <= 3:
		v := hpDangle[0][fl][3-fl][3-rc]
		h, sd, units = h+v[0], sd+v[1], units+0.5
	case rc == 4:
		units += 0.5
	case fl == 0 || fl == 3: // not seen in the data; by symmetry with the 5' end
		h, sd = h+2.2, sd+6.9
	}
	sd += 0.368 * units * math.Log(c.NaEq())
	return 1000*h/(sd+hpGasR*math.Log(c.Oligo*1e-9/4)) - 273.15, h - t37*sd/1000, true
}

// Search returns Tm and dG the way `mfeprimer search` scores a hit of the same length as the query: the perfect
// duplex of q with every stack that holds a mismatch replaced by its table value (read as written), no end or flank
// terms and no truncation of mismatch runs. hit may hold N / ambiguity codes (scored as N).
func Search(q, hit []byte, c Cond) (tm, dg float64, ok bool) {
	n := len(q)
	if n != len(hit) || n < 2 {
		return 0, 0, false
	}
	for i := 0; i < n; i++ {
		if Code(q[i]) < 0 || siteCode(hit[i]) < 0 {
			return 0, 0, false
		}
	}
	h, sd := perfectHS(q, c.NaEq())
	for i := 0; i+1 < n; i++ {
		f1, f2, s1, s2 := Code(q[i]), Code(q[i+1]), siteCode(hit[i]), siteCode(hit[i+1])
		if f1 == s1 && f2 == s2 {
			continue
		}
		mh, ms := tStack(f1, f2, compCode(s1), compCode(s2))
		h += mh - nnH[f1][f2]
		sd += ms - nnS[f1][f2]
	}
	x := 4.0
	if isSelfComp(q) {
		x = 1
	}
	return 1000*h/(sd+gasR*math.Log(c.Oligo*1e-9/x)) - 273.15, h - t37*sd/1000, true
}
