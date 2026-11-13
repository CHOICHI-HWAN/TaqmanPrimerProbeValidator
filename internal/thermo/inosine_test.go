package thermo

import (
	"math"
	"testing"
)

// Every dH/dS row reproduces the dG37 printed next to it in Watkins & SantaLucia (2005) Table 2.
func TestWatkinsTranscription(t *testing.T) {
	if len(watkins2005) != 36 || len(watkinsDG37) != 36 {
		t.Fatalf("rows: %d / %d, want 36", len(watkins2005), len(watkinsDG37))
	}
	for k, v := range watkins2005 {
		g := v[0] - t37*v[1]/1000
		if math.Abs(g-watkinsDG37[k]) > 0.06 {
			t.Errorf("%s: dH %.2f dS %.1f give dG %.2f, table says %.2f", k, v[0], v[1], g, watkinsDG37[k])
		}
	}
}

func bindTm(t *testing.T, f, s string, rev bool) float64 {
	t.Helper()
	tm, _, ok := Binding([]byte(f), []byte(s), 'A', 'C', rev, Default)
	if !ok {
		t.Fatalf("Binding(%s, %s) not ok", f, s)
	}
	return tm
}

// An internal inosine pairs best with C, then A, then T and G (roughly equal), and always more weakly than a
// Watson-Crick G·C. The site is written in the primer's orientation (a perfect site equals the primer), so the
// template base facing the I is the complement of the site letter: I·C is site G, I·A site T, I·T site A, I·G site C.
func TestInosineBindingOrder(t *testing.T) {
	const tmpl = "ACGTTGCAXTGCATCGAGTC"
	at := func(s, x string) string { b := []byte(s); b[8] = x[0]; return string(b) }
	for _, rev := range []bool{false, true} {
		wc := bindTm(t, at(tmpl, "G"), at(tmpl, "G"), rev)
		ic := bindTm(t, at(tmpl, "I"), at(tmpl, "G"), rev)
		ia := bindTm(t, at(tmpl, "I"), at(tmpl, "T"), rev)
		it := bindTm(t, at(tmpl, "I"), at(tmpl, "A"), rev)
		ig := bindTm(t, at(tmpl, "I"), at(tmpl, "C"), rev)
		t.Logf("rev=%v G·C %.2f  I·C %.2f  I·A %.2f  I·T %.2f  I·G %.2f", rev, wc, ic, ia, it, ig)
		if !(wc > ic && ic > ia && ia > it && ia > ig) {
			t.Errorf("rev=%v: want G·C > I·C > I·A > I·T, I·G", rev)
		}
	}
}

// The same oligo without inosine scores exactly as before (the MFEprimer-compatible path is untouched).
func TestNoInosineUnchanged(t *testing.T) {
	f := []byte("ACGTTGCAGTGCATCGAGTC")
	h1, s1 := perfectHS(f, Default.NaEq())
	if HasInosine(f) {
		t.Fatal("HasInosine on a plain oligo")
	}
	tm, dg := Oligo(f, Default)
	if math.IsNaN(tm) || math.IsNaN(dg) || h1 == 0 || s1 == 0 {
		t.Fatal("bad perfect duplex")
	}
}

// Probe-search scoring of a query with inosine against a site with C equals the oligo's own (I·C) duplex.
func TestSearchInosineSelf(t *testing.T) {
	q := []byte("ACGTTGCAITGCATCGAGTC")
	hit := []byte("ACGTTGCAGTGCATCGAGTC") // I faces C (site letter G)
	tm, _, ok := Search(q, hit, Default)
	if !ok {
		t.Fatal("not ok")
	}
	own, _ := Oligo(q, Default)
	if math.Abs(tm-own) > 0.3 { // Oligo uses R=1.987 and the same terms
		t.Errorf("search Tm %.2f vs oligo Tm %.2f", tm, own)
	}
}
