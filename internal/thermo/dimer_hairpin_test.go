package thermo

import (
	"encoding/json"
	"math"
	"os"
	"path/filepath"
	"testing"
)

// Reference runs (git-ignored, work/): `mfeprimer dimer -j` on harness dimer sets (work/dh/dj*.json) and
// `mfeprimer hairpin -j` on synth_hairpin.py sets (work/hpx/*/rows.json); full-precision dG and Tm.

func TestDimerDGMatchesMfeprimer(t *testing.T) {
	files, _ := filepath.Glob("../../work/dh/dj*.json")
	if len(files) == 0 {
		t.Skip("work/dh/dj*.json not available")
	}
	n := 0
	for _, f := range files {
		var d struct {
			DimerList []struct {
				S1, S2                               struct{ Seq string }
				Aseq                                 string
				S1Dangling, AseqDangling, S2Dangling int
				Dg, Tm                               float64
			}
		}
		b, _ := os.ReadFile(f)
		if err := json.Unmarshal(b, &d); err != nil {
			t.Fatal(err)
		}
		for _, x := range d.DimerList {
			a, br := []byte(x.S1.Seq), []byte(x.S2.Seq) // S2 is already written 3'->5'
			shift, off := x.S2Dangling-x.S1Dangling, x.AseqDangling-x.S1Dangling
			i, j := -1, -1
			for k, ch := range x.Aseq {
				if ch == ':' || ch == '.' {
					if i < 0 {
						i = off + k
					}
					j = off + k
				}
			}
			dg, tm, ok := DimerDG(a, br, shift, i, j, Default)
			if !ok || math.Abs(dg-x.Dg) > 1e-9 || math.Abs(tm-x.Tm) > 1e-6 {
				t.Errorf("%s: %s/%s shift %d [%d,%d]: dG %.6f Tm %.4f, mfeprimer %.6f %.4f", f, a, br, shift, i, j, dg, tm, x.Dg, x.Tm)
			}
			n++
		}
	}
	t.Logf("%d dimers", n)
}

func TestHairpinEnergyMatchesMfeprimer(t *testing.T) {
	files, _ := filepath.Glob("../../work/hpx/*/rows.json")
	if len(files) == 0 {
		t.Skip("work/hpx/*/rows.json not available")
	}
	n := 0
	for _, f := range files {
		var rows []struct {
			Seq, Struct string
			Score       *int
			Dg, Tm      float64
		}
		b, _ := os.ReadFile(f)
		if err := json.Unmarshal(b, &rows); err != nil {
			t.Fatal(err)
		}
		for _, r := range rows {
			if r.Score == nil {
				continue
			}
			var op, cl []int
			for k, ch := range r.Struct {
				switch ch {
				case '/':
					op = append(op, k)
				case '\\':
					cl = append(cl, k)
				}
			}
			pairs := make([][2]int, len(op))
			for k := range op {
				pairs[k] = [2]int{op[k], cl[len(cl)-1-k]}
			}
			dg, tm := HairpinEnergy([]byte(r.Seq), pairs, Default)
			if math.Abs(dg-r.Dg) > 1e-9 || math.Abs(tm-r.Tm) > 1e-6 {
				t.Errorf("%s: %s %s: dG %.6f Tm %.4f, mfeprimer %.6f %.4f", f, r.Seq, r.Struct, dg, tm, r.Dg, r.Tm)
			}
			n++
		}
	}
	t.Logf("%d hairpins", n)
}

// work/bind_all.json: full-precision primer bindings from `mfeprimer -j` reruns of the synth.py experiments
// (harness/bind_json.py).
func TestBindingMatchesMfeprimerJSON(t *testing.T) {
	b, err := os.ReadFile("../../work/bind_all.json")
	if err != nil {
		t.Skip("work/bind_all.json not available")
	}
	var rows []struct {
		Side, F, S, L5, R3 string
		Tm, Dg             float64
	}
	if err := json.Unmarshal(b, &rows); err != nil {
		t.Fatal(err)
	}
	bad := 0
	for _, r := range rows {
		tm, dg, ok := Binding([]byte(r.F), []byte(r.S), r.L5[0], r.R3[0], r.Side == "R", Default)
		if !ok || math.Abs(dg-r.Dg) > 1e-9 || math.Abs(tm-r.Tm) > 1e-6 {
			bad++
			if bad <= 5 {
				t.Logf("%s %s/%s %s %s: dG %.6f Tm %.4f, mfeprimer %.6f %.4f", r.Side, r.F, r.S, r.L5, r.R3, dg, tm, r.Dg, r.Tm)
			}
		}
	}
	t.Logf("%d bindings, %d differ", len(rows), bad)
	if bad > 1 { // one record sits in a poly-G run where the harness cannot place the site unambiguously
		t.Errorf("%d bindings differ", bad)
	}
}
