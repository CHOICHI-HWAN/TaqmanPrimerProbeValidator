package design

import (
	"encoding/json"
	"os"
	"testing"
)

func TestPySetOrder(t *testing.T) {
	// expected orders printed by CPython 3.12: list(set(vals))
	cases := []struct{ in, want []byte }{
		{[]byte{84, 71, 65, 67}, []byte{65, 67, 84, 71}},
		{[]byte{71, 84}, []byte{84, 71}},
		{[]byte{78, 45, 65, 84, 71, 67}, []byte{65, 67, 71, 45, 78, 84}},
		{[]byte{84, 84, 71, 65, 82, 89, 75, 77, 83, 87, 67, 78, 45}, []byte{65, 67, 71, 75, 77, 78, 45, 82, 83, 84, 87, 89}},
	}
	for _, c := range cases {
		if got := pySetOrder(c.in); string(got) != string(c.want) {
			t.Errorf("pySetOrder(%v) = %v, want %v", c.in, got, c.want)
		}
	}
}

func TestTmNNMatchesPython(t *testing.T) {
	if got, want := TmNN("ATGCGTACGTTAGCATGCAAT"), 55.13288431164182; got != want {
		t.Fatalf("TmNN = %.17g, want %.17g (bit-identical to probe_designer.tm_nn)", got, want)
	}
}

// TestCoverMatchesPython: designs on a synthetic alignment equal the ones scripts/redesign_coverage.py made
// (testdata/synth_expected.json).
func TestCoverMatchesPython(t *testing.T) {
	ids, rows, err := ReadAligned("testdata/synth.fa")
	if err != nil {
		t.Fatal(err)
	}
	al, err := NewAligned(ids, rows)
	if err != nil {
		t.Fatal(err)
	}
	b, err := os.ReadFile("testdata/synth_expected.json")
	if err != nil {
		t.Fatal(err)
	}
	var exp []struct {
		Params  map[string]json.RawMessage `json:"params"`
		Designs []Design                   `json:"designs"`
	}
	if err := json.Unmarshal(b, &exp); err != nil {
		t.Fatal(err)
	}
	for i, e := range exp {
		p := DefaultParams()
		p.Top = 4
		get := func(k string, v interface{}) {
			if raw, ok := e.Params[k]; ok {
				_ = json.Unmarshal(raw, v)
			}
		}
		get("ambig", &p.Ambig)
		get("probe_ambig", &p.ProbeAmbig)
		get("pass_mm", &p.PassMM)
		get("tp_n", &p.TpN)
		get("amp", &p.Amp)
		get("primer_tm", &p.PrimerTm)
		get("probe_above", &p.ProbeAbove)
		run := DesignFromAlignment(al, p)
		if len(run.Designs) != len(e.Designs) {
			t.Fatalf("case %d: %d designs, want %d", i, len(run.Designs), len(e.Designs))
		}
		for j, d := range run.Designs {
			w := e.Designs[j]
			if d.F != w.F || d.P != w.P || d.R != w.R || d.EstCov != w.EstCov || d.EstAdded != w.EstAdded || d.AmpliconLen != w.AmpliconLen ||
				d.FTm != w.FTm || d.PTm != w.PTm || d.RTm != w.RTm {
				t.Errorf("case %d design %d differs:\n got  %+v\n want %+v", i, j, d, w)
			}
		}
	}
}

func TestChem(t *testing.T) {
	c, _ := ParseChem("mgb", nan(), 0, 0)
	if c.MaxLen != 22 || c.Offset != 10 {
		t.Fatalf("mgb defaults: %+v", c)
	}
	if p := LabelProblems("FAM/BHQ2"); len(p) == 0 {
		t.Fatal("FAM/BHQ2 should be flagged")
	}
	if p := LabelProblems("FAM-ZEN-IBFQ"); len(p) != 0 {
		t.Fatalf("FAM-ZEN-IBFQ flagged: %v", p)
	}
	if p := MultiplexProblems(map[string]string{"S": "FAM/BHQ1", "L": "FAM/BHQ1"}); len(p) == 0 {
		t.Fatal("two FAM sets in one tube should be flagged")
	}
}

func nan() float64 { var z float64; return z / z }
