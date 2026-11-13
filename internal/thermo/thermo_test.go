package thermo

import (
	"bufio"
	"math"
	"os"
	"strconv"
	"strings"
	"testing"
)

// The reference files live in work/ (git-ignored); tests skip when they are absent.

func readTSV(t *testing.T, path string) [][]string {
	fh, err := os.Open(path)
	if err != nil {
		t.Skipf("%s not available: %v", path, err)
	}
	defer fh.Close()
	var rows [][]string
	sc := bufio.NewScanner(fh)
	sc.Buffer(make([]byte, 1<<20), 1<<20)
	for sc.Scan() {
		if l := sc.Text(); l != "" && !strings.HasPrefix(l, "#") {
			rows = append(rows, strings.Split(l, "\t"))
		}
	}
	return rows
}

func TestOligoMatchesMfeprimerThermo(t *testing.T) {
	conds := map[string]Cond{
		"rand.thermo.tsv":                            Default,
		"r_mono50diva0dntp0.thermo.tsv":              {50, 0, 0, 50},
		"r_mono100diva0dntp0.thermo.tsv":             {100, 0, 0, 50},
		"r_mono50diva3dntp0.thermo.tsv":              {50, 3, 0, 50},
		"r_mono50diva1.5dntp0.25oligo200.thermo.tsv": {50, 1.5, 0.25, 200},
	}
	for file, c := range conds {
		rows := readTSV(t, "../../work/"+file)
		var etm, edg float64
		for _, r := range rows {
			tm, dg := Oligo([]byte(r[1]), c)
			wtm, _ := strconv.ParseFloat(r[4], 64)
			wdg, _ := strconv.ParseFloat(r[5], 64)
			etm = math.Max(etm, math.Abs(tm-wtm))
			edg = math.Max(edg, math.Abs(dg-wdg))
		}
		if etm > 0.011 || edg > 0.011 {
			t.Errorf("%s: max |dTm| %.4f |ddG| %.4f", file, etm, edg)
		}
		t.Logf("%s: %d oligos, max |dTm| %.4f |ddG| %.4f", file, len(rows), etm, edg)
	}
}

// joined.tsv columns from harness/synth.py
func TestBindingMatchesSpec(t *testing.T) {
	for _, exp := range []string{"exp_holdout", "exp_R_rand", "exp_tandem", "exp_R_tandem", "exp_t5"} {
		rows := readTSV(t, "../../work/"+exp+"/joined.tsv")
		head := rows[0]
		col := map[string]int{}
		for i, h := range head {
			col[h] = i
		}
		var etm, edg float64
		n := 0
		for _, r := range rows[1:] {
			if r[col["found"]] != "1" {
				continue
			}
			f, s := []byte(r[col["fpSeq"]]), []byte(r[col["site"]])
			rev := r[col["side"]] == "R"
			tm, dg, ok := Binding(f, s, r[col["l5"]][0], r[col["r3"]][0], rev, Default)
			if !ok {
				t.Fatalf("%s: binding not computable for %v", exp, r)
			}
			wtm, _ := strconv.ParseFloat(r[col["fpTm"]], 64)
			wdg, _ := strconv.ParseFloat(r[col["fpDg"]], 64)
			etm = math.Max(etm, math.Abs(tm-wtm))
			edg = math.Max(edg, math.Abs(dg-wdg))
			n++
		}
		if etm > 0.03 || edg > 0.015 {
			t.Errorf("%s: max |dTm| %.4f |ddG| %.4f", exp, etm, edg)
		}
		t.Logf("%s: %d bindings, max |dTm| %.4f |ddG| %.4f", exp, n, etm, edg)
	}
}
