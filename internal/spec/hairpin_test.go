package spec

import (
	"encoding/json"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"chichipes/internal/thermo"
)

// Default-option `mfeprimer hairpin -j` runs on synth_hairpin.py inputs (work/hpx/def_*/, git-ignored): every oligo's
// reported hairpin (or none) must be reproduced, bar the known equal-score tie-breaks (docs/FINDINGS.md §8).
func TestHairpinsMatchMfeprimer(t *testing.T) {
	dirs, _ := filepath.Glob("../../work/hpx/def_*")
	if len(dirs) == 0 {
		t.Skip("work/hpx/def_* not available")
	}
	total, bad := 0, 0
	for _, d := range dirs {
		fa, err := os.ReadFile(filepath.Join(d, "in.fa"))
		if err != nil {
			continue
		}
		var js struct {
			HairpinList []struct {
				Seq   struct{ ID string }
				Aseq  string
				Score int
				Dg    float64
			}
		}
		b, _ := os.ReadFile(filepath.Join(d, "out.json"))
		if err := json.Unmarshal(b, &js); err != nil {
			t.Fatal(err)
		}
		want := map[string]string{}
		for _, h := range js.HairpinList {
			want[h.Seq.ID] = h.Aseq
		}
		lines := strings.Split(strings.TrimSpace(string(fa)), "\n")
		for k := 0; k+1 < len(lines); k += 2 {
			id, seq := strings.TrimPrefix(lines[k], ">"), []byte(lines[k+1])
			got := ""
			if pairs, score, dg, tm, ok := FindHairpin(seq, DefaultHairpinOptions, thermo.Default); ok {
				got = Hairpin{Expanded{id, seq}, pairs, score, dg, tm}.Structure()
			}
			total++
			if got != want[id] {
				bad++
				t.Logf("%s %s: %q, mfeprimer %q", d, id, got, want[id])
			}
		}
	}
	t.Logf("%d oligos, %d differ", total, bad)
	if bad > total/2000 {
		t.Errorf("%d of %d oligos differ (more than the known tie-breaks)", bad, total)
	}
}
