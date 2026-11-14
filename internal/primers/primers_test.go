package primers

import (
	"bufio"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestExpandOrder(t *testing.T) {
	got, err := Expand("ACRTGYAACGTNA")
	if err != nil {
		t.Fatal(err)
	}
	want := []string{"ACGTGTAACGTGA", "ACGTGTAACGTAA", "ACGTGTAACGTTA", "ACGTGTAACGTCA", "ACGTGCAACGTGA"}
	for i, w := range want {
		if got[i] != w {
			t.Fatalf("expansion %d: got %s want %s", i+1, got[i], w)
		}
	}
	if len(got) != 16 {
		t.Fatalf("got %d expansions, want 16", len(got))
	}
}

// The primer table at the top of a real mfeprimer report (every concrete sequence and variant, by name) must be
// reproduced exactly from the FR.fa it was run on. PE_REPORT points at a <set>_mm2.txt next to <set>_FR.fa;
// nothing confidential is copied into the repository.
func TestVariantsMatchReport(t *testing.T) {
	path := os.Getenv("PE_REPORT")
	if path == "" {
		path = filepath.Join(os.Getenv("EDITIONS"), "ebov/output/primers/P001/primer_qc/P001_mm2.txt")
	}
	base := filepath.Base(path)
	fr := filepath.Join(filepath.Dir(path), base[:strings.Index(base, "_mm")]+"_FR.fa")
	oligos, err := ReadFasta(fr)
	if err != nil {
		t.Skip("no FR.fa:", err)
	}
	var want [][2]string
	for _, o := range oligos {
		cs, err := Build(o, 2, 1, 9)
		if err != nil {
			t.Fatal(err)
		}
		for _, c := range cs {
			for _, v := range c.Variants {
				want = append(want, [2]string{c.VariantName(v), v.Seq})
			}
		}
	}
	fh, err := os.Open(path)
	if err != nil {
		t.Skip("no report:", err)
	}
	defer fh.Close()
	var got [][2]string
	sc := bufio.NewScanner(fh)
	sc.Buffer(make([]byte, 1<<20), 1<<20)
	for sc.Scan() {
		l := sc.Text()
		if strings.HasPrefix(l, "Descriptions") {
			break
		}
		f := strings.Fields(l)
		if len(f) < 2 {
			continue
		}
		for _, o := range oligos {
			if f[0] == o.Name || strings.HasPrefix(f[0], o.Name+".") {
				got = append(got, [2]string{f[0], f[1]})
				break
			}
		}
	}
	if len(got) != len(want) {
		t.Fatalf("report lists %d sequences, Build gives %d", len(got), len(want))
	}
	for i := range want {
		if got[i] != want[i] {
			t.Fatalf("row %d: report %s, Build %s", i, got[i][0], want[i][0])
		}
	}
	t.Logf("%d sequences identical by name and order", len(want))
}
