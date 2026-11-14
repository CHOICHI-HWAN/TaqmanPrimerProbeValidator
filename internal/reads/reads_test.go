package reads

import "testing"

func searcher(t *testing.T, seq string, mm int, indel bool) *Searcher {
	t.Helper()
	s, err := NewSearcher([]Oligo{{Name: "S_F", Set: "S", Role: "F", Seq: []byte(seq)}}, mm, 0, indel)
	if err != nil {
		t.Fatal(err)
	}
	return s
}

func TestScanHammingBothStrands(t *testing.T) {
	s := searcher(t, "ACGTTGCAGTGCATCGAGTC", 2, false)
	read := "TTTTT" + "ACGTTGCAGTGAATCGAGTC" + "GGGGG" // one mismatch
	h := s.Scan([]byte(read))
	if len(h) != 1 || h[0].Mismatch != 1 || h[0].Minus || h[0].Start != 5 {
		t.Fatalf("plus strand: %+v", h)
	}
	rc := string(revComp([]byte(read)))
	h = s.Scan([]byte(rc))
	if len(h) != 1 || h[0].Mismatch != 1 || !h[0].Minus || string(h[0].Bases) != "ACGTTGCAGTGAATCGAGTC" {
		t.Fatalf("minus strand: %+v", h)
	}
	if h := s.Scan([]byte("ACGTTGCAGTGAATCGAG")); len(h) != 0 { // oligo not wholly inside the read
		t.Fatalf("partial oligo found: %+v", h)
	}
}

func TestScanInosineAndIUPAC(t *testing.T) {
	s := searcher(t, "ACGTTGCAITGCATCGRGTC", 0, false)
	for _, x := range []string{"A", "C", "G", "T"} {
		for _, y := range []string{"A", "G"} {
			if h := s.Scan([]byte("ACGTTGCA" + x + "TGCATCG" + y + "GTC")); len(h) != 1 || h[0].Mismatch != 0 {
				t.Fatalf("I=%s R=%s: %+v", x, y, h)
			}
		}
	}
	if h := s.Scan([]byte("ACGTTGCAATGCATCGCGTC")); len(h) != 0 { // R does not take C
		t.Fatalf("R matched C: %+v", h)
	}
}

func TestScanIndel(t *testing.T) {
	s := searcher(t, "ACGTTGCAGTGCATCGAGTC", 2, true)
	read := "CCCC" + "ACGTTGCAGTGCAATCGAGTC" + "CCCC" // one inserted base
	h := s.Scan([]byte(read))
	if len(h) != 1 || h[0].Mismatch != 1 {
		t.Fatalf("insertion: %+v", h)
	}
	h = s.Scan([]byte("CCCC" + "ACGTTGCAGTCATCGAGTC" + "CCCC")) // one deleted base
	if len(h) != 1 || h[0].Mismatch != 1 {
		t.Fatalf("deletion: %+v", h)
	}
}
