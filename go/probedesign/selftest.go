package main

import (
	"fmt"
)

// selftest: fixed-value checks against numbers the Python implementation produced (synthetic oligos only).
func selftest() int {
	fail := 0
	check := func(name string, ok bool, detail string) {
		st := "PASS"
		if !ok {
			st, fail = "FAIL", fail+1
		}
		fmt.Printf("%s  %s  %s\n", st, name, detail)
	}
	for _, c := range []struct {
		seq string
		tm  float64
	}{{"ACGTACGTACGTACGT", 48.719697993432135}, {"ATATATATATATATATAT", 23.955395224268273}} {
		t, _ := tmNN(c.seq)
		check("tm_nn "+c.seq, t == c.tm, fmt.Sprintf("go=%.6f py=%.6f", t, c.tm))
	}
	check("longest_run", longestRun("AACCCCGT") == 4, "")
	check("resolutions RY", len(resolutions("ARY", 64)) == 4, "")
	check("resolutions limit", len(resolutions("NNNNNNNNNNNNNNNNNNNNNNNN", 64)) == 0, "no 4^24 expansion")
	check("revcomp", revComp("ACGTRY") == "RYACGT", revComp("ACGTRY"))
	return fail
}
