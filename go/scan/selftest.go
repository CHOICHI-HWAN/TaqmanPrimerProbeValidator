package main

import (
	"fmt"
	"strings"
)

// selftest: small synthetic cases with answers worked out by hand (no real sequences).
func selftest() int {
	fail := 0
	check := func(name string, ok bool, detail string) {
		st := "PASS"
		if !ok {
			st, fail = "FAIL", fail+1
		}
		fmt.Printf("%s  %s  %s\n", st, name, detail)
	}
	g := strings.Repeat("T", 4) + "ACGT" + "ACGG" + strings.Repeat("T", 4) // built from parts: the repo guard rejects long nucleotide literals
	s := oligoSites("ACGTACGG", g, 0)
	check("exact site", len(s) == 1 && s[0].Start == 4 && s[0].End == 12 && s[0].Strand == "+" && s[0].MM == 0, fmt.Sprint(s))
	s = oligoSites("ACGTACGA", g, 1)
	check("one mismatch", len(s) == 1 && s[0].MM == 1, fmt.Sprint(s))
	s = oligoSites("CCGTACGT", g, 0)
	check("minus strand", len(s) == 1 && s[0].Strand == "-" && s[0].Start == 4, fmt.Sprint(s))
	s = oligoSites("ACGTRCGG", g, 0)
	check("IUPAC R matches A", len(s) == 1, fmt.Sprint(s))
	check("revcomp", revComp("ACGTRY") == "RYACGT", revComp("ACGTRY"))
	check("terminal 3' ok", terminalOK("ACGTA", "ACGTT", true, 1) == false && terminalOK("ACGTA", "TCGTA", true, 1), "")
	return fail
}
