package design

import (
	"fmt"
	"math"
	"strings"

	"chichipes/internal/primers"
	"chichipes/internal/spec"
	"chichipes/internal/thermo"
)

// Set is a reported primer/probe set with its coverage and the chichipes check.
type Set struct {
	Name, F, P, R      string
	Change             string // shift: which oligos differ from the starting set
	BaseCov            int    // shift: coverage of the starting set
	Cov, Added         int
	Amplicon           int
	FS, RE             int // alignment columns (F start, R end)
	FTm, PTm, RTm      float64
	PTmLabel           float64 // probe Tm with the format's offset
	FHp, PHp, RHp      float64 // most stable hairpin dG (0 = none)
	Dimer3, DimerProbe float64 // most stable extendable (primer 3' end) dimer; most stable probe-involving dimer
	Penalty            float64
	Flags              []string
	OffTarget          string
	Bits               Bits
}

// SetColumns is the TSV header of reported sets.
var SetColumns = []string{"name", "change", "coverage", "coverage_pct", "base_pct", "added", "amplicon", "F", "P", "R",
	"F_tm", "P_tm", "P_tm_label", "R_tm", "hairpin_dg_F/P/R", "dimer3_dg", "probe_dimer_dg", "penalty", "off_target", "flags"}

// Row renders a set for the TSV.
func (s *Set) Row(n int) []string {
	pct := func(c int) string { return fmt.Sprintf("%.2f", 100*float64(c)/float64(max(n, 1))) }
	base := ""
	if s.Change != "" {
		base = pct(s.BaseCov)
	}
	f1 := func(x float64) string { return fmt.Sprintf("%.1f", x) }
	return []string{s.Name, s.Change, fmt.Sprint(s.Cov), pct(s.Cov), base, fmt.Sprint(s.Added), fmt.Sprint(s.Amplicon), s.F, s.P, s.R,
		f1(s.FTm), f1(s.PTm), f1(s.PTmLabel), f1(s.RTm), f1(s.FHp) + "/" + f1(s.PHp) + "/" + f1(s.RHp), f1(s.Dimer3), f1(s.DimerProbe),
		fmt.Sprintf("%.2f", s.Penalty), s.OffTarget, strings.Join(s.Flags, "; ")}
}

// tm04 is the mean chichipes Tm (SantaLucia & Hicks 2004, spec-search buffer) over the IUPAC resolutions.
func tm04(seq string, c thermo.Cond) float64 {
	res, err := primers.Expand(seq)
	if err != nil || len(res) == 0 {
		return 0
	}
	t := 0.0
	for _, r := range res {
		x, _ := thermo.Oligo([]byte(r), c)
		t += x
	}
	return t / float64(len(res))
}

func worstHairpin(seq string, c thermo.Cond) float64 {
	res, err := primers.Expand(seq)
	if err != nil {
		return 0
	}
	w := 0.0
	for _, r := range res {
		if _, _, dg, _, ok := spec.FindHairpin([]byte(r), spec.DefaultHairpinOptions, c); ok && dg < w {
			w = dg
		}
	}
	return w
}

// QC checks a set with chichipes' thermodynamics and the probe-format rules, and scores a varVAMP-like penalty
// (lower is better): primer Tm mismatch, GC away from 40-60 %, hairpins and dimers below -2 / -5 kcal/mol, IUPAC codes.
// The probe carries a 3' quencher, so it cannot be extended: only dimers at a primer's 3' end count as extendable.
func QC(s *Set, p Params) {
	c := thermo.Default
	s.Flags = nil
	s.FTm, s.PTm, s.RTm = tm04(s.F, c), tm04(s.P, c), tm04(s.R, c)
	s.PTmLabel = s.PTm + p.Probe.Offset
	s.FHp, s.PHp, s.RHp = worstHairpin(s.F, c), worstHairpin(s.P, c), worstHairpin(s.R, c)
	hi := math.Max(s.FTm, s.RTm)
	if s.PTmLabel < hi+2 {
		s.Flags = append(s.Flags, fmt.Sprintf("probe Tm %.1f < primers+2 (%.1f)", s.PTmLabel, hi+2))
	}
	if math.Abs(s.FTm-s.RTm) > 5 {
		s.Flags = append(s.Flags, fmt.Sprintf("primer Tm differ by %.1f", math.Abs(s.FTm-s.RTm)))
	}
	if len(s.P) > 0 && s.P[0] == 'G' {
		s.Flags = append(s.Flags, "probe 5' G (quenches the reporter)")
	}
	if !p.Probe.lengthOK(len(s.P)) {
		s.Flags = append(s.Flags, fmt.Sprintf("probe %d nt outside the %s range %v", len(s.P), p.Probe.Name, p.Probe.ProbeLenRange()))
	}
	if !MoreCThanG(s.P) {
		s.Flags = append(s.Flags, "probe has more G than C (the other strand is preferred)")
	}
	// dimers
	in := []primers.Oligo{{Name: "F", Seq: s.F}, {Name: "R", Seq: s.R}, {Name: "P", Seq: s.P}}
	ex, err := spec.ExpandAll(in)
	if err == nil && len(ex) <= 64 {
		for _, d := range spec.Dimers(ex, 3, 3, 0, c) {
			ra, rb := d.A.Name[:1], d.B.Name[:1]
			ext := (ra != "P" && d.J >= len(d.A.Seq)-1) || (rb != "P" && d.I <= d.Shift && d.Shift <= d.J)
			if (ra == "P" || rb == "P") && d.Dg < s.DimerProbe {
				s.DimerProbe = d.Dg
			}
			if ext && d.Dg < s.Dimer3 {
				s.Dimer3 = d.Dg
			}
		}
	}
	if s.Dimer3 < -6 {
		s.Flags = append(s.Flags, fmt.Sprintf("3' dimer %.1f kcal/mol", s.Dimer3))
	}
	pen := math.Abs(s.FTm - s.RTm)
	for _, x := range []string{s.F, s.R, s.P} {
		g := gcOf(x)
		if g < 40 {
			pen += (40 - g) / 5
		} else if g > 60 {
			pen += (g - 60) / 5
		}
		pen += 0.5 * float64(NAmbig(x))
	}
	for _, h := range []float64{s.FHp, s.RHp, s.PHp} {
		pen += math.Max(0, -h-2)
	}
	pen += 2*math.Max(0, -s.Dimer3-5) + math.Max(0, -s.DimerProbe-8)
	s.Penalty = pen
}

func gcOf(s string) float64 {
	res := Resolutions(s, 1<<12)
	if len(res) == 0 {
		return 50
	}
	t := 0.0
	for _, r := range res {
		t += gcPercent(r)
	}
	return t / float64(len(res))
}

// PreferCStrand puts the probe on the strand with more C than G (the ABI rule for hydrolysis probes) when the other
// strand also has no 5' G and no 4-base run. A perfect-match probe passes the same genomes on either strand, and the
// screening Tm is strand-symmetric, so coverage does not change.
func PreferCStrand(s *Set) {
	if MoreCThanG(s.P) || s.P == "" {
		return
	}
	rc := RevComp(s.P)
	if !MoreCThanG(rc) {
		return
	}
	res := Resolutions(rc, 64)
	if len(res) == 0 {
		return
	}
	for _, r := range res {
		if r[0] == 'G' || longestRun(r) >= 4 {
			return
		}
	}
	s.P = rc
}
