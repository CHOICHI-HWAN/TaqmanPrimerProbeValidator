package design

import (
	"fmt"
	"math"
	"sort"
	"strings"
)

// Chem is the hydrolysis-probe format: it limits the probe length and shifts the probe Tm the design compares
// against the primers.
//
//   - linear: 5' reporter, 3' quencher (e.g. FAM/BHQ1). Quenching depends on the dye-quencher distance, so the probe is
//     kept at <= 30 nt. A terminal dye or quencher changes duplex stability by only about a degree (Moreira et al.
//     2005), so no Tm offset.
//   - zen: double-quenched (an internal quencher ~9 nt from the reporter plus the 3' quencher, e.g. FAM/ZEN/IBFQ).
//     Quenching no longer limits the length, so probes up to 40 nt are allowed; no Tm offset unless given.
//   - mgb: 3' minor-groove binder + non-fluorescent quencher. The MGB raises the duplex Tm strongly (Kutyavin et al.
//     2000), so MGB probes are short (13-22 nt). The default +10 C is a rough figure: calibrate it against the
//     vendor's Tm of a probe you already use (-probe-tm-offset).
//
// Every format: no G at the 5' end (it quenches the reporter), no run of 4+ identical bases.
type Chem struct {
	Name           string
	Offset         float64
	MinLen, MaxLen int
	Dye, Quencher  string
}

// ParseChem builds a Chem; offset NaN and min/max 0 take the format's defaults.
func ParseChem(name string, offset float64, minLen, maxLen int) (Chem, error) {
	c := Chem{Name: strings.ToLower(name)}
	switch c.Name {
	case "", "linear":
		c.Name, c.Offset, c.MinLen, c.MaxLen = "linear", 0, 0, 30
	case "zen", "double", "double-quenched":
		c.Name, c.Offset, c.MinLen, c.MaxLen = "zen", 0, 0, 40
	case "mgb":
		c.Offset, c.MinLen, c.MaxLen = 10, 13, 22
	default:
		return c, fmt.Errorf("unknown probe format %q (linear, zen, mgb)", name)
	}
	if !math.IsNaN(offset) {
		c.Offset = offset
	}
	if minLen > 0 {
		c.MinLen = minLen
	}
	if maxLen > 0 {
		c.MaxLen = maxLen
	}
	return c, nil
}

func (c Chem) lengthOK(n int) bool {
	return (c.MinLen == 0 || n >= c.MinLen) && (c.MaxLen == 0 || n <= c.MaxLen)
}

func (c Chem) tmOffset() float64 { return c.Offset }

// ProbeLenRange is the probe length range a format suggests (design defaults use it when not set explicitly).
func (c Chem) ProbeLenRange() [2]int {
	lo, hi := 20, 30
	if c.MinLen > 0 {
		lo = c.MinLen
	}
	if c.MaxLen > 0 {
		hi = c.MaxLen
	}
	return [2]int{lo, hi}
}

// MoreCThanG is the ABI strand rule for TaqMan probes: prefer the strand with more C than G.
func MoreCThanG(s string) bool {
	c, g := strings.Count(s, "C"), strings.Count(s, "G")
	return c > g
}

// dyeQuencher: quenchers that cover each common reporter's emission (for the multiplex check).
var dyeQuencher = map[string][]string{
	"FAM":      {"BHQ1", "BHQ-1", "IBFQ", "ZEN", "NFQ", "MGB", "TAMRA", "QSY7"},
	"HEX":      {"BHQ1", "BHQ-1", "IBFQ", "ZEN", "NFQ", "MGB", "QSY7"},
	"VIC":      {"BHQ1", "BHQ-1", "NFQ", "MGB", "QSY7"},
	"JOE":      {"BHQ1", "BHQ-1", "NFQ", "MGB"},
	"TET":      {"BHQ1", "BHQ-1", "NFQ", "MGB"},
	"ROX":      {"BHQ2", "BHQ-2", "IBRQ", "TAO", "QSY7", "QSY21"},
	"TEXASRED": {"BHQ2", "BHQ-2", "IBRQ", "TAO"},
	"CY3":      {"BHQ2", "BHQ-2", "IBRQ", "TAO"},
	"CY5":      {"BHQ2", "BHQ-2", "BHQ3", "IBRQ", "TAO", "QSY21"},
	"CY5.5":    {"BHQ3", "BHQ2", "IBRQ"},
}

// normLabel splits a label into reporter and quenchers with vendor spellings unified ("6FAM / MGBNFQ" -> FAM, NFQ).
func normLabel(label string) []string {
	parts := strings.FieldsFunc(strings.ToUpper(label), func(r rune) bool { return r == '/' || r == '-' || r == ' ' || r == ',' })
	for i, x := range parts {
		switch x {
		case "6FAM", "6-FAM", "FAM6":
			x = "FAM"
		case "SUN", "VIC", "HEX", "JOE":
			if i == 0 {
				x = "HEX"
			}
		case "MGBNFQ", "NFQMGB":
			x = "NFQ"
		case "IABKFQ", "IOWABLACKFQ":
			x = "IBFQ"
		case "IABRQSP", "IOWABLACKRQ":
			x = "IBRQ"
		case "CY5", "CY-5":
			x = "CY5"
		}
		parts[i] = x
	}
	return parts
}

// Channel is the detection channel of a label's reporter ("" when unknown).
func Channel(label string) string {
	p := normLabel(label)
	if len(p) == 0 {
		return ""
	}
	return p[0]
}

// LabelProblems checks a label (e.g. "FAM/BHQ1", "FAM-ZEN-IBFQ", "VIC-MGB") for a quencher that does not cover the
// reporter; empty when fine or unknown.
func LabelProblems(label string) []string {
	parts := normLabel(label)
	if len(parts) < 2 {
		return nil
	}
	ok, known := dyeQuencher[parts[0]]
	if !known {
		return nil
	}
	var out []string
	for _, q := range parts[1:] {
		hit := false
		for _, x := range ok {
			if q == x {
				hit = true
			}
		}
		if !hit && q != "BHQ" {
			out = append(out, fmt.Sprintf("%s does not quench %s well", q, parts[0]))
		}
	}
	return out
}

// MultiplexProblems: oligo sets meant to run in one tube must use different reporters (one per detection channel).
func MultiplexProblems(labels map[string]string) []string {
	seen := map[string]string{}
	var out []string
	sets := make([]string, 0, len(labels))
	for s := range labels {
		sets = append(sets, s)
	}
	sort.Strings(sets)
	for _, set := range sets {
		l := labels[set]
		d := Channel(l)
		if d == "" {
			continue
		}
		if o, ok := seen[d]; ok {
			out = append(out, fmt.Sprintf("%s and %s share the %s channel", o, set, d))
		}
		seen[d] = set
		out = append(out, LabelProblems(l)...)
	}
	return out
}
