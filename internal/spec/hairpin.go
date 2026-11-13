package spec

import (
	"bufio"
	"fmt"
	"os"
	"sort"
	"strings"
	"time"

	"chichipes/internal/thermo"
)

// HairpinOptions are `mfeprimer hairpin`'s settings (its defaults in DefaultHairpinOptions).
type HairpinOptions struct {
	Score      int     // -s: minimum score
	Dg         float64 // -d: maximum dG
	Tm         float64 // -t: minimum Tm
	MinDH      int     // --mindh (accepted; no effect found)
	MinL, MaxL int     // --minl, --maxl (only the default maxl is verified)
	MaxGap     int     // --maxgap
}

var DefaultHairpinOptions = HairpinOptions{Score: 5, Dg: 0, Tm: 0, MinDH: 2, MinL: 3, MaxL: 12, MaxGap: 2}

// Hairpin is the one hairpin reported for an oligo.
type Hairpin struct {
	Oligo  Expanded
	Pairs  [][2]int // outermost first
	Score  int
	Dg, Tm float64
}

type hpCand struct {
	pairs [][2]int
	score int
	at, k int // column at which it was found, diagonal
}

// hpWalk walks diagonal K (column c pairs s[c] with s[K-c]) from its outer end towards the loop while at least minl
// bases would stay unpaired. Mismatched columns are gaps; when a segment collects more than maxGap of them it is
// reported at that column with its end gaps dropped and scored by its pairs alone, and a new segment starts after
// it. The segment that reaches the loop keeps its leading gaps, drops the loop-side pieces (gaps plus the helix
// after them) that do not raise its score, and scores pairs - gaps.
func hpWalk(s []byte, K int, o HairpinOptions, out []hpCand) []hpCand {
	n := len(s)
	c := 0
	if K-n+1 > 0 {
		c = K - n + 1
	}
	emit := func(start int, cols []bool, overflow bool, at int) {
		for len(cols) > 0 && !cols[len(cols)-1] {
			cols = cols[:len(cols)-1]
		}
		if overflow {
			for len(cols) > 0 && !cols[0] {
				cols = cols[1:]
				start++
			}
		} else {
			for {
				k := len(cols)
				for k > 0 && cols[k-1] {
					k--
				}
				if k == 0 {
					break
				}
				g := k
				for g > 0 && !cols[g-1] {
					g--
				}
				if (len(cols)-k)-(k-g) > 0 {
					break
				}
				cols = cols[:g]
			}
		}
		var pairs [][2]int
		gaps := 0
		for x, ok := range cols {
			if ok {
				pairs = append(pairs, [2]int{start + x, K - start - x})
			} else {
				gaps++
			}
		}
		if len(pairs) == 0 {
			return
		}
		score := len(pairs)
		if !overflow {
			score -= gaps
		}
		out = append(out, hpCand{pairs, score, at, K})
	}
	start := c
	var cols []bool
	gaps := 0
	for (K-c)-c-1 >= o.MinL {
		ok := pairOK[s[c]][s[K-c]]
		cols = append(cols, ok)
		if !ok {
			gaps++
		}
		if gaps > o.MaxGap {
			emit(start, cols[:len(cols)-1], true, c)
			start, cols, gaps = c+1, nil, 0
		}
		c++
	}
	emit(start, cols, false, c-1)
	return out
}

// FindHairpin returns the hairpin mfeprimer reports for s (ok false when none passes the cutoffs): the first
// candidate with the highest score, candidates taken in the order they are found (by column, then diagonal); the
// cutoffs are applied to that one only.
func FindHairpin(s []byte, o HairpinOptions, c thermo.Cond) (pairs [][2]int, score int, dg, tm float64, ok bool) {
	var cands []hpCand
	for K := 0; K <= 2*len(s)-2; K++ {
		cands = hpWalk(s, K, o, cands)
	}
	sort.SliceStable(cands, func(a, b int) bool {
		if cands[a].at != cands[b].at {
			return cands[a].at < cands[b].at
		}
		return cands[a].k < cands[b].k
	})
	best := -1
	for x := range cands {
		if best < 0 || cands[x].score > cands[best].score {
			best = x
		}
	}
	if best < 0 || cands[best].score < o.Score {
		return nil, 0, 0, 0, false
	}
	b := cands[best]
	dg, tm = thermo.HairpinEnergy(s, b.pairs, c)
	if dg > o.Dg || tm < o.Tm {
		return nil, 0, 0, 0, false
	}
	return b.pairs, b.score, dg, tm, true
}

// Hairpins finds the hairpin of every expanded oligo; sorted by dG (most stable first).
func Hairpins(ex []Expanded, o HairpinOptions, c thermo.Cond) []Hairpin {
	var out []Hairpin
	for _, e := range ex {
		if pairs, score, dg, tm, ok := FindHairpin(e.Seq, o, c); ok {
			out = append(out, Hairpin{e, pairs, score, dg, tm})
		}
	}
	sort.SliceStable(out, func(a, b int) bool { return out[a].Dg < out[b].Dg })
	return out
}

// Structure is the hairpin drawn as mfeprimer does: '/' for the 5' arm, '\' for the 3' arm, '-' unpaired.
func (h Hairpin) Structure() string {
	b := []byte(strings.Repeat("-", len(h.Oligo.Seq)))
	for _, p := range h.Pairs {
		b[p[0]], b[p[1]] = '/', '\\'
	}
	return string(b)
}

// WriteHairpins writes a `mfeprimer hairpin`-style report.
func WriteHairpins(path, input string, ex []Expanded, hs []Hairpin, o HairpinOptions, c thermo.Cond) error {
	fh, err := os.Create(path)
	if err != nil {
		return err
	}
	defer fh.Close()
	w := bufio.NewWriter(fh)
	fmt.Fprintf(w, "MFEprimer-3.0 Hairpin Reports (%s)\n\n", time.Now().Format("2006-01-02 15:04:05"))
	writeOligoTable(w, ex, c)
	fmt.Fprintf(w, "\n\n\nHairpin List (%d)\n\n", len(hs))
	if len(hs) == 0 {
		fmt.Fprint(w, "No hairpins found.\n\n\n")
	} else {
		fmt.Fprint(w, "\n")
		for k, h := range hs {
			if k > 0 {
				fmt.Fprint(w, "\n\n\n")
			}
			fmt.Fprintf(w, "Hairpin %d: %s\n\n    Score: %d, Tm = %.2f °C, Delta G = %.2f kcal/mol\n    \n    %s\n    %s",
				k+1, h.Oligo.Name, h.Score, h.Tm, h.Dg, h.Structure(), h.Oligo.Seq)
		}
		fmt.Fprint(w, "\n\n\n\n\n")
	}
	fmt.Fprintf(w, "Parameters\n\n%30s: %s\n%30s: %d\n%30s: %.2f\n%30s: %.2f\n%30s: %d\n%30s: %d\n%30s: %d\n%30s: %d\n",
		"Primer file", input, "Score cutoff", o.Score, "Delta G cutoff", o.Dg, "Tm cutoff (°C)", o.Tm,
		"Min double helix size", o.MinDH, "Min loop size", o.MinL, "Max loop size", o.MaxL, "Max gap count", o.MaxGap)
	fmt.Fprintf(w, "%30s: %.2f\n%30s: %.2f\n%30s: %.2f\n%30s: %.2f\n\nchichipes (mfeprimer-compatible hairpin report)\n",
		"Monovalent cations [mM]", c.Mono, "Divalent cations [mM]", c.Diva, "dNTPs [mM]", c.DNTP, "Annealing oligos [nM]", c.Oligo)
	return w.Flush()
}
