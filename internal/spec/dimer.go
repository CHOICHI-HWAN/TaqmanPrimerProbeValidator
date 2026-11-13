package spec

import (
	"bufio"
	"fmt"
	"os"
	"sort"
	"strings"
	"time"

	"chichipes/internal/primers"
	"chichipes/internal/thermo"
)

// Dimer is one oligo-oligo duplex as `mfeprimer dimer` reports it. A is the later oligo of the (expanded) input,
// B the earlier one (or the same one); B is written reversed under A, shifted by Shift columns.
type Dimer struct {
	A, B       Expanded
	Shift      int
	I, J       int // paired span on A's columns
	Score      int // pairs - mismatches inside the span
	Dg         float64
	order      int
	marksCache string
}

// Expanded is one concrete oligo of the input (degenerate ones expanded, named like the primer search).
type Expanded struct {
	Name string
	Seq  []byte
}

// ExpandAll expands the input oligos in order.
func ExpandAll(oligos []primers.Oligo) ([]Expanded, error) {
	var out []Expanded
	for _, o := range oligos {
		seqs, err := primers.Expand(o.Seq)
		if err != nil {
			return nil, fmt.Errorf("%s: %w", o.Name, err)
		}
		for i, s := range seqs {
			name := o.Name
			if len(seqs) > 1 {
				name = fmt.Sprintf("%s.%d", o.Name, i+1)
			}
			out = append(out, Expanded{name, []byte(s)})
		}
	}
	return out, nil
}

var pairOK = func() (t [256][256]bool) {
	for _, p := range []string{"AT", "TA", "CG", "GC"} {
		t[p[0]][p[1]] = true
	}
	return
}()

// Dimers finds, for every pair of expanded oligos (and every oligo with itself) and every relative shift, the span
// from the first to the last complementary column; it is a dimer when it has at most maxMM mismatches, a score
// (pairs - mismatches) of at least minScore and dG <= maxDg. Sorted by dG.
func Dimers(ex []Expanded, minScore, maxMM int, maxDg float64, c thermo.Cond) []Dimer {
	var out []Dimer
	for ai := range ex {
		for bi := 0; bi <= ai; bi++ {
			a, b := ex[ai].Seq, ex[bi].Seq
			br := make([]byte, len(b))
			for k := range b {
				br[k] = b[len(b)-1-k]
			}
			for sh := -len(br) + 1; sh < len(a); sh++ {
				lo, hi := max(0, sh), min(len(a), sh+len(br))
				first, last, pairs := -1, -1, 0
				for col := lo; col < hi; col++ {
					if pairOK[a[col]][br[col-sh]] {
						if first < 0 {
							first = col
						}
						last = col
						pairs++
					}
				}
				if first < 0 {
					continue
				}
				mm := last - first + 1 - pairs
				if mm > maxMM || pairs-mm < minScore {
					continue
				}
				dg, _, ok := thermo.DimerDG(a, br, sh, first, last, c)
				if !ok || dg > maxDg {
					continue
				}
				out = append(out, Dimer{A: ex[ai], B: ex[bi], Shift: sh, I: first, J: last, Score: pairs - mm, Dg: dg,
					order: len(out)})
			}
		}
	}
	sort.SliceStable(out, func(x, y int) bool { return out[x].Dg < out[y].Dg })
	return out
}

// Align is the three-line picture mfeprimer prints for a dimer.
func (d Dimer) Align() string {
	a, b := d.A.Seq, d.B.Seq
	br := make([]byte, len(b))
	for k := range b {
		br[k] = b[len(b)-1-k]
	}
	alead, blead := 0, d.Shift
	if d.Shift < 0 {
		alead, blead = -d.Shift, 0
	}
	marks := []byte(strings.Repeat(" ", alead+3+d.J+1))
	for col := d.I; col <= d.J; col++ {
		if pairOK[a[col]][br[col-d.Shift]] {
			marks[alead+3+col] = ':'
		} else {
			marks[alead+3+col] = '.'
		}
	}
	return fmt.Sprintf("%s5' %s 3'\n%s\n%s3' %s 5'", strings.Repeat(" ", alead), a, marks, strings.Repeat(" ", blead), br)
}

func writeOligoTable(w *bufio.Writer, ex []Expanded, c thermo.Cond) {
	fmt.Fprintf(w, "Primer ID                      Sequence (5'-->3')                    Length     GC      Tm      Dg     \n")
	fmt.Fprintf(w, "                                                                      (bp)      (%%)    (°C)  (kcal/mol)\n\n")
	for _, e := range ex {
		tm, dg := thermo.OligoReport(e.Seq, c)
		fmt.Fprintf(w, "%-31s%-36s%8d%8.2f%8.2f%10.2f\n", e.Name, e.Seq, len(e.Seq), thermo.GC(e.Seq), tm, dg)
	}
}

// WriteDimers writes a `mfeprimer dimer`-style report.
func WriteDimers(path, input string, ex []Expanded, ds []Dimer, minScore, maxMM int, maxDg float64, c thermo.Cond) error {
	fh, err := os.Create(path)
	if err != nil {
		return err
	}
	defer fh.Close()
	w := bufio.NewWriter(fh)
	fmt.Fprintf(w, "MFEprimer-3.0 Dimer Reports (%s)\n\n", time.Now().Format("2006-01-02 15:04:05"))
	writeOligoTable(w, ex, c)
	fmt.Fprintf(w, "\n\n\nDimer List (%d)\n\n", len(ds))
	if len(ds) == 0 {
		fmt.Fprint(w, "No dimer found.\n\n\n\n")
	} else {
		fmt.Fprint(w, "\n")
	}
	for k, d := range ds {
		if k > 0 {
			fmt.Fprint(w, "\n\n\n")
		}
		fmt.Fprintf(w, "Dimer %d: %s x %s\n\n    Score: %d, Delta G = %.2f kcal/mol\n\n%s", k+1, d.A.Name, d.B.Name, d.Score, d.Dg,
			d.Align())
	}
	if len(ds) > 0 {
		fmt.Fprint(w, "\n\n\n\n\n\n")
	}
	fmt.Fprintf(w, "Parameters\n\n%30s: %s\n%30s: %d\n%30s: %.2f\n%30s: %d\n%30s: false\n", "Primer file", input,
		"Score cutoff", minScore, "Delta G cutoff", maxDg, "MisMatch", maxMM, "Primer mode")
	fmt.Fprintf(w, "%30s: %.2f\n%30s: %.2f\n%30s: %.2f\n%30s: %.2f\n\nchichipes (mfeprimer-compatible dimer report)\n",
		"Monovalent cations [mM]", c.Mono, "Divalent cations [mM]", c.Diva, "dNTPs [mM]", c.DNTP, "Annealing oligos [nM]", c.Oligo)
	return w.Flush()
}
