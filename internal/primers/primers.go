// Package primers reads oligos and expands them the way MFEprimer names them: degenerate bases into concrete
// sequences (NAME.1, NAME.2, ...) and each concrete sequence into its 3'-window mismatch variants (NAME.1.1, ...).
package primers

import (
	"bufio"
	"fmt"
	"os"
	"strings"
)

// expansion order of each IUPAC code, as `mfeprimer degen` writes it (docs/FINDINGS.md)
var degen = map[byte]string{
	'A': "A", 'C': "C", 'G': "G", 'T': "T",
	'R': "GA", 'Y': "TC", 'M': "AC", 'K': "GT", 'S': "GC", 'W': "TA",
	'B': "GTC", 'D': "GTA", 'H': "TAC", 'V': "GAC", 'N': "GATC",
}

// substitution order for a mismatch at one position, per original base (read off real reports)
var subOrder = map[byte]string{'A': "TCG", 'C': "TAG", 'G': "TCA", 'T': "ACG"}

// Oligo is one input sequence.
type Oligo struct {
	Name string
	Seq  string
}

// ReadFasta reads oligos from a FASTA file (name = first word of the header), or from "name seq" lines.
func ReadFasta(path string) ([]Oligo, error) {
	fh, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer fh.Close()
	var out []Oligo
	sc := bufio.NewScanner(fh)
	sc.Buffer(make([]byte, 1<<20), 1<<24)
	for sc.Scan() {
		l := strings.TrimSpace(sc.Text())
		if l == "" {
			continue
		}
		if strings.HasPrefix(l, ">") {
			name := strings.Fields(l[1:])
			if len(name) == 0 {
				return nil, fmt.Errorf("%s: empty FASTA header", path)
			}
			out = append(out, Oligo{Name: name[0]})
			continue
		}
		if len(out) > 0 && !strings.HasPrefix(l, ">") && len(strings.Fields(l)) == 1 && out[len(out)-1].Seq == "" ||
			len(out) > 0 && len(strings.Fields(l)) == 1 {
			out[len(out)-1].Seq += strings.ToUpper(l)
			continue
		}
		f := strings.Fields(l)
		if len(f) >= 2 {
			out = append(out, Oligo{Name: f[0], Seq: strings.ToUpper(f[1])})
			continue
		}
		return nil, fmt.Errorf("%s: cannot parse line %q", path, l)
	}
	return out, sc.Err()
}

// Expand returns the concrete sequences of a (possibly degenerate) oligo in MFEprimer's order.
func Expand(seq string) ([]string, error) {
	out := []string{""}
	for i := 0; i < len(seq); i++ {
		alts, ok := degen[seq[i]]
		if !ok {
			return nil, fmt.Errorf("unsupported base %q in %s", seq[i], seq)
		}
		next := make([]string, 0, len(out)*len(alts))
		for _, p := range out {
			for j := 0; j < len(alts); j++ {
				next = append(next, p+string(alts[j]))
			}
		}
		out = next
	}
	return out, nil
}

// Variant is one 3'-window variant of a concrete primer: index 0 is the primer itself (name NAME.k), index i > 0 is
// named NAME.k.i.
type Variant struct {
	Index int
	Seq   string
}

// Variants lists the concrete primer and every sequence with 1..maxMM substitutions at positions misStart..misEnd
// counted from the 3' end (1 = last base), in MFEprimer's numbering order: by number of mismatches, then by the
// 3'-nearest position, then the next position, with the later position's substitution varying fastest.
func Variants(seq string, maxMM, misStart, misEnd int) []Variant {
	n := len(seq)
	if misEnd > n {
		misEnd = n
	}
	var pos []int // string indices, 3'-nearest first
	for d := misStart; d <= misEnd; d++ {
		pos = append(pos, n-d)
	}
	out := []Variant{{0, seq}}
	idx := 1
	for k := 1; k <= maxMM; k++ {
		// position combinations first (3'-nearest first), then the substitutions of that combination
		forCombos(len(pos), k, func(combo []int) {
			cur := []byte(seq)
			var subs func(i int)
			subs = func(i int) {
				if i == len(combo) {
					out = append(out, Variant{idx, string(cur)})
					idx++
					return
				}
				p := pos[combo[i]]
				alts := subOrder[seq[p]]
				for j := 0; j < len(alts); j++ {
					cur[p] = alts[j]
					subs(i + 1)
				}
				cur[p] = seq[p]
			}
			subs(0)
		})
	}
	return out
}

// forCombos calls fn with every k-subset of 0..n-1 in lexicographic order.
func forCombos(n, k int, fn func([]int)) {
	c := make([]int, k)
	var rec func(i, start int)
	rec = func(i, start int) {
		if i == k {
			fn(c)
			return
		}
		for v := start; v <= n-(k-i); v++ {
			c[i] = v
			rec(i+1, v+1)
		}
	}
	rec(0, 0)
}

// Concrete is one concrete (non-degenerate) sequence of an oligo with the 3'-window variants searched for it.
type Concrete struct {
	Name     string    // NAME.k, or NAME when the oligo is not degenerate
	Seq      string    // the concrete sequence
	Variants []Variant // [0] is Seq itself; the rest are named Name.i
}

// VariantName is the name MFEprimer gives variant v of c.
func (c Concrete) VariantName(v Variant) string {
	if v.Index == 0 {
		return c.Name
	}
	return fmt.Sprintf("%s.%d", c.Name, v.Index)
}

// Build expands an oligo into the concrete sequences and variants MFEprimer searches. The expansions of one oligo
// share a set of sequences already generated: an expansion that is itself already generated (a variant of an earlier
// expansion) gets no variants, and the others skip variants already generated; the kept variants are numbered
// 1, 2, ... in generation order.
func Build(o Oligo, maxMM, misStart, misEnd int) ([]Concrete, error) {
	seqs, err := Expand(o.Seq)
	if err != nil {
		return nil, fmt.Errorf("%s: %w", o.Name, err)
	}
	seen := map[string]bool{}
	out := make([]Concrete, 0, len(seqs))
	for i, s := range seqs {
		name := o.Name
		if len(seqs) > 1 {
			name = fmt.Sprintf("%s.%d", o.Name, i+1)
		}
		c := Concrete{Name: name, Seq: s, Variants: []Variant{{0, s}}}
		if !seen[s] {
			seen[s] = true
			idx := 1
			for _, v := range Variants(s, maxMM, misStart, misEnd)[1:] {
				if !seen[v.Seq] {
					seen[v.Seq] = true
					c.Variants = append(c.Variants, Variant{idx, v.Seq})
					idx++
				}
			}
		}
		out = append(out, c)
	}
	return out, nil
}
