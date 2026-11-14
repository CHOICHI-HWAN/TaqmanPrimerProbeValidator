package design

import (
	"fmt"
	"io"
	"sort"
	"strings"

	"chichipes/internal/reads"
	"chichipes/internal/seqio"
)

// Exclusion holds non-target sequences (other viruses, the host genome, ...). A set is cross-reactive on a record when
// two of its primers bind it (<= mm mismatches, 3' bases intact) facing each other within the amplicon limit; the
// probe binding inside that product makes it a false signal. The search is done in overlapping chunks, keeping the
// best hit of each oligo per chunk.
type Exclusion struct {
	recs   []seqio.Record
	mm, tp int
	ampMax int
}

// LoadExclusion reads the non-target FASTA (gzip ok).
func LoadExclusion(path string, p Params) (*Exclusion, error) {
	r, err := seqio.Open(path)
	if err != nil {
		return nil, err
	}
	defer r.Close()
	x := &Exclusion{mm: p.PassMM, tp: p.TpN, ampMax: max(p.Amp[1], 400)}
	for {
		rec, err := r.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
		rec.Seq = []byte(strings.ToUpper(string(rec.Seq)))
		x.recs = append(x.recs, rec)
	}
	if len(x.recs) == 0 {
		return nil, fmt.Errorf("%s: no sequences", path)
	}
	return x, nil
}

func iupacOK(o, b byte) bool { return strings.IndexByte(iupac[o], b) >= 0 }

type xhit struct {
	role       byte
	fwd        bool // extends towards higher coordinates
	start, end int
}

// Check sets s.OffTarget, e.g. "2 products (1 with probe) on 2 records; oligo hits F3 R1 P0", and flags the set.
func (x *Exclusion) Check(s *Set) {
	olig := []reads.Oligo{{Name: "F", Role: "F", Seq: []byte(s.F)}, {Name: "R", Role: "R", Seq: []byte(s.R)}, {Name: "P", Role: "P", Seq: []byte(s.P)}}
	srch, err := reads.NewSearcher(olig, x.mm, 0, false)
	if err != nil {
		s.OffTarget = "error: " + err.Error()
		return
	}
	hitsBy := map[byte]int{}
	products, withProbe := 0, 0
	recsHit := map[string]bool{}
	chunk := 4000
	step := chunk - x.ampMax - 60
	for _, rec := range x.recs {
		seen := map[[3]int]bool{}
		var hs []xhit
		for off := 0; off < len(rec.Seq); off += step {
			end := min(len(rec.Seq), off+chunk)
			for _, h := range srch.Scan(rec.Seq[off:end]) {
				o := olig[h.Oligo]
				role := o.Role[0]
				if role != 'P' {
					ok := true
					for j := len(o.Seq) - x.tp; j < len(o.Seq); j++ {
						if !iupacOK(o.Seq[j], h.Bases[j]) {
							ok = false
						}
					}
					if !ok {
						continue
					}
				}
				k := [3]int{int(role), off + h.Start, boolI(h.Minus)}
				if seen[k] {
					continue
				}
				seen[k] = true
				hs = append(hs, xhit{role, !h.Minus, off + h.Start, off + h.End})
			}
			if end == len(rec.Seq) {
				break
			}
		}
		for _, h := range hs {
			hitsBy[h.role]++
		}
		sort.Slice(hs, func(a, b int) bool { return hs[a].start < hs[b].start })
		for _, a := range hs {
			if a.role == 'P' || !a.fwd {
				continue
			}
			for _, b := range hs {
				if b.role == 'P' || b.fwd || b.end <= a.start || b.end-a.start > x.ampMax {
					continue
				}
				products++
				recsHit[rec.ID] = true
				for _, p := range hs {
					if p.role == 'P' && p.start >= a.start && p.end <= b.end {
						withProbe++
						break
					}
				}
			}
		}
	}
	s.OffTarget = fmt.Sprintf("%d products (%d with probe) on %d records; oligo hits F%d R%d P%d", products, withProbe, len(recsHit), hitsBy['F'], hitsBy['R'], hitsBy['P'])
	if withProbe > 0 {
		s.Flags = append(s.Flags, "CROSS-REACTIVE (non-target product with probe)")
	} else if products > 0 {
		s.Flags = append(s.Flags, "non-target product (no probe)")
	}
	s.Penalty += 10*float64(withProbe) + 2*float64(products) + 0.2*float64(hitsBy['F']+hitsBy['R']+hitsBy['P'])
}

func boolI(b bool) int {
	if b {
		return 1
	}
	return 0
}

// Rejects reports whether the set would give a non-target product with the probe inside (for Params.Reject).
func (x *Exclusion) Rejects(F, P, R string) bool {
	s := Set{F: F, P: P, R: R}
	x.Check(&s)
	for _, f := range s.Flags {
		if strings.HasPrefix(f, "CROSS-REACTIVE") {
			return true
		}
	}
	return false
}
