package reads

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"sort"
	"strings"
	"sync"

	"chichipes/internal/seqio"
)

// Stats are the counts over all reads (single-end) or read pairs (paired-end: a pair counts once, as a fragment,
// and carries an oligo when either mate does).
type Stats struct {
	Units    int64
	Bases    int64
	MM       [][]int64    // [oligo][mismatches] units whose best hit has that many mismatches
	Plus     []int64      // units where the oligo reads on the read's plus strand
	Pos      [][][6]int64 // [oligo][position] A, C, G, T, other, deleted
	Sets     map[string]*SetCount
	setOrder []string
}

// SetCount counts units carrying the oligos of one primer set.
type SetCount struct{ F, R, P, FR, FRP, Any int64 }

func newStats(s *Searcher) *Stats {
	st := &Stats{MM: make([][]int64, len(s.Oligos)), Plus: make([]int64, len(s.Oligos)),
		Pos: make([][][6]int64, len(s.Oligos)), Sets: map[string]*SetCount{}}
	for i, o := range s.Oligos {
		st.MM[i] = make([]int64, s.MaxMM+1)
		st.Pos[i] = make([][6]int64, len(o.Seq))
		if o.Set != "" && st.Sets[o.Set] == nil {
			st.Sets[o.Set] = &SetCount{}
			st.setOrder = append(st.setOrder, o.Set)
		}
	}
	return st
}

func (st *Stats) add(o *Stats) {
	st.Units += o.Units
	st.Bases += o.Bases
	for i := range st.MM {
		for k := range st.MM[i] {
			st.MM[i][k] += o.MM[i][k]
		}
		st.Plus[i] += o.Plus[i]
		for p := range st.Pos[i] {
			for b := 0; b < 6; b++ {
				st.Pos[i][p][b] += o.Pos[i][p][b]
			}
		}
	}
	for k, v := range o.Sets {
		x := st.Sets[k]
		x.F, x.R, x.P, x.FR, x.FRP, x.Any = x.F+v.F, x.R+v.R, x.P+v.P, x.FR+v.FR, x.FRP+v.FRP, x.Any+v.Any
	}
}

// Found reports whether oligo i is in at least one unit.
func (st *Stats) Found(i int) int64 {
	var n int64
	for _, v := range st.MM[i] {
		n += v
	}
	return n
}

type unit struct {
	ids  []string
	seqs [][]byte
}

var baseIdx = func() (t [256]int) {
	for i := range t {
		t[i] = 4
	}
	for i, b := range []byte("ACGT") {
		t[b], t[b+32] = i, i
	}
	t['-'] = 5
	return
}()

// count adds one unit's hits (all mates together) to st, and writes hit lines when hits is not nil.
func (s *Searcher) count(st *Stats, u unit, hits *strings.Builder) {
	st.Units++
	best := map[int]Hit{}
	for m, seq := range u.seqs {
		st.Bases += int64(len(seq))
		for _, h := range s.Scan(seq) {
			if hits != nil {
				strand := '+'
				if h.Minus {
					strand = '-'
				}
				fmt.Fprintf(hits, "%s\t%d\t%s\t%c\t%d\t%d\t%d\t%s\n", u.ids[m], m+1, s.Oligos[h.Oligo].Name, strand, h.Start+1, h.End,
					h.Mismatch, h.Bases)
			}
			if old, ok := best[h.Oligo]; !ok || h.Mismatch < old.Mismatch {
				best[h.Oligo] = h
			}
		}
	}
	type roles struct{ f, r, p bool }
	sets := map[string]*roles{}
	for oi, h := range best {
		st.MM[oi][h.Mismatch]++
		if !h.Minus {
			st.Plus[oi]++
		}
		for p, b := range h.Bases {
			st.Pos[oi][p][baseIdx[b]]++
		}
		o := s.Oligos[oi]
		if o.Set == "" {
			continue
		}
		x := sets[o.Set]
		if x == nil {
			x = &roles{}
			sets[o.Set] = x
		}
		switch o.Role {
		case "F":
			x.f = true
		case "R":
			x.r = true
		case "P":
			x.p = true
		}
	}
	for name, x := range sets {
		c := st.Sets[name]
		c.Any++
		if x.f {
			c.F++
		}
		if x.r {
			c.R++
		}
		if x.p {
			c.P++
		}
		if x.f && x.r {
			c.FR++
			if x.p {
				c.FRP++
			}
		}
	}
}

// Run scans single-end files (each read is a unit) or, with r2 set, paired files read in lockstep (each pair is a
// unit). Hit lines go to hitsOut when it is not nil.
func (s *Searcher) Run(r1 []string, r2 []string, cpu int, hitsOut io.Writer) (*Stats, error) {
	if len(r2) > 0 && len(r2) != len(r1) {
		return nil, fmt.Errorf("%d R1 files but %d R2 files", len(r1), len(r2))
	}
	if cpu < 1 {
		cpu = 1
	}
	batches := make(chan []unit, cpu*2)
	results := make(chan *Stats, cpu)
	hitLines := make(chan string, cpu*2)
	var wg sync.WaitGroup
	for w := 0; w < cpu; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			st := newStats(s)
			for b := range batches {
				var hb *strings.Builder
				if hitsOut != nil {
					hb = &strings.Builder{}
				}
				for _, u := range b {
					s.count(st, u, hb)
				}
				if hb != nil && hb.Len() > 0 {
					hitLines <- hb.String()
				}
			}
			results <- st
		}()
	}
	var hwg sync.WaitGroup
	if hitsOut != nil {
		hwg.Add(1)
		go func() {
			defer hwg.Done()
			bw := bufio.NewWriterSize(hitsOut, 1<<20)
			fmt.Fprintln(bw, "#read_id\tmate\toligo\tstrand\tstart\tend\tmismatch\tread_bases")
			for l := range hitLines {
				bw.WriteString(l)
			}
			bw.Flush()
		}()
	}
	readErr := s.feed(r1, r2, batches)
	close(batches)
	wg.Wait()
	close(results)
	close(hitLines)
	hwg.Wait()
	total := newStats(s)
	for st := range results {
		total.add(st)
	}
	return total, readErr
}

func (s *Searcher) feed(r1, r2 []string, out chan<- []unit) error {
	const batch = 4096
	for fi := range r1 {
		a, err := seqio.Open(r1[fi])
		if err != nil {
			return err
		}
		var b *seqio.Reader
		if len(r2) > 0 {
			if b, err = seqio.Open(r2[fi]); err != nil {
				a.Close()
				return err
			}
		}
		buf := make([]unit, 0, batch)
		for {
			x, err := a.Next()
			if err == io.EOF {
				if b != nil {
					if _, err2 := b.Next(); err2 != io.EOF {
						err = fmt.Errorf("%s has more reads than %s", r2[fi], r1[fi])
						a.Close()
						b.Close()
						return err
					}
				}
				break
			}
			if err != nil {
				return err
			}
			u := unit{ids: []string{x.ID}, seqs: [][]byte{x.Seq}}
			if b != nil {
				y, err := b.Next()
				if err != nil {
					return fmt.Errorf("%s ends before %s (%v)", r2[fi], r1[fi], err)
				}
				if mateID(x.ID) != mateID(y.ID) {
					return fmt.Errorf("read pair out of order: %s / %s", x.ID, y.ID)
				}
				u.ids, u.seqs = append(u.ids, y.ID), append(u.seqs, y.Seq)
			}
			buf = append(buf, u)
			if len(buf) == batch {
				out <- buf
				buf = make([]unit, 0, batch)
			}
		}
		if len(buf) > 0 {
			out <- buf
		}
		a.Close()
		if b != nil {
			b.Close()
		}
	}
	return nil
}

// mateID strips a /1 or /2 suffix so the two mates of a pair compare equal.
func mateID(id string) string {
	if n := len(id); n > 2 && id[n-2] == '/' && (id[n-1] == '1' || id[n-1] == '2') {
		return id[:n-2]
	}
	return id
}

// WriteReports writes <prefix>.summary.tsv, .sets.tsv (when the oligos are named SET_F/_R/_P) and .positions.tsv,
// and returns a short human-readable report.
func WriteReports(prefix string, s *Searcher, st *Stats, paired bool) (string, error) {
	unitName := "reads"
	if paired {
		unitName = "read_pairs"
	}
	pct := func(n int64) float64 {
		if st.Units == 0 {
			return 0
		}
		return 100 * float64(n) / float64(st.Units)
	}
	var sum strings.Builder
	fmt.Fprintf(&sum, "#oligo\tset\trole\tlength\t%s\twith_oligo\tpct", unitName)
	for k := 0; k <= s.MaxMM; k++ {
		fmt.Fprintf(&sum, "\tmm%d", k)
	}
	sum.WriteString("\tplus_strand\tminus_strand\n")
	var rep strings.Builder
	mode := "substitutions"
	if s.Indel {
		mode = "edits (indels allowed)"
	}
	fmt.Fprintf(&rep, "chichipes reads: %d %s, %d bases, up to %d %s per oligo\n\n", st.Units, strings.ReplaceAll(unitName, "_", " "), st.Bases, s.MaxMM, mode)
	fmt.Fprintf(&rep, "%-32s %12s %9s\n", "oligo", "with oligo", "%")
	for i, o := range s.Oligos {
		n := st.Found(i)
		fmt.Fprintf(&sum, "%s\t%s\t%s\t%d\t%d\t%d\t%.4f", o.Name, o.Set, o.Role, len(o.Seq), st.Units, n, pct(n))
		for _, v := range st.MM[i] {
			fmt.Fprintf(&sum, "\t%d", v)
		}
		fmt.Fprintf(&sum, "\t%d\t%d\n", st.Plus[i], n-st.Plus[i])
		fmt.Fprintf(&rep, "%-32s %12d %9.4f\n", o.Name, n, pct(n))
	}
	if err := os.WriteFile(prefix+".summary.tsv", []byte(sum.String()), 0o644); err != nil {
		return "", err
	}
	if len(st.setOrder) > 0 {
		var b strings.Builder
		fmt.Fprintf(&b, "#set\t%s\tF\tR\tP\tF_and_R\tF_R_and_P\tany\tpct_F_and_R\tpct_F_R_and_P\n", unitName)
		fmt.Fprintf(&rep, "\n%-24s %10s %10s %10s %10s %10s\n", "set", "F", "R", "P", "F+R", "F+R+P")
		names := append([]string(nil), st.setOrder...)
		sort.Strings(names)
		for _, n := range names {
			c := st.Sets[n]
			fmt.Fprintf(&b, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\n", n, st.Units, c.F, c.R, c.P, c.FR, c.FRP, c.Any, pct(c.FR), pct(c.FRP))
			fmt.Fprintf(&rep, "%-24s %10d %10d %10d %10d %10d\n", n, c.F, c.R, c.P, c.FR, c.FRP)
		}
		if err := os.WriteFile(prefix+".sets.tsv", []byte(b.String()), 0o644); err != nil {
			return "", err
		}
	}
	var pb strings.Builder
	pb.WriteString("#oligo\tpos\toligo_base\tdepth\tA\tC\tG\tT\tother\tdeleted\tmatch_pct\n")
	for i, o := range s.Oligos {
		for p, c := range st.Pos[i] {
			depth := c[0] + c[1] + c[2] + c[3] + c[4] + c[5]
			var match int64
			for b := 0; b < 4; b++ {
				if mask[o.Seq[p]]&(1<<b) != 0 {
					match += c[b]
				}
			}
			mp := 0.0
			if depth > 0 {
				mp = 100 * float64(match) / float64(depth)
			}
			fmt.Fprintf(&pb, "%s\t%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\n", o.Name, p+1, o.Seq[p], depth, c[0], c[1], c[2], c[3], c[4], c[5], mp)
		}
	}
	if err := os.WriteFile(prefix+".positions.tsv", []byte(pb.String()), 0o644); err != nil {
		return "", err
	}
	if err := os.WriteFile(prefix+".txt", []byte(rep.String()), 0o644); err != nil {
		return "", err
	}
	return rep.String(), nil
}
