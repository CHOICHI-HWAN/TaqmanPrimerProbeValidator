package spec

import (
	"bufio"
	"fmt"
	"os"
	"sort"
	"sync"

	"chichipes/internal/primers"
	"chichipes/internal/seqio"
	"chichipes/internal/thermo"
)

// SearchHit is one hit of `mfeprimer search`: the query's 5' and 3' k-mers both match the genome exactly, and the
// span between them is within MisMatch edits (Levenshtein) of the query.
type SearchHit struct {
	Rec      int
	Start    int // 0-based, plus strand
	End      int // exclusive
	Name     string
	Strand   byte // '+': the query reads on the plus strand
	Mismatch int
	Tm, Dg   float64 // 0 for a gapped hit (span length differs from the query)
	Query    string
	Hit      string // the span in the query's orientation
}

type query struct {
	name string
	seq  []byte
}

// Search runs `mfeprimer search` semantics over db: every oligo (degenerate ones expanded, named NAME.k) is anchored
// by its 5' and 3' k-mers on either strand.
func Search(db []seqio.Record, oligos []primers.Oligo, k, maxMM, cpu int, c thermo.Cond) ([]SearchHit, error) {
	var qs []query
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
			if len(s) >= k {
				qs = append(qs, query{name, []byte(s)})
			}
		}
	}
	per := make([][]SearchHit, len(db))
	jobs := make(chan int, len(db))
	for i := range db {
		jobs <- i
	}
	close(jobs)
	if cpu < 1 {
		cpu = 1
	}
	var wg sync.WaitGroup
	for w := 0; w < cpu; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for i := range jobs {
				per[i] = searchOne(i, db[i].Seq, qs, k, maxMM, c)
			}
		}()
	}
	wg.Wait()
	var out []SearchHit
	for _, h := range per {
		out = append(out, h...)
	}
	return out, nil
}

// kmerPositions lists the start of every exact occurrence of each wanted k-mer (A/C/G/T only).
func kmerPositions(g []byte, k int, want map[uint32]bool) map[uint32][]int {
	pos := map[uint32][]int{}
	mask := uint32(1)<<(2*uint(k)) - 1
	var w uint32
	valid := 0
	for i := 0; i < len(g); i++ {
		c := thermo.Code(g[i])
		if c < 0 {
			valid, w = 0, 0
			continue
		}
		w = (w<<2 | uint32(c)) & mask
		valid++
		if valid >= k && want[w] {
			pos[w] = append(pos[w], i+1-k)
		}
	}
	return pos
}

func searchOne(rec int, g []byte, qs []query, k, maxMM int, c thermo.Cond) []SearchHit {
	want := map[uint32]bool{}
	type anchors struct{ a, b, ra, rb []uint32 }
	an := make([]anchors, len(qs))
	for i, q := range qs {
		n := len(q.seq)
		rcq := revComp(q.seq)
		// on the plus strand the reverse complement starts with rc(3' k-mer); an inosine gives one key per base
		an[i] = anchors{encodeAll(q.seq[:k]), encodeAll(q.seq[n-k:]), encodeAll(rcq[:k]), encodeAll(rcq[n-k:])}
		for _, ks := range [][]uint32{an[i].a, an[i].b, an[i].ra, an[i].rb} {
			for _, w := range ks {
				want[w] = true
			}
		}
	}
	pos := kmerPositions(g, k, want)
	var hits []SearchHit
	up := make([]byte, len(g))
	for i, x := range g {
		up[i] = upper(x)
	}
	for i, q := range qs {
		n := len(q.seq)
		for _, strand := range []byte{'+', '-'} {
			lead, tail := positions(pos, an[i].a), positions(pos, an[i].b)
			if strand == '-' {
				lead, tail = positions(pos, an[i].ra), positions(pos, an[i].rb)
			}
			for _, s := range lead {
				for _, t := range tail {
					end := t + k
					span := end - s
					if span < n-maxMM || span > n+maxMM || span < k {
						continue
					}
					hit := up[s:end]
					if strand == '-' {
						hit = revComp(hit)
					}
					mm := levenshtein(q.seq, hit)
					if mm > maxMM {
						continue
					}
					h := SearchHit{Rec: rec, Start: s, End: end, Name: q.name, Strand: strand, Mismatch: mm,
						Query: string(q.seq), Hit: string(hit)}
					if span == n {
						if tm, dg, ok := thermo.Search(q.seq, hit, c); ok {
							h.Tm, h.Dg = tm, dg
						}
					}
					hits = append(hits, h)
				}
			}
		}
	}
	sort.SliceStable(hits, func(a, b int) bool { return hits[a].Start < hits[b].Start })
	return hits
}

// positions merges the occurrences of several anchor keys (one key unless the anchor holds inosine).
func positions(pos map[uint32][]int, keys []uint32) []int {
	if len(keys) == 1 {
		return pos[keys[0]]
	}
	var out []int
	for _, w := range keys {
		out = append(out, pos[w]...)
	}
	sort.Ints(out)
	return out
}

func levenshtein(a, b []byte) int {
	prev := make([]int, len(b)+1)
	cur := make([]int, len(b)+1)
	for j := range prev {
		prev[j] = j
	}
	for i := 1; i <= len(a); i++ {
		cur[0] = i
		for j := 1; j <= len(b); j++ {
			d := prev[j-1]
			if a[i-1] != b[j-1] && !(a[i-1] == 'I' && thermo.Code(b[j-1]) >= 0) { // inosine faces any base
				d++
			}
			if prev[j]+1 < d {
				d = prev[j] + 1
			}
			if cur[j-1]+1 < d {
				d = cur[j-1] + 1
			}
			cur[j] = d
		}
		prev, cur = cur, prev
	}
	return prev[len(b)]
}

// WriteSearch writes hits in `mfeprimer search` format.
func WriteSearch(path string, db []seqio.Record, hits []SearchHit, maxMM int) error {
	fh, err := os.Create(path)
	if err != nil {
		return err
	}
	defer fh.Close()
	w := bufio.NewWriter(fh)
	fmt.Fprintf(w, "# mfeprimer search\n# max allowed mismatch: %d\n\n", maxMM)
	fmt.Fprintln(w, "#Chrom\tStart(0-based)\tEnd\tQueryName\tSize\tStrand\tMisMatch\tDeltaG\tTm\tQuerySeq\tHitSeq")
	for _, h := range hits {
		fmt.Fprintf(w, "%s\t%d\t%d\t%s\t%d\t%c\t%d\t%.2f\t%.2f\t%s\t%s\n", db[h.Rec].ID, h.Start, h.End, h.Name,
			h.End-h.Start, h.Strand, h.Mismatch, h.Dg, h.Tm, h.Query, h.Hit)
	}
	return w.Flush()
}
