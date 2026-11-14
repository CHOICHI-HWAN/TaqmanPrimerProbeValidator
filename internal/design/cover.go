package design

import (
	"container/heap"
	"fmt"
	"math"
	"runtime"
	"sort"
	"sync"
)

// Params of one coverage-first design run (the Python engine's design_from_alignment arguments).
type Params struct {
	Ambig      int        `json:"ambig"`
	ProbeAmbig int        `json:"probe_ambig"`
	PassMM     int        `json:"pass_mm"`
	TpN        int        `json:"tp_n"`
	Amp        [2]int     `json:"amp"`
	PrimerTm   [2]float64 `json:"primer_tm"`
	GC         [2]float64 `json:"-"`
	ProbeAbove [2]float64 `json:"probe_above"`
	MaxDTm     float64    `json:"-"`
	LenPrimer  [2]int     `json:"-"` // inclusive
	LenProbe   [2]int     `json:"-"`
	PoolBound  int        `json:"-"`
	PoolKeep   int        `json:"-"`
	Top        int        `json:"-"`
	Probe      Chem       `json:"-"`
	Workers    int        `json:"-"`
	// Reject (optional) vetoes a triple in the greedy pick, e.g. one that amplifies a non-target with its probe;
	// the next best triple is taken instead.
	Reject func(F, P, R string) bool `json:"-"`
}

// DefaultParams are the Python engine's defaults.
func DefaultParams() Params {
	return Params{Ambig: 2, ProbeAmbig: 2, PassMM: 2, TpN: 2, Amp: [2]int{70, 200}, PrimerTm: [2]float64{52, 66},
		GC: [2]float64{30, 75}, ProbeAbove: [2]float64{2, 14}, MaxDTm: 5, LenPrimer: [2]int{18, 25}, LenProbe: [2]int{20, 30},
		PoolBound: 4000, PoolKeep: 1200, Top: 3, Probe: Chem{Name: "linear"}, Workers: runtime.NumCPU()}
}

// Cand is one evaluated oligo window.
type Cand struct {
	Cnt, S, E int
	Seq       string // 5'->3' as ordered
	Strand    byte   // '+' or '-' (probe); primers '+'
	Tm, GC    float64
	Bits      Bits
	Amb       int
	Kind      Kind
}

// Design is one chosen set (JSON like the Python engine's designs).
type Design struct {
	Name        string  `json:"name"`
	F           string  `json:"F"`
	P           string  `json:"P"`
	R           string  `json:"R"`
	FTm         float64 `json:"F_tm"`
	PTm         float64 `json:"P_tm"`
	RTm         float64 `json:"R_tm"`
	FGC         float64 `json:"F_gc"`
	PGC         float64 `json:"P_gc"`
	RGC         float64 `json:"R_gc"`
	AmpliconLen int     `json:"amplicon_len"`
	StdEst      *bool   `json:"std_est"`
	EstCov      int     `json:"est_cov"`
	EstCovPct   float64 `json:"est_cov_pct"`
	EstAdded    int     `json:"est_added"`
	PStrand     string  `json:"P_strand"`
	FS          int     `json:"-"`
	RE          int     `json:"-"`
	Bits        Bits    `json:"-"`
}

func round1(x float64) float64 { return math.Round(x*10) / 10 }

type poolItem struct{ ub, s, k int }

// parallel runs fn(i) for i in [0, n) on w workers.
func parallel(n, w int, fn func(i int)) {
	if w < 1 {
		w = 1
	}
	var wg sync.WaitGroup
	next := make(chan int, 256)
	for x := 0; x < w; x++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for i := range next {
				fn(i)
			}
		}()
	}
	for i := 0; i < n; i++ {
		next <- i
	}
	close(next)
	wg.Wait()
}

// orient turns a genome-sense oligo into the strands the kind is ordered on.
func orientations(kind Kind, oligo string) [][2]string {
	switch kind {
	case 'P':
		return [][2]string{{"+", oligo}, {"-", RevComp(oligo)}}
	case 'R':
		return [][2]string{{"+", RevComp(oligo)}}
	}
	return [][2]string{{"+", oligo}}
}

// screen applies the per-oligo rules of the Python engine (and the probe chemistry); ok false = rejected.
func screen(kind Kind, seq string, p *Params) (tm, gc float64, ok bool) {
	res, tm, gc, ok := stats(seq)
	if !ok {
		return 0, 0, false
	}
	if kind == 'P' {
		for _, r := range res {
			if r[0] == 'G' || longestRun(r) >= 4 {
				return 0, 0, false
			}
		}
		if !p.Probe.lengthOK(len(seq)) {
			return 0, 0, false
		}
		return tm + p.Probe.tmOffset(), gc, true
	}
	if tm < p.PrimerTm[0] || tm > p.PrimerTm[1] || gc < p.GC[0] || gc > p.GC[1] {
		return 0, 0, false
	}
	for _, r := range res {
		if longestRun(r) >= 5 {
			return 0, 0, false
		}
	}
	return tm, gc, true
}

// Candidates evaluates the most promising windows of one oligo kind (the Python engine's candidates()).
func Candidates(al *Aligned, kind Kind, p *Params) []Cand {
	amb := p.Ambig
	lr := p.LenPrimer
	thr := p.PassMM + amb
	if kind == 'P' {
		amb, lr, thr = p.ProbeAmbig, p.LenProbe, p.ProbeAmbig
	}
	var pool []poolItem
	for k := lr[0]; k <= lr[1]; k++ {
		for _, s := range al.Windows(k) {
			pool = append(pool, poolItem{al.upperBound(s, k, thr), s, k})
		}
	}
	sort.SliceStable(pool, func(a, b int) bool {
		if pool[a].ub != pool[b].ub {
			return pool[a].ub > pool[b].ub
		}
		return pool[a].s < pool[b].s
	})
	if len(pool) > p.PoolBound {
		pool = pool[:p.PoolBound]
	}
	per := make([][]Cand, len(pool))
	parallel(len(pool), p.Workers, func(i int) {
		it := pool[i]
		codes := al.evalCodes(it.s, it.k, kind, amb, p.PassMM, p.TpN)
		oligo := al.oligoWith(it.s, it.k, codes)
		bits := al.PassBits(oligo, it.s, kind, p.PassMM, p.TpN)
		cnt := bits.Count()
		for _, o := range orientations(kind, oligo) {
			tm, gc, ok := screen(kind, o[1], p)
			if !ok {
				continue
			}
			per[i] = append(per[i], Cand{Cnt: cnt, S: it.s, E: it.s + it.k, Seq: o[1], Strand: o[0][0], Tm: tm, GC: gc, Bits: bits, Amb: NAmbig(o[1]), Kind: kind})
		}
	})
	var out []Cand
	for _, c := range per {
		out = append(out, c...)
	}
	sortCands(out)
	if len(out) > p.PoolKeep {
		out = out[:p.PoolKeep]
	}
	return out
}

func sortCands(out []Cand) {
	sort.SliceStable(out, func(a, b int) bool {
		if out[a].Cnt != out[b].Cnt {
			return out[a].Cnt > out[b].Cnt
		}
		if out[a].Amb != out[b].Amb {
			return out[a].Amb < out[b].Amb
		}
		return out[a].S < out[b].S
	})
}

// Triple is a ranked F/R/P combination.
type Triple struct {
	New, Total int
	F, R, P    *Cand
	key        [5]float64
}

type tripleHeap []Triple

func (h tripleHeap) Len() int            { return len(h) }
func (h tripleHeap) Less(i, j int) bool  { return lessKey(h[i].key, h[j].key) }
func (h tripleHeap) Swap(i, j int)       { h[i], h[j] = h[j], h[i] }
func (h *tripleHeap) Push(x interface{}) { *h = append(*h, x.(Triple)) }
func (h *tripleHeap) Pop() interface{}   { o := *h; x := o[len(o)-1]; *h = o[:len(o)-1]; return x }

func lessKey(a, b [5]float64) bool {
	for i := range a {
		if a[i] != b[i] {
			return a[i] < b[i]
		}
	}
	return false
}

// BestTriples ranks F/R/P combinations by the genomes they pass that `covered` does not (then all they pass, fewer
// IUPAC codes, closer primer Tm, found first); the Python engine's best_triples() with the same pruning.
func BestTriples(F, R, P []Cand, covered Bits, p *Params, want int) []Triple {
	h := &tripleHeap{}
	floor, seq := 0, 0
	fp := make(Bits, len(covered))
	for pi := range P {
		pr := &P[pi]
		if pr.Cnt <= floor {
			break
		}
		for fi := range F {
			f := &F[fi]
			if !(f.E <= pr.S && pr.E-f.S <= p.Amp[1]) {
				continue
			}
			if f.Cnt <= floor {
				break
			}
			nfp := 0
			for i := range fp {
				fp[i] = f.Bits[i] & pr.Bits[i] &^ covered[i]
				nfp += popcount(fp[i])
			}
			if nfp <= floor {
				continue
			}
			for ri := range R {
				r := &R[ri]
				if !(r.S >= pr.E && r.E-pr.S <= p.Amp[1]) {
					continue
				}
				if r.Cnt <= floor {
					break
				}
				if !(p.Amp[0] <= r.E-f.S && r.E-f.S <= p.Amp[1]) || math.Abs(f.Tm-r.Tm) > p.MaxDTm {
					continue
				}
				hi := math.Max(f.Tm, r.Tm)
				if !(hi+p.ProbeAbove[0] <= pr.Tm && pr.Tm <= hi+p.ProbeAbove[1]) {
					continue
				}
				nw := 0
				for i := range fp {
					nw += popcount(fp[i] & r.Bits[i])
				}
				if nw <= floor {
					continue
				}
				tot := and3Count(f.Bits, pr.Bits, r.Bits, nil)
				seq++
				t := Triple{New: nw, Total: tot, F: f, R: r, P: pr,
					key: [5]float64{float64(nw), float64(tot), -float64(f.Amb + r.Amb + pr.Amb), -math.Abs(f.Tm - r.Tm), -float64(seq)}}
				if h.Len() < want {
					heap.Push(h, t)
				} else if lessKey((*h)[0].key, t.key) {
					(*h)[0] = t
					heap.Fix(h, 0)
				}
				if h.Len() >= want {
					if m := int((*h)[0].key[0]); m > 0 {
						floor = m - 1
					} else {
						floor = 0
					}
				}
			}
		}
	}
	out := append([]Triple(nil), (*h)...)
	sort.Slice(out, func(a, b int) bool { return lessKey(out[b].key, out[a].key) })
	return out
}

// Run is the result of one design run.
type Run struct {
	N       int
	Status  string
	Designs []Design
	Message string
	F, R, P []Cand
}

// DesignFromAlignment is the coverage-first greedy panel: each round takes the best triple for the genomes the earlier
// rounds miss (a new triple must not repeat an earlier amplicon within 6 bases at both ends).
func DesignFromAlignment(al *Aligned, p Params) Run {
	F := Candidates(al, 'F', &p)
	R := Candidates(al, 'R', &p)
	P := Candidates(al, 'P', &p)
	run := Run{N: al.G, F: F, R: R, P: P}
	if len(F) == 0 || len(R) == 0 || len(P) == 0 {
		run.Status, run.Message = "none", "not enough F/R/probe candidates for the Tm, GC and IUPAC limits"
		return run
	}
	covered := newBits(al.G)
	for i := 0; i < p.Top; i++ {
		cand := BestTriples(F, R, P, covered, &p, 12)
		var pick *Triple
		for ci := range cand {
			c := &cand[ci]
			dup := false
			for _, d := range run.Designs {
				if abs(c.F.S-d.FS) < 6 && abs(c.R.E-d.RE) < 6 {
					dup = true
					break
				}
			}
			if !dup && (p.Reject == nil || !p.Reject(c.F.Seq, c.P.Seq, c.R.Seq)) {
				pick = c
				break
			}
		}
		if pick == nil {
			break
		}
		f, r, pr := pick.F, pick.R, pick.P
		b := f.Bits.Clone()
		for w := range b {
			b[w] &= r.Bits[w] & pr.Bits[w]
		}
		covered.Or(b)
		run.Designs = append(run.Designs, Design{Name: fmt.Sprintf("COV_%d", i), F: f.Seq, P: pr.Seq, R: r.Seq,
			FTm: round1(f.Tm), PTm: round1(pr.Tm), RTm: round1(r.Tm), FGC: round1(f.GC), PGC: round1(pr.GC), RGC: round1(r.GC),
			AmpliconLen: al.AmpliconLen(f.S, r.E), EstCov: pick.Total, EstCovPct: round1(100 * float64(pick.Total) / float64(al.G)),
			EstAdded: pick.New, PStrand: string(pr.Strand), FS: f.S, RE: r.E, Bits: b})
	}
	if len(run.Designs) == 0 {
		run.Status, run.Message = "none", "no F/R/probe combination fits (check amplicon length and Tm limits, or every candidate amplifies a non-target)"
		return run
	}
	run.Status = "ok"
	return run
}

func abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}
