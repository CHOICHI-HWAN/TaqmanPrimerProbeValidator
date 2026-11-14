package design

import (
	"fmt"
	"math"
	"sort"
)

// variants evaluates, for every conserved window start in starts and every length in lens, the consensus with 0..amb
// greedy IUPAC codes (each level is one variant) and returns the variants that pass the per-oligo rules.
// tmOK (may be nil) is an extra Tm test on the ordered oligo.
func variants(al *Aligned, kind Kind, starts []int, lens [2]int, amb int, p *Params, tmOK func(float64) bool) []Cand {
	type job struct{ s, k int }
	var jobs []job
	for _, s := range starts {
		for k := lens[0]; k <= lens[1]; k++ {
			if al.Conserved(s, k) {
				jobs = append(jobs, job{s, k})
			}
		}
	}
	per := make([][]Cand, len(jobs))
	q := *p
	q.PrimerTm = [2]float64{-1e9, 1e9} // shift: Tm is judged against the original oligo (tmOK), GC not at all
	q.GC = [2]float64{0, 100}
	parallel(len(jobs), p.Workers, func(i int) {
		j := jobs[i]
		codes := al.evalCodes(j.s, j.k, kind, amb, p.PassMM, p.TpN)
		seen := map[string]bool{}
		for lv := 0; lv <= len(codes); lv++ {
			oligo := al.oligoWith(j.s, j.k, codes[:lv])
			if seen[oligo] {
				continue
			}
			seen[oligo] = true
			bits := al.PassBits(oligo, j.s, kind, p.PassMM, p.TpN)
			cnt := bits.Count()
			for _, o := range orientations(kind, oligo) {
				pp := p
				if tmOK != nil {
					pp = &q
				}
				tm, gc, ok := screen(kind, o[1], pp)
				if !ok || (tmOK != nil && kind != 'P' && !tmOK(tm)) {
					continue
				}
				per[i] = append(per[i], Cand{Cnt: cnt, S: j.s, E: j.s + j.k, Seq: o[1], Strand: o[0][0], Tm: tm, GC: gc, Bits: bits, Amb: NAmbig(o[1]), Kind: kind})
			}
		}
	})
	var out []Cand
	for _, c := range per {
		out = append(out, c...)
	}
	return out
}

func rangeStep(lo, hi, step int) []int {
	var out []int
	for s := lo; s <= hi; s += max(1, step) {
		out = append(out, s)
	}
	return out
}

// Seed is a set to optimise, as ordered (F and P in any orientation; R as ordered).
type Seed struct{ Name, F, P, R string }

// ShiftOptions of the exhaustive search around a set.
type ShiftOptions struct {
	Shift, Step, Keep   int
	PrimerAmb, ProbeAmb int
	DTm                 float64 // primer Tm may move this much from the original
}

type located struct {
	s   int
	seq string // genome-sense
	n   int
}

func locateEither(al *Aligned, o string) (plus, minus located) {
	if s, n, ok := al.Locate(o); ok {
		plus = located{s, o, n}
	}
	if s, n, ok := al.Locate(RevComp(o)); ok {
		minus = located{s, RevComp(o), n}
	}
	return
}

// Shift searches every variant of the set's F, R and probe: each start moved up to opt.Shift bases either way in
// opt.Step steps, every length of the length ranges, 0..amb greedy IUPAC codes; keeps the opt.Keep best of each and
// combines all of them. Returns the best set for each change class (F, R, P, FR, FP, RP, FRP) that beats the original.
func Shift(al *Aligned, seed Seed, p Params, opt ShiftOptions) ([]Set, Set, error) {
	fp, fm := locateEither(al, seed.F)
	rp, rm := locateEither(al, seed.R)
	// the primer found in genome sense is the forward one (manifests sometimes swap the columns)
	// one primer must sit in genome sense (forward) and the other on the minus strand; manifests sometimes swap the
	// columns, so both assignments are tried and the one more genomes match wins
	var fl, rl located
	s1, s2 := 0, 0
	if fp.n > 0 && rm.n > 0 {
		s1 = fp.n + rm.n
	}
	if rp.n > 0 && fm.n > 0 {
		s2 = rp.n + fm.n
	}
	switch {
	case s1 == 0 && s2 == 0:
		return nil, Set{}, fmt.Errorf("%s: primers not located in the alignment", seed.Name)
	case s1 >= s2:
		fl, rl = fp, rm
	default:
		fl, rl = rp, fm
		seed.F, seed.R = seed.R, seed.F
	}
	if rl.s < fl.s {
		return nil, Set{}, fmt.Errorf("%s: primers located but not facing each other", seed.Name)
	}
	pp, pm := locateEither(al, seed.P)
	pl := pp
	if pm.n > pp.n {
		pl = pm
	}
	if pl.n == 0 {
		return nil, Set{}, fmt.Errorf("%s: probe not located in the alignment", seed.Name)
	}
	orig := func(kind Kind, l located, ordered string) Cand {
		b := al.PassBits(l.seq, l.s, kind, p.PassMM, p.TpN)
		_, tm, gc, _ := stats(ordered)
		return Cand{Cnt: b.Count(), S: l.s, E: l.s + len(l.seq), Seq: ordered, Strand: '+', Tm: tm + map[bool]float64{true: p.Probe.Offset}[kind == 'P'], GC: gc, Bits: b, Amb: NAmbig(ordered), Kind: kind}
	}
	F0, R0, P0 := orig('F', fl, seed.F), orig('R', rl, seed.R), orig('P', pl, seed.P)
	base := Set{Name: seed.Name, F: seed.F, P: seed.P, R: seed.R, Cov: and3Count(F0.Bits, R0.Bits, P0.Bits, nil), FS: F0.S, RE: R0.E, Amplicon: al.AmpliconLen(F0.S, R0.E)}
	near := func(t0 float64) func(float64) bool {
		return func(t float64) bool { return math.Abs(t-t0) <= opt.DTm }
	}
	Fs := variants(al, 'F', rangeStep(F0.S-opt.Shift, F0.S+opt.Shift, opt.Step), p.LenPrimer, opt.PrimerAmb, &p, near(F0.Tm))
	// R: the window's right end is the 5' end; shift the window start
	Rs := variants(al, 'R', rangeStep(R0.S-opt.Shift, R0.S+opt.Shift, opt.Step), p.LenPrimer, opt.PrimerAmb, &p, near(R0.Tm))
	Ps := variants(al, 'P', rangeStep(F0.S-opt.Shift, R0.E+opt.Shift, opt.Step), p.LenProbe, opt.ProbeAmb, &p, nil)
	prune := func(cs []Cand, o Cand) []Cand {
		sort.SliceStable(cs, func(a, b int) bool {
			if cs[a].Cnt != cs[b].Cnt {
				return cs[a].Cnt > cs[b].Cnt
			}
			if cs[a].Amb != cs[b].Amb {
				return cs[a].Amb < cs[b].Amb
			}
			return abs(cs[a].S-o.S) < abs(cs[b].S-o.S)
		})
		if len(cs) > opt.Keep {
			cs = cs[:opt.Keep]
		}
		return append(cs, o)
	}
	Fs, Rs, Ps = prune(Fs, F0), prune(Rs, R0), prune(Ps, P0)
	same := func(a, b *Cand) bool { return a.S == b.S && a.E == b.E && a.Seq == b.Seq }
	type best struct {
		key     [3]float64
		f, r, q *Cand
		cnt     int
	}
	classes := map[string]*best{}
	fr := newBits(al.G)
	for fi := range Fs {
		f := &Fs[fi]
		for ri := range Rs {
			r := &Rs[ri]
			amp := r.E - f.S
			if amp < p.Amp[0] || amp > p.Amp[1] || r.S < f.E {
				continue
			}
			fc, rc := !same(f, &F0), !same(r, &R0)
			nfr := 0
			for w := range fr {
				fr[w] = f.Bits[w] & r.Bits[w]
				nfr += popcount(fr[w])
			}
			if nfr <= base.Cov {
				continue
			}
			hi := math.Max(f.Tm, r.Tm)
			for qi := range Ps {
				q := &Ps[qi]
				if q.Cnt <= base.Cov || q.S < f.E || q.E > r.S {
					continue
				}
				qc := !same(q, &P0)
				if !fc && !rc && !qc {
					continue
				}
				// as optimize_user_sets.py: each primer stays within DTm of its original, the probe at least
				// ProbeAbove[0] above the higher primer (no upper limit)
				if q.Tm < hi+p.ProbeAbove[0] {
					continue
				}
				cnt := 0
				for w := range fr {
					cnt += popcount(fr[w] & q.Bits[w])
				}
				if cnt <= base.Cov {
					continue
				}
				cls := ""
				if fc {
					cls += "F"
				}
				if rc {
					cls += "R"
				}
				if qc {
					cls += "P"
				}
				dist := abs(f.S-F0.S) + abs(r.E-R0.E) + abs(q.S-P0.S)
				key := [3]float64{float64(cnt), -float64(f.Amb + r.Amb + q.Amb), -float64(dist)}
				if b := classes[cls]; b == nil || greater3(key, b.key) {
					classes[cls] = &best{key, f, r, q, cnt}
				}
			}
		}
	}
	var out []Set
	for _, cls := range []string{"P", "F", "R", "FR", "FP", "RP", "FRP"} {
		b := classes[cls]
		if b == nil {
			continue
		}
		bits := b.f.Bits.Clone()
		for w := range bits {
			bits[w] &= b.r.Bits[w] & b.q.Bits[w]
		}
		out = append(out, Set{Name: seed.Name, Change: cls, BaseCov: base.Cov, Cov: b.cnt, Added: b.cnt - base.Cov, F: b.f.Seq, P: b.q.Seq, R: b.r.Seq,
			FS: b.f.S, RE: b.r.E, Amplicon: al.AmpliconLen(b.f.S, b.r.E), Bits: bits})
	}
	return out, base, nil
}

func greater3(a, b [3]float64) bool {
	for i := range a {
		if a[i] != b[i] {
			return a[i] > b[i]
		}
	}
	return false
}

// TileRow is the best set found in one amplicon-sized window of the alignment.
type TileRow struct {
	Start int
	Set   Set
	Found bool
}

// Tile evaluates EVERY conserved window of every length (0..amb IUPAC codes each), then for each region
// [x, x+amp_max) with x stepping by step finds the best triple inside it. The result is a coverage map along the genome.
func Tile(al *Aligned, p Params, step, keep int) []TileRow {
	all := func(kind Kind, lens [2]int, amb int) []Cand {
		cs := variants(al, kind, rangeStep(0, al.L-lens[0], 1), lens, amb, &p, nil)
		sortCands(cs)
		return cs
	}
	F := all('F', p.LenPrimer, p.Ambig)
	R := all('R', p.LenPrimer, p.Ambig)
	P := all('P', p.LenProbe, p.ProbeAmb())
	var starts []int
	for x := 0; x < al.L; x += max(1, step) {
		starts = append(starts, x)
	}
	rows := make([]TileRow, len(starts))
	empty := newBits(al.G)
	parallel(len(starts), p.Workers, func(i int) {
		x := starts[i]
		lim := x + p.Amp[1]
		pick := func(cs []Cand) []Cand {
			var o []Cand
			for _, c := range cs {
				if c.S >= x && c.E <= lim {
					o = append(o, c)
					if len(o) >= keep {
						break
					}
				}
			}
			return o
		}
		f, r, q := pick(F), pick(R), pick(P)
		rows[i].Start = x
		if len(f) == 0 || len(r) == 0 || len(q) == 0 {
			return
		}
		t := BestTriples(f, r, q, empty, &p, 1)
		if len(t) == 0 {
			return
		}
		b := t[0]
		bits := b.F.Bits.Clone()
		for w := range bits {
			bits[w] &= b.R.Bits[w] & b.P.Bits[w]
		}
		rows[i] = TileRow{Start: x, Found: true, Set: Set{Name: fmt.Sprintf("TILE_%d", x), F: b.F.Seq, P: b.P.Seq, R: b.R.Seq, Cov: b.Total,
			FS: b.F.S, RE: b.R.E, Amplicon: al.AmpliconLen(b.F.S, b.R.E), Bits: bits}}
	})
	return rows
}

// ProbeAmb is the probe IUPAC limit.
func (p Params) ProbeAmb() int { return p.ProbeAmbig }

// GreedyPanel picks, from tile results, sets that together cover the most genomes (each next one adds the most).
func GreedyPanel(rows []TileRow, top int, n int) []Set {
	covered := newBits(n)
	var out []Set
	used := map[[2]int]bool{}
	for len(out) < top {
		bi, bn := -1, 0
		for i, r := range rows {
			if !r.Found || used[[2]int{r.Set.FS, r.Set.RE}] {
				continue
			}
			add := 0
			for w := range covered {
				add += popcount(r.Set.Bits[w] &^ covered[w])
			}
			if add > bn || (add == bn && bi >= 0 && add > 0 && r.Set.Cov > rows[bi].Set.Cov) {
				bi, bn = i, add
			}
		}
		if bi < 0 || bn == 0 {
			break
		}
		s := rows[bi].Set
		s.Added = bn
		s.Name = fmt.Sprintf("TILE_%d", len(out))
		used[[2]int{s.FS, s.RE}] = true
		covered.Or(s.Bits)
		out = append(out, s)
	}
	return out
}
