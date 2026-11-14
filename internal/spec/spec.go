// Package spec finds primer binding sites and the amplicons they form across a genome database — the clean-room
// counterpart of `mfeprimer --virus`.
package spec

import (
	"math"
	"runtime"
	"sort"
	"sync"

	"chichipes/internal/primers"
	"chichipes/internal/seqio"
	"chichipes/internal/thermo"
)

// Options mirror the mfeprimer flags the pipeline uses.
type Options struct {
	K        int     // 3'-window length (-k, 9)
	MisMatch int     // --misMatch
	MisStart int     // --misStart (1)
	MisEnd   int     // --misEnd (9)
	MinSize  int     // -s
	MaxSize  int     // -S (2000)
	MinTm    float64 // -t (30)
	MaxTm    float64 // -T (100)
	CPU      int
	Cond     thermo.Cond
}

// DefaultOptions are mfeprimer's defaults.
func DefaultOptions() Options {
	return Options{K: 9, MisStart: 1, MisEnd: 9, MaxSize: 2000, MinTm: 30, MaxTm: 100, CPU: runtime.NumCPU(),
		Cond: thermo.Default}
}

// probe is one searchable sequence: a variant of a concrete primer.
type probe struct {
	name string
	seq  []byte // 5'->3'
	rcs  []byte // reverse complement
}

// Site is one primer binding site on a genome, in plus-strand coordinates (0-based start, length N).
type Site struct {
	Pos, N int
	Rev    bool // primer binds the plus strand: its reverse complement is read at Pos
	P      *probe
	Tm, Dg float64
}

// Name and Seq of the oligo that binds (the variant), Seq in the primer's own 5'->3' orientation.
func (s Site) Name() string { return s.P.name }
func (s Site) Seq() []byte  { return s.P.seq }

// BorderBand is how close (°C) a binding Tm may lie to the -t cut-off to be listed as a borderline case. The Tm now
// equals mfeprimer's (docs/FINDINGS.md §3), so the list only flags decisions that hinge on the last decimals.
const BorderBand = 0.02

// Border is one binding whose Tm lies within BorderBand of the cut-off.
type Border struct {
	Site
	Kept bool
}

// Amplicon is one product: forward site F (binds the minus strand) and reverse site R (binds the plus strand).
type Amplicon struct {
	Rec  int // index of the genome record
	F, R Site
}

// Start (0-based) and End (exclusive) of the product on the genome.
func (a Amplicon) Start() int { return a.F.Pos }
func (a Amplicon) End() int   { return a.R.Pos + a.R.N }
func (a Amplicon) Size() int  { return a.End() - a.Start() }

// Engine holds the probes indexed by their 3' window.
type Engine struct {
	opt    Options
	probes []*probe
	fwd    map[uint32][]*probe // 3' window (as written 5'->3')        -> probes
	rev    map[uint32][]*probe // reverse complement of the 3' window -> probes
	maxLen int
}

var comp [256]byte

func init() {
	for i := range comp {
		comp[i] = 'N'
	}
	// IUPAC complements, both cases
	for _, p := range []string{"AT", "CG", "GC", "TA", "RY", "YR", "KM", "MK", "SS", "WW", "BV", "VB", "DH", "HD", "NN", "II"} {
		comp[p[0]], comp[p[0]+32] = p[1], p[1]+32
	}
}

func revComp(s []byte) []byte {
	out := make([]byte, len(s))
	for i, b := range s {
		out[len(s)-1-i] = comp[b]
	}
	return out
}

func upper(b byte) byte {
	if b >= 'a' && b <= 'z' {
		return b - 32
	}
	return b
}

// encode packs an A/C/G/T k-mer into 2 bits per base; ok is false when it holds another letter.
func encode(s []byte) (uint32, bool) {
	var v uint32
	for _, b := range s {
		c := thermo.Code(b)
		if c < 0 {
			return 0, false
		}
		v = v<<2 | uint32(c)
	}
	return v, true
}

// encodeAll is encode for an oligo window that may hold inosine: one key per base an I can face (4^m keys for m
// inosines), so the window is found whatever the genome has there.
func encodeAll(s []byte) []uint32 {
	keys := []uint32{0}
	for _, b := range s {
		if b == 'I' || b == 'i' {
			next := make([]uint32, 0, 4*len(keys))
			for _, v := range keys {
				for c := uint32(0); c < 4; c++ {
					next = append(next, v<<2|c)
				}
			}
			keys = next
			continue
		}
		c := thermo.Code(b)
		if c < 0 {
			return nil
		}
		for i := range keys {
			keys[i] = keys[i]<<2 | uint32(c)
		}
	}
	return keys
}

// New builds the search set from the input oligos.
func New(oligos []primers.Oligo, opt Options) (*Engine, error) {
	e := &Engine{opt: opt, fwd: map[uint32][]*probe{}, rev: map[uint32][]*probe{}}
	for _, o := range oligos {
		cs, err := primers.Build(o, opt.MisMatch, opt.MisStart, opt.MisEnd)
		if err != nil {
			return nil, err
		}
		seen := map[string]bool{} // an expansion equal to an earlier variant is searched once, under the first name
		for _, c := range cs {
			for _, v := range c.Variants {
				if seen[v.Seq] {
					continue
				}
				seen[v.Seq] = true
				p := &probe{name: c.VariantName(v), seq: []byte(v.Seq), rcs: revComp([]byte(v.Seq))}
				e.probes = append(e.probes, p)
				if len(p.seq) < opt.K {
					continue
				}
				for _, w := range encodeAll(p.seq[len(p.seq)-opt.K:]) {
					e.fwd[w] = append(e.fwd[w], p)
				}
				for _, w := range encodeAll(p.rcs[:opt.K]) {
					e.rev[w] = append(e.rev[w], p)
				}
				if len(p.seq) > e.maxLen {
					e.maxLen = len(p.seq)
				}
			}
		}
	}
	return e, nil
}

// Sites lists every binding site on one genome that passes the Tm window.
func (e *Engine) Sites(g []byte) (fs, rs []Site, border []Border) {
	k := e.opt.K
	n := len(g)
	mask := uint32(1)<<(2*uint(k)) - 1
	var w uint32
	valid := 0 // number of consecutive A/C/G/T bases ending at the current position
	for i := 0; i < n; i++ {
		c := thermo.Code(g[i])
		if c < 0 {
			valid = 0
			w = 0
			continue
		}
		w = (w<<2 | uint32(c)) & mask
		valid++
		if valid < k {
			continue
		}
		end := i + 1 // the window is g[end-k:end]
		for _, p := range e.fwd[w] {
			pl := len(p.seq)
			pos := end - pl
			if pos < 0 {
				continue
			}
			l5 := byte('-') // the site starts the genome: no base 5' of it
			if pos > 0 {
				l5 = upper(g[pos-1])
			}
			s := make([]byte, pl)
			for j := 0; j < pl; j++ {
				s[j] = upper(g[pos+j])
			}
			r3 := byte('A')
			if end < n {
				r3 = upper(g[end])
			}
			tm, dg, ok := thermo.Binding(p.seq, s, l5, r3, false, e.opt.Cond)
			site := Site{Pos: pos, N: pl, P: p, Tm: tm, Dg: dg}
			keep := ok && tm >= e.opt.MinTm && tm <= e.opt.MaxTm
			if keep {
				fs = append(fs, site)
			}
			if ok && math.Abs(tm-e.opt.MinTm) <= BorderBand {
				border = append(border, Border{site, keep})
			}
		}
		start := end - k // the reverse complement of the primer's 3' window starts here
		for _, p := range e.rev[w] {
			pl := len(p.seq)
			if start+pl > n || start < 1 { // (a site at the very start cannot be the right end of an amplicon)
				continue
			}
			s := make([]byte, pl) // the site read on the minus strand, 5'->3' like the primer
			for j := 0; j < pl; j++ {
				s[j] = comp[upper(g[start+pl-1-j])]
			}
			l5 := byte('-') // the site ends the genome
			if start+pl < n {
				l5 = comp[upper(g[start+pl])]
			}
			r3 := comp[upper(g[start-1])]
			tm, dg, ok := thermo.Binding(p.seq, s, l5, r3, true, e.opt.Cond)
			site := Site{Pos: start, N: pl, Rev: true, P: p, Tm: tm, Dg: dg}
			keep := ok && tm >= e.opt.MinTm && tm <= e.opt.MaxTm
			if keep {
				rs = append(rs, site)
			}
			if ok && math.Abs(tm-e.opt.MinTm) <= BorderBand {
				border = append(border, Border{site, keep})
			}
		}
	}
	return fs, rs, border
}

// Amplicons pairs the sites of one genome into products within the size limits.
func (e *Engine) Amplicons(rec int, fs, rs []Site) []Amplicon {
	sort.Slice(rs, func(a, b int) bool { return rs[a].Pos < rs[b].Pos })
	var out []Amplicon
	for _, f := range fs {
		for _, r := range rs {
			a := Amplicon{Rec: rec, F: f, R: r}
			if r.Pos < f.Pos {
				continue
			}
			sz := a.Size()
			if sz > e.opt.MaxSize {
				break
			}
			// the two primers may overlap by fewer than K bases (larger overlaps are never reported)
			if sz < e.opt.MinSize || f.Pos+f.N-r.Pos >= e.opt.K {
				continue
			}
			out = append(out, a)
		}
	}
	return out
}

// Result is the outcome over a whole database.
type Result struct {
	Amplicons []Amplicon
	Border    [][]Border // per record: bindings near the Tm cut-off
	Failed    []int      // records without any amplicon, in database order
	Total     int
}

// Run searches every record in parallel and returns the amplicons in database order.
func (e *Engine) Run(db []seqio.Record) Result {
	per := make([][]Amplicon, len(db))
	border := make([][]Border, len(db))
	jobs := make(chan int, len(db))
	for i := range db {
		jobs <- i
	}
	close(jobs)
	cpu := e.opt.CPU
	if cpu < 1 {
		cpu = 1
	}
	var wg sync.WaitGroup
	for w := 0; w < cpu; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for i := range jobs {
				fs, rs, b := e.Sites(db[i].Seq)
				per[i] = e.Amplicons(i, fs, rs)
				border[i] = b
			}
		}()
	}
	wg.Wait()
	res := Result{Total: len(db), Border: border}
	for i, a := range per {
		if len(a) == 0 {
			res.Failed = append(res.Failed, i)
		}
		res.Amplicons = append(res.Amplicons, a...)
	}
	return res
}
