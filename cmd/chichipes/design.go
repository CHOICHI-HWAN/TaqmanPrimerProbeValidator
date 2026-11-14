package main

import (
	"bufio"
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"
	"sync"
	"time"

	"chichipes/internal/design"
)

const designUsage = `chichipes design: primer / probe design from an alignment

  chichipes design cover -a aligned.fa [options]            coverage-first panel (greedy set cover)
  chichipes design sweep -a aligned.fa --grid mega -o results.jsonl [--jobs N]
                                                            many parameter combinations (design_sweep.py compatible)
  chichipes design shift -a aligned.fa --sets sets.tsv [--shift 50 --step 1]
                                                            exhaustive shifted / resized / degenerate variants of sets
  chichipes design tile  -a aligned.fa [--step 10]          best set in every amplicon-sized window along the genome
  chichipes design labels --sets primers.tsv [--groups set_groups.tsv]
                                                            reporter/quencher pairs; channel clashes inside a group

Common options (cover, shift, tile):
  --ambig N --probe-ambig N   IUPAC codes allowed in a primer / the probe (2 / 2)
  --mm N --tp N               primer mismatches allowed / 3' bases that must match (2 / 2)
  --amp 70,200                amplicon length range
  --primer-tm 52,66           primer Tm window (screening Tm, SantaLucia 1998, 250 nM / 50 mM Na+)
  --probe-above 2,14          probe Tm above the higher primer Tm
  --primer-len 18,25 --probe-len 20,30
  --probe-chem linear|zen|mgb probe format (length limits, Tm offset; see README)
  --probe-tm-offset X         override the format's probe Tm offset
  --exclude nontarget.fa      reject sets that would amplify a non-target sequence; penalise oligos that bind one
  --top N                     designs to report (3)
  --cpu N                     worker threads
`

type designFlags struct {
	fs                                    *flag.FlagSet
	aln, out, exclude, sets, grid, chem   string
	amp, primerTm, probeAbove, pLen, qLen string
	ambig, probeAmbig, mm, tp, top, cpu   int
	shift, step, keep, jobs, maxAmb       int
	chemOffset                            float64
	qc                                    bool
	groups                                string
}

func newDesignFlags(name string) *designFlags {
	d := &designFlags{fs: flag.NewFlagSet("design "+name, flag.ExitOnError)}
	f := d.fs
	f.StringVar(&d.aln, "a", "", "aligned FASTA")
	f.StringVar(&d.out, "o", "", "output file (default: stdout / <name>.tsv)")
	f.StringVar(&d.exclude, "exclude", "", "non-target FASTA")
	f.StringVar(&d.sets, "sets", "", "sets to optimise: TSV with name, F, R, P columns (or primers.tsv)")
	f.StringVar(&d.grid, "grid", "wide", "sweep grid: wide, hightm, midtm, longoligo, mega, mega2, mega3")
	f.StringVar(&d.chem, "probe-chem", "linear", "probe format: linear, zen, mgb")
	f.Float64Var(&d.chemOffset, "probe-tm-offset", math.NaN(), "probe Tm offset (C)")
	f.StringVar(&d.amp, "amp", "70,200", "amplicon length range")
	f.StringVar(&d.primerTm, "primer-tm", "52,66", "primer Tm range")
	f.StringVar(&d.probeAbove, "probe-above", "2,14", "probe Tm above the higher primer")
	f.StringVar(&d.pLen, "primer-len", "18,25", "primer length range")
	f.StringVar(&d.qLen, "probe-len", "", "probe length range (default 20,30 or the format's)")
	f.IntVar(&d.ambig, "ambig", 2, "IUPAC codes per primer")
	f.IntVar(&d.probeAmbig, "probe-ambig", 2, "IUPAC codes in the probe")
	f.IntVar(&d.mm, "mm", 2, "primer mismatches allowed")
	f.IntVar(&d.tp, "tp", 2, "3' bases that must match")
	f.IntVar(&d.top, "top", 3, "designs to report")
	f.IntVar(&d.cpu, "cpu", runtime.NumCPU(), "worker threads")
	f.IntVar(&d.shift, "shift", 50, "shift: bases each oligo may move either way")
	f.IntVar(&d.step, "step", 1, "shift / tile: step in bases")
	f.IntVar(&d.keep, "keep", 300, "shift / tile: variants kept per oligo for the combination step")
	f.IntVar(&d.maxAmb, "max-amb", -1, "shift: IUPAC codes a variant may have (default --ambig / --probe-ambig)")
	f.IntVar(&d.jobs, "jobs", max(1, runtime.NumCPU()/4), "sweep: runs in parallel")
	f.String("segment", "", "shift: only sets of this segment (primers.tsv)")
	f.StringVar(&d.groups, "groups", "", "labels: set groups (group_id, label, member_ids) run in one tube")
	f.BoolVar(&d.qc, "qc", true, "check the chosen sets with chichipes thermodynamics (Tm 2004, hairpin, dimer)")
	f.Usage = func() { fmt.Fprint(os.Stderr, designUsage) }
	return d
}

func pair(s string) ([2]float64, error) {
	p := strings.Split(s, ",")
	if len(p) != 2 {
		return [2]float64{}, fmt.Errorf("%q: expected lo,hi", s)
	}
	a, err1 := strconv.ParseFloat(strings.TrimSpace(p[0]), 64)
	b, err2 := strconv.ParseFloat(strings.TrimSpace(p[1]), 64)
	if err1 != nil || err2 != nil {
		return [2]float64{}, fmt.Errorf("%q: expected lo,hi", s)
	}
	return [2]float64{a, b}, nil
}

func ipair(s string) ([2]int, error) {
	p, err := pair(s)
	return [2]int{int(p[0]), int(p[1])}, err
}

func (d *designFlags) params() (design.Params, error) {
	p := design.DefaultParams()
	var err error
	p.Ambig, p.ProbeAmbig, p.PassMM, p.TpN, p.Top, p.Workers = d.ambig, d.probeAmbig, d.mm, d.tp, d.top, d.cpu
	if p.Amp, err = ipair(d.amp); err != nil {
		return p, err
	}
	if p.PrimerTm, err = pair(d.primerTm); err != nil {
		return p, err
	}
	if p.ProbeAbove, err = pair(d.probeAbove); err != nil {
		return p, err
	}
	if p.LenPrimer, err = ipair(d.pLen); err != nil {
		return p, err
	}
	if p.Probe, err = design.ParseChem(d.chem, d.chemOffset, 0, 0); err != nil {
		return p, err
	}
	if d.qLen != "" {
		if p.LenProbe, err = ipair(d.qLen); err != nil {
			return p, err
		}
	} else if p.Probe.Name != "linear" {
		p.LenProbe = p.Probe.ProbeLenRange()
	}
	return p, nil
}

func loadAligned(path string) (*design.Aligned, error) {
	if path == "" {
		return nil, fmt.Errorf("-a aligned.fa is required")
	}
	ids, rows, err := design.ReadAligned(path)
	if err != nil {
		return nil, err
	}
	return design.NewAligned(ids, rows)
}

func runDesign(args []string) error {
	if len(args) == 0 || args[0] == "-h" || args[0] == "--help" {
		fmt.Print(designUsage)
		return nil
	}
	d := newDesignFlags(args[0])
	if err := d.fs.Parse(args[1:]); err != nil {
		return err
	}
	switch args[0] {
	case "cover":
		return designCover(d)
	case "sweep":
		return designSweep(d)
	case "shift":
		return designShift(d)
	case "tile":
		return designTile(d)
	case "labels":
		return designLabels(d)
	}
	return fmt.Errorf("unknown design command %q (cover, sweep, shift, tile)", args[0])
}

func outWriter(path string) (*bufio.Writer, func() error, error) {
	if path == "" || path == "-" {
		w := bufio.NewWriter(os.Stdout)
		return w, w.Flush, nil
	}
	f, err := os.Create(path)
	if err != nil {
		return nil, nil, err
	}
	w := bufio.NewWriter(f)
	return w, func() error {
		if err := w.Flush(); err != nil {
			return err
		}
		return f.Close()
	}, nil
}

func designCover(d *designFlags) error {
	al, err := loadAligned(d.aln)
	if err != nil {
		return err
	}
	p, err := d.params()
	if err != nil {
		return err
	}
	var xs *design.Exclusion
	if d.exclude != "" {
		if xs, err = design.LoadExclusion(d.exclude, p); err != nil {
			return err
		}
		p.Reject = xs.Rejects
	}
	t0 := time.Now()
	run := design.DesignFromAlignment(al, p)
	w, done, err := outWriter(d.out)
	if err != nil {
		return err
	}
	fmt.Fprintf(os.Stderr, "%d genomes; candidates F %d R %d P %d; %d designs in %.2fs\n", al.G, len(run.F), len(run.R), len(run.P), len(run.Designs), time.Since(t0).Seconds())
	if run.Status != "ok" {
		fmt.Fprintln(os.Stderr, run.Message)
	}
	var sets []design.Set
	for _, x := range run.Designs {
		sets = append(sets, design.Set{Name: x.Name, F: x.F, P: x.P, R: x.R, Cov: x.EstCov, Added: x.EstAdded, Amplicon: x.AmpliconLen})
	}
	if err := writeSets(w, al, sets, p, xs, d.qc); err != nil {
		return err
	}
	return done()
}

func writeSets(w *bufio.Writer, al *design.Aligned, sets []design.Set, p design.Params, xs *design.Exclusion, qc bool) error {
	fmt.Fprintln(w, strings.Join(design.SetColumns, "\t"))
	for i := range sets {
		s := &sets[i]
		design.PreferCStrand(s)
		if qc {
			design.QC(s, p)
		}
		if xs != nil {
			xs.Check(s)
		}
		fmt.Fprintln(w, strings.Join(s.Row(al.G), "\t"))
	}
	return nil
}

// ---- sweep ----

type gridSpec struct {
	keys []string
	vals [][]interface{}
}

func grid(name string) (gridSpec, error) {
	f := func(xs ...float64) []interface{} {
		o := make([]interface{}, len(xs))
		for i, x := range xs {
			o[i] = x
		}
		return o
	}
	ii := func(xs ...int) []interface{} {
		o := make([]interface{}, len(xs))
		for i, x := range xs {
			o[i] = x
		}
		return o
	}
	pr := func(xs ...[2]float64) []interface{} {
		o := make([]interface{}, len(xs))
		for i, x := range xs {
			o[i] = x
		}
		return o
	}
	ap := func(xs ...[2]int) []interface{} {
		o := make([]interface{}, len(xs))
		for i, x := range xs {
			o[i] = x
		}
		return o
	}
	_ = f
	g := gridSpec{}
	add := func(k string, v []interface{}) { g.keys = append(g.keys, k); g.vals = append(g.vals, v) }
	switch name {
	case "hightm":
		add("ambig", ii(2, 3, 4, 5, 6))
		add("probe_ambig", ii(2, 3, 4))
		add("pass_mm", ii(1))
		add("tp_n", ii(2))
		add("amp", ap([2]int{70, 200}, [2]int{70, 250}))
		add("primer_tm", pr([2]float64{56, 66}, [2]float64{58, 68}))
		add("probe_above", pr([2]float64{6, 14}, [2]float64{8, 16}))
	case "longoligo":
		add("ambig", ii(3, 4, 5, 6))
		add("probe_ambig", ii(2, 3, 4))
		add("pass_mm", ii(1))
		add("tp_n", ii(2))
		add("amp", ap([2]int{70, 200}, [2]int{70, 250}))
		add("len_primer_max", ii(30))
		add("len_probe_max", ii(36))
		add("primer_tm", pr([2]float64{54, 68}, [2]float64{56, 68}))
		add("probe_above", pr([2]float64{3, 16}, [2]float64{6, 16}))
	case "midtm":
		add("ambig", ii(3, 4, 5, 6))
		add("probe_ambig", ii(2, 3, 4))
		add("pass_mm", ii(1))
		add("tp_n", ii(2))
		add("amp", ap([2]int{70, 200}, [2]int{70, 250}))
		add("primer_tm", pr([2]float64{54, 66}, [2]float64{55, 67}))
		add("probe_above", pr([2]float64{3, 14}, [2]float64{5, 14}))
	case "mega":
		add("ambig", ii(1, 2, 3, 4, 5, 6))
		add("probe_ambig", ii(1, 2, 3, 4, 5))
		add("pass_mm", ii(1, 2, 3))
		add("tp_n", ii(1, 2))
		add("amp", ap([2]int{60, 150}, [2]int{70, 200}, [2]int{70, 250}, [2]int{70, 300}, [2]int{100, 250}, [2]int{150, 300}))
		add("primer_tm", pr([2]float64{50, 64}, [2]float64{52, 66}, [2]float64{54, 68}, [2]float64{56, 70}))
		add("probe_above", pr([2]float64{2, 12}, [2]float64{4, 14}, [2]float64{6, 16}))
		add("len_primer_max", ii(24, 30))
		add("len_probe_max", ii(28, 36))
	case "mega2":
		add("ambig", ii(1, 2, 3, 4, 5, 6))
		add("probe_ambig", ii(1, 2, 3, 4, 5))
		add("pass_mm", ii(1, 2, 3))
		add("tp_n", ii(1, 2))
		add("amp", ap([2]int{70, 180}, [2]int{90, 220}, [2]int{110, 260}, [2]int{130, 280}, [2]int{150, 320}, [2]int{180, 350}))
		add("primer_tm", pr([2]float64{51, 63}, [2]float64{53, 65}, [2]float64{55, 67}, [2]float64{57, 69}, [2]float64{59, 71}))
		add("probe_above", pr([2]float64{3, 13}, [2]float64{5, 15}, [2]float64{7, 17}))
		add("len_primer_max", ii(28, 32))
		add("len_probe_max", ii(32, 40))
	case "mega3":
		add("ambig", ii(1, 2, 3, 4, 5, 6))
		add("probe_ambig", ii(1, 2, 3, 4, 5))
		add("pass_mm", ii(1, 2, 3))
		add("tp_n", ii(1, 2, 3))
		add("amp", ap([2]int{60, 160}, [2]int{80, 210}, [2]int{100, 240}, [2]int{120, 270}, [2]int{140, 300}, [2]int{160, 330}, [2]int{200, 360}))
		add("primer_tm", pr([2]float64{48, 62}, [2]float64{60, 72}, [2]float64{53, 69}))
		add("probe_above", pr([2]float64{2, 10}, [2]float64{4, 16}, [2]float64{8, 18}))
		add("len_primer_max", ii(26, 32))
		add("len_probe_max", ii(30, 42))
	case "wide":
		add("ambig", ii(1, 2, 3, 4, 5))
		add("probe_ambig", ii(1, 2, 3, 4))
		add("pass_mm", ii(1, 2))
		add("tp_n", ii(2))
		add("amp", ap([2]int{70, 200}, [2]int{70, 300}))
		add("primer_tm", pr([2]float64{52, 66}, [2]float64{50, 68}))
	default:
		return g, fmt.Errorf("unknown grid %q", name)
	}
	return g, nil
}

type orderedParams struct {
	keys []string
	vals []interface{}
}

func (o orderedParams) MarshalJSON() ([]byte, error) {
	var b strings.Builder
	b.WriteByte('{')
	for i, k := range o.keys {
		if i > 0 {
			b.WriteString(", ")
		}
		v := o.vals[i]
		var js []byte
		switch x := v.(type) {
		case [2]float64:
			js = []byte(fmt.Sprintf("[%s, %s]", pyFloat(x[0]), pyFloat(x[1])))
		default:
			js, _ = json.Marshal(v)
		}
		fmt.Fprintf(&b, "%q: %s", k, js)
	}
	b.WriteByte('}')
	return []byte(b.String()), nil
}

func pyFloat(x float64) string {
	s := strconv.FormatFloat(x, 'f', -1, 64)
	if !strings.Contains(s, ".") {
		s += ".0"
	}
	return s
}

func (o orderedParams) apply(p *design.Params) {
	for i, k := range o.keys {
		switch v := o.vals[i].(type) {
		case int:
			switch k {
			case "ambig":
				p.Ambig = v
			case "probe_ambig":
				p.ProbeAmbig = v
			case "pass_mm":
				p.PassMM = v
			case "tp_n":
				p.TpN = v
			case "len_primer_max":
				p.LenPrimer[1] = v
			case "len_probe_max":
				p.LenProbe[1] = v
			}
		case [2]int:
			p.Amp = v
		case [2]float64:
			if k == "primer_tm" {
				p.PrimerTm = v
			} else {
				p.ProbeAbove = v
			}
		}
	}
}

// canonical key (sorted keys, as design_sweep.py's resume check).
func canonKey(m map[string]interface{}) string {
	b, _ := json.Marshal(m)
	return string(b)
}

func designSweep(d *designFlags) error {
	al, err := loadAligned(d.aln)
	if err != nil {
		return err
	}
	g, err := grid(d.grid)
	if err != nil {
		return err
	}
	base, err := d.params()
	if err != nil {
		return err
	}
	if d.top == 3 {
		base.Top = 10
	}
	var todo []orderedParams
	idx := make([]int, len(g.keys))
	for {
		o := orderedParams{keys: g.keys}
		for i := range g.keys {
			o.vals = append(o.vals, g.vals[i][idx[i]])
		}
		todo = append(todo, o)
		i := len(idx) - 1
		for i >= 0 {
			idx[i]++
			if idx[i] < len(g.vals[i]) {
				break
			}
			idx[i] = 0
			i--
		}
		if i < 0 {
			break
		}
	}
	if d.out == "" {
		return fmt.Errorf("-o results.jsonl is required")
	}
	done := map[string]bool{}
	if b, err := os.ReadFile(d.out); err == nil {
		for _, ln := range strings.Split(string(b), "\n") {
			var r struct {
				Params map[string]interface{} `json:"params"`
			}
			if json.Unmarshal([]byte(ln), &r) == nil && r.Params != nil {
				done[canonKey(r.Params)] = true
			}
		}
	}
	var left []orderedParams
	for _, o := range todo {
		js, _ := o.MarshalJSON()
		var m map[string]interface{}
		_ = json.Unmarshal(js, &m)
		if !done[canonKey(m)] {
			left = append(left, o)
		}
	}
	segment := filepath.Base(filepath.Dir(d.aln))
	fmt.Fprintf(os.Stderr, "%d genomes; %d runs to do (%d done), %d in parallel\n", al.G, len(left), len(todo)-len(left), d.jobs)
	f, err := os.OpenFile(d.out, os.O_CREATE|os.O_APPEND|os.O_WRONLY, 0o644)
	if err != nil {
		return err
	}
	defer f.Close()
	var mu sync.Mutex
	t0 := time.Now()
	n := 0
	per := max(1, d.cpu/max(1, d.jobs))
	jobs := make(chan orderedParams)
	var wg sync.WaitGroup
	for j := 0; j < d.jobs; j++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for o := range jobs {
				p := base
				p.Workers = per
				o.apply(&p)
				t := time.Now()
				run := design.DesignFromAlignment(al, p)
				designs := run.Designs
				if designs == nil {
					designs = []design.Design{}
				}
				rec := struct {
					Segment string          `json:"segment"`
					Params  orderedParams   `json:"params"`
					N       int             `json:"n"`
					Status  string          `json:"status"`
					Designs []design.Design `json:"designs"`
					Seconds float64         `json:"seconds"`
					Engine  string          `json:"engine"`
				}{segment, o, run.N, run.Status, designs, math.Round(time.Since(t).Seconds()*100) / 100, "chichipes"}
				js, _ := json.Marshal(rec)
				mu.Lock()
				f.Write(append(js, '\n'))
				n++
				if n%50 == 0 || n == len(left) {
					fmt.Fprintf(os.Stderr, "[%d/%d] %.0fs\n", n, len(left), time.Since(t0).Seconds())
				}
				mu.Unlock()
			}
		}()
	}
	for _, o := range left {
		jobs <- o
	}
	close(jobs)
	wg.Wait()
	return nil
}

// ---- shift / tile ----

func readSeeds(path, segment string) ([]design.Seed, error) {
	b, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	lines := strings.Split(strings.TrimRight(string(b), "\n"), "\n")
	if len(lines) < 2 {
		return nil, fmt.Errorf("%s: no sets", path)
	}
	hdr := strings.Split(lines[0], "\t")
	col := func(names ...string) int {
		for _, n := range names {
			for i, h := range hdr {
				if strings.EqualFold(strings.TrimSpace(h), n) {
					return i
				}
			}
		}
		return -1
	}
	cn, cf, cr, cp := col("name", "primer_name", "primer_id"), col("F", "forward_primer"), col("R", "reverse_primer"), col("P", "probe")
	cs := col("segment")
	if cf < 0 || cr < 0 || cp < 0 {
		return nil, fmt.Errorf("%s: needs F, R and P (or forward_primer, reverse_primer, probe) columns", path)
	}
	var out []design.Seed
	for i, ln := range lines[1:] {
		f := strings.Split(ln, "\t")
		get := func(c int) string {
			if c < 0 || c >= len(f) {
				return ""
			}
			return strings.ToUpper(strings.TrimSpace(f[c]))
		}
		if segment != "" && cs >= 0 && !strings.EqualFold(get(cs), segment) {
			continue
		}
		name := fmt.Sprintf("set%d", i+1)
		if cn >= 0 && cn < len(f) {
			name = strings.TrimSpace(f[cn])
		}
		if get(cf) == "" || get(cr) == "" || get(cp) == "" {
			continue
		}
		out = append(out, design.Seed{Name: name, F: get(cf), R: get(cr), P: get(cp)})
	}
	return out, nil
}

func designShift(d *designFlags) error {
	al, err := loadAligned(d.aln)
	if err != nil {
		return err
	}
	p, err := d.params()
	if err != nil {
		return err
	}
	if d.sets == "" {
		return fmt.Errorf("--sets sets.tsv is required")
	}
	seg := ""
	if d.fs.Lookup("segment") != nil {
		seg = d.fs.Lookup("segment").Value.String()
	}
	seeds, err := readSeeds(d.sets, seg)
	if err != nil {
		return err
	}
	var xs *design.Exclusion
	if d.exclude != "" {
		if xs, err = design.LoadExclusion(d.exclude, p); err != nil {
			return err
		}
	}
	opt := design.ShiftOptions{Shift: d.shift, Step: d.step, Keep: d.keep, PrimerAmb: d.ambig, ProbeAmb: max(d.probeAmbig, 3), DTm: 3}
	if d.maxAmb >= 0 {
		opt.PrimerAmb, opt.ProbeAmb = d.maxAmb, d.maxAmb
	}
	w, done, err := outWriter(d.out)
	if err != nil {
		return err
	}
	fmt.Fprintln(w, strings.Join(design.SetColumns, "\t"))
	for _, s := range seeds {
		t0 := time.Now()
		res, base, err := design.Shift(al, s, p, opt)
		if err != nil {
			fmt.Fprintln(os.Stderr, err)
			continue
		}
		fmt.Fprintf(os.Stderr, "%s: now %d/%d (%.1f%%); %d improved classes (%.1fs)\n", s.Name, base.Cov, al.G, 100*float64(base.Cov)/float64(al.G), len(res), time.Since(t0).Seconds())
		base.Change = "-"
		base.BaseCov = base.Cov
		all := append([]design.Set{base}, res...)
		for i := range all {
			x := &all[i]
			if x.Change != "-" {
				design.PreferCStrand(x)
			}
			if d.qc {
				design.QC(x, p)
			}
			if xs != nil {
				xs.Check(x)
			}
			fmt.Fprintln(w, strings.Join(x.Row(al.G), "\t"))
		}
	}
	return done()
}

func designTile(d *designFlags) error {
	al, err := loadAligned(d.aln)
	if err != nil {
		return err
	}
	p, err := d.params()
	if err != nil {
		return err
	}
	t0 := time.Now()
	rows := design.Tile(al, p, d.step, d.keep)
	fmt.Fprintf(os.Stderr, "%d genomes, %d windows (step %d) in %.1fs\n", al.G, len(rows), d.step, time.Since(t0).Seconds())
	out := d.out
	if out == "" {
		out = "tile.tsv"
	}
	w, done, err := outWriter(out)
	if err != nil {
		return err
	}
	fmt.Fprintln(w, "window_start\twindow_end\tcoverage\tcoverage_pct\tf_start\tr_end\tamplicon\tF\tP\tR")
	for _, r := range rows {
		if !r.Found {
			fmt.Fprintf(w, "%d\t%d\t0\t0.00\t\t\t\t\t\t\n", r.Start, r.Start+p.Amp[1])
			continue
		}
		s := r.Set
		fmt.Fprintf(w, "%d\t%d\t%d\t%.2f\t%d\t%d\t%d\t%s\t%s\t%s\n", r.Start, r.Start+p.Amp[1], s.Cov, 100*float64(s.Cov)/float64(al.G), s.FS, s.RE, s.Amplicon, s.F, s.P, s.R)
	}
	if err := done(); err != nil {
		return err
	}
	var xs *design.Exclusion
	if d.exclude != "" {
		if xs, err = design.LoadExclusion(d.exclude, p); err != nil {
			return err
		}
	}
	panel := design.GreedyPanel(rows, d.top, al.G)
	pw, pdone, err := outWriter(strings.TrimSuffix(out, ".tsv") + ".panel.tsv")
	if err != nil {
		return err
	}
	if err := writeSets(pw, al, panel, p, xs, d.qc); err != nil {
		return err
	}
	return pdone()
}

func readTSV(path string) ([]map[string]string, error) {
	b, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	lines := strings.Split(strings.TrimRight(string(b), "\n"), "\n")
	hdr := strings.Split(lines[0], "\t")
	var out []map[string]string
	for _, ln := range lines[1:] {
		f := strings.Split(ln, "\t")
		m := map[string]string{}
		for i, h := range hdr {
			if i < len(f) {
				m[strings.TrimSpace(h)] = strings.TrimSpace(f[i])
			}
		}
		out = append(out, m)
	}
	return out, nil
}

func designLabels(d *designFlags) error {
	if d.sets == "" {
		return fmt.Errorf("--sets primers.tsv is required")
	}
	rows, err := readTSV(d.sets)
	if err != nil {
		return err
	}
	label := map[string]string{}
	name := map[string]string{}
	bad := 0
	for _, r := range rows {
		id := r["primer_id"]
		if id == "" {
			id = r["name"]
		}
		l := r["fluorophore"]
		label[id], name[id] = l, r["primer_name"]
		if l == "" {
			continue
		}
		for _, p := range design.LabelProblems(l) {
			fmt.Printf("%s %s: %s\n", id, r["primer_name"], p)
			bad++
		}
	}
	if d.groups != "" {
		gs, err := readTSV(d.groups)
		if err != nil {
			return err
		}
		for _, g := range gs {
			m := map[string]string{}
			for _, id := range strings.Split(g["member_ids"], ";") {
				if l := label[strings.TrimSpace(id)]; l != "" {
					m[id] = l
				}
			}
			for _, p := range design.MultiplexProblems(m) {
				if strings.Contains(p, "share") {
					fmt.Printf("group %s %s (if run in one tube): %s\n", g["group_id"], g["label"], p)
					bad++
				}
			}
		}
	}
	fmt.Fprintf(os.Stderr, "%d label problems\n", bad)
	return nil
}
