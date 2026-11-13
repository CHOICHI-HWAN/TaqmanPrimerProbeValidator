// chichipes: clean-room primer engine with an mfeprimer-compatible command line (the subset the pipeline uses).
//
//	chichipes -d db.fa -i primers.fa -o out.txt [--fasta] [--virus] [--misMatch 2] [--cpu 8] [-t 30] [-T 100] [-s 0] [-S 2000]
//	chichipes index -i db.fa        (writes the .fai/.json the pipeline expects; the search needs no index)
//	chichipes dimer -i oligos.fa -o dimer.txt / chichipes hairpin -i oligos.fa -o hairpin.txt
//	chichipes version
package main

import (
	"flag"
	"fmt"
	"os"
	"runtime"
	"time"

	"chichipes/internal/primers"
	"chichipes/internal/seqio"
	"chichipes/internal/spec"
	"chichipes/internal/thermo"
)

// Version is set at release build time with -ldflags "-X main.Version=v0.1.0".
var Version = "dev"

var version = "chichipes " + Version + " (mfeprimer-4.2.4 compatible subset)"

func main() {
	if len(os.Args) > 1 {
		switch os.Args[1] {
		case "version", "-v", "--version", "-version":
			fmt.Println(version)
			return
		case "index":
			runIndex(os.Args[2:])
			return
		case "search":
			runSearch(os.Args[2:])
			return
		case "dimer":
			runDimer(os.Args[2:])
			return
		case "hairpin":
			runHairpin(os.Args[2:])
			return
		}
	}
	opt := spec.DefaultOptions()
	fs := flag.NewFlagSet("chichipes", flag.ExitOnError)
	db := fs.String("d", "", "database FASTA")
	in := fs.String("i", "", "primer FASTA")
	out := fs.String("o", "", "output report")
	fs.Bool("fasta", false, "write amplicons in FASTA (always written)")
	fs.Bool("virus", false, "virus coverage mode (always on)")
	fs.IntVar(&opt.MisMatch, "misMatch", 0, "max mismatches in the 3' window")
	fs.IntVar(&opt.MisStart, "misStart", 1, "")
	fs.IntVar(&opt.MisEnd, "misEnd", 9, "")
	fs.IntVar(&opt.K, "k", 9, "")
	fs.IntVar(&opt.CPU, "cpu", opt.CPU, "")
	fs.IntVar(&opt.CPU, "c", opt.CPU, "")
	fs.Float64Var(&opt.MinTm, "t", 30, "")
	fs.Float64Var(&opt.MaxTm, "T", 100, "")
	fs.IntVar(&opt.MinSize, "s", 0, "")
	fs.IntVar(&opt.MaxSize, "S", 2000, "")
	fs.Float64Var(&opt.Cond.Mono, "mono", 50, "")
	fs.Float64Var(&opt.Cond.Diva, "diva", 1.5, "")
	fs.Float64Var(&opt.Cond.DNTP, "dntp", 0.25, "")
	fs.Float64Var(&opt.Cond.Oligo, "oligo", 50, "")
	maxAmp := fs.Int("M", 10000, "max amplicons written to the FASTA")
	args := os.Args[1:]
	if len(args) > 0 && args[0] == "spec" {
		args = args[1:]
	}
	if err := fs.Parse(normalise(args)); err != nil {
		os.Exit(2)
	}
	if *db == "" || *in == "" || *out == "" {
		fmt.Fprintln(os.Stderr, "need -d, -i and -o")
		os.Exit(2)
	}
	t0 := time.Now()
	oligos, err := primers.ReadFasta(*in)
	check(err)
	e, err := spec.New(oligos, opt)
	check(err)
	recs, err := seqio.ReadAll(*db)
	check(err)
	res := e.Run(recs)
	check(spec.WriteOutputs(*out, *db, recs, res, *maxAmp))
	fmt.Fprintf(os.Stderr, "%d records, %d amplicons, %d missed, %.1fs\n", res.Total, len(res.Amplicons), len(res.Failed),
		time.Since(t0).Seconds())
}

// normalise turns GNU-style "--flag" into Go's "-flag" (both work with package flag already) and keeps the rest.
func normalise(a []string) []string { return a }

func check(err error) {
	if err != nil {
		fmt.Fprintln(os.Stderr, "error:", err)
		os.Exit(1)
	}
}

// runSearch: chichipes search db.fa query.fa out.txt [-m 2] [-k 9] [--cpu N]
func runSearch(args []string) {
	fs := flag.NewFlagSet("search", flag.ExitOnError)
	mm := fs.Int("m", 2, "max mismatches (edit distance)")
	fs.IntVar(mm, "mismatch", 2, "")
	k := fs.Int("k", 9, "")
	cpu := fs.Int("cpu", runtime.NumCPU(), "")
	c := thermo.Default
	fs.Float64Var(&c.Mono, "mono", 50, "")
	fs.Float64Var(&c.Diva, "diva", 1.5, "")
	fs.Float64Var(&c.DNTP, "dntp", 0.25, "")
	fs.Float64Var(&c.Oligo, "oligo", 50, "")
	var pos []string
	for len(args) > 0 { // positional arguments may come before or after the flags
		if err := fs.Parse(args); err != nil {
			os.Exit(2)
		}
		args = fs.Args()
		if len(args) > 0 {
			pos = append(pos, args[0])
			args = args[1:]
		}
	}
	if len(pos) != 3 {
		fmt.Fprintln(os.Stderr, "usage: chichipes search db.fa query.fa out.txt [-m N]")
		os.Exit(2)
	}
	t0 := time.Now()
	recs, err := seqio.ReadAll(pos[0])
	check(err)
	oligos, err := primers.ReadFasta(pos[1])
	check(err)
	hits, err := spec.Search(recs, oligos, *k, *mm, *cpu, c)
	check(err)
	check(spec.WriteSearch(pos[2], recs, hits, *mm))
	fmt.Fprintf(os.Stderr, "%d records, %d hits, %.1fs\n", len(recs), len(hits), time.Since(t0).Seconds())
}

// runDimer: chichipes dimer -i oligos.fa -o dimer.txt [-s 5] [-m 2] [-d 0] [--cpu N]
func runDimer(args []string) {
	fs := flag.NewFlagSet("dimer", flag.ExitOnError)
	in := fs.String("i", "", "oligo FASTA")
	out := fs.String("o", "", "report")
	score := fs.Int("s", 5, "score cutoff")
	mm := fs.Int("m", 2, "max mismatches")
	dg := fs.Float64("d", 0, "dG cutoff")
	fs.Int("cpu", 1, "")
	fs.Int("c", 1, "")
	fs.Bool("p", false, "only dimers with 3' end binding (not implemented)")
	c := thermo.Default
	fs.Float64Var(&c.Mono, "mono", 50, "")
	fs.Float64Var(&c.Diva, "diva", 1.5, "")
	fs.Float64Var(&c.DNTP, "dntp", 0.25, "")
	fs.Float64Var(&c.Oligo, "oligo", 50, "")
	check(fs.Parse(args))
	oligos, err := primers.ReadFasta(*in)
	check(err)
	ex, err := spec.ExpandAll(oligos)
	check(err)
	ds := spec.Dimers(ex, *score, *mm, *dg, c)
	check(spec.WriteDimers(*out, *in, ex, ds, *score, *mm, *dg, c))
}

// runHairpin: chichipes hairpin -i oligos.fa -o hairpin.txt [-s 5] [-d 0] [-t 0] [--mindh 2] [--minl 3] [--maxl 12] [--maxgap 2]
func runHairpin(args []string) {
	fs := flag.NewFlagSet("hairpin", flag.ExitOnError)
	in := fs.String("i", "", "oligo FASTA")
	out := fs.String("o", "", "report")
	o := spec.DefaultHairpinOptions
	fs.IntVar(&o.Score, "s", o.Score, "score cutoff")
	fs.Float64Var(&o.Dg, "d", o.Dg, "dG cutoff")
	fs.Float64Var(&o.Tm, "t", o.Tm, "Tm cutoff")
	fs.IntVar(&o.MinDH, "mindh", o.MinDH, "")
	fs.IntVar(&o.MinL, "minl", o.MinL, "")
	fs.IntVar(&o.MaxL, "maxl", o.MaxL, "")
	fs.IntVar(&o.MaxGap, "maxgap", o.MaxGap, "")
	fs.Int("cpu", 1, "")
	fs.Int("c", 1, "")
	c := thermo.Default
	fs.Float64Var(&c.Mono, "mono", 50, "")
	fs.Float64Var(&c.Diva, "diva", 1.5, "")
	fs.Float64Var(&c.DNTP, "dntp", 0.25, "")
	fs.Float64Var(&c.Oligo, "oligo", 50, "")
	check(fs.Parse(args))
	oligos, err := primers.ReadFasta(*in)
	check(err)
	ex, err := spec.ExpandAll(oligos)
	check(err)
	check(spec.WriteHairpins(*out, *in, ex, spec.Hairpins(ex, o, c), o, c))
}
