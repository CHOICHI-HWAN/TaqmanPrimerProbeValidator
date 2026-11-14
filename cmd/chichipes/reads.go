package main

import (
	"flag"
	"fmt"
	"os"
	"runtime"
	"strings"
	"time"

	"chichipes/internal/primers"
	"chichipes/internal/reads"
)

type fileList []string

func (f *fileList) String() string { return strings.Join(*f, ",") }
func (f *fileList) Set(v string) error {
	for _, p := range strings.Split(v, ",") {
		if p = strings.TrimSpace(p); p != "" {
			*f = append(*f, p)
		}
	}
	return nil
}

const readsUsage = `usage:
  chichipes reads -i oligos.fa -r reads.fastq.gz [-r more.fq.gz] -o prefix [-m 2] [--indel] [--cpu N] [--hits]
  chichipes reads -i oligos.fa -1 R1.fastq.gz -2 R2.fastq.gz -o prefix   (paired-end: a pair counts once)

Checks which reads carry each oligo (FASTQ or FASTA, gzip or not). Name oligos SET_F / SET_R / SET_P to also get
per-set counts (reads with F and R, with F, R and P). -m is mismatches per oligo (substitutions; with --indel, edits,
for error-prone long reads). Writes <prefix>.summary.tsv, .sets.tsv, .positions.tsv (bases the reads have at every
oligo position), .txt and, with --hits, .hits.tsv (one line per read and oligo found).
`

// runReads: chichipes reads ...
func runReads(args []string) {
	fs := flag.NewFlagSet("reads", flag.ExitOnError)
	fs.Usage = func() { fmt.Fprint(os.Stderr, readsUsage); fs.PrintDefaults() }
	in := fs.String("i", "", "oligo FASTA (primers / probes)")
	var se, r1, r2 fileList
	fs.Var(&se, "r", "single-end or long-read file(s), repeatable or comma-separated")
	fs.Var(&r1, "1", "paired-end R1 file(s)")
	fs.Var(&r2, "2", "paired-end R2 file(s), same order as -1")
	out := fs.String("o", "", "output prefix")
	mm := fs.Int("m", 2, "max mismatches (or edits with --indel) per oligo")
	indel := fs.Bool("indel", false, "allow insertions/deletions (ONT, PacBio CLR)")
	seed := fs.Int("seed", 0, "seed length (0 = automatic)")
	cpu := fs.Int("cpu", runtime.NumCPU(), "threads")
	hits := fs.Bool("hits", false, "also write <prefix>.hits.tsv")
	check(fs.Parse(args))
	paired := len(r1) > 0
	if *in == "" || *out == "" || (len(se) == 0) == !paired || (paired && len(r2) == 0) {
		fs.Usage()
		os.Exit(2)
	}
	t0 := time.Now()
	ol, err := primers.ReadFasta(*in)
	check(err)
	oligos, err := reads.FromPrimers(ol)
	check(err)
	s, err := reads.NewSearcher(oligos, *mm, *seed, *indel)
	check(err)
	var hf *os.File
	if *hits {
		hf, err = os.Create(*out + ".hits.tsv")
		check(err)
		defer hf.Close()
	}
	var st *reads.Stats
	if paired {
		st, err = runWith(s, r1, r2, *cpu, hf)
	} else {
		st, err = runWith(s, se, nil, *cpu, hf)
	}
	check(err)
	rep, err := reads.WriteReports(*out, s, st, paired)
	check(err)
	fmt.Print(rep)
	fmt.Fprintf(os.Stderr, "%.1fs\n", time.Since(t0).Seconds())
}

func runWith(s *reads.Searcher, a, b []string, cpu int, hf *os.File) (*reads.Stats, error) {
	if hf != nil {
		return s.Run(a, b, cpu, hf)
	}
	return s.Run(a, b, cpu, nil)
}
