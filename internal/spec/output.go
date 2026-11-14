package spec

import (
	"bufio"
	"fmt"
	"os"
	"strings"

	"chichipes/internal/seqio"
	"chichipes/internal/thermo"
)

// WriteOutputs writes the files the pipeline reads, in mfeprimer's formats:
// <out> (report with the coverage line), <out>.spec.tsv, <out>.amp.fasta (first maxAmp products) and
// <out>.virus.Failed.txt.
func WriteOutputs(out, dbPath string, db []seqio.Record, res Result, maxAmp int) error {
	os.Remove(out + ".spec.tsv")
	specPath := out + ".spec.tsv"
	if len(res.Amplicons) == 0 { // mfeprimer writes no spec.tsv at all when nothing amplifies
		specPath = os.DevNull
	}
	spec, err := os.Create(specPath)
	if err != nil {
		return err
	}
	sw := bufio.NewWriterSize(spec, 1<<20)
	fmt.Fprintln(sw, "#1-based coordinate, should be idential to blat output if align primer sequences to genome")
	fmt.Fprintln(sw, "#name\tchrom\tampStart\tampEnd\tampGC\tampSize\tfpName\tfpStart\tfpEnd\tfpSeq\tfpTm\tfpGC\tfpDg\t"+
		"rpName\trpEnd\trpStart\trpSeq\trpTm\trpGC\trpDg\tnote")
	amp, err := os.Create(out + ".amp.fasta")
	if err != nil {
		return err
	}
	aw := bufio.NewWriterSize(amp, 1<<20)
	for i, a := range res.Amplicons {
		g := db[a.Rec].Seq
		fseq := a.F.Seq()
		rplus := a.R.P.rcs // the reverse primer as it reads on the plus strand
		var mid, prod []byte
		overlap := a.R.Pos < a.F.Pos+a.F.N
		if overlap { // the two primers overlap: no middle; GC from the genome region (docs/FINDINGS.md, known difference)
			prod = g[a.Start():a.End()]
		} else {
			mid = g[a.F.Pos+a.F.N : a.R.Pos]
			prod = make([]byte, 0, a.Size()) // the product as amplified: primer sequences at both ends
			prod = append(append(append(prod, fseq...), mid...), rplus...)
		}
		name := fmt.Sprintf("Amp_%d", i+1)
		fmt.Fprintf(sw, "%s\t%s\t%d\t%d\t%.2f\t%d\t%s\t%d\t%d\t%s\t%.2f\t%.2f\t%.2f\t%s\t%d\t%d\t%s\t%.2f\t%.2f\t%.2f\t\n",
			name, db[a.Rec].ID, a.Start()+1, a.End(), thermo.GC(prod), a.Size(),
			a.F.Name(), a.F.Pos+1, a.F.Pos+a.F.N, fseq, a.F.Tm, thermo.GC(fseq), a.F.Dg,
			a.R.Name(), a.R.Pos+1, a.R.Pos+a.R.N, rplus, a.R.Tm, thermo.GC(a.R.Seq()), a.R.Dg)
		if i < maxAmp && overlap {
			fmt.Fprintf(aw, ">%s %s+%s==>%s:%d-%d %s\n%s\n", name, a.F.Name(), a.R.Name(), db[a.Rec].ID,
				a.Start(), a.End(), dbPath, strings.ToUpper(string(prod)))
		} else if i < maxAmp {
			fmt.Fprintf(aw, ">%s %s+%s==>%s:%d-%d %s\n%s%s%s\n", name, a.F.Name(), a.R.Name(), db[a.Rec].ID,
				a.Start(), a.End(), dbPath, fseq, strings.ToLower(string(mid)), rplus)
		}
	}
	if err := sw.Flush(); err != nil {
		return err
	}
	spec.Close()
	if err := aw.Flush(); err != nil {
		return err
	}
	amp.Close()

	failed, err := os.Create(out + ".virus.Failed.txt")
	if err != nil {
		return err
	}
	fw := bufio.NewWriter(failed)
	for _, i := range res.Failed {
		fmt.Fprintf(fw, "%s \n", db[i].ID)
	}
	fw.Flush()
	failed.Close()

	if err := writeBorder(out+".borderline.tsv", db, res); err != nil {
		return err
	}

	rep, err := os.Create(out)
	if err != nil {
		return err
	}
	defer rep.Close()
	cov := 0.0
	if res.Total > 0 {
		cov = 100 * float64(res.Total-len(res.Failed)) / float64(res.Total)
	}
	fmt.Fprintf(rep, "chichipes report (mfeprimer-compatible outputs)\n\n")
	fmt.Fprintf(rep, "Primer Coverate Rate (virus mode, e.g., to check whether these primers can amplify all the variants of SARS-CoV-2 genomes)\n")
	fmt.Fprintf(rep, "Coverage %.2f%%, total records in database: %d, amplicon number: %d, missed number: %d.\n",
		cov, res.Total, len(res.Amplicons), len(res.Failed))
	return nil
}

// writeBorder lists the bindings whose Tm lies within BorderBand of the -t cut-off (engine-only side file): there the
// keep/drop decision may differ from mfeprimer's, because the fitted parameters match its Tm to about ±0.01 °C.
func writeBorder(path string, db []seqio.Record, res Result) error {
	fh, err := os.Create(path)
	if err != nil {
		return err
	}
	defer fh.Close()
	w := bufio.NewWriter(fh)
	fmt.Fprintf(w, "# bindings with Tm within %.2f C of the cut-off; decision may differ from mfeprimer\n", BorderBand)
	fmt.Fprintln(w, "#chrom\tstart\tend\tstrand\tname\tseq\ttm\tdg\tdecision")
	for i, bs := range res.Border {
		for _, b := range bs {
			strand, dec := "+", "dropped"
			if b.Rev {
				strand = "-"
			}
			if b.Kept {
				dec = "kept"
			}
			fmt.Fprintf(w, "%s\t%d\t%d\t%s\t%s\t%s\t%.3f\t%.3f\t%s\n", db[i].ID, b.Pos+1, b.Pos+b.N, strand, b.Name(),
				b.Seq(), b.Tm, b.Dg, dec)
		}
	}
	return w.Flush()
}
