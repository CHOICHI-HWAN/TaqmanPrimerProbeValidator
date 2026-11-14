//! Primer / probe design from an alignment: coverage-first panels (a port of the panels' Python engine, same designs),
//! exhaustive shift / tile searches, non-target exclusion, probe formats and reporter / quencher checks.

pub mod align;
pub mod chem;
pub mod cli;
pub mod cover;
pub mod exclude;
pub mod oligo;
pub mod qc;
pub mod shift;

#[cfg(test)]
mod tests {
    use super::align::{py_set_order, read_aligned, Aligned};
    use super::chem::{label_problems, multiplex_problems, parse_chem};
    use super::cover::{design_from_alignment, Params};
    use super::oligo::tm_nn;
    use crate::json;

    #[test]
    fn set_order_matches_cpython() {
        // list(set(vals)) printed by CPython 3.12
        let cases: [(&[u8], &[u8]); 4] = [
            (&[84, 71, 65, 67], &[65, 67, 84, 71]),
            (&[71, 84], &[84, 71]),
            (&[78, 45, 65, 84, 71, 67], &[65, 67, 71, 45, 78, 84]),
            (&[84, 84, 71, 65, 82, 89, 75, 77, 83, 87, 67, 78, 45], &[65, 67, 71, 75, 77, 78, 45, 82, 83, 84, 87, 89]),
        ];
        for (inp, want) in cases {
            assert_eq!(py_set_order(inp), want.to_vec(), "{inp:?}");
        }
    }

    #[test]
    fn round1_matches_python() {
        use super::cover::round1;
        // round(x, 1) in CPython
        for (x, want) in [(56.25, 56.2), (56.35, 56.4), (0.05, 0.1), (2.675, 2.7), (-1.25, -1.2), (83.04580152671755, 83.0)] {
            assert_eq!(round1(x), want, "{x}");
        }
    }

    #[test]
    fn tm_nn_matches_python() {
        assert_eq!(tm_nn(b"ATGCGTACGTTAGCATGCAAT"), 55.13288431164182);
    }

    /// Designs on a synthetic alignment equal the ones scripts/redesign_coverage.py made.
    #[test]
    fn cover_matches_python() {
        let dir = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/");
        let (ids, rows) = read_aligned(&format!("{dir}design_synth.fa")).unwrap();
        let al = Aligned::new(ids, rows).unwrap();
        let exp = json::parse(&std::fs::read_to_string(format!("{dir}design_synth_expected.json")).unwrap()).unwrap();
        for (ci, e) in exp.as_array().iter().enumerate() {
            let mut p = Params { top: 4, workers: 2, ..Default::default() };
            let pr = e.get("params");
            let num = |k: &str| pr.get(k).as_f64();
            let arr = |k: &str| (pr.get(k).as_array()[0].as_f64(), pr.get(k).as_array()[1].as_f64());
            if !pr.as_object().is_empty() {
                p.ambig = num("ambig") as usize;
                p.probe_ambig = num("probe_ambig") as usize;
                p.pass_mm = num("pass_mm") as usize;
                p.tp_n = num("tp_n") as usize;
                let a = arr("amp");
                p.amp = (a.0 as usize, a.1 as usize);
                p.primer_tm = arr("primer_tm");
                p.probe_above = arr("probe_above");
            }
            let run = design_from_alignment(&al, &p);
            let want = e.get("designs").as_array();
            assert_eq!(run.designs.len(), want.len(), "case {ci}");
            for (d, w) in run.designs.iter().zip(want) {
                assert_eq!(String::from_utf8_lossy(&d.f), w.get("F").as_str(), "case {ci}");
                assert_eq!(String::from_utf8_lossy(&d.p), w.get("P").as_str(), "case {ci}");
                assert_eq!(String::from_utf8_lossy(&d.r), w.get("R").as_str(), "case {ci}");
                assert_eq!(d.est_cov as f64, w.get("est_cov").as_f64());
                assert_eq!(d.est_added as f64, w.get("est_added").as_f64());
                assert_eq!(d.amplicon_len as f64, w.get("amplicon_len").as_f64());
                assert_eq!((d.f_tm, d.p_tm, d.r_tm), (w.get("F_tm").as_f64(), w.get("P_tm").as_f64(), w.get("R_tm").as_f64()));
            }
        }
    }

    #[test]
    fn chem_and_labels() {
        let c = parse_chem("mgb", f64::NAN, 0, 0).unwrap();
        assert_eq!((c.max_len, c.offset), (22, 10.0));
        assert!(!label_problems("FAM/BHQ2").is_empty());
        assert!(label_problems("FAM-ZEN-IBFQ").is_empty());
        assert!(label_problems("6FAM / MGBNFQ").is_empty());
        let m = vec![("S".to_string(), "FAM/BHQ1".to_string()), ("L".to_string(), "FAM/BHQ1".to_string())];
        assert!(!multiplex_problems(&m).is_empty());
    }
}
