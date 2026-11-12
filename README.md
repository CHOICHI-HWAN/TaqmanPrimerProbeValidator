# Pipeline list page (branch `entry`)

A small stand-alone front door for the per-pathogen primer/probe validation pipelines of this repository
(one branch per pathogen: `main` = LASV, `nipah`, `ebov`, `sudv`, `bdbv`, `restv`, `tafv`, `marv`, `cchfv`, `mpxv`, `rvfv`, `flua`).

It shows every pipeline as a card with live status, a few numbers read from the project's own result files, and can start / stop a
pipeline's web app. Pure Python standard library.

```bash
./start.sh            # http://localhost:8760   (the dashboards' "전체 파이프라인 목록" link points here)
```

`pipelines.json` lists the pipelines: `path` is the folder of the pipeline **relative to the parent of this folder**
(check out each pathogen branch into that folder layout), `port` is the port of its web app.
No sequences, vendor names or credentials are in this branch.
