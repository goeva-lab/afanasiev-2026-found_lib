# publication-related scripts for `found` article

## [`./bench`](./bench): benchmarking infrastructure

- [`bench.py`](./bench/bench.py): benchmarking script, provides a CLI interface
- [`fetch.py`](./bench/fetch.py): utilities to fetch/download input datasets
- [`perf.py`](./bench/perf.py): script to aggregate/collect performance metrics following benchmarking

## [`./figs`](./figs): figure generation following benchmarking

- [`analysis.R`](./figs/analysis.R): analyzes benchmarking outputs and generates `ggplot2::ggplot` objects for figures 2, s1, s2, and s3
- [`il15_vignette.R`](./figs/il15_vignette.R): R script demonstrating R `found` package use, generates `ggplot2::ggplot` objects for figure 1e
- [`il15_vignette.py`](./figs/il15_vignette.py): python/jupytext script demonstrating python `found` package use, generates `altair.Chart` objects used in figure 1d

note: above scripts require presence of properly laid-out `./cache` (used as cache directory for [`bench.py`](./bench/bench.py) script, via `-c` flag) and `./.cache` (used to hold miscellaneous analysis & vignette data/output files) subdirectories.

## [`./manuscript`](./manuscript): manuscript text related documents

- [`article.typ`](./manuscript/article.typ): [typst](https://github.com/typst/typst?tab=readme-ov-file#installation) file for manuscript text \
  note: building [`article.typ`](./manuscript/article.typ) assumes the presence of a `./manuscript/figs` subdirectory containing the following files: `fig_1.pdf`, `fig_2.pdf`, `fig_s1.svg`, `fig_s2.svg`, `fig_s3.svg` (these are not currently present in this repository)
- [`citations.yml`](./manuscript/citations.yml): [hayagriva](https://github.com/typst/hayagriva/blob/main/docs/file-format.md) file for manuscript citations
