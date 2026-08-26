# Multi-ancestry fine-mapping simulation example

Simulates a 50-variant locus shared by two ancestries (EUR, AFR) with one
true causal variant, and shows that fine-mapping each ancestry alone
(`susieR`) yields a longer credible set than joint multi-ancestry
fine-mapping (`sushie`).

Design doc: `../docs/superpowers/specs/2026-08-26-fine-mapping-simulation-design.md`.

## Setup

```bash
cd example
uv sync
Rscript -e 'if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "https://cloud.r-project.org")'
Rscript -e 'remotes::install_version("susieR", version = "0.14.2", repos = "https://cloud.r-project.org")'
```

## Reproduce everything

```bash
cd example
./scripts/run_finemapping.sh
uv run python scripts/plot_results.py
```

This regenerates `data/*.tsv`, `results/*.cs.tsv` / `results/*.weights.tsv`,
`results/versions.txt`, and one combined figure:

- `results/fig_finemapping_summary.png` — 2 rows x 3 columns (plotnine,
  transparent background). Top row: GWAS locus zoom (EUR vs AFR) and both
  ancestries' LD matrices. Bottom row: PIP/credible-set panels for AFR
  susieR, EUR susieR, and joint sushie, each credible-set SNP labeled by
  ID.

`scripts/check_cs_sizes.py` (run automatically at the end of
`run_finemapping.sh`) asserts that both single-ancestry credible sets are
longer than the joint sushie credible set, and that the true causal SNP
(`snp25`) is included in all three. With the committed parameters this
gives EUR CS size 3, AFR CS size 4, and sushie CS size 1.

All randomness is controlled by a single fixed seed in `scripts/simulate.py`
(`SEED = 20260826`), and the fine-mapping tools are pinned to
`sushie==0.20` (PyPI) and `susieR==0.14.2` (CRAN) — the released packages,
not local development forks.
