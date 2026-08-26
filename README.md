# fine-mapping-simulation-example

A small, reproducible demonstration that multi-ancestry fine-mapping
(`sushie`) can resolve a causal variant more precisely than single-ancestry
fine-mapping (`susieR`) alone, by exploiting the fact that different
ancestries have different linkage disequilibrium (LD) structure around
the same causal SNP.

The example: a simulated 50-variant locus shared by two ancestries (EUR,
AFR) with one true causal variant. Fine-mapping EUR or AFR alone yields a
credible set with several LD-tagging variants that can't be told apart;
fine-mapping both jointly with `sushie` collapses that credible set down
to essentially the one true causal SNP.

## Contents

- **[`example/`](example/)** — the runnable pipeline: simulation, fine-mapping
  (pinned `susieR` 0.14.2 from CRAN, `sushie` 0.20 from PyPI), and the
  result figure. See [`example/README.md`](example/README.md) for setup
  and reproduction instructions.
- **[`docs/superpowers/specs/`](docs/superpowers/specs/)** — the design
  doc explaining the simulation approach and why it demonstrates the
  intended effect.
- **[`docs/superpowers/plans/`](docs/superpowers/plans/)** — the
  task-by-task implementation plan.

## Quick start

```bash
cd example
uv sync
Rscript -e 'if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "https://cloud.r-project.org")'
Rscript -e 'remotes::install_version("susieR", version = "0.14.2", repos = "https://cloud.r-project.org")'
./scripts/run_finemapping.sh
uv run python scripts/plot_results.py
```

This produces `example/results/fig_finemapping_summary.png` — a single
2 x 3 figure with the GWAS locus and both ancestries' LD matrices on top,
and the PIP/credible-set comparison across AFR susieR, EUR susieR, and
joint sushie on the bottom.
