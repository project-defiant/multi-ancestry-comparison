# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

A single self-contained example (`example/`) that simulates a 50-variant
GWAS locus shared by two ancestries (EUR, AFR) with one true causal
variant, then demonstrates that single-ancestry fine-mapping (`susieR`)
produces a longer credible set than joint multi-ancestry fine-mapping
(`sushie`). The design rationale and full parameter/tuning history live in
`docs/superpowers/specs/2026-08-26-fine-mapping-simulation-design.md` and
the implementation plan in `docs/superpowers/plans/2026-08-26-fine-mapping-simulation.md`.

## Commands

All commands run from `example/`.

```bash
# one-time setup
uv sync
Rscript -e 'if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "https://cloud.r-project.org")'
Rscript -e 'remotes::install_version("susieR", version = "0.14.2", repos = "https://cloud.r-project.org")'

# regenerate everything (testdata, fine-mapping outputs, both figures)
./scripts/run_finemapping.sh
uv run python scripts/plot_results.py

# individual stages
uv run python scripts/simulate.py                                          # writes data/*.gwas.tsv, data/*.ld.tsv
Rscript scripts/finemap_susieR.R data/EUR.gwas.tsv data/EUR.ld.tsv 50000 EUR results/EUR.susieR.cs.tsv
uv run sushie finemap --summary --gwas data/EUR.gwas.tsv data/AFR.gwas.tsv \
  --ld data/EUR.ld.tsv data/AFR.ld.tsv --sample-size 50000 20000 --L 1 --min-snps 50 \
  --trait locus1 --output results/locus1
uv run python scripts/check_cs_sizes.py                                    # asserts the CS-size pattern
```

There is no unit test suite by design — this is a demonstration pipeline,
not a library. `scripts/check_cs_sizes.py` (run automatically at the end
of `run_finemapping.sh`) is the correctness check: it asserts both
single-ancestry credible sets are longer than the joint sushie credible
set, and that the true causal SNP (`snp25`) is in all three.

## Architecture

The pipeline is three isolated stages that communicate only through flat
TSV files under `example/data/` and `example/results/` — no code imports
across the Python/R boundary:

1. **`scripts/simulate.py`** — generates ancestry-specific LD via a
   block-factor liability model (`build_blocks` / `simulate_genotypes`):
   each ancestry gets its own, differently-bounded partition of the 50
   SNPs into correlated blocks, so the "plausible causal SNP" set differs
   by ancestry even though both share the same true causal variant. From
   the simulated genotypes it derives phenotypes at a target local
   heritability, runs per-SNP marginal regression, and writes
   `EUR.gwas.tsv` / `AFR.gwas.tsv` (columns: `chrom, snp, pos, a1, a0,
   beta, se, pval, z`, matching sushie's default `--gwas-header`) and
   `EUR.ld.tsv` / `AFR.ld.tsv` (tab-separated correlation matrix, SNP IDs
   as the header row).
2. **`scripts/finemap_susieR.R`** — single-ancestry fine-mapping via
   `susieR::susie_rss(z, R, n, L=1)`; run once per ancestry, writing
   `snp, pos, z, pip, cs_id` to a results tsv.
3. **`uv run sushie finemap --summary ...`** — joint fine-mapping via the
   pinned `sushie==0.20` CLI. Its output filenames get a `.sushie.`
   infix: `{output}.sushie.weights.tsv` (per-SNP PIP across all 50 SNPs,
   columns include `sushie_pip_all`, `sushie_cs_index`) and
   `{output}.sushie.cs.tsv` (only the SNPs inside a credible set).
4. **`scripts/plot_results.py`** — renders `fig1_locus_ld.png` (GWAS
   locus zoom + both LD heatmaps) and `fig2_finemapping.png` (PIP panels
   with credible-set size in the title, for EUR susieR, AFR susieR, and
   joint sushie).

`scripts/run_finemapping.sh` runs simulate → both susieR calls → sushie →
`check_cs_sizes.py`, and records the exact tool versions used in
`results/versions.txt`.

### Reproducibility and package provenance

Everything is keyed off one fixed seed (`SEED = 20260826` in
`simulate.py`); re-running the pipeline reproduces bit-identical outputs.

Fine-mapping deliberately uses the **released** packages — `sushie==0.20`
from PyPI (a normal dependency of `example/`'s own `uv` project) and
`susieR` `0.14.2` from CRAN (installed via `remotes::install_version` into
a personal R library) — never the `../sushie` / `../susier` checkouts that
happen to sit in the parent `OpenTargets/` directory, which are personal
forks, not the canonical packages.

Two non-obvious CLI details worth knowing before touching the sushie
invocation: `sushie finemap` defaults `--min-snps` to 100 and will reject
this 50-SNP locus unless it's explicitly overridden (`--min-snps 50`);
and its output paths get the `.sushie.` infix noted above rather than the
bare `{output}.cs.tsv` / `{output}.weights.tsv` the tool's own `--help`
text might suggest.

### Simulation tuning

The block-factor model's causal-block correlation (`RHO_CAUSAL_BLOCK` in
`simulate.py`) is tuned to `0.999` — much higher than a naive reading of
"strong LD" would suggest. At more moderate values (0.85-0.99) susieR
still pinned the PIP almost entirely to the true causal SNP even with
z-scores of 6-9 on flanking SNPs; only at 0.999 does each ancestry's
credible set genuinely include multiple indistinguishable tag SNPs while
sushie's joint result still collapses to size 1. If re-tuning, re-run
`./scripts/run_finemapping.sh` and check `check_cs_sizes.py`'s PASS/FAIL
lines rather than assuming a given rho/heritability value works.
