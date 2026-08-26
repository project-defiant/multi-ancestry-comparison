#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

mkdir -p results

echo "== simulate =="
uv run python scripts/simulate.py

echo "== susieR: EUR =="
Rscript scripts/finemap_susieR.R data/EUR.gwas.tsv data/EUR.ld.tsv 50000 EUR results/EUR.susieR.cs.tsv

echo "== susieR: AFR =="
Rscript scripts/finemap_susieR.R data/AFR.gwas.tsv data/AFR.ld.tsv 20000 AFR results/AFR.susieR.cs.tsv

echo "== sushie: EUR+AFR joint =="
uv run sushie finemap \
  --summary \
  --gwas data/EUR.gwas.tsv data/AFR.gwas.tsv \
  --ld data/EUR.ld.tsv data/AFR.ld.tsv \
  --sample-size 50000 20000 \
  --L 1 \
  --min-snps 50 \
  --trait locus1 \
  --output results/locus1

{
  echo "R: $(Rscript -e 'cat(R.version.string)')"
  echo "susieR: $(Rscript -e 'cat(as.character(packageVersion("susieR")))')"
  echo "python: $(uv run python --version)"
  echo "sushie: $(uv run python -c 'import importlib.metadata as m; print(m.version("sushie"))')"
} > results/versions.txt

echo "== Credible set size comparison =="
uv run python scripts/check_cs_sizes.py
