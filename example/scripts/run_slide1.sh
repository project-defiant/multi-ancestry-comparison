#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

mkdir -p results

echo "== simulate two-locus scan =="
uv run python scripts/simulate_locus_breaker.py

echo "== susieR: locus A =="
Rscript scripts/finemap_susieR.R data/locus_scan.locusA.gwas.tsv data/locus_scan.locusA.ld.tsv 50000 locusA results/locusA.susieR.cs.tsv

echo "== plot slide 1 figures =="
uv run python scripts/plot_slide1.py
