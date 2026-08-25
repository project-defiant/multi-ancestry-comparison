# Multi-ancestry fine-mapping simulation example — design

## Goal

Produce a small, self-contained, reproducible example under `example/`
that:

1. Simulates individual-level genotype + phenotype data for two
   ancestries (EUR, AFR) sharing the exact same 50-variant locus.
2. Derives GWAS summary statistics (β, SE, p, z) and empirical LD
   matrices per ancestry from that simulated data.
3. Runs single-ancestry fine-mapping (`susieR::susie_rss`) separately
   on EUR and AFR, and joint multi-ancestry fine-mapping (`sushie`) on
   both together.
4. Demonstrates that both single-ancestry credible sets (CS) are
   longer than the credible set sushie produces after combining
   ancestries, because population-specific LD drags in different tag
   SNPs alone but disagrees across ancestries everywhere except the
   true causal SNP.
5. Renders two multi-panel PNG figures visualizing the locus, the LD
   structure, and the fine-mapping results (PIP + CS size) for each
   method.

Non-goals: this is not a production pipeline, not a general-purpose
simulator, and not testing sushie/susieR's own correctness — it is a
demonstration/testdata generator with parameters tuned so the effect
is visible.

## Locus & causal variant

- 50 SNPs, `chrom=1`, positions `1,000,000` to `1,098,000` in steps of
  2,000 bp (index `0..49`).
- One true causal SNP at index 24 (0-indexed), i.e. the 25th variant,
  position `1,048,000`.
- Both ancestries share identical SNP IDs/positions/alleles (full
  overlap by construction — no allele-flipping logic needed).

## Genotype simulation (per ancestry)

A block-factor liability model generates genuinely-PSD, ancestry-specific
LD without hand-specifying a correlation matrix:

- Each ancestry gets its own partition of the 50 SNP indices into
  disjoint blocks. Two haplotypes per individual; for SNP `j` in block
  `b` with within-block loading `ρ_b`:
  `liability = sqrt(ρ_b) · F_b + sqrt(1-ρ_b) · E_j`, `F_b ~ N(0,1)`
  shared per block per haplotype draw, `E_j ~ N(0,1)` idiosyncratic.
  Genotype dosage = sum of the two haplotypes' calls, each haplotype
  called by thresholding its liability at the SNP's per-ancestry MAF
  quantile.
- **EUR partition**: causal SNP (idx 24) sits in a wide block spanning
  indices **17–31** (15 SNPs), `ρ_within = 0.85` — mirrors the long LD
  tracts typical of the EUR bottleneck.
- **AFR partition**: causal SNP sits in a narrower, differently-bounded
  block spanning indices **20–28** (9 SNPs), `ρ_within = 0.85` —
  mirrors faster AFR LD decay.
- All other SNPs (outside each ancestry's causal block) are assigned
  to small filler blocks of 5 consecutive SNPs each, `ρ_within = 0.3`,
  just to avoid an unrealistic all-independent background. These never
  compete with the causal block for PIP mass.
- Per-SNP MAF drawn independently per ancestry, `Uniform(0.05, 0.5)`.
- Sample sizes: **EUR N=50,000**, **AFR N=20,000** (realistic
  imbalance).

The key property: the EUR-plausible tag-SNP set (17–31) and the
AFR-plausible tag-SNP set (20–28) overlap only near the causal SNP —
so each ancestry alone can't distinguish within its own block, but the
two blocks disagree with each other everywhere except the causal SNP.

## Phenotype & GWAS summary stats

- Phenotype = standardized causal-SNP dosage × true effect + Gaussian
  noise, scaled to hit a target local heritability. Starting point:
  `h2_local ≈ 0.002`, which given the sample sizes above targets a
  causal-SNP z-score of roughly 6 (AFR) to 10 (EUR) — strong but not
  so strong that LD-driven ambiguity disappears.
- **This parameter (and the block `ρ_within` values above) are
  starting points, not fixed constants.** During implementation, after
  the first simulate → fine-map pass, check the realized CS sizes; if
  the qualitative pattern (EUR CS long, AFR CS long, sushie CS short)
  doesn't hold, adjust `h2_local` / `ρ_within` / block widths and
  re-run. The random seed stays fixed once the pattern holds, so the
  final example is reproducible.
- Marginal GWAS regression (β, SE, p-value, z) computed per SNP per
  ancestry via simple OLS, same approach as sushie's own
  `data/make_example.py`.
- Empirical in-sample LD computed as the Pearson correlation matrix of
  standardized dosages, per ancestry (matches sushie's own
  `_compute_ld`).
- Fixed random seed for full reproducibility.

## Output testdata files (`example/data/`)

- `EUR.gwas.tsv`, `AFR.gwas.tsv` — columns
  `chrom, snp, pos, a1, a0, beta, se, pval, z` (matches sushie's
  default `--gwas-header ['chrom','snp','pos','a1','a0','z']`, extra
  columns ignored by sushie/susieR but useful for the plots).
- `EUR.ld.tsv`, `AFR.ld.tsv` — tab-separated correlation matrix, SNP
  IDs as the header row, same SNP order as the GWAS file (matches
  sushie's LD file format and is trivially loadable in R).

## Fine-mapping execution

Both tools are installed from their **official released packages**
(CRAN / PyPI), pinned to an exact version — not from the local
`../susier` / `../sushie` forks in the parent directory.

- **susieR (single-ancestry, ×2)**: install from CRAN, pinned to the
  current release **0.14.2** via
  `remotes::install_version("susieR", version = "0.14.2")` (falls back
  to `install.packages("susieR")` + recorded `packageVersion()` if
  pinning is unavailable). A small R script
  `example/scripts/finemap_susieR.R` takes a GWAS tsv + LD tsv +
  sample size, calls `susieR::susie_rss(z=z, R=R, n=n, L=1)`, and
  writes a tsv of `snp, pos, z, pip, cs_id` (`cs_id` is `NA`/0 for SNPs
  outside any credible set) plus the CS size.
- **sushie (joint)**: install from PyPI, pinned to the current release
  **`sushie==0.20`**, as a regular dependency of `example/`'s own
  `uv`-managed venv (`uv add sushie==0.20`; requires Python ≥3.11,
  satisfied by the system's 3.12). Confirmed the released package's
  console entry point is `sushie` (pointing at the classic
  `sushie.cli:run_cli`, the same summary-stats interface described
  below) — not `sushie-legacy`, which only exists in the fork. Run:
  `uv run sushie finemap --summary --gwas EUR.gwas.tsv AFR.gwas.tsv --ld EUR.ld.tsv AFR.ld.tsv --sample-size 50000 20000 --L 1 --trait locus1 --output <path>`,
  reading back its `.cs.tsv` / PIP output.
- Both use `L=1` (matches the single simulated causal variant and
  keeps credible-set output unambiguous).

## Reproducibility

- **Fixed random seed** for the entire simulation (genotypes,
  phenotype noise, block-factor draws) — once the qualitative CS
  pattern holds, the seed is frozen and never changed again.
- **Pinned, released tool versions** rather than local forks:
  `sushie==0.20` (PyPI) and `susieR` `0.14.2` (CRAN), both recorded
  explicitly (see below) so re-running months later reproduces the
  same fine-mapping behavior even if newer releases change defaults.
- `example/pyproject.toml` + committed `uv.lock` pin the full Python
  dependency graph (Python ≥3.11, `sushie==0.20`, numpy/scipy/
  matplotlib for simulation and plotting).
- The R side has no project-level lockfile manager available here, so
  `example/scripts/finemap_susieR.R` pins the package version at
  install time (`remotes::install_version`) and the orchestrator
  writes the resolved `R.version.string` + `packageVersion("susieR")`
  into `example/results/versions.txt` for the record.
- `example/README.md` documents the exact commands to reproduce every
  artifact from scratch (simulate → fine-map → plot), so the whole
  pipeline is a handful of copy-pasteable commands.

## Plots (`example/results/`, matplotlib, saved as PNG)

**Figure 1 — locus & LD** (2×2 grid, one tile unused):
- upper-left: GWAS locus-zoom, EUR and AFR overlaid (position vs.
  −log10 p or |z|, distinct colors/markers per ancestry).
- upper-right: AFR LD heatmap.
- lower-left: EUR LD heatmap.
- lower-right: blank.

**Figure 2 — fine-mapping results** (2×2 grid, one tile unused):
- upper-right: AFR susieR result — locus-zoom of PIP vs. position,
  points colored/shaped by CS membership, CS size annotated in the
  title.
- lower-left: EUR susieR result — same, for EUR.
- lower-right: sushie joint result — same, for the merged run.
- upper-left: blank.

## File layout

```
example/
  pyproject.toml            # pins sushie==0.20, python>=3.11, numpy/scipy/matplotlib
  uv.lock                    # committed, full resolved dependency graph
  data/
    EUR.gwas.tsv
    AFR.gwas.tsv
    EUR.ld.tsv
    AFR.ld.tsv
  scripts/
    simulate.py            # genotype/phenotype/GWAS/LD simulation (fixed seed)
    finemap_susieR.R        # single-ancestry susie_rss runner (pinned susieR 0.14.2)
    run_finemapping.sh       # orchestrates susieR (x2) + sushie, writes versions.txt
    plot_results.py         # builds both figures
  results/
    EUR.susieR.cs.tsv
    AFR.susieR.cs.tsv
    locus1.sushie.cs.tsv (+ whatever else `sushie finemap` emits)
    versions.txt              # R/Python/susieR/sushie versions actually used
    fig1_locus_ld.png
    fig2_finemapping.png
  README.md                  # how to reproduce, one command each step
```

## Testing / verification

- No unit tests — this is a demonstration script, not library code.
- Verification is the demonstration itself: after running the full
  pipeline, assert (in `run_finemapping.sh` or by eye) that
  `len(CS_EUR) > len(CS_sushie)` and `len(CS_AFR) > len(CS_sushie)`,
  and that the true causal SNP is in all three credible sets. If this
  fails, that's the signal to go back and retune the simulation
  parameters (see "Phenotype & GWAS summary stats" above).
