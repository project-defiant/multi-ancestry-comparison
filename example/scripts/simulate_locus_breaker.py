"""Two-signal toy GWAS scan for the slide-1 narrative (marginal scan -> locus
definition -> fine-mapping). Independent of simulate.py: this is a
presentation aid, not part of the EUR/AFR/sushie pipeline, and touches none
of its pinned seed/rho/heritability invariants.

Locus A is a tight-LD causal block (same block-factor design and rho as the
main pipeline's EUR/AFR loci) -- this is the locus carried through all three
slide-1 figures and the one susieR fine-maps. Locus B is a second,
genome-wide-significant signal placed far enough away that a distance-based
locus-breaker genuinely has something to split from Locus A; it is never
fine-mapped.

Writes:
    data/locus_scan.gwas.tsv       both loci, one combined marginal scan
    data/locus_scan.locusA.gwas.tsv  Locus A only, same columns, for susieR
    data/locus_scan.locusA.ld.tsv    Locus A only, SNP-ordered correlation matrix
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import norm

SEED = 20260907
CHROM = 1
STEP = 2_000
N = 50_000

LOCUS_A_START = 5_000_000
LOCUS_A_N_VARIANTS = 30
LOCUS_A_CAUSAL_BLOCK = range(10, 25)  # 15 SNPs, tight LD, mirrors EUR/AFR design
LOCUS_A_CAUSAL_IDX = 17  # 0-indexed within locus A -> "locusA_snp18"
RHO_A_CAUSAL = 0.999
RHO_A_FILLER = 0.3
H2_A = 0.002

GAP = 150_000  # exceeds the 100kb distance_cutoff used in the slide-1 figures
LOCUS_B_N_VARIANTS = 15
LOCUS_B_CAUSAL_BLOCK = range(5, 10)  # 5 SNPs, looser LD -- a simpler second signal
LOCUS_B_CAUSAL_IDX = 7
RHO_B_CAUSAL = 0.6
RHO_B_FILLER = 0.2
H2_B = 0.0015

FILLER_BLOCK_SIZE = 5

DATA_DIR = Path(__file__).resolve().parent.parent / "data"


def build_blocks(n_variants, causal_block):
    causal_set = set(causal_block)
    blocks = [sorted(causal_set)]
    remaining = [i for i in range(n_variants) if i not in causal_set]
    for start in range(0, len(remaining), FILLER_BLOCK_SIZE):
        blocks.append(remaining[start:start + FILLER_BLOCK_SIZE])
    return blocks


def simulate_genotypes(n_individuals, blocks, block_rhos, mafs, rng):
    n_variants = sum(len(b) for b in blocks)
    dosage = np.zeros((n_individuals, n_variants))
    thresholds = norm.ppf(1 - mafs)
    for _ in range(2):
        liability = np.empty((n_individuals, n_variants))
        for block, rho in zip(blocks, block_rhos):
            if not block:
                continue
            factor = rng.standard_normal(n_individuals)
            idio = rng.standard_normal((n_individuals, len(block)))
            liability[:, block] = np.sqrt(rho) * factor[:, None] + np.sqrt(1 - rho) * idio
        dosage += (liability > thresholds).astype(float)
    return dosage


def simulate_locus_genotype(n_variants, causal_block, rho_causal, rho_filler, rng):
    blocks = build_blocks(n_variants, causal_block)
    rhos = [rho_causal] + [rho_filler] * (len(blocks) - 1)
    mafs = rng.uniform(0.05, 0.5, size=n_variants)
    return simulate_genotypes(N, blocks, rhos, mafs, rng)


def marginal_regression(genotype, y):
    n_variants = genotype.shape[1]
    betas = np.empty(n_variants)
    ses = np.empty(n_variants)
    pvals = np.empty(n_variants)
    zs = np.empty(n_variants)
    for j in range(n_variants):
        res = stats.linregress(genotype[:, j], y)
        betas[j] = res.slope
        ses[j] = res.stderr
        pvals[j] = res.pvalue
        zs[j] = res.slope / res.stderr
    return betas, ses, pvals, zs


def compute_ld(genotype):
    x = (genotype - genotype.mean(axis=0)) / genotype.std(axis=0)
    return (x.T @ x) / x.shape[0]


def main():
    DATA_DIR.mkdir(exist_ok=True)
    rng = np.random.default_rng(SEED)

    geno_a = simulate_locus_genotype(
        LOCUS_A_N_VARIANTS, LOCUS_A_CAUSAL_BLOCK, RHO_A_CAUSAL, RHO_A_FILLER, rng)
    geno_b = simulate_locus_genotype(
        LOCUS_B_N_VARIANTS, LOCUS_B_CAUSAL_BLOCK, RHO_B_CAUSAL, RHO_B_FILLER, rng)

    def standardize_col(g, idx):
        col = g[:, idx]
        return (col - col.mean()) / col.std()

    g_a = standardize_col(geno_a, LOCUS_A_CAUSAL_IDX)
    g_b = standardize_col(geno_b, LOCUS_B_CAUSAL_IDX)
    noise = rng.standard_normal(N) * np.sqrt(max(1e-6, 1 - H2_A - H2_B))
    y = np.sqrt(H2_A) * g_a + np.sqrt(H2_B) * g_b + noise

    beta_a, se_a, pval_a, z_a = marginal_regression(geno_a, y)
    beta_b, se_b, pval_b, z_b = marginal_regression(geno_b, y)

    ids_a = [f"locusA_snp{i + 1}" for i in range(LOCUS_A_N_VARIANTS)]
    ids_b = [f"locusB_snp{i + 1}" for i in range(LOCUS_B_N_VARIANTS)]
    pos_a = LOCUS_A_START + np.arange(LOCUS_A_N_VARIANTS) * STEP
    locus_b_start = pos_a[-1] + GAP
    pos_b = locus_b_start + np.arange(LOCUS_B_N_VARIANTS) * STEP

    df_a = pd.DataFrame({
        "chrom": CHROM, "snp": ids_a, "pos": pos_a, "a1": "A", "a0": "G",
        "beta": beta_a, "se": se_a, "pval": pval_a, "z": z_a, "locus": "A",
    })
    df_b = pd.DataFrame({
        "chrom": CHROM, "snp": ids_b, "pos": pos_b, "a1": "A", "a0": "G",
        "beta": beta_b, "se": se_b, "pval": pval_b, "z": z_b, "locus": "B",
    })
    combined = pd.concat([df_a, df_b], ignore_index=True)
    combined.to_csv(DATA_DIR / "locus_scan.gwas.tsv", sep="\t", index=False)

    df_a.drop(columns="locus").to_csv(DATA_DIR / "locus_scan.locusA.gwas.tsv", sep="\t", index=False)
    ld_a = compute_ld(geno_a)
    pd.DataFrame(ld_a, columns=ids_a).round(6).to_csv(
        DATA_DIR / "locus_scan.locusA.ld.tsv", sep="\t", index=False)

    print(f"gap between locus A and locus B: {locus_b_start - pos_a[-1]:,} bp")
    print(f"locus A causal ({ids_a[LOCUS_A_CAUSAL_IDX]}) z={z_a[LOCUS_A_CAUSAL_IDX]:.2f} "
          f"p={pval_a[LOCUS_A_CAUSAL_IDX]:.2e}")
    print(f"locus B causal ({ids_b[LOCUS_B_CAUSAL_IDX]}) z={z_b[LOCUS_B_CAUSAL_IDX]:.2f} "
          f"p={pval_b[LOCUS_B_CAUSAL_IDX]:.2e}")
    print(f"locus A min p={pval_a.min():.2e}, snps at p<5e-8: {(pval_a < 5e-8).sum()}")
    print(f"locus B min p={pval_b.min():.2e}, snps at p<5e-8: {(pval_b < 5e-8).sum()}")


if __name__ == "__main__":
    main()
