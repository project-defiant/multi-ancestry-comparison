from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import norm

SEED = 20260826
N_VARIANTS = 50
CHROM = 1
START_POS = 1_000_000
STEP = 2_000
CAUSAL_IDX = 24  # 0-indexed -> the 25th variant, "snp25"

SNP_IDS = [f"snp{i + 1}" for i in range(N_VARIANTS)]
POSITIONS = START_POS + np.arange(N_VARIANTS) * STEP

N_EUR = 50_000
N_AFR = 20_000

H2_EUR = 0.002
H2_AFR = 0.002

EUR_CAUSAL_BLOCK = range(17, 32)  # 15 SNPs: 17..31 inclusive
AFR_CAUSAL_BLOCK = range(20, 29)  # 9 SNPs: 20..28 inclusive
RHO_CAUSAL_BLOCK = 0.85
RHO_FILLER_BLOCK = 0.3
FILLER_BLOCK_SIZE = 5

DATA_DIR = Path(__file__).resolve().parent.parent / "data"


def build_blocks(n_variants, causal_block):
    """Partition variant indices into the causal block plus filler blocks of 5."""
    causal_set = set(causal_block)
    blocks = [sorted(causal_set)]
    remaining = [i for i in range(n_variants) if i not in causal_set]
    for start in range(0, len(remaining), FILLER_BLOCK_SIZE):
        blocks.append(remaining[start:start + FILLER_BLOCK_SIZE])
    return blocks


def simulate_genotypes(n_individuals, blocks, block_rhos, mafs, rng):
    """Block-factor liability model: 2 haplotypes per individual, each SNP's
    liability is a weighted mix of its block's shared factor and idiosyncratic
    noise, thresholded at the SNP's MAF quantile."""
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


def simulate_phenotype(genotype, causal_idx, h2, rng):
    x = (genotype - genotype.mean(axis=0)) / genotype.std(axis=0)
    g = x[:, causal_idx]
    s2g = np.var(g)
    s2e = s2g * (1 / h2 - 1)
    noise = rng.standard_normal(genotype.shape[0]) * np.sqrt(s2e)
    return g + noise


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


def simulate_ancestry(n_individuals, blocks, block_rhos, h2, rng):
    mafs = rng.uniform(0.05, 0.5, size=N_VARIANTS)
    genotype = simulate_genotypes(n_individuals, blocks, block_rhos, mafs, rng)
    y = simulate_phenotype(genotype, CAUSAL_IDX, h2, rng)
    betas, ses, pvals, zs = marginal_regression(genotype, y)
    ld = compute_ld(genotype)
    gwas_df = pd.DataFrame({
        "chrom": CHROM,
        "snp": SNP_IDS,
        "pos": POSITIONS,
        "a1": "A",
        "a0": "G",
        "beta": betas,
        "se": ses,
        "pval": pvals,
        "z": zs,
    })
    return gwas_df, ld


def write_gwas(df, path):
    df.to_csv(path, sep="\t", index=False)


def write_ld(ld, snp_ids, path):
    pd.DataFrame(ld, columns=snp_ids).round(6).to_csv(path, sep="\t", index=False)


def main():
    DATA_DIR.mkdir(exist_ok=True)

    rng_eur = np.random.default_rng(SEED)
    rng_afr = np.random.default_rng(SEED + 1)

    eur_blocks = build_blocks(N_VARIANTS, EUR_CAUSAL_BLOCK)
    afr_blocks = build_blocks(N_VARIANTS, AFR_CAUSAL_BLOCK)
    eur_rhos = [RHO_CAUSAL_BLOCK] + [RHO_FILLER_BLOCK] * (len(eur_blocks) - 1)
    afr_rhos = [RHO_CAUSAL_BLOCK] + [RHO_FILLER_BLOCK] * (len(afr_blocks) - 1)

    eur_gwas, eur_ld = simulate_ancestry(N_EUR, eur_blocks, eur_rhos, H2_EUR, rng_eur)
    afr_gwas, afr_ld = simulate_ancestry(N_AFR, afr_blocks, afr_rhos, H2_AFR, rng_afr)

    write_gwas(eur_gwas, DATA_DIR / "EUR.gwas.tsv")
    write_gwas(afr_gwas, DATA_DIR / "AFR.gwas.tsv")
    write_ld(eur_ld, SNP_IDS, DATA_DIR / "EUR.ld.tsv")
    write_ld(afr_ld, SNP_IDS, DATA_DIR / "AFR.ld.tsv")

    print(f"EUR causal SNP ({SNP_IDS[CAUSAL_IDX]}) z: {eur_gwas.loc[CAUSAL_IDX, 'z']:.2f}")
    print(f"AFR causal SNP ({SNP_IDS[CAUSAL_IDX]}) z: {afr_gwas.loc[CAUSAL_IDX, 'z']:.2f}")


if __name__ == "__main__":
    main()
