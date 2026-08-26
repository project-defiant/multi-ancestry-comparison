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
