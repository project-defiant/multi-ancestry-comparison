from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
DATA_DIR = Path(__file__).resolve().parent.parent / "data"


def load_gwas(path):
    return pd.read_csv(path, sep="\t")


def load_ld(path):
    return pd.read_csv(path, sep="\t").values


def plot_locus_and_ld():
    eur_gwas = load_gwas(DATA_DIR / "EUR.gwas.tsv")
    afr_gwas = load_gwas(DATA_DIR / "AFR.gwas.tsv")
    eur_ld = load_ld(DATA_DIR / "EUR.ld.tsv")
    afr_ld = load_ld(DATA_DIR / "AFR.ld.tsv")

    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    ax_locus, ax_afr_ld, ax_eur_ld, ax_blank = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

    ax_locus.scatter(eur_gwas["pos"], -np.log10(eur_gwas["pval"]), label="EUR", color="#4C72B0")
    ax_locus.scatter(afr_gwas["pos"], -np.log10(afr_gwas["pval"]), label="AFR", color="#DD8452")
    ax_locus.set_xlabel("Position (bp)")
    ax_locus.set_ylabel("-log10(p)")
    ax_locus.set_title("GWAS locus zoom")
    ax_locus.legend()

    im_afr = ax_afr_ld.imshow(afr_ld, cmap="RdBu_r", vmin=-1, vmax=1)
    ax_afr_ld.set_title("AFR LD matrix")
    fig.colorbar(im_afr, ax=ax_afr_ld, fraction=0.046, pad=0.04)

    im_eur = ax_eur_ld.imshow(eur_ld, cmap="RdBu_r", vmin=-1, vmax=1)
    ax_eur_ld.set_title("EUR LD matrix")
    fig.colorbar(im_eur, ax=ax_eur_ld, fraction=0.046, pad=0.04)

    ax_blank.axis("off")

    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "fig1_locus_ld.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    plot_locus_and_ld()
