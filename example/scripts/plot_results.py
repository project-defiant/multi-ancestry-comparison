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


def load_susieR_cs(path):
    return pd.read_csv(path, sep="\t")


def load_sushie_weights(path):
    df = pd.read_csv(path, sep="\t")
    df["pip"] = df["sushie_pip_all"]
    df["in_cs"] = df["sushie_cs_index"] != "No CS"
    return df


def plot_pip_panel(ax, pos, pip, in_cs, title):
    colors = np.where(in_cs, "#C44E52", "#8C8C8C")
    ax.scatter(pos, pip, c=colors)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("Position (bp)")
    ax.set_ylabel("PIP")
    cs_size = int(np.sum(in_cs))
    ax.set_title(f"{title} (CS size = {cs_size})")


def plot_finemapping_results():
    eur = load_susieR_cs(RESULTS_DIR / "EUR.susieR.cs.tsv")
    afr = load_susieR_cs(RESULTS_DIR / "AFR.susieR.cs.tsv")
    sushie = load_sushie_weights(RESULTS_DIR / "locus1.sushie.weights.tsv")

    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    ax_blank, ax_afr, ax_eur, ax_sushie = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

    ax_blank.axis("off")
    plot_pip_panel(ax_afr, afr["pos"], afr["pip"], afr["cs_id"] == 1, "AFR susieR")
    plot_pip_panel(ax_eur, eur["pos"], eur["pip"], eur["cs_id"] == 1, "EUR susieR")
    plot_pip_panel(ax_sushie, sushie["pos"], sushie["pip"], sushie["in_cs"], "sushie (EUR+AFR)")

    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "fig2_finemapping.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    plot_locus_and_ld()
    plot_finemapping_results()
