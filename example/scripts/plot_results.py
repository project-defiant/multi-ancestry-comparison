import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image
from plotnine import (
    aes,
    coord_equal,
    element_blank,
    element_rect,
    element_text,
    geom_point,
    geom_text,
    geom_tile,
    ggplot,
    labs,
    scale_color_manual,
    scale_fill_gradient2,
    scale_y_reverse,
    theme,
    theme_minimal,
    ylim,
)

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
DATA_DIR = Path(__file__).resolve().parent.parent / "data"

PANEL_WIDTH = 6.5
PANEL_HEIGHT = 5.5
DPI = 200

RED = "#C44E52"
GREY = "#8C8C8C"


def theme_presentation(base_size=18, legend_position="none"):
    return (
        theme_minimal(base_size=base_size)
        + theme(
            legend_position=legend_position,
            plot_background=element_rect(fill="none", color="none"),
            panel_background=element_rect(fill="none", color="none"),
            legend_background=element_rect(fill="none", color="none"),
            legend_key=element_rect(fill="none", color="none"),
            panel_grid_minor=element_blank(),
            axis_text=element_text(size=base_size * 0.7),
            axis_title=element_text(size=base_size * 0.85),
            plot_title=element_text(size=base_size, weight="bold"),
            legend_text=element_text(size=base_size * 0.7),
            legend_title=element_text(size=base_size * 0.8),
        )
    )


def load_gwas(path):
    return pd.read_csv(path, sep="\t")


def load_ld(path):
    return pd.read_csv(path, sep="\t").values


def load_susieR_cs(path):
    df = pd.read_csv(path, sep="\t")
    df["in_cs"] = df["cs_id"] == 1
    return df


def load_sushie_weights(path):
    df = pd.read_csv(path, sep="\t")
    df["pip"] = df["sushie_pip_all"]
    df["in_cs"] = df["sushie_cs_index"] != "No CS"
    return df


def render_panel(plot, path):
    plot.save(
        path,
        width=PANEL_WIDTH,
        height=PANEL_HEIGHT,
        dpi=DPI,
        transparent=True,
        verbose=False,
    )


def compose_grid(panels, output_path):
    """panels: length-4 list in [top-left, top-right, bottom-left, bottom-right]
    order; a None entry leaves that tile transparent/blank."""
    px_w, px_h = int(PANEL_WIDTH * DPI), int(PANEL_HEIGHT * DPI)
    canvas = Image.new("RGBA", (px_w * 2, px_h * 2), (0, 0, 0, 0))
    with tempfile.TemporaryDirectory() as tmp:
        for idx, plot in enumerate(panels):
            if plot is None:
                continue
            panel_path = Path(tmp) / f"panel_{idx}.png"
            render_panel(plot, panel_path)
            img = Image.open(panel_path).convert("RGBA")
            if img.size != (px_w, px_h):
                img = img.resize((px_w, px_h), Image.LANCZOS)
            row, col = divmod(idx, 2)
            canvas.paste(img, (col * px_w, row * px_h), img)
    canvas.save(output_path)


def build_locus_zoom_panel(eur_gwas, afr_gwas):
    eur = eur_gwas.assign(ancestry="EUR", neglog10p=-np.log10(eur_gwas["pval"]))
    afr = afr_gwas.assign(ancestry="AFR", neglog10p=-np.log10(afr_gwas["pval"]))
    long_df = pd.concat([eur, afr], ignore_index=True)
    return (
        ggplot(long_df, aes(x="pos", y="neglog10p", color="ancestry"))
        + geom_point(size=4, alpha=0.85)
        + scale_color_manual(values={"EUR": "#4C72B0", "AFR": "#DD8452"}, name="Ancestry")
        + labs(x="Position (bp)", y="-log10(p)", title="GWAS locus zoom")
        + theme_presentation(legend_position="right")
    )


def build_ld_panel(ld_matrix, title):
    n = ld_matrix.shape[0]
    long_df = pd.DataFrame(
        [(i, j, ld_matrix[i, j]) for i in range(n) for j in range(n)],
        columns=["i", "j", "r"],
    )
    return (
        ggplot(long_df, aes(x="j", y="i", fill="r"))
        + geom_tile()
        + scale_fill_gradient2(low="#2166AC", mid="#F7F7F7", high="#B2182B", midpoint=0, limits=(-1, 1), name="r")
        + scale_y_reverse()
        + coord_equal()
        + labs(x="SNP index", y="SNP index", title=title)
        + theme_presentation(legend_position="right")
    )


def plot_locus_and_ld():
    eur_gwas = load_gwas(DATA_DIR / "EUR.gwas.tsv")
    afr_gwas = load_gwas(DATA_DIR / "AFR.gwas.tsv")
    eur_ld = load_ld(DATA_DIR / "EUR.ld.tsv")
    afr_ld = load_ld(DATA_DIR / "AFR.ld.tsv")

    locus_panel = build_locus_zoom_panel(eur_gwas, afr_gwas)
    afr_ld_panel = build_ld_panel(afr_ld, "AFR LD matrix")
    eur_ld_panel = build_ld_panel(eur_ld, "EUR LD matrix")

    compose_grid(
        [locus_panel, afr_ld_panel, eur_ld_panel, None],
        RESULTS_DIR / "fig1_locus_ld.png",
    )


def build_pip_panel(df, title):
    cs_size = int(df["in_cs"].sum())
    labels = df[df["in_cs"]].sort_values("pos").reset_index(drop=True)
    # Stack labels in tiers above the panel's highest point (rather than each
    # point's own PIP) so closely-spaced tag SNPs don't overlap each other.
    tier_gaps = [0.06, 0.18, 0.30]
    top = df["pip"].max()
    labels["label_y"] = [top + tier_gaps[i % len(tier_gaps)] for i in range(len(labels))]
    return (
        ggplot(df, aes(x="pos", y="pip", color="in_cs"))
        + geom_point(size=5)
        + geom_text(
            data=labels,
            mapping=aes(y="label_y", label="snp"),
            size=13,
            color=RED,
            fontweight="bold",
        )
        + scale_color_manual(values={True: RED, False: GREY})
        + ylim(-0.05, 1.4)
        + labs(x="Position (bp)", y="PIP", title=f"{title} (CS size = {cs_size})")
        + theme_presentation()
    )


def plot_finemapping_results():
    eur = load_susieR_cs(RESULTS_DIR / "EUR.susieR.cs.tsv")
    afr = load_susieR_cs(RESULTS_DIR / "AFR.susieR.cs.tsv")
    sushie = load_sushie_weights(RESULTS_DIR / "locus1.sushie.weights.tsv")

    afr_panel = build_pip_panel(afr, "AFR susieR")
    eur_panel = build_pip_panel(eur, "EUR susieR")
    sushie_panel = build_pip_panel(sushie, "sushie (EUR+AFR)")

    compose_grid(
        [None, afr_panel, eur_panel, sushie_panel],
        RESULTS_DIR / "fig2_finemapping.png",
    )


if __name__ == "__main__":
    plot_locus_and_ld()
    plot_finemapping_results()
