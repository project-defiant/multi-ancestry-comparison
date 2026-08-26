import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from adjustText import adjust_text
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
GRID_ROWS = 2
GRID_COLS = 3

RED = "#C44E52"
GREY = "#8C8C8C"


def theme_presentation(base_size=18, legend_position="none"):
    return (
        theme_minimal(base_size=base_size)
        + theme(
            legend_position=legend_position,
            figure_size=(PANEL_WIDTH, PANEL_HEIGHT),
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


def render_panel(plot, label_anchors=None):
    """Draw a plotnine plot to a matplotlib Figure. If label_anchors (a
    DataFrame with 'pos'/'pip' columns, one row per geom_text label already
    on the plot, in the same order) is given, nudge the labels apart with
    adjustText so they stay as close as possible to their own point while
    never overlapping, drawing a thin leader line back to the point when a
    label has to move away from it."""
    fig = plot.draw()
    ax = fig.axes[0]
    if label_anchors is not None and len(label_anchors) > 0:
        texts = list(ax.texts)
        adjust_text(
            texts,
            x=list(label_anchors["pos"]),
            y=list(label_anchors["pip"]),
            ax=ax,
            expand_axes=True,
            force_text=(0.4, 0.8),
            expand=(1.3, 1.6),
            arrowprops=dict(arrowstyle="-", color=RED, lw=1.2, alpha=0.7),
        )
    return fig


def panel_to_image(fig):
    px_w, px_h = int(PANEL_WIDTH * DPI), int(PANEL_HEIGHT * DPI)
    with tempfile.TemporaryDirectory() as tmp:
        panel_path = Path(tmp) / "panel.png"
        fig.savefig(panel_path, dpi=DPI, transparent=True)
        img = Image.open(panel_path).convert("RGBA")
    if img.size != (px_w, px_h):
        img = img.resize((px_w, px_h), Image.LANCZOS)
    return img


def compose_grid(panels, output_path, rows=GRID_ROWS, cols=GRID_COLS):
    """panels: row-major list of (plot, label_anchors) tuples, or None for a
    blank tile. label_anchors is passed through to render_panel."""
    px_w, px_h = int(PANEL_WIDTH * DPI), int(PANEL_HEIGHT * DPI)
    canvas = Image.new("RGBA", (px_w * cols, px_h * rows), (0, 0, 0, 0))
    for idx, panel in enumerate(panels):
        if panel is None:
            continue
        plot, label_anchors = panel
        fig = render_panel(plot, label_anchors)
        img = panel_to_image(fig)
        row, col = divmod(idx, cols)
        canvas.paste(img, (col * px_w, row * px_h), img)
    canvas.save(output_path)


def build_locus_zoom_panel(eur_gwas, afr_gwas):
    eur = eur_gwas.assign(ancestry="EUR", neglog10p=-np.log10(eur_gwas["pval"]))
    afr = afr_gwas.assign(ancestry="AFR", neglog10p=-np.log10(afr_gwas["pval"]))
    long_df = pd.concat([eur, afr], ignore_index=True)
    plot = (
        ggplot(long_df, aes(x="pos", y="neglog10p", color="ancestry"))
        + geom_point(size=4, alpha=0.85)
        + scale_color_manual(values={"EUR": "#4C72B0", "AFR": "#DD8452"}, name="Ancestry")
        + labs(x="Position (bp)", y="-log10(p)", title="GWAS locus zoom")
        + theme_presentation(legend_position="right")
    )
    return plot, None


def build_ld_panel(ld_matrix, title):
    n = ld_matrix.shape[0]
    long_df = pd.DataFrame(
        [(i, j, ld_matrix[i, j]) for i in range(n) for j in range(n)],
        columns=["i", "j", "r"],
    )
    plot = (
        ggplot(long_df, aes(x="j", y="i", fill="r"))
        + geom_tile()
        + scale_fill_gradient2(low="#2166AC", mid="#F7F7F7", high="#B2182B", midpoint=0, limits=(-1, 1), name="r")
        + scale_y_reverse()
        + coord_equal()
        + labs(x="SNP index", y="SNP index", title=title)
        + theme_presentation(legend_position="right")
    )
    return plot, None


def build_pip_panel(df, title):
    cs_size = int(df["in_cs"].sum())
    labels = df[df["in_cs"]].sort_values("pos").reset_index(drop=True)
    plot = (
        ggplot(df, aes(x="pos", y="pip", color="in_cs"))
        + geom_point(size=5)
        + geom_text(
            data=labels,
            mapping=aes(label="snp"),
            nudge_y=0.035,
            size=13,
            color=RED,
            fontweight="bold",
        )
        + scale_color_manual(values={True: RED, False: GREY})
        + ylim(-0.05, 1.15)
        + labs(x="Position (bp)", y="PIP", title=f"{title} (CS size = {cs_size})")
        + theme_presentation()
    )
    return plot, labels


def plot_summary_figure():
    eur_gwas = load_gwas(DATA_DIR / "EUR.gwas.tsv")
    afr_gwas = load_gwas(DATA_DIR / "AFR.gwas.tsv")
    eur_ld = load_ld(DATA_DIR / "EUR.ld.tsv")
    afr_ld = load_ld(DATA_DIR / "AFR.ld.tsv")

    eur_cs = load_susieR_cs(RESULTS_DIR / "EUR.susieR.cs.tsv")
    afr_cs = load_susieR_cs(RESULTS_DIR / "AFR.susieR.cs.tsv")
    sushie = load_sushie_weights(RESULTS_DIR / "locus1.sushie.weights.tsv")

    panels = [
        build_locus_zoom_panel(eur_gwas, afr_gwas),
        build_ld_panel(afr_ld, "AFR LD matrix"),
        build_ld_panel(eur_ld, "EUR LD matrix"),
        build_pip_panel(afr_cs, "AFR susieR"),
        build_pip_panel(eur_cs, "EUR susieR"),
        build_pip_panel(sushie, "sushie (EUR+AFR)"),
    ]

    compose_grid(panels, RESULTS_DIR / "fig_finemapping_summary.png")


if __name__ == "__main__":
    plot_summary_figure()
