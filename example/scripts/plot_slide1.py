"""Three standalone figures for the "why fine-mapping" slide-1 narrative:
marginal scan (panel 1) -> locus definition (panel 2) -> fine-mapping (panel
3). Each panel is its own PNG, paired with its own markdown+katex file under
slides/ (gwas.md, locus_definition.md, finemapping.md) -- the equations live
in those markdown files, not in these figures.

Styled to match plot_results.py (the susieR-vs-sushie summary figure):
plotnine, theme_presentation(), transparent background, the same RED/GREY
significance palette and adjustText-nudged SNP labels, rather than the
Open-Targets-navy matplotlib look used in plot_why_finemapping.py.

Inputs (from simulate_locus_breaker.py + finemap_susieR.R; no re-run needed
unless those change):
    data/locus_scan.gwas.tsv          both loci, combined marginal scan
    data/locus_scan.locusA.gwas.tsv   locus A only
    results/locusA.susieR.cs.tsv      susieR PIP/CS for locus A

Locus-breaker parameters below are tuned to this toy window's scale (a
150kb gap between two ~30-60kb loci), not gentropy's production defaults
(250kb distance / 100kb flanking on real chromosome-scale data) -- the
mechanics are identical, just rescaled so the split is visible on a
230kb-wide figure.

Run: uv run python scripts/plot_slide1.py    (from example/)
"""

from pathlib import Path

import numpy as np
import pandas as pd
from adjustText import adjust_text
from plotnine import (
    aes,
    element_blank,
    element_rect,
    element_text,
    geom_hline,
    geom_point,
    geom_rect,
    geom_text,
    geom_vline,
    ggplot,
    labs,
    scale_color_manual,
    scale_fill_manual,
    theme,
    theme_minimal,
    ylim,
)

HERE = Path(__file__).resolve().parent
EXAMPLE = HERE.parent
DATA = EXAMPLE / "data"
RESULTS = EXAMPLE / "results"

PANEL_WIDTH = 8.0
PANEL_HEIGHT = 5.0
DPI = 200

RED = "#C44E52"
GREY = "#8C8C8C"
LOCUS_FILL = {"1": "#4C72B0", "2": "#DD8452"}

GW_SIG = 5e-8
LEAD_PVALUE = 5e-8
BASELINE_PVALUE = 1e-4
DISTANCE_CUTOFF = 100_000
FLANKING_DISTANCE = 20_000

CAUSAL_A = "locusA_snp18"
CAUSAL_B = "locusB_snp8"


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
            plot_subtitle=element_text(size=base_size * 0.6, color=GREY),
            legend_text=element_text(size=base_size * 0.7),
            legend_title=element_text(size=base_size * 0.8),
        )
    )


def render_panel(plot, label_anchors=None):
    """Draw a plotnine plot to a matplotlib Figure, nudging any geom_text
    labels apart with adjustText so they stay close to their own point while
    never overlapping (mirrors plot_results.py's render_panel)."""
    fig = plot.draw()
    ax = fig.axes[0]
    if label_anchors is not None and len(label_anchors) > 0:
        texts = list(ax.texts)
        adjust_text(
            texts,
            x=list(label_anchors["x"]),
            y=list(label_anchors["y"]),
            ax=ax,
            expand_axes=True,
            force_text=(0.4, 0.8),
            expand=(1.3, 1.6),
            arrowprops=dict(arrowstyle="-", color=RED, lw=1.2, alpha=0.7),
        )
    return fig


def save_panel(plot, label_anchors, path):
    fig = render_panel(plot, label_anchors)
    fig.savefig(path, dpi=DPI, transparent=True)


def load_scan():
    df = pd.read_csv(DATA / "locus_scan.gwas.tsv", sep="\t")
    df["kb"] = (df.pos - df.pos.min()) / 1000.0
    df["neglog10p"] = -np.log10(df.pval)
    df["sig"] = df.pval < GW_SIG
    return df


def locus_breaker(df, baseline_pvalue, distance_cutoff, pvalue_cutoff, flanking_distance):
    """Plain pandas re-implementation of gentropy's LocusBreakerClumping,
    minus Spark: drop points below the baseline, split on gaps between the
    survivors that exceed distance_cutoff, keep only windows whose top point
    clears the (stricter) lead pvalue_cutoff, then flank the boundaries."""
    surv = df[df.pval <= baseline_pvalue].sort_values("pos").reset_index(drop=True)
    gap = surv.pos.diff()
    new_locus = gap.isna() | (gap > distance_cutoff)
    surv["locus_id"] = new_locus.cumsum()

    loci = []
    for locus_id, group in surv.groupby("locus_id"):
        top = group.loc[group.pval.idxmin()]
        if top.pval <= pvalue_cutoff:
            loci.append({
                "locus_id": locus_id,
                "start": group.pos.min() - flanking_distance,
                "end": group.pos.max() + flanking_distance,
                "lead_snp": top.snp,
                "lead_pval": top.pval,
            })
    return pd.DataFrame(loci)


def causal_labels(df, snps):
    labels = df[df.snp.isin(snps)].copy()
    labels["label"] = labels.snp + " (causal)"
    return labels


def build_gwas_manhattan_panel(df):
    n_sig_a = int((df.sig & (df.locus == "A")).sum())
    labels = causal_labels(df, [CAUSAL_A, CAUSAL_B])

    plot = (
        ggplot(df, aes(x="kb", y="neglog10p", color="sig"))
        + geom_point(size=4, alpha=0.85)
        + geom_hline(yintercept=-np.log10(GW_SIG), linetype="dashed", color=GREY, size=0.8)
        + geom_text(
            data=labels,
            mapping=aes(x="kb", y="neglog10p", label="label"),
            inherit_aes=False,
            nudge_y=df.neglog10p.max() * 0.02,
            size=13,
            color=RED,
            fontweight="bold",
        )
        + scale_color_manual(values={True: RED, False: GREY})
        + labs(
            x="Position (kb)",
            y="-log10(p)",
            title="1. Single-variant association",
            subtitle=f"{n_sig_a} SNPs in locus A cross genome-wide significance —\n"
                     "the marginal test alone can't tell which is causal",
        )
        + theme_presentation()
    )
    label_anchors = labels[["kb", "neglog10p"]].rename(columns={"kb": "x", "neglog10p": "y"})
    return plot, label_anchors


def build_locus_definition_panel(df):
    loci = locus_breaker(df, BASELINE_PVALUE, DISTANCE_CUTOFF, LEAD_PVALUE, FLANKING_DISTANCE)
    x0 = df.pos.min()
    rects = pd.DataFrame({
        "xmin": (loci.start - x0) / 1000,
        "xmax": (loci.end - x0) / 1000,
        "ymin": -1.0,
        "ymax": df.neglog10p.max() * 1.15,
        "locus": [str(i + 1) for i in range(len(loci))],
    })
    a_end_kb = (df[df.locus == "A"].pos.max() - x0) / 1000
    b_start_kb = (df[df.locus == "B"].pos.min() - x0) / 1000
    gap_label = pd.DataFrame({
        "x": [(a_end_kb + b_start_kb) / 2],
        "y": [df.neglog10p.max() * 1.05],
        "label": [f"gap {DISTANCE_CUTOFF / 1000:.0f}kb+ → split into two loci"],
    })
    locus_labels = pd.DataFrame({
        "x": ((loci.start + loci.end) / 2 - x0) / 1000,
        "y": -0.6,
        "label": [f"locus {i + 1}\nlead: {row.lead_snp}" for i, row in loci.iterrows()],
    })

    plot = (
        ggplot(df, aes(x="kb", y="neglog10p"))
        + geom_rect(
            data=rects,
            mapping=aes(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax", fill="locus"),
            inherit_aes=False, alpha=0.15,
        )
        + scale_fill_manual(values=LOCUS_FILL, guide=None)
        + geom_point(mapping=aes(color="sig"), size=4, alpha=0.85)
        + scale_color_manual(values={True: RED, False: GREY})
        + geom_hline(yintercept=-np.log10(GW_SIG), linetype="dashed", color=GREY, size=0.8)
        + geom_hline(yintercept=-np.log10(BASELINE_PVALUE), linetype="dotted", color=GREY, size=0.8)
        + geom_vline(xintercept=a_end_kb, linetype="dotted", color=GREY, size=0.6)
        + geom_vline(xintercept=b_start_kb, linetype="dotted", color=GREY, size=0.6)
        + geom_text(data=gap_label, mapping=aes(x="x", y="y", label="label"),
                    inherit_aes=False, size=11, color=GREY)
        + geom_text(data=locus_labels, mapping=aes(x="x", y="y", label="label"),
                    inherit_aes=False, size=10, color="black", lineheight=1.0)
        + labs(
            x="Position (kb)",
            y="-log10(p)",
            title="2. Locus definition (locus-breaker)",
            subtitle=f"baseline p ≤ {BASELINE_PVALUE:.0e} keeps candidates; a gap > {DISTANCE_CUTOFF // 1000}kb splits the locus;\n"
                     f"only windows with a lead p ≤ {LEAD_PVALUE:.0e} are kept",
        )
        + theme_presentation()
    )
    return plot, None, loci


def build_finemapping_pip_panel():
    gwas = pd.read_csv(DATA / "locus_scan.locusA.gwas.tsv", sep="\t")
    cs = pd.read_csv(RESULTS / "locusA.susieR.cs.tsv", sep="\t")
    df = gwas.merge(cs[["snp", "pip", "cs_id"]], on="snp", validate="one_to_one")
    df["kb"] = (df.pos - df.pos.min()) / 1000.0
    df["in_cs"] = df.cs_id == 1
    n_cs, cs_mass = int(df.in_cs.sum()), float(df.pip[df.in_cs].sum())

    labels = df[df.in_cs].sort_values("pos").copy()
    labels["label"] = labels.snp + np.where(labels.snp == CAUSAL_A, " (causal)", "")

    plot = (
        ggplot(df, aes(x="kb", y="pip", color="in_cs"))
        + geom_point(size=5)
        + geom_text(
            data=labels,
            mapping=aes(x="kb", y="pip", label="label"),
            inherit_aes=False,
            nudge_y=0.05,
            size=13,
            color=RED,
            fontweight="bold",
        )
        + scale_color_manual(values={True: RED, False: GREY})
        + ylim(-0.05, 1.15)
        + labs(
            x="Position in locus A (kb)",
            y="PIP",
            title="3. Fine-mapping (SuSiE)",
            subtitle=f"one joint regression over the whole locus — 95% credible set: "
                     f"{n_cs} of {len(df)} SNPs, PIP {cs_mass:.2f}",
        )
        + theme_presentation()
    )
    label_anchors = labels[["kb", "pip"]].rename(columns={"kb": "x", "pip": "y"})
    return plot, label_anchors


def main():
    df = load_scan()

    plot, anchors = build_gwas_manhattan_panel(df)
    save_panel(plot, anchors, RESULTS / "fig_gwas_manhattan.png")
    print(f"wrote {RESULTS / 'fig_gwas_manhattan.png'}")

    plot, anchors, loci = build_locus_definition_panel(df)
    save_panel(plot, anchors, RESULTS / "fig_locus_breaker.png")
    print(f"wrote {RESULTS / 'fig_locus_breaker.png'}: {len(loci)} loci -> "
          f"{loci.lead_snp.tolist()}")

    plot, anchors = build_finemapping_pip_panel()
    save_panel(plot, anchors, RESULTS / "fig_finemapping_pip.png")
    print(f"wrote {RESULTS / 'fig_finemapping_pip.png'}")


if __name__ == "__main__":
    main()
