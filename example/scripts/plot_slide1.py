"""Three standalone figures for the "why fine-mapping" slide-1 narrative:
marginal scan (panel 1) -> locus definition (panel 2) -> fine-mapping (panel
3). Each panel is its own PNG, paired with its own markdown+katex file under
slides/ (gwas.md, locus_definition.md, finemapping.md) -- the equations live
in those markdown files, not in these figures.

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

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parent
EXAMPLE = HERE.parent
DATA = EXAMPLE / "data"
RESULTS = EXAMPLE / "results"

# Open Targets palette, matching plot_why_finemapping.py so every slide-1
# figure sits inside the same deck.
NAVY = "#163A5F"
BLUE_MID = "#4A93CB"
GREY = "#58595B"
MUTE = "#B9BEC4"
FAINT = "#CFDCE8"
BAND = "#EDF3F9"
ACCENT = "#B4530A"          # reserved for the causal variant and nothing else
LOCUS_A_COLOR = "#2C77B5"
LOCUS_B_COLOR = "#8FB8DE"

GW_SIG = 5e-8
LEAD_PVALUE = 5e-8
BASELINE_PVALUE = 1e-4
DISTANCE_CUTOFF = 100_000
FLANKING_DISTANCE = 20_000

CAUSAL_A = "locusA_snp18"
CAUSAL_B = "locusB_snp8"


def style_axes(ax):
    ax.grid(axis="y", color=BAND, lw=1.0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_edgecolor(FAINT)
    ax.tick_params(colors=GREY, labelsize=9.5)


def head(ax, title, subtitle):
    ax.text(0, 1.16, title, transform=ax.transAxes, fontsize=14,
             fontweight="bold", color=NAVY, va="baseline")
    ax.text(0, 1.07, subtitle, transform=ax.transAxes, fontsize=10.5,
             color=GREY, va="baseline")


def load_scan():
    df = pd.read_csv(DATA / "locus_scan.gwas.tsv", sep="\t")
    df["kb"] = (df.pos - df.pos.min()) / 1000.0
    df["logp"] = -np.log10(df.pval)
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


def scatter_scan(ax, df):
    sig = df.pval < GW_SIG
    colors = np.where(df.locus == "A", LOCUS_A_COLOR, LOCUS_B_COLOR)
    ax.scatter(df.kb[~sig], df.logp[~sig], s=50, c=MUTE, edgecolor="white",
               linewidth=0.6, zorder=2)
    ax.scatter(df.kb[sig], df.logp[sig], s=95, c=list(colors[sig.values]),
               edgecolor="white", linewidth=0.9, zorder=4)
    ax.axhline(-np.log10(GW_SIG), color=GREY, ls=(0, (5, 4)), lw=1.1, zorder=1)
    ax.text(df.kb.max(), -np.log10(GW_SIG) + 0.25, "p = 5×10⁻⁸",
            ha="right", va="bottom", fontsize=9.5, color=GREY)
    ax.set_ylabel("−log₁₀ p", fontsize=11, color=GREY)
    ax.set_xlabel("position (kb)", fontsize=10.5, color=GREY)
    ax.set_ylim(-0.5, df.logp.max() * 1.25)
    return sig


def annotate_causal(ax, df, snp, dy=0.0):
    row = df[df.snp == snp].iloc[0]
    y = row.logp + df.logp.max() * 0.08 + dy
    ax.plot([row.kb, row.kb], [row.logp + 0.15, y - 0.12], color=ACCENT, lw=1.0, zorder=5)
    ax.scatter([row.kb], [y], marker="*", s=260, c=ACCENT, edgecolor="white",
               linewidth=0.7, zorder=6)
    near_right_edge = row.kb > df.kb.max() * 0.85
    ha, xoffset = ("right", -10) if near_right_edge else ("left", 10)
    ax.annotate(f"{row.snp} (causal)", (row.kb, y), textcoords="offset points",
                xytext=(xoffset, 0), ha=ha, va="center", fontsize=10, color=ACCENT)


def plot_gwas_manhattan(df):
    fig, ax = plt.subplots(figsize=(10.5, 5.6), dpi=200)
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=0.09, right=0.97, top=0.80, bottom=0.15)
    style_axes(ax)

    sig = scatter_scan(ax, df)
    annotate_causal(ax, df, CAUSAL_A)
    annotate_causal(ax, df, CAUSAL_B, dy=df.logp.max() * 0.18)

    n_sig_a = int(((df.locus == "A") & sig).sum())
    ax.text(0.015, 0.965,
            f"Locus A: {n_sig_a} variants cross genome-wide significance — "
            "the marginal test alone cannot tell which is causal",
            transform=ax.transAxes, ha="left", va="top", fontsize=10.5, color=NAVY)

    head(ax, "1. Single-variant (marginal) association",
         "one regression per SNP, testing each in isolation")
    fig.savefig(RESULTS / "fig_gwas_manhattan.png", facecolor="white")
    plt.close(fig)
    print(f"wrote {RESULTS / 'fig_gwas_manhattan.png'}: {n_sig_a} significant in locus A")


def plot_locus_definition(df):
    loci = locus_breaker(df, BASELINE_PVALUE, DISTANCE_CUTOFF, LEAD_PVALUE, FLANKING_DISTANCE)

    fig, ax = plt.subplots(figsize=(10.5, 5.6), dpi=200)
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=0.09, right=0.97, top=0.80, bottom=0.15)
    style_axes(ax)

    x0 = df.pos.min()
    for _, locus in loci.iterrows():
        ax.axvspan((locus.start - x0) / 1000, (locus.end - x0) / 1000,
                   color=LOCUS_A_COLOR if locus.locus_id == loci.locus_id.min() else LOCUS_B_COLOR,
                   alpha=0.12, zorder=0)

    scatter_scan(ax, df)
    ax.axhline(-np.log10(BASELINE_PVALUE), color=BLUE_MID, ls=(0, (2, 3)), lw=1.1, zorder=1)
    ax.text(df.kb.max(), -np.log10(BASELINE_PVALUE) + 0.25,
            f"baseline p = {BASELINE_PVALUE:.0e} (corridor floor)",
            ha="right", va="bottom", fontsize=9.5, color=BLUE_MID)

    a_end_kb = (df[df.locus == "A"].pos.max() - x0) / 1000
    b_start_kb = (df[df.locus == "B"].pos.min() - x0) / 1000
    y_bracket = df.logp.max() * 1.05
    ax.annotate("", xy=(b_start_kb, y_bracket), xytext=(a_end_kb, y_bracket),
                arrowprops=dict(arrowstyle="<->", color=GREY, lw=1.1))
    ax.text((a_end_kb + b_start_kb) / 2, y_bracket + df.logp.max() * 0.03,
            f"gap {DISTANCE_CUTOFF / 1000:.0f}kb+ → split into two loci",
            ha="center", va="bottom", fontsize=9.5, color=GREY)
    ax.set_ylim(-0.5, df.logp.max() * 1.30)

    for _, locus in loci.iterrows():
        mid = ((locus.start + locus.end) / 2 - x0) / 1000
        ax.text(mid, -0.35, f"locus {int(locus.locus_id)}\nlead: {locus.lead_snp}",
                ha="center", va="top", fontsize=9, color=NAVY)

    head(ax, "2. Locus definition (locus-breaker)",
         f"baseline p ≤ {BASELINE_PVALUE:.0e} keeps candidates; a gap > "
         f"{DISTANCE_CUTOFF // 1000}kb between them splits the locus; only "
         f"windows with a lead p ≤ {LEAD_PVALUE:.0e} are kept")
    fig.savefig(RESULTS / "fig_locus_breaker.png", facecolor="white")
    plt.close(fig)
    print(f"wrote {RESULTS / 'fig_locus_breaker.png'}: {len(loci)} loci -> "
          f"{loci.lead_snp.tolist()}")


def plot_finemapping_pip():
    gwas = pd.read_csv(DATA / "locus_scan.locusA.gwas.tsv", sep="\t")
    cs = pd.read_csv(RESULTS / "locusA.susieR.cs.tsv", sep="\t")
    df = gwas.merge(cs[["snp", "pip", "cs_id"]], on="snp", validate="one_to_one")
    df["kb"] = (df.pos - df.pos.min()) / 1000.0
    inset = df.cs_id == 1
    n_cs, cs_mass = int(inset.sum()), float(df.pip[inset].sum())

    fig, ax = plt.subplots(figsize=(10.5, 5.6), dpi=200)
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=0.09, right=0.97, top=0.80, bottom=0.15)
    style_axes(ax)

    ax.vlines(df.kb[~inset], 0, df.pip[~inset], color=MUTE, lw=2.0, zorder=2)
    ax.scatter(df.kb[~inset], df.pip[~inset], s=45, c=MUTE, edgecolor="white",
               linewidth=0.6, zorder=3)
    ax.vlines(df.kb[inset], 0, df.pip[inset], color=NAVY, lw=2.8, zorder=4)
    ax.scatter(df.kb[inset], df.pip[inset], s=100, c=NAVY, edgecolor="white",
               linewidth=0.9, zorder=5)

    ymax = max(0.5, df.pip.max() * 1.35)
    row = df[df.snp == CAUSAL_A].iloc[0]
    star_y = row.pip + ymax * 0.08
    ax.plot([row.kb, row.kb], [row.pip + 0.012, star_y - 0.012], color=ACCENT, lw=1.0, zorder=6)
    ax.scatter([row.kb], [star_y], marker="*", s=260, c=ACCENT, edgecolor="white",
               linewidth=0.7, zorder=7)
    ax.annotate(f"{row.snp} (causal)", (row.kb, star_y), textcoords="offset points",
                xytext=(10, 0), ha="left", va="center", fontsize=10, color=ACCENT)

    ax.set_ylabel("posterior inclusion probability", fontsize=11, color=GREY)
    ax.set_xlabel("position in locus A (kb)", fontsize=10.5, color=GREY)
    ax.set_ylim(0, ymax)
    ax.text(0.985, 0.965, f"95% credible set: {n_cs} of {len(df)} SNPs, PIP {cs_mass:.2f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=10.5, color=NAVY)

    head(ax, "3. Fine-mapping (SuSiE)",
         "one joint regression over the whole locus, using its LD structure")
    fig.savefig(RESULTS / "fig_finemapping_pip.png", facecolor="white")
    plt.close(fig)
    print(f"wrote {RESULTS / 'fig_finemapping_pip.png'}: CS size {n_cs} ({', '.join(df.snp[inset])})")


def main():
    df = load_scan()
    plot_gwas_manhattan(df)
    plot_locus_definition(df)
    plot_finemapping_pip()


if __name__ == "__main__":
    main()
