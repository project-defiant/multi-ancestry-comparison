"""Single-ancestry teaching facet: GWAS manhattan : LD matrix : SuSiE fine-mapping.

One row of three panels, per ancestry (EUR, AFR), reusing plot_results.py's
LD-matrix and PIP panels as-is -- both already draw every SNP in the 95%
credible set (not just the top-PIP one), so the fine-mapping panel here is
genuinely ambiguous: a reader sees several tied dots, not a single "answer"
variant. Only the manhattan panel (single ancestry, not the EUR+AFR overlay
used in plot_results.py's locus-zoom panel) is new.

Run: uv run python scripts/plot_finemapping_facet.py    (from example/)
"""

import numpy as np
from plotnine import (
    aes,
    geom_hline,
    geom_point,
    geom_text,
    ggplot,
    labs,
    scale_color_manual,
)

from plot_results import (
    DATA_DIR,
    GREY,
    RED,
    RESULTS_DIR,
    build_ld_panel,
    build_pip_panel,
    compose_grid,
    load_gwas,
    load_ld,
    load_susieR_cs,
    theme_presentation,
)

GW_SIG = 5e-8
CAUSAL_SNP = "snp25"  # known only because this is a simulation


def build_manhattan_panel(gwas_df, ancestry, causal_snp=CAUSAL_SNP, gw_sig=GW_SIG):
    df = gwas_df.assign(neglog10p=-np.log10(gwas_df["pval"]))
    df["sig"] = df["pval"] < gw_sig
    labels = df.loc[df["snp"] == causal_snp].assign(label=lambda d: d["snp"] + " (causal)")

    plot = (
        ggplot(df, aes(x="pos", y="neglog10p", color="sig"))
        + geom_point(size=4, alpha=0.85)
        + geom_hline(yintercept=-np.log10(gw_sig), linetype="dashed", color=GREY, size=0.8)
        + geom_text(
            data=labels,
            mapping=aes(label="label"),
            nudge_y=df["neglog10p"].max() * 0.04,
            size=13,
            color=RED,
            fontweight="bold",
        )
        + scale_color_manual(values={True: RED, False: GREY})
        + labs(x="Position (bp)", y="-log10(p)", title=f"{ancestry} GWAS marginal scan")
        + theme_presentation()
    )
    # render_panel's adjustText call expects columns named "pos"/"pip" -- reused
    # here purely as generic (x, y) anchor coordinates, not literal PIP values.
    label_anchors = labels[["pos", "neglog10p"]].rename(columns={"neglog10p": "pip"})
    return plot, label_anchors


def plot_ancestry_facet(ancestry):
    gwas = load_gwas(DATA_DIR / f"{ancestry}.gwas.tsv")
    ld = load_ld(DATA_DIR / f"{ancestry}.ld.tsv")
    cs = load_susieR_cs(RESULTS_DIR / f"{ancestry}.susieR.cs.tsv")

    panels = [
        build_manhattan_panel(gwas, ancestry),
        build_ld_panel(ld, f"{ancestry} LD matrix"),
        build_pip_panel(cs, f"{ancestry} susieR fine-mapping"),
    ]
    out_path = RESULTS_DIR / f"fig_finemapping_facet_{ancestry}.png"
    compose_grid(panels, out_path, rows=1, cols=3)
    print(f"wrote {out_path.relative_to(RESULTS_DIR.parent)}")


def main():
    for ancestry in ("EUR", "AFR"):
        plot_ancestry_facet(ancestry)


if __name__ == "__main__":
    main()
