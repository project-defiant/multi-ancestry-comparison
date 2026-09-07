"""Single-variant association against fine-mapping, on the AFR arm of the locus.

Builds ``results/fig_why_finemapping.png``: three panels, each carrying the
genotype regression it corresponds to, written out elementwise.

    1. single-variant association   y = x_j beta_j + eps,  j = 1..p
                                    E[beta_hat_j] = beta_j + sum_(k!=j) r_jk beta_k
    2. linkage disequilibrium       R = X'X / n,  r_jk = x_j' x_k / n
    3. fine-mapping                 y = X beta + eps over the locus
                                    beta = sum_l gamma_l b_l           (SuSiE)

Panels 1 and 3 are the same regression of phenotype on genotype, with the same y
and the same error vector. The only difference is how many genotype columns are
on the right-hand side: one at a time in panel 1, all of the locus at once in
panel 3. That is the whole argument, and writing the design out elementwise is
what makes it visible -- panel 1's second line names the bias explicitly, as
beta_j plus a sum of every other true effect weighted by its correlation.

The error term is drawn in both. It is part of the model, and a regression shown
without it reads as an identity.

WHAT "OVER THE LOCUS" MEANS
---------------------------
Panel 3's denominator is the locus, not the summary statistics as a whole:
clumping fixes a region first, and the multivariable fit runs over the variants
inside that region. Earlier wording ("one regression over all variants") implied
a genome-wide joint fit, which is not what is done and is not feasible.

FOUR EARLIER VERSIONS OF THIS NOTATION WERE WRONG
-------------------------------------------------
* beta_hat_j = (x_j' x_j)^-1 x_j' y -- the scalar estimator for one variant.
  Correct, but it cannot express "the estimate absorbs its neighbours".
* bold-symbol E[beta_marg] = R beta -- states the relation without showing it.
* the same expanded as a matrix product of R and beta -- shows the mixing, but
  starts from the estimator rather than the model, so the regression being fitted
  never appears and neither does the error term.
* "joint fine-mapping" as panel 3's title -- fine-mapping is the joint fit; the
  qualifier belongs to the multi-ancestry step, which this figure excludes.

SCOPE: ONE ANCESTRY
-------------------
The figure stops before the multi-ancestry step, so it motivates fine-mapping
without conflating it with the ancestry question. ``plot_results.py`` makes the
multi-ancestry comparison.

AFR rather than EUR because its credible set is the larger of the two (4 variants
against 3), so the figure understates rather than flatters.

WHAT IT DOES NOT CLAIM
----------------------
The multivariable fit goes from 50 candidates to 4. It does not identify the causal
variant: snp25 (causal, known only because this is a simulation) carries PIP
0.237, below snp21's 0.334. Both variants are labelled in panels 1 and 3 so the
reader can see it rather than being told.

Every annotated number is computed from the input files at run time.

INPUTS (already in the repo; no re-run needed)
    data/AFR.gwas.tsv          marginal beta, se, p, z per variant
    data/AFR.ld.tsv            50 x 50 correlation matrix
    results/AFR.susieR.cs.tsv  PIP and credible-set membership per variant

PLOTTING LIBRARY, AND HOW THE MATHS IS SET
------------------------------------------
``plot_results.py`` uses plotnine composited with Pillow. This uses matplotlib:
it needs a fixed-aspect heatmap with a colourbar between two annotated scatter
panels, and one equation strip on a common baseline beneath all three. plotnine
has no native form for that layout, which is why the existing script composites
rasters instead.

The equations are set by real LaTeX (``standalone`` + ``amsmath`` + ``bm``),
rendered to PDF and rasterised, then placed into the figure with
``fig.figimage``. matplotlib's own mathtext was tried first and is not adequate
here: it has no ``\operatorname``, its ``\bm`` support is partial, and bold
upright matrix symbols came out inconsistent with the surrounding text.

KaTeX would have been the other option -- the Marp deck this figure feeds
declares ``math: katex``, and the same source strings render there natively --
but there is no npm access in this environment, and KaTeX is a browser
reimplementation of TeX maths, so rendering with TeX itself gives the same
notation from the same markup.

Run: uv run python scripts/plot_why_finemapping.py    (from example/)
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
EXAMPLE = HERE.parent
OUT = EXAMPLE / "results" / "fig_why_finemapping.png"

# Open Targets palette, so the figure sits inside the deck it serves.
NAVY = "#163A5F"
BLUE = "#2C77B5"
BLUE_MID = "#4A93CB"
LIGHT = "#8FB8DE"
GREY = "#58595B"
MUTE = "#B9BEC4"
FAINT = "#CFDCE8"
BAND = "#EDF3F9"
ACCENT = "#B4530A"          # reserved for the causal variant and nothing else

# r2 is a magnitude on [0, 1] with a meaningful zero, so a single-hue sequential
# ramp. The categorical hues above are never used for it; the accent never
# encodes a magnitude.
LD_CMAP = LinearSegmentedColormap.from_list("ot_seq", ["#FFFFFF", LIGHT, BLUE, NAVY])

GW_SIG = 5e-8
TRUE_CAUSAL = "snp25"       # known only because this is a simulation
R2_TIGHT = 0.8

# Equation strips, on a common baseline under all three panels. Kept as data so
# the layout code stays free of literals.
# Bottom of the expanded matrix equations (figure fraction). Panel 2's shorter
# lines are centred against that block at run time rather than pinned here,
# because their offset depends on the rendered height of the tall ones.
EQ_BOTTOM = 0.055


def latex_png(tex, colour=NAVY, pt=12, dpi=220, supersample=4):
    """Render one inline-maths string with LaTeX and return an RGB array.

    `standalone` with `preview` crops to the formula's own bounding box, so the
    returned array carries no padding to align around.

    Rendered at `supersample` x the figure's dpi and then reduced, so the result
    is `pt` points tall at `dpi` with antialiasing. Rendering directly at the
    target dpi gives visibly rough glyph edges; rendering high and pasting 1:1
    puts the equation on screen at four times its intended size.
    """
    from PIL import Image
    import shutil
    import subprocess
    import tempfile

    if shutil.which("pdflatex") is None:
        raise SystemExit(
            "pdflatex not found. The equation strips are set with real LaTeX; "
            "install a TeX distribution, or drop the strips from the figure.")

    doc = (
        f"\\documentclass[preview,border=1pt,{pt}pt]{{standalone}}\n"
        "\\usepackage{amsmath,amssymb,bm,xcolor}\n"
        f"\\definecolor{{fg}}{{HTML}}{{{colour.lstrip('#')}}}\n"
        "\\begin{document}%\n"
        f"\\color{{fg}}{tex}%\n"
        "\\end{document}\n"
    )
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        (tmp / "eq.tex").write_text(doc)
        run = subprocess.run(
            ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "eq.tex"],
            cwd=tmp, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if not (tmp / "eq.pdf").exists():
            raise SystemExit("pdflatex failed on:\n" + tex + "\n" +
                             run.stdout.decode()[-1500:])
        subprocess.run(["pdftoppm", "-r", str(dpi * supersample), "-png",
                        "-singlefile", "eq.pdf", "eq"], cwd=tmp, check=True)
        img = Image.open(tmp / "eq.png").convert("RGB")
        img = img.resize((max(1, img.width // supersample),
                          max(1, img.height // supersample)), Image.LANCZOS)
        return np.asarray(img) / 255.0


def head(ax, index, title, subtitle):
    """Panel number, title and a plain description of what is plotted.

    Drawn in axes coordinates rather than via set_title(loc="left"), which is not
    width-constrained and silently runs across the neighbouring panel.
    """
    ax.text(0, 1.235, f"{index}", transform=ax.transAxes, fontsize=13,
            fontweight="bold", color=BLUE_MID, va="baseline")
    ax.annotate(title, xy=(0, 1.235), xycoords=ax.transAxes,
                textcoords="offset points", xytext=(15, 0), va="baseline",
                fontsize=13.5, fontweight="bold", color=NAVY)
    ax.text(0, 1.135, subtitle, transform=ax.transAxes, fontsize=10.5,
            color=GREY, va="baseline")


def load():
    gwas = pd.read_csv(EXAMPLE / "data" / "AFR.gwas.tsv", sep="\t")
    ld = pd.read_csv(EXAMPLE / "data" / "AFR.ld.tsv", sep="\t")
    cs = pd.read_csv(EXAMPLE / "results" / "AFR.susieR.cs.tsv", sep="\t")
    if list(ld.columns) != list(gwas.snp):
        raise SystemExit("LD column order does not match the GWAS variant order")
    df = gwas.merge(cs[["snp", "pip", "cs_id"]], on="snp", validate="one_to_one")
    df["kb"] = (df.pos - df.pos.min()) / 1000.0
    df["logp"] = -np.log10(df.pval)
    return df, ld.values.astype(float)


def panel_gwas(ax, df, lead, n_sig):
    sig = df.pval < GW_SIG
    ax.scatter(df.kb[~sig], df.logp[~sig], s=52, c=MUTE, edgecolor="white",
               linewidth=0.6, zorder=2)

    ld_row = df[df.snp == lead].iloc[0]
    ax.scatter([ld_row.kb], [ld_row.logp], s=460, facecolor="none",
               edgecolor=NAVY, linewidth=2.2, zorder=3)
    # significant dots drawn last, so the ring cannot hide them
    ax.scatter(df.kb[sig], df.logp[sig], s=95, c=BLUE_MID, edgecolor="white",
               linewidth=0.9, zorder=7)
    ax.annotate(f"{lead} (lead)", (ld_row.kb, ld_row.logp),
                textcoords="offset points", xytext=(-16, 16), ha="right",
                fontsize=10.5, color=NAVY)

    tc = df[df.snp == TRUE_CAUSAL].iloc[0]
    star_y = tc.logp + df.logp.max() * 0.10
    ax.plot([tc.kb, tc.kb], [tc.logp + 0.14, star_y - 0.12], color=ACCENT,
            lw=1.0, zorder=4)
    ax.scatter([tc.kb], [star_y], marker="*", s=300, c=ACCENT,
               edgecolor="white", linewidth=0.7, zorder=6)
    ax.annotate(f"{TRUE_CAUSAL} (causal)", (tc.kb, star_y),
                textcoords="offset points", xytext=(13, 0), ha="left",
                va="center", fontsize=10.5, color=ACCENT)

    ax.axhline(-np.log10(GW_SIG), color=GREY, ls=(0, (5, 4)), lw=1.1, zorder=1)
    ax.text(0, -np.log10(GW_SIG) + 0.20, "p = 5×10⁻⁸", ha="left", va="bottom",
            fontsize=9.5, color=GREY)

    ax.set_ylabel("−log₁₀ p", fontsize=11, color=GREY)
    ax.set_xlabel("position in locus (kb)", fontsize=10.5, color=GREY)
    ax.set_ylim(-0.5, df.logp.max() * 1.42)
    head(ax, "1", "Single-variant association",
         "one regression per variant, genome-wide")
    ax.text(0.985, 0.965, f"{n_sig} variants at p < 5×10⁻⁸",
            transform=ax.transAxes, ha="right", va="top", fontsize=10.5,
            color=NAVY)


def panel_ld(ax, fig, r2, n_tight, lead_i):
    im = ax.imshow(r2, cmap=LD_CMAP, vmin=0, vmax=1, origin="upper",
                   interpolation="nearest")
    ax.set_xticks([0, 24, 49])
    ax.set_yticks([0, 24, 49])
    ax.set_xticklabels(["1", "25", "50"], fontsize=9.5, color=GREY)
    ax.set_yticklabels(["1", "25", "50"], fontsize=9.5, color=GREY)
    ax.set_xlabel("variant", fontsize=10.5, color=GREY)
    ax.set_ylabel("variant", fontsize=10.5, color=GREY)
    for side in ax.spines.values():
        side.set_edgecolor(FAINT)

    ax.axhline(lead_i, color=ACCENT, lw=0.9, alpha=0.8)
    ax.axvline(lead_i, color=ACCENT, lw=0.9, alpha=0.8)

    head(ax, "2", "Linkage disequilibrium", "correlation among those variants")

    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03, ticks=[0, 0.5, 1.0])
    cb.set_label("r²", fontsize=10.5, color=GREY)
    cb.ax.tick_params(labelsize=9, colors=GREY)
    cb.outline.set_edgecolor(FAINT)
    return n_tight


def panel_pip(ax, df, n_cs, cs_mass):
    inset = df.cs_id == 1
    ax.vlines(df.kb[~inset], 0, df.pip[~inset], color=MUTE, lw=2.0, zorder=2)
    ax.scatter(df.kb[~inset], df.pip[~inset], s=40, c=MUTE, edgecolor="white",
               linewidth=0.6, zorder=3)
    ax.vlines(df.kb[inset], 0, df.pip[inset], color=NAVY, lw=2.8, zorder=4)
    ax.scatter(df.kb[inset], df.pip[inset], s=92, c=NAVY, edgecolor="white",
               linewidth=0.9, zorder=5)

    ymax = max(0.46, df.pip.max() * 1.50)
    tc = df[df.snp == TRUE_CAUSAL].iloc[0]
    star_y = tc.pip + ymax * 0.085
    ax.plot([tc.kb, tc.kb], [tc.pip + 0.010, star_y - 0.010], color=ACCENT,
            lw=1.0, zorder=6)
    ax.scatter([tc.kb], [star_y], marker="*", s=300, c=ACCENT,
               edgecolor="white", linewidth=0.7, zorder=7)
    ax.annotate(f"{TRUE_CAUSAL} (causal)", (tc.kb, star_y),
                textcoords="offset points", xytext=(12, 0), ha="left",
                va="center", fontsize=10.5, color=ACCENT)

    top = df.loc[df.pip.idxmax()]
    ax.annotate(f"{top.snp}", (top.kb, top.pip), textcoords="offset points",
                xytext=(-9, 4), ha="right", va="bottom", fontsize=10.5,
                color=NAVY)

    ax.set_ylabel("posterior inclusion probability", fontsize=11, color=GREY)
    ax.set_xlabel("position in locus (kb)", fontsize=10.5, color=GREY)
    ax.set_ylim(0, ymax)
    head(ax, "3", "Fine-mapping",
         "one regression over the variants in the locus")
    ax.text(0.985, 0.965, f"95% credible set: {n_cs} of {len(df)}, "
                          f"PIP {cs_mass:.2f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=10.5,
            color=NAVY)


def main():
    df, R = load()
    r2 = R ** 2
    snps = list(df.snp)
    lead = df.loc[df.pval.idxmin(), "snp"]
    lead_i = snps.index(lead)
    n_sig = int((df.pval < GW_SIG).sum())
    n_tight = int((r2[lead_i] >= R2_TIGHT).sum())
    cs = df[df.cs_id == 1]
    n_cs, cs_mass = len(cs), float(cs.pip.sum())

    fig = plt.figure(figsize=(15.0, 6.7), dpi=220)
    fig.patch.set_facecolor("white")
    gs = fig.add_gridspec(1, 3, width_ratios=[1.18, 0.80, 1.18],
                          left=0.048, right=0.985, top=0.805, bottom=0.350,
                          wspace=0.30)
    ax1, ax2, ax3 = (fig.add_subplot(gs[0, i]) for i in range(3))
    for ax in (ax1, ax3):
        ax.grid(axis="y", color=BAND, lw=1.0)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_edgecolor(FAINT)
        ax.tick_params(colors=GREY, labelsize=9.5)

    panel_gwas(ax1, df, lead, n_sig)
    panel_ld(ax2, fig, r2, n_tight, lead_i)
    panel_pip(ax3, df, n_cs, cs_mass)

    # Equation strips. Panel 2 holds a fixed aspect, so matplotlib adjusts its
    # axes box at draw time; querying positions after a draw is what keeps the
    # three strips on one baseline instead of one per axes height.
    fig.canvas.draw()
    # The genotype regression, written out. Panels 1 and 3 differ only in how
    # many genotype columns sit on the right-hand side.
    tall = [
        (ax1, r"$\begin{bmatrix}y_1\\ y_2\\ \vdots\\ y_n\end{bmatrix} = \begin{bmatrix}x_{1j}\\ x_{2j}\\ \vdots\\ x_{nj}\end{bmatrix}\beta_j + \begin{bmatrix}\varepsilon_1\\ \varepsilon_2\\ \vdots\\ \varepsilon_n\end{bmatrix}\quad j = 1,\dots,p$"),
        (ax3, r"$\begin{bmatrix}y_1\\ y_2\\ \vdots\\ y_n\end{bmatrix} = \begin{bmatrix}x_{11} & x_{12} & \cdots & x_{1p}\\x_{21} & x_{22} & \cdots & x_{2p}\\\vdots & \vdots & \ddots & \vdots\\x_{n1} & x_{n2} & \cdots & x_{np}\end{bmatrix}\begin{bmatrix}\beta_1\\ \beta_2\\ \vdots\\ \beta_p\end{bmatrix} + \begin{bmatrix}\varepsilon_1\\ \varepsilon_2\\ \vdots\\ \varepsilon_n\end{bmatrix}$"),
    ]
    fw, fh = fig.get_size_inches() * fig.dpi

    # the consequence line sits at the bottom, the model above it, so the strip
    # reads downward as "this is the fit" then "this is what it gives you"
    second = [
        (ax1, r"$\mathbb{E}[\hat\beta_j] = \beta_j"
              r" + \sum_{k \neq j} r_{jk}\,\beta_k$"),
        (ax3, r"$\bm\beta = \sum_{l=1}^{L} \gamma_l b_l,"
              r"\quad \gamma_l \sim \mathrm{Mult}(1, \bm\pi)$"),
    ]
    second_h = 0.0
    for ax, tex in second:
        img = latex_png(tex, dpi=int(fig.dpi), pt=11)
        fig.figimage(img, xo=ax.get_position().x0 * fw, yo=EQ_BOTTOM * fh,
                     zorder=5)
        second_h = max(second_h, img.shape[0] / fh)

    model_y = EQ_BOTTOM + second_h + 0.022
    tall_h = 0.0
    for ax, tex in tall:
        box = ax.get_position()
        img = latex_png(tex, dpi=int(fig.dpi))
        avail = (box.x1 - box.x0) * fw * 1.20
        if img.shape[1] > avail:
            raise SystemExit(
                f"equation is {img.shape[1]}px wide but its column allows "
                f"{avail:.0f}px; shorten it or drop a matrix column:\n  {tex}")
        fig.figimage(img, xo=box.x0 * fw, yo=model_y * fh, zorder=5)
        tall_h = max(tall_h, img.shape[0] / fh)

    # panel 2's lines are centred against the tall block, not baseline-aligned
    mid = model_y + tall_h / 2
    lines = [r"$\bm R = \tfrac{1}{n}\bm X^{\!\top}\bm X"
             r" \in \mathbb{R}^{p \times p}$",
             r"$r_{jk} = \tfrac{1}{n}\bm x_j^{\!\top}\bm x_k$"]
    imgs = [latex_png(t, dpi=int(fig.dpi)) for t in lines]
    gap = 0.030
    total = sum(i.shape[0] / fh for i in imgs) + gap
    y = mid + total / 2
    for img in imgs:
        y -= img.shape[0] / fh
        fig.figimage(img, xo=ax2.get_position().x0 * fw, yo=y * fh, zorder=5)
        y -= gap
    fig.text(ax2.get_position().x0, mid - total / 2 - 0.042,
             f"{n_tight} variants at r\u00b2 \u2265 {R2_TIGHT} with the lead",
             fontsize=10, color=GREY, va="baseline")

    fig.legend(handles=[
        Line2D([], [], marker="o", ls="", markerfacecolor=MUTE,
               markeredgecolor="white", markersize=8, label="variant"),
        Line2D([], [], marker="o", ls="", markerfacecolor=BLUE_MID,
               markeredgecolor="white", markersize=8, label="p < 5×10⁻⁸"),
        Line2D([], [], marker="o", ls="", markerfacecolor=NAVY,
               markeredgecolor="white", markersize=8, label="in 95% credible set"),
        Line2D([], [], marker="*", ls="", markerfacecolor=ACCENT,
               markeredgecolor="white", markersize=13,
               label="causal variant (simulation)"),
    ], loc="upper right", ncol=4, frameon=False, fontsize=10,
        bbox_to_anchor=(0.988, 1.010), handletextpad=0.4, columnspacing=1.8,
        labelcolor=GREY)

    fig.savefig(OUT, facecolor="white")
    print(f"wrote {OUT.relative_to(EXAMPLE)}")
    print(f"  {len(df)} variants over {df.kb.max():.0f} kb, AFR arm")
    print(f"  lead {lead} (p = {df.pval.min():.2e}); causal {TRUE_CAUSAL} "
          f"PIP {df.loc[df.snp == TRUE_CAUSAL, 'pip'].iat[0]:.3f}")
    print(f"  {n_sig} variants at p < 5e-8; {n_tight} at r2 >= {R2_TIGHT} "
          f"with the lead")
    print(f"  credible set {n_cs} variants ({', '.join(cs.snp)}), "
          f"PIP sum {cs_mass:.3f}")
    print(f"  highest PIP {df.loc[df.pip.idxmax(), 'snp']} "
          f"({df.pip.max():.3f}); causal variant is not the top-ranked one")


if __name__ == "__main__":
    main()
