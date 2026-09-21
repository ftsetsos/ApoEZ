#!/usr/bin/env python3
"""Generate cohort-level APOE genotype and magnitude summary plots.

Input: the per-sample TSV produced by apoe_genotyper.py.
Output: publication-ready PNG/SVG figures, a multipage PDF, and summary TSVs.
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import PercentFormatter


APOE_ORDER = [
    "ε1/ε1", "ε1/ε2", "ε1/ε3", "ε1/ε4",
    "ε2/ε2", "ε2/ε3", "ε2/ε4",
    "ε3/ε3", "ε3/ε4", "ε4/ε4",
    "ε1/ε3 OR ε2/ε4", "UNKNOWN",
]

MARKERS = [
    "rs449647", "rs267606664", "rs121918393", "rs387906567",
    "rs121918394", "rs199768005", "rs4420638",
]

COLORS = {
    "navy": "#26332E",
    "green": "#6E7C5B",
    "terracotta": "#B66F4A",
    "gold": "#D7BE91",
    "brown": "#9B8155",
    "pale_green": "#EEF1EA",
    "pale_orange": "#F8F1EA",
    "light": "#F7F4EC",
    "gray": "#A8AEA9",
}


def set_style() -> None:
    plt.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.edgecolor": COLORS["navy"],
        "axes.labelcolor": COLORS["navy"],
        "axes.titlecolor": COLORS["navy"],
        "axes.titleweight": "bold",
        "font.family": "sans-serif",
        "font.sans-serif": ["Aptos", "Arial", "DejaVu Sans"],
        "font.size": 10,
        "xtick.color": COLORS["navy"],
        "ytick.color": COLORS["navy"],
        "grid.color": "#D9DDD8",
        "grid.linewidth": 0.7,
        "grid.alpha": 0.65,
        "legend.frameon": False,
        "savefig.facecolor": "white",
        "savefig.bbox": "tight",
    })


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series.replace(".", np.nan), errors="coerce")


def ordered_genotypes(values: pd.Series) -> list[str]:
    observed = [str(x) for x in values.dropna().unique()]
    return [x for x in APOE_ORDER if x in observed] + sorted(set(observed) - set(APOE_ORDER))


def clean_axes(ax: plt.Axes, grid_axis: str | None = None) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid_axis:
        ax.grid(axis=grid_axis)
        ax.set_axisbelow(True)


def save_figure(fig: plt.Figure, outdir: Path, stem: str, dpi: int, pdf: PdfPages) -> None:
    fig.savefig(outdir / f"{stem}.png", dpi=dpi)
    fig.savefig(outdir / f"{stem}.svg")
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def genotype_prevalence(df: pd.DataFrame) -> plt.Figure:
    order = ordered_genotypes(df["APOE_genotype"])
    counts = df["APOE_genotype"].value_counts().reindex(order, fill_value=0)
    pct = counts / max(counts.sum(), 1) * 100
    fig, ax = plt.subplots(figsize=(10, 5.6))
    bars = ax.bar(counts.index, counts.values, color=COLORS["green"], width=0.72)
    for bar, n, p in zip(bars, counts, pct):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f"{n:,}\n{p:.1f}%",
                ha="center", va="bottom", fontsize=8)
    ax.set(
        xlabel="APOE genotype",
        ylabel="Samples"
    )
    ax.set_title("APOE genotype distribution", pad=24)
    ax.tick_params(axis="x", rotation=35)
    clean_axes(ax, "y")
    fig.tight_layout()
    return fig


def allele_frequencies(df: pd.DataFrame) -> tuple[plt.Figure, pd.DataFrame]:
    rows = []
    for allele in ("ε1", "ε2", "ε3", "ε4"):
        col = f"{allele}_dose"
        values = numeric(df[col]) if col in df else pd.Series(dtype=float)
        rows.append({"allele": allele, "copies": values.sum(min_count=1), "called_samples": values.notna().sum()})
    out = pd.DataFrame(rows)
    out["denominator_alleles"] = out["called_samples"] * 2
    out["frequency"] = out["copies"] / out["denominator_alleles"].replace(0, np.nan)
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    bars = ax.bar(out["allele"], out["frequency"] * 100,
                  color=[COLORS["gold"], COLORS["green"], COLORS["brown"], COLORS["terracotta"]])
    for bar, value in zip(bars, out["frequency"]):
        if pd.notna(value):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(), f"{value*100:.2f}%",
                    ha="center", va="bottom")
    ax.set(title="APOE allele frequencies", xlabel="Allele", ylabel="Allele frequency")
    ax.yaxis.set_major_formatter(PercentFormatter(100))
    clean_axes(ax, "y")
    fig.tight_layout()
    return fig, out


def overall_distribution(df: pd.DataFrame) -> plt.Figure:
    scores = numeric(df["overall_magnitude_score"]).dropna()
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    bins = min(35, max(8, int(math.sqrt(max(len(scores), 1)))))
    ax.hist(scores, bins=bins, color=COLORS["terracotta"], edgecolor="white", linewidth=0.7)
    if not scores.empty:
        median = scores.median()
        ax.axvline(median, color=COLORS["navy"], linestyle="--", linewidth=1.6,
                   label=f"Median = {median:.2f}")
        ax.legend()
    ax.set(title="Overall magnitude distribution", xlabel="Overall magnitude score", ylabel="Samples")
    clean_axes(ax, "y")
    fig.tight_layout()
    return fig


def score_by_genotype(df: pd.DataFrame) -> plt.Figure:
    work = df[["APOE_genotype", "overall_magnitude_score"]].copy()
    work["overall_magnitude_score"] = numeric(work["overall_magnitude_score"])
    work = work.dropna()
    order = ordered_genotypes(work["APOE_genotype"])
    groups = [work.loc[work["APOE_genotype"] == genotype, "overall_magnitude_score"].values for genotype in order]
    fig, ax = plt.subplots(figsize=(10.5, 5.8))
    if groups:
        bp = ax.boxplot(groups, labels=order, patch_artist=True, showfliers=False, widths=0.62)
        for patch in bp["boxes"]:
            patch.set_facecolor(COLORS["pale_orange"])
            patch.set_edgecolor(COLORS["terracotta"])
        for median in bp["medians"]:
            median.set_color(COLORS["navy"])
            median.set_linewidth(1.5)
        for i, group in enumerate(groups, 1):
            ax.text(i, ax.get_ylim()[1], f"n={len(group):,}", ha="center", va="bottom", fontsize=7)
    ax.set(
        xlabel="APOE genotype",
        ylabel="Overall magnitude score"
    )
    ax.set_title("Overall magnitude by APOE genotype", pad=24)
    ax.tick_params(axis="x", rotation=35)
    clean_axes(ax, "y")
    fig.tight_layout()
    return fig


def marker_genotype_distribution(df: pd.DataFrame) -> tuple[plt.Figure, pd.DataFrame]:
    rows = []
    for marker in MARKERS:
        col = f"{marker}_genotype"
        if col not in df:
            continue
        values = df[col].astype(str)
        values = values[~values.isin([".", "nan", "None"])]
        total = len(values)
        for genotype, count in values.value_counts().items():
            rows.append({"marker": marker, "genotype": genotype, "count": count,
                         "percent": count / total * 100 if total else np.nan})
    summary = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(10.5, 5.8))
    if summary.empty:
        ax.text(0.5, 0.5, "No auxiliary-marker genotypes available", ha="center", va="center")
        ax.axis("off")
        return fig, summary
    pivot = summary.pivot(index="marker", columns="genotype", values="percent").fillna(0).reindex(MARKERS).dropna(how="all")
    palette = [COLORS["green"], COLORS["terracotta"], COLORS["gold"], COLORS["brown"], "#8796A5", "#B59CA8"]
    left = np.zeros(len(pivot))
    for i, genotype in enumerate(pivot.columns):
        vals = pivot[genotype].values
        ax.barh(pivot.index, vals, left=left, label=genotype, color=palette[i % len(palette)], height=0.7)
        left += vals
    ax.set(title="Genotype composition at scored auxiliary markers", xlabel="Samples with an observed genotype", ylabel="")
    ax.xaxis.set_major_formatter(PercentFormatter(100))
    ax.set_xlim(0, 100)
    ax.legend(title="Genotype", bbox_to_anchor=(1.02, 1), loc="upper left")
    clean_axes(ax, "x")
    fig.tight_layout()
    return fig, summary


def marker_contributions(df: pd.DataFrame) -> tuple[plt.Figure, pd.DataFrame]:
    rows = []
    n = len(df)
    for marker in MARKERS:
        col = f"{marker}_magnitude"
        if col not in df:
            continue
        values = numeric(df[col]).fillna(0)
        rows.append({
            "marker": marker,
            "nonzero_samples": int((values > 0).sum()),
            "nonzero_percent": (values > 0).mean() * 100 if n else np.nan,
            "total_magnitude": values.sum(),
            "mean_magnitude": values.mean() if n else np.nan,
        })
    summary = pd.DataFrame(rows).sort_values("nonzero_percent", ascending=True)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.8), sharey=True)
    if summary.empty:
        axes[0].text(0.5, 0.5, "No marker magnitudes available", ha="center", va="center")
        axes[0].axis("off")
        axes[1].axis("off")
        return fig, summary
    axes[0].barh(summary["marker"], summary["nonzero_percent"], color=COLORS["green"])
    axes[0].set(title="Nonzero-magnitude prevalence", xlabel="Samples", ylabel="")
    axes[0].xaxis.set_major_formatter(PercentFormatter(100))
    axes[1].barh(summary["marker"], summary["total_magnitude"], color=COLORS["terracotta"])
    axes[1].set(title="Cumulative magnitude contribution", xlabel="Sum across cohort", ylabel="")
    clean_axes(axes[0], "x")
    clean_axes(axes[1], "x")
    fig.suptitle("Auxiliary-marker contribution overview", fontsize=14, fontweight="bold", color=COLORS["navy"])
    fig.tight_layout()
    return fig, summary


def score_qc(df: pd.DataFrame) -> tuple[plt.Figure, pd.DataFrame]:
    status = df.get("magnitude_score_status", pd.Series("UNKNOWN", index=df.index)).fillna("UNKNOWN").value_counts()
    unresolved = []
    if "magnitude_score_unresolved" in df:
        for value in df["magnitude_score_unresolved"].dropna().astype(str):
            if value == ".":
                continue
            for item in value.split(";"):
                unresolved.append(re.sub(r":.*$", "", item))
    unresolved_counts = pd.Series(unresolved).value_counts().reindex(MARKERS, fill_value=0)
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))
    axes[0].bar(status.index, status.values, color=[COLORS["green"], COLORS["terracotta"], COLORS["gray"]][:len(status)])
    axes[0].set(title="Magnitude-score status", xlabel="", ylabel="Samples")
    axes[0].tick_params(axis="x", rotation=25)
    clean_axes(axes[0], "y")
    axes[1].barh(unresolved_counts.index, unresolved_counts.values, color=COLORS["gold"])
    axes[1].set(title="Samples assigned zero by marker", xlabel="Samples", ylabel="")
    clean_axes(axes[1], "x")
    fig.tight_layout()
    qc = pd.DataFrame({"status": status.index, "count": status.values})
    return fig, qc


def dashboard(df: pd.DataFrame) -> plt.Figure:
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    order = ordered_genotypes(df["APOE_genotype"])
    counts = df["APOE_genotype"].value_counts().reindex(order, fill_value=0)
    axes[0, 0].bar(counts.index, counts.values, color=COLORS["green"])
    axes[0, 0].set_title("APOE genotypes")
    axes[0, 0].tick_params(axis="x", rotation=35, labelsize=8)
    clean_axes(axes[0, 0], "y")

    scores = numeric(df["overall_magnitude_score"]).dropna()
    axes[0, 1].hist(scores, bins=min(30, max(8, int(math.sqrt(max(len(scores), 1))))),
                    color=COLORS["terracotta"], edgecolor="white")
    axes[0, 1].set_title("Overall magnitude")
    axes[0, 1].set_xlabel("Score")
    clean_axes(axes[0, 1], "y")

    allele_values = []
    for allele in ("ε1", "ε2", "ε3", "ε4"):
        values = numeric(df[f"{allele}_dose"])
        allele_values.append(values.sum() / (2 * values.notna().sum()) * 100 if values.notna().sum() else 0)
    axes[1, 0].bar(["ε1", "ε2", "ε3", "ε4"], allele_values,
                   color=[COLORS["gold"], COLORS["green"], COLORS["brown"], COLORS["terracotta"]])
    axes[1, 0].set_title("Allele frequencies")
    axes[1, 0].yaxis.set_major_formatter(PercentFormatter(100))
    clean_axes(axes[1, 0], "y")

    prevalences = []
    for marker in MARKERS:
        col = f"{marker}_magnitude"
        values = numeric(df[col]).fillna(0) if col in df else pd.Series(0, index=df.index)
        prevalences.append((values > 0).mean() * 100)
    axes[1, 1].barh(MARKERS, prevalences, color=COLORS["gold"])
    axes[1, 1].set_title("Nonzero auxiliary magnitudes")
    axes[1, 1].xaxis.set_major_formatter(PercentFormatter(100))
    clean_axes(axes[1, 1], "x")

    fig.suptitle(f"APOE cohort overview · n={len(df):,}", fontsize=18, fontweight="bold", color=COLORS["navy"])
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return fig


def write_cohort_summary(df: pd.DataFrame, outdir: Path, prefix: str) -> None:
    scores = numeric(df["overall_magnitude_score"])
    apoe_called = ~df["APOE_genotype"].isin(["UNKNOWN", "."])
    ch = df.get("APOE_Christchurch_carrier", pd.Series("UNKNOWN", index=df.index)).astype(str)
    metrics = [
        ("samples", len(df)),
        ("APOE_called", int(apoe_called.sum())),
        ("APOE_call_rate", apoe_called.mean()),
        ("overall_score_available", int(scores.notna().sum())),
        ("overall_score_mean", scores.mean()),
        ("overall_score_median", scores.median()),
        ("overall_score_min", scores.min()),
        ("overall_score_max", scores.max()),
        ("Christchurch_YES", int((ch == "YES").sum())),
        ("Christchurch_NO_VARIANT_RECORD", int((ch == "NO_VARIANT_RECORD").sum())),
    ]
    pd.DataFrame(metrics, columns=["metric", "value"]).to_csv(
        outdir / f"{prefix}.cohort_summary.tsv", sep="\t", index=False
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="Per-sample TSV from apoe_genotyper.py")
    parser.add_argument("-o", "--output-dir", default="apoe_plots", help="Output directory")
    parser.add_argument("--prefix", default="apoe", help="Output filename prefix")
    parser.add_argument("--dpi", type=int, default=300, help="PNG resolution")
    args = parser.parse_args()

    set_style()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.input, sep="\t", dtype=str, keep_default_na=False)
    required = {"sample", "APOE_genotype", "overall_magnitude_score"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise SystemExit("Missing required column(s): " + ", ".join(missing))

    pdf_path = outdir / f"{args.prefix}.cohort_overview.pdf"
    with PdfPages(pdf_path) as pdf:
        save_figure(dashboard(df), outdir, f"{args.prefix}.01_dashboard", args.dpi, pdf)
        save_figure(genotype_prevalence(df), outdir, f"{args.prefix}.02_genotype_distribution", args.dpi, pdf)
        fig, allele_summary = allele_frequencies(df)
        save_figure(fig, outdir, f"{args.prefix}.03_allele_frequencies", args.dpi, pdf)
        save_figure(overall_distribution(df), outdir, f"{args.prefix}.04_overall_magnitude", args.dpi, pdf)
        save_figure(score_by_genotype(df), outdir, f"{args.prefix}.05_score_by_genotype", args.dpi, pdf)
        fig, genotype_summary = marker_genotype_distribution(df)
        save_figure(fig, outdir, f"{args.prefix}.06_marker_genotypes", args.dpi, pdf)
        fig, marker_summary = marker_contributions(df)
        save_figure(fig, outdir, f"{args.prefix}.07_marker_contributions", args.dpi, pdf)
        fig, qc_summary = score_qc(df)
        save_figure(fig, outdir, f"{args.prefix}.08_score_qc", args.dpi, pdf)

    allele_summary.to_csv(outdir / f"{args.prefix}.allele_summary.tsv", sep="\t", index=False)
    genotype_summary.to_csv(outdir / f"{args.prefix}.marker_genotype_summary.tsv", sep="\t", index=False)
    marker_summary.to_csv(outdir / f"{args.prefix}.marker_magnitude_summary.tsv", sep="\t", index=False)
    qc_summary.to_csv(outdir / f"{args.prefix}.score_qc_summary.tsv", sep="\t", index=False)
    df["APOE_genotype"].value_counts(dropna=False).rename_axis("APOE_genotype").reset_index(name="count").to_csv(
        outdir / f"{args.prefix}.APOE_genotype_counts.tsv", sep="\t", index=False
    )
    write_cohort_summary(df, outdir, args.prefix)
    print(f"Wrote APOE plots and summaries to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
