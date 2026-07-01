# Argiments to set in the bash command line for the launch of the code:
"""
python motif_disco_vs_ortho.py \
    --run-all \   # remembre that this is necessary to run the whole script
    --disco-dir 05_aggregate/05_totalscore \
    --ortho-dir 05_aggregate/05_totalscore_orthofinder \
    --outdir risultati_ortho_vs_disco
"""
# Script used to verify the quantity of information lost by the DISCO pipeline compared to the Orthofinder output, for each orthogroup and for each motif. It produces a summary table with statistics and graphs for each motif, and a batch summary table aggregating the results of all motifs.
# This code requires as imput only the files produced by the DISCO pipeline and the Orthofinder output (totalscore_MAXXXX_XX.tsv), in order to compare the scores of the two methods for each orthogroup and classify them into categories based on the modifications they underwent.

import argparse
import os
import sys
import warnings
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy import stats
from matplotlib.patches import Patch

warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)


# ─────────────────────────────────────────────────────────────
# 1. Loading and preparation
# ─────────────────────────────────────────────────────────────

# Loadings of tables file and aggregation by species name
def load_tables(disco_path: str, ortho_path: str):
    disco = pd.read_csv(disco_path, sep="\t", index_col=0)
    ortho = pd.read_csv(ortho_path, sep="\t", index_col=0)
    common_species = [c for c in disco.columns if c in ortho.columns]
    disco = disco[common_species]
    ortho = ortho[common_species]
    return disco, ortho


# Extraction of basename of OG
def parse_og_base(index_name: str) -> str:
    return index_name.rsplit("_", 1)[0]


def classify_orthogroups(disco: pd.DataFrame, ortho: pd.DataFrame):
    # Extrct the OG base name
    og_base_series = disco.index.to_series().apply(parse_og_base)
    disco_copy = disco.copy()
    disco_copy["_og_base"] = og_base_series.values

    species_cols = [c for c in disco.columns]

    # Number of forms per OG
    n_forms = disco_copy.groupby("_og_base").size().rename("n_forms_disco")

    # Aggregation of forms by Og base
    disco_agg = disco_copy.groupby("_og_base")[species_cols].sum(min_count=1)

    # Base OG present in DISCO
    disco_og_bases = set(n_forms.index)
    # Base OG present in ORTHO
    ortho_og_bases = set(ortho.index)

    # Common OG present in both DISCO and ORHO
    common = disco_og_bases & ortho_og_bases
    # OG presents only in ORTHO
    ortho_only = ortho_og_bases - disco_og_bases

    # Creaztion of categories: single_complete, split and ortho_only
    categories = {}
    for og in common:
        nf = n_forms[og]
        if nf == 1:
            categories[og] = "single_complete" # OG which didn't got split in DISCO
        else:
            categories[og] = "split" # OG that has more than one form in DISCO
    for og in ortho_only:
        categories[og] = "ortho_only" # OG absent in DISCO file

    og_meta = pd.DataFrame({
        "og_base": list(categories.keys()),
        "n_forms_disco": [n_forms.get(og, 0) for og in categories.keys()],
        "category": list(categories.values()),
    }).set_index("og_base")

    return og_meta, disco_agg


# ─────────────────────────────────────────────────────────────
# 2. Comparison of DISCO and ORTHO metrics per OG
# ─────────────────────────────────────────────────────────────

# Calculation of statistical value per OG
def summarize_row(row: pd.Series) -> dict:
    vals = row.dropna()
    if len(vals) == 0:
        return {"n_species": 0, "sum": np.nan, "mean": np.nan,
                "median": np.nan, "n_nonzero": 0, "frac_nonzero": np.nan}
    return {
        "n_species":   len(vals),
        "sum":         vals.sum(),
        "mean":        vals.mean(),
        "median":      vals.median(),
        "n_nonzero":   (vals > 0).sum(),
        "frac_nonzero": (vals > 0).mean(),
    }

# Construction of a dataframe with comparison of DISCO and ORTHO metrics per OG
def build_per_og_comparison(disco_agg: pd.DataFrame,
                             ortho: pd.DataFrame,
                             og_meta: pd.DataFrame) -> pd.DataFrame:
    common_ogs = disco_agg.index.intersection(ortho.index)
    records = []
    for og in common_ogs:
        d_stats = summarize_row(disco_agg.loc[og])
        o_stats = summarize_row(ortho.loc[og])
        cat = og_meta.loc[og, "category"] if og in og_meta.index else "unknown"
        nf  = og_meta.loc[og, "n_forms_disco"] if og in og_meta.index else 1

        ratio = (d_stats["sum"] / o_stats["sum"]
                 if (o_stats["sum"] and not np.isnan(o_stats["sum"]) and o_stats["sum"] > 0)
                 else np.nan)

        records.append({
            "og_base":            og,
            "category":           cat,
            "n_forms_disco":      nf,
            "disco_sum":          d_stats["sum"],
            "disco_mean":         d_stats["mean"],
            "disco_frac_nonzero": d_stats["frac_nonzero"],
            "ortho_sum":          o_stats["sum"],
            "ortho_mean":         o_stats["mean"],
            "ortho_frac_nonzero": o_stats["frac_nonzero"],
            "ratio_sum":          ratio,
            "delta_mean":         (o_stats["mean"] - d_stats["mean"])
                                  if not (np.isnan(o_stats["mean"]) or np.isnan(d_stats["mean"]))
                                  else np.nan,
        })

    df = pd.DataFrame(records).set_index("og_base")
    return df

# Compute summary statistics for each category of OG
def compute_summary_stats(per_og: pd.DataFrame, og_meta: pd.DataFrame,
                           ortho: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cat in ["single_complete", "split"]:
        sub = per_og[per_og["category"] == cat]
        if sub.empty:
            continue
        rows.append({
            "category":                cat,
            "n_og":                    len(sub),
            "median_ratio_sum":        sub["ratio_sum"].median(),
            "mean_ratio_sum":          sub["ratio_sum"].mean(),
            "frac_ratio_gt_0.5":       (sub["ratio_sum"] > 0.5).mean(),
            "median_disco_sum":        sub["disco_sum"].median(),
            "median_ortho_sum":        sub["ortho_sum"].median(),
            "median_delta_mean":       sub["delta_mean"].median(),
            "median_disco_frac_nz":    sub["disco_frac_nonzero"].median(),
            "median_ortho_frac_nz":    sub["ortho_frac_nonzero"].median(),
        })

    # DataFrame only for ortho_only categories
    ortho_only_ogs = og_meta[og_meta["category"] == "ortho_only"].index
    ortho_only_ogs = [og for og in ortho_only_ogs if og in ortho.index]
    if ortho_only_ogs:
        sums = ortho.loc[ortho_only_ogs].sum(axis=1, skipna=True)
        frac_nz = (ortho.loc[ortho_only_ogs] > 0).mean(axis=1)
        rows.append({
            "category":                "ortho_only",
            "n_og":                    len(ortho_only_ogs),
            "median_ratio_sum":        np.nan,   # non confrontabile
            "mean_ratio_sum":          np.nan,
            "frac_ratio_gt_0.5":       np.nan,
            "median_disco_sum":        np.nan,
            "median_ortho_sum":        sums.median(),
            "median_delta_mean":       np.nan,
            "median_disco_frac_nz":    np.nan,
            "median_ortho_frac_nz":    frac_nz.median(),
        })

    return pd.DataFrame(rows).set_index("category")


# ─────────────────────────────────────────────────────────────
# 3. Visualization and Graph
# ─────────────────────────────────────────────────────────────

COLORS = {
    "single_complete": "#4C8FBF",
    "split":           "#E8703A",
    "ortho_only":      "#6AAF6A",
}

CAT_LABELS = {
    "single_complete": "Single-complete",
    "split":           "Split (≥2 forme)",
    "ortho_only":      "Solo ORTHO",
}

# Rappresentation of ration DISCO/ORTHO
# When < 1 it means that DISCO has lonst scores compared to ORTHO
def plot_ratio_distribution(per_og: pd.DataFrame, outdir: Path, motif_name: str):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(f"Motivo: {motif_name} — Ratio punteggio DISCO/ORTHO", fontsize=13)

    cats = ["single_complete", "split"]
    data = [per_og[per_og["category"] == c]["ratio_sum"].dropna() for c in cats]
    labels = [CAT_LABELS[c] for c in cats]
    colors = [COLORS[c] for c in cats]

    # Boxplot
    ax = axes[0]
    bp = ax.boxplot(data, patch_artist=True, labels=labels,
                    medianprops=dict(color="black", linewidth=2))
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    ax.axhline(1.0, linestyle="--", color="gray", linewidth=1, label="ratio = 1 (nessuna perdita)")
    ax.axhline(0.5, linestyle=":", color="red", linewidth=1, label="ratio = 0.5")
    ax.set_ylabel("ratio (DISCO_sum / ORTHO_sum)")
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=8)
    ax.set_title("Boxplot per categoria")

    # Overlapped histograms
    ax2 = axes[1]
    for vals, color, label in zip(data, colors, labels):
        ax2.hist(vals.clip(0, 2), bins=40, alpha=0.5, color=color, label=label, density=True)
    ax2.axvline(1.0, linestyle="--", color="gray", linewidth=1)
    ax2.set_xlabel("ratio (clip a 2)")
    ax2.set_ylabel("densità")
    ax2.legend()
    ax2.set_title("Distribuzione ratio")

    plt.tight_layout()
    fig.savefig(outdir / "plot_ratio_distribution.png", dpi=150)
    plt.close(fig)

# Scatter plot of DISCO_sum vs ORTHO_sum, colored by category
def plot_scatter_disco_vs_ortho(per_og: pd.DataFrame, outdir: Path, motif_name: str):
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.set_title(f"Motivo: {motif_name} — DISCO vs ORTHO (sum)", fontsize=12)

    for cat in ["single_complete", "split"]:
        sub = per_og[per_og["category"] == cat]
        ax.scatter(sub["ortho_sum"], sub["disco_sum"],
                   c=COLORS[cat], label=CAT_LABELS[cat], alpha=0.4, s=15, linewidths=0)

    lim = max(per_og["ortho_sum"].max(), per_og["disco_sum"].max()) * 1.05
    ax.plot([0, lim], [0, lim], "k--", linewidth=1, label="y = x (parità)")
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel("ORTHO sum punteggio")
    ax.set_ylabel("DISCO sum punteggio")
    ax.legend(fontsize=9)
    plt.tight_layout()
    fig.savefig(outdir / "plot_scatter_disco_vs_ortho.png", dpi=150)
    plt.close(fig)

# Scatter plot which shows the species that has a score > 0 in DISCO and ORTHO
# I THINK THAT THIS MAY BE REMOVED FROM THE SCRIPT, ALSO BECAUSE IT HAS TO BE PRODUCED FOR EVERY MOTIFS
def plot_frac_nonzero(per_og: pd.DataFrame, outdir: Path, motif_name: str):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(f"Motivo: {motif_name} — Frazione specie con punteggio > 0", fontsize=13)

    for ax, cat in zip(axes, ["single_complete", "split"]):
        sub = per_og[per_og["category"] == cat].dropna(
            subset=["disco_frac_nonzero", "ortho_frac_nonzero"])
        ax.scatter(sub["ortho_frac_nonzero"], sub["disco_frac_nonzero"],
                   c=COLORS[cat], alpha=0.4, s=15, linewidths=0)
        ax.plot([0, 1], [0, 1], "k--", linewidth=1)
        ax.set_xlim(0, 1.05)
        ax.set_ylim(0, 1.05)
        ax.set_xlabel("ORTHO frac specie > 0")
        ax.set_ylabel("DISCO frac specie > 0")
        ax.set_title(CAT_LABELS[cat])
        # Percentuale punti sotto la diagonale (ORTHO > DISCO)
        below = (sub["ortho_frac_nonzero"] > sub["disco_frac_nonzero"]).mean()
        ax.text(0.05, 0.92, f"ORTHO > DISCO: {below:.1%}",
                transform=ax.transAxes, fontsize=9, color="darkred")

    plt.tight_layout()
    fig.savefig(outdir / "plot_frac_nonzero.png", dpi=150)
    plt.close(fig)

# Scatter plot for the distribution of the scores of the OG absent in DISCO
def plot_ortho_only_distribution(ortho: pd.DataFrame, og_meta: pd.DataFrame,
                                  outdir: Path, motif_name: str):
    ortho_only_ogs = og_meta[og_meta["category"] == "ortho_only"].index
    ortho_only_ogs = [og for og in ortho_only_ogs if og in ortho.index]
    if not ortho_only_ogs:
        return

    sums = ortho.loc[ortho_only_ogs].sum(axis=1, skipna=True)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(sums[sums > 0], bins=50, color=COLORS["ortho_only"], alpha=0.8)
    ax.set_xlabel("Somma punteggio ORTHO (tra specie)")
    ax.set_ylabel("Numero OG")
    ax.set_title(f"Motivo: {motif_name} — OG solo in ORTHO (n={len(ortho_only_ogs)})\n"
                 f"Distribuzione punteggi")
    ax.text(0.97, 0.95,
            f"OG con punteggio > 0: {(sums > 0).sum()} / {len(sums)}",
            transform=ax.transAxes, ha="right", va="top", fontsize=9)
    plt.tight_layout()
    fig.savefig(outdir / "plot_ortho_only_distribution.png", dpi=150)
    plt.close(fig)


# ─────────────────────────────────────────────────────────────
# 4. Processing single OG
# ─────────────────────────────────────────────────────────────

# Process a single motif and save results in outdir
def process_motif(disco_path: str, ortho_path: str, outdir: str, save_plots: bool = True):

    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    # Retrive the motif name from the disco_path
    motif_name = Path(disco_path).stem.replace("totalscore_", "")

    # 1. Loading datas
    print(f"[{motif_name}] Caricamento dati...")
    disco, ortho = load_tables(disco_path, ortho_path)
    print(f"  DISCO: {disco.shape[0]} righe × {disco.shape[1]} specie")
    print(f"  ORTHO: {ortho.shape[0]} OG × {ortho.shape[1]} specie")

    # 2. Classification of OG
    print(f"[{motif_name}] Classificazione ortogruppi...")
    og_meta, disco_agg = classify_orthogroups(disco, ortho)
    cat_counts = og_meta["category"].value_counts()
    for cat, n in cat_counts.items():
        print(f"  {CAT_LABELS.get(cat, cat)}: {n} OG")

    # 3. Comparative metrics per OG
    print(f"[{motif_name}] Calcolo metriche comparative...")
    per_og = build_per_og_comparison(disco_agg, ortho, og_meta)

    # Ratio single_complete/split: Mann-Whitney U test
    ratio_sc = per_og[per_og["category"] == "single_complete"]["ratio_sum"].dropna()
    ratio_sp = per_og[per_og["category"] == "split"]["ratio_sum"].dropna()
    if len(ratio_sc) > 5 and len(ratio_sp) > 5:
        stat, pval = stats.mannwhitneyu(ratio_sc, ratio_sp, alternative="two-sided")
        print(f"  Mann-Whitney ratio (single vs split): U={stat:.0f}, p={pval:.4f}")
    else:
        pval = np.nan
        print("  Campione troppo piccolo per test statistico")

    # Aggregation of summary statistics per category
    summary = compute_summary_stats(per_og, og_meta, ortho)

    # 5. Saving files
    summary.to_csv(out / "summary_stats.tsv", sep="\t", float_format="%.4f") # It contains the summary statistics for each category of OG
    per_og.to_csv(out / "per_og_comparison.tsv", sep="\t", float_format="%.4f") # It contains the comparison of DISCO and ORTHO metrics for each OG
    og_meta.to_csv(out / "og_classification.tsv", sep="\t") # It contains the classification of OG into categories
    print(f"[{motif_name}] Tabelle salvate in: {out}")

    # 6. Saving plots only if save_plots is True
    if save_plots:
        plots_dir = out / "plots"
        plots_dir.mkdir(exist_ok=True)
        print(f"[{motif_name}] Generazione grafici...")
        plot_ratio_distribution(per_og, plots_dir, motif_name)
        plot_scatter_disco_vs_ortho(per_og, plots_dir, motif_name)
        plot_frac_nonzero(per_og, plots_dir, motif_name)
        plot_ortho_only_distribution(ortho, og_meta, plots_dir, motif_name)
        print(f"[{motif_name}] Grafici salvati in: {plots_dir}")

    return {"summary_stats": summary, "per_og_comparison": per_og, "og_meta": og_meta}


# ─────────────────────────────────────────────────────────────
# 5. Comparison between every single motif in one run
# ─────────────────────────────────────────────────────────────

# Aggregation of several results from every motif processed with process_motif() into a single summary table for cross-motif comparison
def run_all(disco_dir: str, ortho_dir: str, outdir: str, save_per_motif_plots: bool = False):

    disco_dir = Path(disco_dir)
    ortho_dir = Path(ortho_dir)
    out_root  = Path(outdir)
    out_root.mkdir(parents=True, exist_ok=True)

    # Searches for available motif file
    disco_files = sorted(disco_dir.glob("totalscore_*.tsv"))
    if not disco_files:
        print(f"ERRORE: nessun file totalscore_*.tsv trovato in {disco_dir}")
        sys.exit(1)

    print(f"Trovati {len(disco_files)} motivi in {disco_dir}")

    skipped = []
    processed = []

    for i, disco_path in enumerate(disco_files, 1):
        motif_name = disco_path.stem.replace("totalscore_", "")
        ortho_path = ortho_dir / f"totalscore_{motif_name}.tsv"

        if not ortho_path.exists():
            print(f"[{i}/{len(disco_files)}] SKIP {motif_name} — file ORTHO mancante: {ortho_path}")
            skipped.append(motif_name)
            continue

        motif_outdir = out_root / motif_name
        print(f"\n[{i}/{len(disco_files)}] ── {motif_name} ──────────────────")
        try:
            process_motif(str(disco_path), str(ortho_path),
                          str(motif_outdir), save_plots=save_per_motif_plots)
            processed.append(motif_name)
        except Exception as e:
            print(f"  ERRORE durante il processing di {motif_name}: {e}")
            skipped.append(motif_name)

    print(f"\n{'='*60}")
    print(f"Processati: {len(processed)}  |  Saltati/Errori: {len(skipped)}")
    if skipped:
        print(f"Saltati: {skipped}")

    # Generate a summary table for all motifs processed
    print("\nGenerazione batch_summary.tsv...")
    batch_tsv = out_root / "batch_summary.tsv"
    df_batch = batch_summary(str(out_root), str(batch_tsv))

    if df_batch is not None:
        print("Generazione batch_overview.png...")
        plot_batch_summary(str(batch_tsv), str(out_root / "batch_overview.png"))

    print(f"\nDone. Output globale in: {out_root}")


# ─────────────────────────────────────────────────────────────
# 6. Setting command line for bash launch
# ─────────────────────────────────────────────────────────────

# Aggregation of every summary.tsv of all motifs already processed into a single table for cross-motif comparison
def batch_summary(results_root: str, output_file: str = "batch_summary.tsv"):
    """
    Aggrega i summary_stats.tsv di tutti i motivi già processati
    in un'unica tabella per confronto cross-motivo.
    """
    root = Path(results_root)
    records = []
    for motif_dir in sorted(root.iterdir()):
        if not motif_dir.is_dir():
            continue
        summary_file = motif_dir / "summary_stats.tsv"
        if not summary_file.exists():
            continue
        df = pd.read_csv(summary_file, sep="\t", index_col=0)
        df["motif"] = motif_dir.name
        records.append(df.reset_index())

    if not records:
        print("Nessun risultato trovato in", results_root)
        return None

    combined = pd.concat(records, ignore_index=True)
    combined.to_csv(output_file, sep="\t", index=False, float_format="%.4f")
    print(f"Batch summary salvato in: {output_file} ({len(records)} motivi)")
    return combined

#─────────────────────────────────────────────────────────────
# 7. Summary chart
#─────────────────────────────────────────────────────────────

# Generats a plot divided into four different parts representing every motifs and sorting them by decreasing ratio
def plot_batch_summary(batch_tsv: str, output_png: str = "batch_overview.png"):

    df = pd.read_csv(batch_tsv, sep="\t")

    def get_cat(df, cat):
        return df[df["category"] == cat].copy()

    sc = get_cat(df, "single_complete").set_index("motif")
    sp = get_cat(df, "split").set_index("motif")

    if sc.empty:
        raise ValueError("Nessun dato 'single_complete' nella batch summary.")

    # Sorting of motif based on their ratio
    order = sc["median_ratio_sum"].sort_values().index.tolist()
    extra = [m for m in sp.index if m not in order]
    order = order + extra

    n = len(order)
    x = np.arange(n)

    C_SC  = COLORS["single_complete"]
    C_SP  = COLORS["split"]
    ALPHA = 0.8

    fig, axes = plt.subplots(2, 2, figsize=(18, 10))
    fig.suptitle(
        "Confronto DISCO vs Orthofinder — tutti i motivi\n"
        "(ordinati per perdita crescente di punteggio: sinistra = piu penalizzato)",
        fontsize=14, fontweight="bold", y=1.01
    )

    # First plot: Ratio single_complete
    ax = axes[0, 0]
    vals_sc = [sc.loc[m, "median_ratio_sum"] if m in sc.index else np.nan for m in order]
    ax.bar(x, vals_sc, color=C_SC, alpha=ALPHA, width=0.8)
    ax.axhline(1.0, color="gray", linestyle="--", linewidth=1, label="ratio = 1 (nessuna perdita)")
    ax.axhline(0.5, color="red",  linestyle=":",  linewidth=1, label="ratio = 0.5")
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(0, max(1.2, np.nanmax(vals_sc) * 1.05))
    ax.set_ylabel("Ratio mediano DISCO/ORTHO")
    ax.set_title(f"[A] Single-complete  (OG mediana: {int(sc['n_og'].median())})",
                 color=C_SC, fontweight="bold")
    ax.set_xticks([])
    ax.legend(fontsize=8)
    n_loss = sum(v < 0.9 for v in vals_sc if not np.isnan(v))
    ax.text(0.02, 0.05, f"Motivi con ratio < 0.9: {n_loss}/{n}",
            transform=ax.transAxes, fontsize=9, color="darkred")

    # Second plot: Ratio split
    ax = axes[0, 1]
    vals_sp = [sp.loc[m, "median_ratio_sum"] if m in sp.index else np.nan for m in order]
    ax.bar(x, vals_sp, color=C_SP, alpha=ALPHA, width=0.8)
    ax.axhline(1.0, color="gray", linestyle="--", linewidth=1)
    ax.axhline(0.5, color="red",  linestyle=":",  linewidth=1)
    ax.set_xlim(-0.5, n - 0.5)
    valid_sp = [v for v in vals_sp if v is not None and not np.isnan(v)]
    ax.set_ylim(0, max(1.2, (max(valid_sp) * 1.05) if valid_sp else 1.2))
    ax.set_ylabel("Ratio mediano DISCO/ORTHO")
    n_sp_og_med = int(sp["n_og"].median()) if not sp.empty else 0
    ax.set_title(f"[B] Split >=2 forme  (OG mediana: {n_sp_og_med})",
                 color=C_SP, fontweight="bold")
    ax.set_xticks([])
    n_loss_sp = sum(v < 0.9 for v in vals_sp if not np.isnan(v))
    ax.text(0.02, 0.05, f"Motivi con ratio < 0.9: {n_loss_sp}/{len(valid_sp)}",
            transform=ax.transAxes, fontsize=9, color="darkred")

    # Third plot: Delta frac_nonzero single_complete
    ax = axes[1, 0]
    delta_sc = [
        (sc.loc[m, "median_ortho_frac_nz"] - sc.loc[m, "median_disco_frac_nz"])
        if m in sc.index else np.nan
        for m in order
    ]
    bar_colors_sc = [
        C_SC if (v is not None and not np.isnan(v) and v >= 0) else "#cc3333"
        for v in delta_sc
    ]
    ax.bar(x, delta_sc, color=bar_colors_sc, alpha=ALPHA, width=0.8)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylabel("Delta frac specie > 0  (ORTHO - DISCO)")
    ax.set_title("[C] Single-complete — perdita copertura specie",
                 color=C_SC, fontweight="bold")
    # I have commented out these lines because the labels under the bars where too short
    # ax.set_xticks(x)
    # ax.set_xticklabels(order, rotation=90, fontsize=5)
    # ax.set_xlabel("Motivo")
    n_loss_frac_sc = sum(v > 0 for v in delta_sc if not np.isnan(v))
    ax.text(0.02, 0.95, f"Motivi dove ORTHO > DISCO: {n_loss_frac_sc}/{n}",
            transform=ax.transAxes, fontsize=9, color="darkred", va="top")

    # Fourth plot: Delta frac_nonzero split
    ax = axes[1, 1]
    delta_sp = [
        (sp.loc[m, "median_ortho_frac_nz"] - sp.loc[m, "median_disco_frac_nz"])
        if m in sp.index else np.nan
        for m in order
    ]
    bar_colors_sp = [
        C_SP if (v is not None and not np.isnan(v) and v >= 0) else "#cc3333"
        for v in delta_sp
    ]
    ax.bar(x, delta_sp, color=bar_colors_sp, alpha=ALPHA, width=0.8)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylabel("Delta frac specie > 0  (ORTHO - DISCO)")
    ax.set_title("[D] Split >=2 forme — perdita copertura specie",
                 color=C_SP, fontweight="bold")
    # I have commented out these lines because the labels under the bars where too short
    # ax.set_xticks(x)
    # ax.set_xticklabels(order, rotation=90, fontsize=5)
    # ax.set_xlabel("Motivo")
    n_loss_frac_sp = sum(v > 0 for v in delta_sp if not np.isnan(v))
    n_sp_valid = sum(1 for v in delta_sp if not np.isnan(v))
    ax.text(0.02, 0.95, f"Motivi dove ORTHO > DISCO: {n_loss_frac_sp}/{n_sp_valid}",
            transform=ax.transAxes, fontsize=9, color="darkred", va="top")

    # Legend for the four plots
    legend_els = [
        Patch(facecolor=C_SC,       alpha=0.8, label="Single-complete"),
        Patch(facecolor=C_SP,       alpha=0.8, label="Split >=2 forme"),
        Patch(facecolor="#cc3333",  alpha=0.8, label="Perdita (DISCO < ORTHO)"),
    ]
    fig.legend(handles=legend_els, loc="lower center", ncol=3,
               fontsize=10, bbox_to_anchor=(0.5, -0.02))

    plt.tight_layout()
    fig.savefig(output_png, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Grafico riassuntivo salvato in: {output_png}")

#─────────────────────────────────────────────────────────────
# 8. Command line to set the arguments
#─────────────────────────────────────────────────────────────
# Command line interface for 
def parse_args():
    p = argparse.ArgumentParser(
        description="Confronto distribuzione punteggi dei motivi DISCO vs Orthofinder.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    # Defining the output directory
    p.add_argument("--outdir", help="Cartella di output")

    # Arguments to run this script for every motif in a single run
    p.add_argument("--run-all", action="store_true",
                   help="Processa tutti i motivi trovati in --disco-dir e --ortho-dir")
    p.add_argument("--disco-dir", default="05_aggregate/05_totalscore",
                   help="Cartella DISCO  (default: 05_aggregate/05_totalscore)")
    p.add_argument("--ortho-dir", default="05_aggregate/05_totalscore_orthofinder",
                   help="Cartella ORTHO  (default: 05_aggregate/05_totalscore_orthofinder)")
    p.add_argument("--save-per-motif-plots", action="store_true",
                   help="In --run-all, genera anche i grafici per singolo motivo "
                        "(più lento; di default solo le tabelle + grafico globale)")

    # Aggregation
    p.add_argument("--batch-summary", metavar="ROOT",
                   help="Aggrega i summary_stats.tsv in ROOT → batch_summary.tsv in --outdir")

    # Graphic overview
    p.add_argument("--plot-batch", metavar="BATCH_TSV",
                   help="Genera batch_overview.png da batch_summary.tsv")
    p.add_argument("--output-png", default="batch_overview.png",
                   help="Nome PNG output per --plot-batch (default: batch_overview.png)")

    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()

    if args.run_all:
        if not args.outdir:
            print("ERRORE: --outdir e obbligatorio con --run-all")
            sys.exit(1)
        run_all(args.disco_dir, args.ortho_dir, args.outdir,
                save_per_motif_plots=args.save_per_motif_plots)

    elif args.plot_batch:
        plot_batch_summary(args.plot_batch, args.output_png)

    elif args.batch_summary:
        out_tsv = Path(args.outdir) / "batch_summary.tsv" if args.outdir else "batch_summary.tsv"
        batch_summary(args.batch_summary, str(out_tsv))

    else:
        print(__doc__)