# Arguments to set in the bash command line for the launch of the code:
# --run-all is required to process all motifs in a single run.
#
# python motif_disco_vs_ortho.py \
#     --run-all \
#     --disco-dir 02_totalscore \
#     --ortho-dir 05_totalscore_orthofinder \
#     --outdir risultati_ortho_vs_disco
# Script used to verify the quantity of information lost by the DISCO pipeline compared to the Orthofinder output, for each orthogroup and for each motif. It produces a summary table with statistics and graphs for each motif, and a batch summary table aggregating the results of all motifs.
# This code requires as imput only the files produced by the DISCO pipeline and the Orthofinder output (totalscore_MAXXXX_XX.tsv), in order to compare the scores of the two methods for each orthogroup and classify them into categories based on the modifications they underwent.

import argparse
import os
import re
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
def load_tables(disco_path: str, ortho_path: str, min_species_ortho: int = 90):
    """
    Carica le tabelle DISCO e ORTHO e applica uno screening preliminare:
    gli OG di ORTHO con un numero di specie <= min_species_ortho vengono
    scartati prima della classificazione perché troppo poco informativi.
    Il filtro conta le specie con valore non-NaN nella riga ORTHO
    (ovvero specie in cui l'OG possiede almeno un gene), non il punteggio.

    Returns
    -------
    disco, ortho, n_og_raw, n_og_filtered
        n_og_raw      : OG in ORTHO prima del filtro
        n_og_filtered : OG eliminati dal filtro (<=min_species_ortho specie)
    """
    disco = pd.read_csv(disco_path, sep="\t", index_col=0)
    ortho = pd.read_csv(ortho_path, sep="\t", index_col=0)
    common_species = [c for c in disco.columns if c in ortho.columns]
    disco = disco[common_species]
    ortho = ortho[common_species]

    n_og_raw = ortho.shape[0]
    ortho_species_count = ortho.notna().sum(axis=1)
    ortho = ortho[ortho_species_count > min_species_ortho]
    n_og_filtered = n_og_raw - ortho.shape[0]
    if n_og_filtered > 0:
        print(f"  Filtro min-species-ortho (> {min_species_ortho} specie): "
              f"rimossi {n_og_filtered}/{n_og_raw} OG ORTHO "
              f"({ortho.shape[0]} rimasti)")

    return disco, ortho, n_og_raw, n_og_filtered


# Extraction of basename of OG.
# Matches a trailing underscore followed by one or more digits (e.g. _00, _01, _002).
# Falls back to the original name if the pattern is not found, so OG names
# without a numeric suffix are passed through unchanged rather than silently
# truncated at the last underscore.
_OG_SUFFIX_RE = re.compile(r'^(.+)_\d+$')

def parse_og_base(index_name: str) -> str:
    m = _OG_SUFFIX_RE.match(index_name)
    return m.group(1) if m else index_name


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

    # Keep individual forms (one row per OG_XX) for per-form comparison of split OGs.
    # For each split OG we need to compare ORTHO against each form separately,
    # and then compute coverage on the union of all forms — not on the sum.
    disco_forms = disco_copy.drop(columns=["_og_base"])  # original rows, indexed by OG_XX

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

    return og_meta, disco_agg, disco_forms


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
                             og_meta: pd.DataFrame,
                             disco_forms: pd.DataFrame,
                             good_og_species_frac: float = 0.6) -> pd.DataFrame:
    """
    Confronto RISTRETTO: vengono analizzati SOLO gli OG di ORTHO che hanno
    una corrispondenza in DISCO (single_complete o split).
    Gli OG ortho_only sono esclusi dal confronto metrico ma restano
    tracciati in og_meta per le statistiche di riepilogo.

    Report OG stampato a schermo per ogni motivo:
      - OG ORTHO post-filtro
      - OG ORTHO con corrispondenza in DISCO (base confrontato)
      - OG ORTHO senza corrispondenza (ortho_only, esclusi dal confronto)
      - OG DISCO totali (righe originali, incluse forme _00/_01/...)
      - OG DISCO distinti per base (dopo aggregazione)
      - Forme aggiuntive generate dallo split
    """
    # Solo OG presenti in entrambi (esclude ortho_only dal confronto)
    common_ogs = disco_agg.index.intersection(ortho.index)
    total_species = ortho.shape[1]

    # ── Report conteggi OG ────────────────────────────────────────────────
    n_ortho_total   = ortho.shape[0]
    n_common        = len(common_ogs)
    n_ortho_only    = n_ortho_total - n_common

    # Forme DISCO: tutte le righe originali (OG_00, OG_01, ...)
    n_disco_forms   = disco_forms.shape[0]
    # OG DISCO distinti per base (dopo aggregazione)
    n_disco_bases   = disco_agg.shape[0]
    # Forme aggiuntive = righe in più rispetto alle basi (ogni split aggiunge forme)
    n_extra_forms   = n_disco_forms - n_disco_bases

    print(f"  ── Report OG ──────────────────────────────────────")
    print(f"  OG ORTHO (post-filtro specie)          : {n_ortho_total:>7}")
    print(f"  OG ORTHO con match in DISCO (confronto): {n_common:>7}")
    print(f"  OG ORTHO senza match (ortho_only, skip): {n_ortho_only:>7}")
    print(f"  OG DISCO (righe totali, incl. forme)   : {n_disco_forms:>7}")
    print(f"  OG DISCO distinti per base             : {n_disco_bases:>7}")
    print(f"  Forme aggiuntive da split (_01, _02...): {n_extra_forms:>7}")
    print(f"  ───────────────────────────────────────────────────")

    records = []
    for og in common_ogs:
        o_stats = summarize_row(ortho.loc[og])
        cat = og_meta.loc[og, "category"] if og in og_meta.index else "unknown"
        nf  = og_meta.loc[og, "n_forms_disco"] if og in og_meta.index else 1

        if cat == "single_complete":
            d_stats = summarize_row(disco_agg.loc[og])
            disco_frac_nz = d_stats["frac_nonzero"]
            ratio = (d_stats["sum"] / o_stats["sum"]
                     if (o_stats["sum"] and not np.isnan(o_stats["sum"]) and o_stats["sum"] > 0)
                     else np.nan)
            delta_mean = (o_stats["mean"] - d_stats["mean"]
                          if not (np.isnan(o_stats["mean"]) or np.isnan(d_stats["mean"]))
                          else np.nan)
            disco_sum  = d_stats["sum"]
            disco_mean = d_stats["mean"]

        else:
            # Split: OR tra le forme per frac_nonzero
            form_rows = disco_forms[disco_forms.index.str.startswith(og + "_")]
            union_nonzero = (form_rows > 0).any(axis=0)
            species_present = form_rows.notna().any(axis=0)
            disco_frac_nz = (union_nonzero[species_present].mean()
                             if species_present.sum() > 0 else np.nan)

            d_stats = summarize_row(disco_agg.loc[og])
            ratio = (d_stats["sum"] / o_stats["sum"]
                     if (o_stats["sum"] and not np.isnan(o_stats["sum"]) and o_stats["sum"] > 0)
                     else np.nan)
            delta_mean = (o_stats["mean"] - d_stats["mean"]
                          if not (np.isnan(o_stats["mean"]) or np.isnan(d_stats["mean"]))
                          else np.nan)
            disco_sum  = d_stats["sum"]
            disco_mean = d_stats["mean"]

        records.append({
            "og_base":            og,
            "category":           cat,
            "n_forms_disco":      nf,
            "disco_sum":          disco_sum,
            "disco_mean":         disco_mean,
            "disco_frac_nonzero": disco_frac_nz,
            "ortho_sum":          o_stats["sum"],
            "ortho_mean":         o_stats["mean"],
            "ortho_frac_nonzero": o_stats["frac_nonzero"],
            "ratio_sum":          ratio,
            "delta_mean":         delta_mean,
        })

    df = pd.DataFrame(records).set_index("og_base")
    return df

# Compute summary statistics for each category of OG
def compute_summary_stats(per_og: pd.DataFrame, og_meta: pd.DataFrame,
                           ortho: pd.DataFrame,
                           n_og_raw: int = 0,
                           n_og_filtered: int = 0,
                           good_species_frac: float = 0.60,
                           mw_pval: float = np.nan,
                           n_disco_forms: int = 0,
                           n_disco_bases: int = 0,
                           n_extra_forms: int = 0) -> pd.DataFrame:
    """
    Parameters added vs original:
      n_og_raw        : OG totali in ORTHO prima del filtro specie
      n_og_filtered   : OG eliminati dal filtro (<=min_species_ortho specie)
      good_species_frac: soglia per candidati paralog-dropped (default 0.60)

    Adds columns to the 'ortho_only' row:
      n_candidate_paralog_dropped : OG ortho_only con copertura >= 60% specie
    And a global row 'filter_info' with n_og_raw / n_og_filtered.
    """
    total_species = ortho.shape[1]
    rows = []
    for cat in ["single_complete", "split"]:
        sub = per_og[per_og["category"] == cat]
        if sub.empty:
            continue
        rows.append({
            "category":                    cat,
            "n_og":                        len(sub),
            "n_og_raw_ortho":              n_og_raw,
            "n_og_filtered":               n_og_filtered,
            "n_og_ortho_only":             og_meta[og_meta["category"] == "ortho_only"].shape[0],
            "n_og_disco_forms":            n_disco_forms,
            "n_og_disco_bases":            n_disco_bases,
            "n_og_extra_forms":            n_extra_forms,
            "n_candidate_paralog_dropped": np.nan,
            "mw_pval_sc_vs_sp":            mw_pval,   # Mann-Whitney p-value
            "median_ratio_sum":            sub["ratio_sum"].median(),
            "mean_ratio_sum":              sub["ratio_sum"].mean(),
            "frac_ratio_gt_0.5":           (sub["ratio_sum"] > 0.5).mean(),
            "median_disco_sum":            sub["disco_sum"].median(),
            "median_ortho_sum":            sub["ortho_sum"].median(),
            "median_delta_mean":           sub["delta_mean"].median(),
            "median_disco_frac_nz":        sub["disco_frac_nonzero"].median(),
            "median_ortho_frac_nz":        sub["ortho_frac_nonzero"].median(),
        })

    # OG present only in ORTHO
    ortho_only_ogs = og_meta[og_meta["category"] == "ortho_only"].index
    ortho_only_ogs = [og for og in ortho_only_ogs if og in ortho.index]
    if ortho_only_ogs:
        sums = ortho.loc[ortho_only_ogs].sum(axis=1, skipna=True)
        frac_nz = (ortho.loc[ortho_only_ogs] > 0).mean(axis=1)
        # Count candidates: OG with species coverage >= good_species_frac
        n_species_per_og = ortho.loc[ortho_only_ogs].notna().sum(axis=1)
        species_frac_per_og = n_species_per_og / total_species if total_species > 0 else 0
        n_candidates = int((species_frac_per_og >= good_species_frac).sum())
        rows.append({
            "category":                    "ortho_only",
            "n_og":                        len(ortho_only_ogs),
            "n_og_raw_ortho":              n_og_raw,
            "n_og_filtered":               n_og_filtered,
            "n_candidate_paralog_dropped": n_candidates,
            "median_ratio_sum":            np.nan,
            "mean_ratio_sum":              np.nan,
            "frac_ratio_gt_0.5":           np.nan,
            "median_disco_sum":            np.nan,
            "median_ortho_sum":            sums.median(),
            "median_delta_mean":           np.nan,
            "median_disco_frac_nz":        np.nan,
            "median_ortho_frac_nz":        frac_nz.median(),
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

# Scatter plot which shows the species that has a score > 0 in DISCO and ORTHO
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
    disco, ortho, n_og_raw, n_og_filtered = load_tables(disco_path, ortho_path)
    print(f"  DISCO: {disco.shape[0]} righe × {disco.shape[1]} specie")
    print(f"  ORTHO: {ortho.shape[0]} OG × {ortho.shape[1]} specie "
          f"(raw={n_og_raw}, filtrati={n_og_filtered})")

    # 2. Classification of OG
    print(f"[{motif_name}] Classificazione ortogruppi...")
    og_meta, disco_agg, disco_forms = classify_orthogroups(disco, ortho)
    cat_counts = og_meta["category"].value_counts()
    for cat, n in cat_counts.items():
        print(f"  {CAT_LABELS.get(cat, cat)}: {n} OG")

    # 3. Comparative metrics per OG
    print(f"[{motif_name}] Calcolo metriche comparative...")
    per_og = build_per_og_comparison(disco_agg, ortho, og_meta, disco_forms)

    # Ratio single_complete/split: Mann-Whitney U test
    ratio_sc = per_og[per_og["category"] == "single_complete"]["ratio_sum"].dropna()
    ratio_sp = per_og[per_og["category"] == "split"]["ratio_sum"].dropna()
    if len(ratio_sc) > 5 and len(ratio_sp) > 5:
        stat, pval = stats.mannwhitneyu(ratio_sc, ratio_sp, alternative="two-sided")
        print(f"  Mann-Whitney ratio (single vs split): U={stat:.0f}, p={pval:.4f}")
    else:
        pval = np.nan
        print("  Campione troppo piccolo per test statistico")

    # Conteggi OG DISCO per il report (calcolati dalla stessa sorgente
    # usata in build_per_og_comparison, ma qui ricalcolati per passarli al summary)
    n_disco_forms_cnt = disco_forms.shape[0]
    n_disco_bases_cnt = disco_agg.shape[0]
    n_extra_forms_cnt = n_disco_forms_cnt - n_disco_bases_cnt

    # Aggregation of summary statistics per category
    summary = compute_summary_stats(per_og, og_meta, ortho,
                                    n_og_raw=n_og_raw,
                                    n_og_filtered=n_og_filtered,
                                    mw_pval=pval,
                                    n_disco_forms=n_disco_forms_cnt,
                                    n_disco_bases=n_disco_bases_cnt,
                                    n_extra_forms=n_extra_forms_cnt)

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
        # Named *_quick to distinguish from the higher-quality ggplot2 version
        # produced by the companion R script (prova_rappresentazione_grafica.R).
        print("Generazione batch_overview_quick.png...")
        plot_batch_summary(str(batch_tsv), str(out_root / "batch_overview_quick.png"))

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
# 7. Summary chart — layout rivisto
#─────────────────────────────────────────────────────────────

def plot_batch_summary(batch_tsv: str, output_png: str = "batch_overview.png",
                       good_species_frac: float = 0.60):
    # Struttura batch_overview.png — 3 pannelli:
    #   RIGA SUPERIORE (35% altezza) — due pannelli piccoli e uguali:
    #     [A] sx: riepilogo globale aggregato + motivi MW p<0.05
    #     [B] dx: barplot OG candidati paralog-dropped per motivo
    #   RIGA INFERIORE (65% altezza) — pannello largo:
    #     [C]:    Delta frac-nonzero split, grande e leggibile

    df = pd.read_csv(batch_tsv, sep="\t")

    def get_cat(df_, cat):
        return df_[df_["category"] == cat].copy()

    sc = get_cat(df, "single_complete").set_index("motif")
    sp = get_cat(df, "split").set_index("motif")
    oo = get_cat(df, "ortho_only").set_index("motif")

    if sc.empty:
        raise ValueError("Nessun dato 'single_complete' nella batch summary.")

    all_motifs = sorted(sc.index.union(sp.index).union(oo.index))
    n = len(all_motifs)
    x = np.arange(n)

    C_SC   = COLORS["single_complete"]
    C_SP   = COLORS["split"]
    C_CAND = "#9B59B6"
    ALPHA  = 0.85

    # ── Layout ───────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(max(20, n * 0.07 + 6), 18))
    gs = fig.add_gridspec(
        2, 2,
        height_ratios=[0.55, 1.0],    # riga alta ~35%, bassa ~65%
        hspace=0.38, wspace=0.30,
        left=0.06, right=0.97, top=0.94, bottom=0.07
    )
    ax_info  = fig.add_subplot(gs[0, 0])   # [A] info globale
    ax_cand  = fig.add_subplot(gs[0, 1])   # [B] candidati paralog-drop
    ax_delta = fig.add_subplot(gs[1, :])   # [C] delta frac-nz, largo

    fig.suptitle(
        "Confronto DISCO vs Orthofinder — riepilogo globale tutti i motivi",
        fontsize=14, fontweight="bold"
    )

    # ─────────────────────────────────────────────────────────────────────────
    # [A] Riepilogo globale aggregato + motivi con MW p<0.05
    # ─────────────────────────────────────────────────────────────────────────
    ax_info.axis("off")

    def uniform_or_median(ser, col):
        vals = ser[col].dropna() if col in ser.columns else pd.Series(dtype=float)
        if vals.empty:
            return "n.d."
        uvals = vals.unique()
        if len(uvals) == 1:
            return str(int(uvals[0]))
        return f"{int(vals.median())} (mediana)"

    n_raw_val  = uniform_or_median(sc, "n_og_raw_ortho")
    n_filt_val = uniform_or_median(sc, "n_og_filtered")
    n_sc_val   = uniform_or_median(sc, "n_og")
    n_sp_val   = uniform_or_median(sp, "n_og")
    n_oo_val   = uniform_or_median(oo, "n_og")

    n_cand_tot    = int(oo["n_candidate_paralog_dropped"].sum()) if "n_candidate_paralog_dropped" in oo.columns else 0
    n_motifs_cand = int((oo["n_candidate_paralog_dropped"] > 0).sum()) if "n_candidate_paralog_dropped" in oo.columns else 0

    # Motivi con MW p<0.05
    mw_col = "mw_pval_sc_vs_sp"
    mw_sig_list = []
    if mw_col in sc.columns:
        sig_mask = sc[mw_col].notna() & (sc[mw_col] < 0.05)
        mw_sig_list = sorted(sc.index[sig_mask].tolist())

    lines = [
        ("N. motivi analizzati",             str(n)),
        ("OG ORTHO per motivo (raw)",        n_raw_val),
        ("OG filtrati (<=2 spp)",            n_filt_val),
        ("Single-complete per motivo",       n_sc_val),
        ("Split >=2 forme per motivo",       n_sp_val),
        ("Ortho-only per motivo",            n_oo_val),
        (f"Candidati paralog-drop totale (>={int(good_species_frac*100)}% spp)", str(n_cand_tot)),
        ("  Motivi con almeno 1 candidato",  f"{n_motifs_cand}/{n}"),
    ]

    tbl = ax_info.table(
        cellText  = [[k, v] for k, v in lines],
        colLabels = ["Metrica", "Valore"],
        loc       = "upper center",
        cellLoc   = "left",
        bbox      = [0.0, 0.42, 1.0, 0.56],
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    for (r, c), cell in tbl.get_celld().items():
        if r == 0:
            cell.set_facecolor("#2C3E50")
            cell.set_text_props(color="white", fontweight="bold")
        elif r % 2 == 0:
            cell.set_facecolor("#F2F2F2")
        else:
            cell.set_facecolor("white")
        cell.set_edgecolor("#CCCCCC")

    if mw_sig_list:
        ax_info.text(0.02, 0.38,
                     "Motivi Mann-Whitney p < 0.05 (single vs split ratio):",
                     transform=ax_info.transAxes, fontsize=8.5, fontweight="bold",
                     va="top", color="#c0392b")
        mw_text = ",  ".join(mw_sig_list[:50])
        if len(mw_sig_list) > 50:
            mw_text += f"  ... (+{len(mw_sig_list)-50} altri)"
        ax_info.text(0.02, 0.30, mw_text,
                     transform=ax_info.transAxes, fontsize=7.5, va="top",
                     color="#333333",
                     bbox=dict(boxstyle="round,pad=0.3", facecolor="#FDECEA", alpha=0.7))
    else:
        ax_info.text(0.02, 0.35,
                     ("Colonna 'mw_pval_sc_vs_sp' non trovata nel batch_summary.\n"
                      "Aggiorna compute_summary_stats per salvare il p-value MW\n"
                      "e vedere qui i motivi significativi (p < 0.05)."),
                     transform=ax_info.transAxes, fontsize=8, va="top",
                     color="#777777", style="italic",
                     bbox=dict(boxstyle="round,pad=0.3", facecolor="#F5F5F5", alpha=0.7))

    ax_info.set_title("[A] Statistiche globali (valori uniformi tra motivi)",
                      fontweight="bold", fontsize=10, loc="left", pad=4)

    # ─────────────────────────────────────────────────────────────────────────
    # [B] Barplot candidati paralog-dropped per motivo
    # ─────────────────────────────────────────────────────────────────────────
    if "n_candidate_paralog_dropped" in oo.columns:
        cand_vals = [oo.loc[m, "n_candidate_paralog_dropped"]
                     if m in oo.index else np.nan for m in all_motifs]
    else:
        cand_vals = [np.nan] * n

    bar_colors_cand = [C_CAND if (not np.isnan(v) and v > 0) else "#DDDDDD"
                       for v in cand_vals]
    cand_plot = [v if not np.isnan(v) else 0 for v in cand_vals]

    ax_cand.bar(x, cand_plot, color=bar_colors_cand, alpha=ALPHA, width=0.8)
    ax_cand.axhline(0, color="black", linewidth=0.6)
    ax_cand.set_xlim(-0.5, n - 0.5)
    ax_cand.set_ylabel("N. OG candidati", fontsize=9)
    ax_cand.set_title(
        (f"[B] OG ortho-only con copertura >={int(good_species_frac*100)}% spp  "
         "(assenti in DISCO — probabili paralog non risolti)"),
        fontweight="bold", fontsize=10, color=C_CAND
    )
    ax_cand.set_xticks([])
    n_motifs_with_cand = sum(1 for v in cand_vals if not np.isnan(v) and v > 0)
    ax_cand.text(0.02, 0.97,
                 (f"Totale OG candidati: {int(sum(cand_plot))}  |  "
                  f"Motivi coinvolti: {n_motifs_with_cand}/{n}"),
                 transform=ax_cand.transAxes, fontsize=8.5, color="darkred", va="top")
    ax_cand.legend(handles=[
        Patch(facecolor=C_CAND,    alpha=0.85, label=f">={int(good_species_frac*100)}% spp"),
        Patch(facecolor="#DDDDDD", alpha=0.85, label="Nessun candidato"),
    ], fontsize=8, loc="upper right")

    # ─────────────────────────────────────────────────────────────────────────
    # [C] Delta frac-nonzero split — pannello grande in basso
    # ─────────────────────────────────────────────────────────────────────────
    delta_sp = [
        (sp.loc[m, "median_ortho_frac_nz"] - sp.loc[m, "median_disco_frac_nz"])
        if m in sp.index else np.nan
        for m in all_motifs
    ]
    bar_colors_sp = [
        C_SP if (not np.isnan(v) and v >= 0) else "#cc3333"
        for v in delta_sp
    ]
    ax_delta.bar(x, [v if not np.isnan(v) else 0 for v in delta_sp],
                 color=bar_colors_sp, alpha=ALPHA, width=0.85)
    ax_delta.axhline(0, color="black", linewidth=0.9)
    ax_delta.set_xlim(-0.5, n - 0.5)
    ax_delta.set_ylabel("Delta frac specie > 0  (ORTHO - DISCO)", fontsize=11)
    ax_delta.set_xlabel("Motivi (ordinati alfabeticamente)", fontsize=10)
    ax_delta.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    ax_delta.set_title(
        ("[C] Split >=2 forme — perdita copertura specie per motivo  "
         "(barre rosse: DISCO copre meno specie di ORTHO;  arancio: nessuna perdita)"),
        fontweight="bold", fontsize=12, color=C_SP
    )
    n_loss  = sum(v > 0  for v in delta_sp if not np.isnan(v))
    n_gain  = sum(v <= 0 for v in delta_sp if not np.isnan(v))
    n_valid = n_loss + n_gain
    ax_delta.text(0.01, 0.97,
                  (f"Motivi con perdita (ORTHO > DISCO): {n_loss}/{n_valid}   |   "
                   f"Senza perdita (DISCO >= ORTHO): {n_gain}/{n_valid}"),
                  transform=ax_delta.transAxes, fontsize=10, color="darkred", va="top")
    ax_delta.legend(handles=[
        Patch(facecolor=C_SP,      alpha=0.85, label="DISCO >= ORTHO (ok)"),
        Patch(facecolor="#cc3333", alpha=0.85, label="DISCO < ORTHO (perdita)"),
    ], fontsize=10, loc="upper right")

    fig.savefig(output_png, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Grafico riassuntivo salvato in: {output_png}")



#─────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(
        description="Confronto distribuzione punteggi dei motivi DISCO vs Orthofinder.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    # Defining the output directory
    p.add_argument("--outdir", help="Cartella di output")

    # Arguments to run this script for every motif in a single run
    p.add_argument("--run-all", action="store_true",
                   help="Processa tutti i motivi trovati in --disco-dir e --ortho-dir")
    p.add_argument("--disco-dir", default="02_totalscore",
                   help="Cartella DISCO  (default: 02_totalscore)")
    p.add_argument("--ortho-dir", default="05_totalscore_orthofinder",
                   help="Cartella ORTHO  (default: 05_totalscore_orthofinder)")
    p.add_argument("--save-per-motif-plots", action="store_true",
                   help="In --run-all, genera anche i grafici per singolo motivo "
                        "(più lento; di default solo le tabelle + grafico globale)")

    # Aggregation
    p.add_argument("--batch-summary", metavar="ROOT",
                   help="Aggrega i summary_stats.tsv in ROOT → batch_summary.tsv in --outdir")

    # Graphic overview
    p.add_argument("--plot-batch", metavar="BATCH_TSV",
                   help="Genera batch_overview.png da batch_summary.tsv")
    p.add_argument("--output-png", default="batch_overview_quick.png",
                   help="Nome PNG output per --plot-batch (default: batch_overview_quick.png). "
                        "Questo è il grafico rapido matplotlib; il grafico definitivo è prodotto "
                        "dallo script R (prova_rappresentazione_grafica.R).")

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