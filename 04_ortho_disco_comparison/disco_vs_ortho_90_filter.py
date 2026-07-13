#!/usr/bin/env python3
# Arguments to set in the bash command line for the launch of the code:
# --run-all is required to process all motifs in a single run.
#
# python motif_disco_vs_ortho.py \
#     --run-all \
#     --disco-dir 02_totalscore \
#     --ortho-dir 05_totalscore_orthofinder \
#     --outdir risultati_ortho_vs_disco
#
# Script used to verify the quantity of information lost by the DISCO pipeline compared to the Orthofinder output, 
# for each orthogroup and for each motif. It produces summary data tables with statistics for each motif, 
# and a batch summary table aggregating the results of all motifs.

import argparse
import os
import re
import sys
import warnings
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)


# ─────────────────────────────────────────────────────────────
# 1. Loading and preparation
# ─────────────────────────────────────────────────────────────

def load_tables(disco_path: str, ortho_path: str, min_species_ortho: int = 90):
    """
    Carica le tabelle DISCO e ORTHO e applica uno screening preliminare:
    gli OG di ORTHO con un numero di specie <= min_species_ortho vengono
    scartati prima della classificazione perché troppo poco informativi.
    Il filtro conta le specie con valore non-NaN nella riga ORTHO
    (ovvero specie in cui l'OG possiede almeno un gene), non il punteggio.
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


_OG_SUFFIX_RE = re.compile(r'^(.+)_\d+$')

def parse_og_base(index_name: str) -> str:
    m = _OG_SUFFIX_RE.match(index_name)
    return m.group(1) if m else index_name


def classify_orthogroups(disco: pd.DataFrame, ortho: pd.DataFrame):
    og_base_series = disco.index.to_series().apply(parse_og_base)
    disco_copy = disco.copy()
    disco_copy["_og_base"] = og_base_series.values

    species_cols = [c for c in disco.columns]

    n_forms = disco_copy.groupby("_og_base").size().rename("n_forms_disco")
    disco_agg = disco_copy.groupby("_og_base")[species_cols].sum(min_count=1)
    disco_forms = disco_copy.drop(columns=["_og_base"])

    disco_og_bases = set(n_forms.index)
    ortho_og_bases = set(ortho.index)

    common = disco_og_bases & ortho_og_bases
    ortho_only = ortho_og_bases - disco_og_bases

    categories = {}
    for og in common:
        nf = n_forms[og]
        if nf == 1:
            categories[og] = "single_complete"
        else:
            categories[og] = "split"
    for og in ortho_only:
        categories[og] = "ortho_only"

    og_meta = pd.DataFrame({
        "og_base": list(categories.keys()),
        "n_forms_disco": [n_forms.get(og, 0) for og in categories.keys()],
        "category": list(categories.values()),
    }).set_index("og_base")

    return og_meta, disco_agg, disco_forms


# ─────────────────────────────────────────────────────────────
# 2. Comparison of DISCO and ORTHO metrics per OG
# ─────────────────────────────────────────────────────────────

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


def build_per_og_comparison(disco_agg: pd.DataFrame,
                             ortho: pd.DataFrame,
                             og_meta: pd.DataFrame,
                             disco_forms: pd.DataFrame,
                             good_og_species_frac: float = 0.6) -> pd.DataFrame:
    """
    Confronto RISTRETTO: vengono analizzati SOLO gli OG di ORTHO che hanno
    una corrispondenza in DISCO (single_complete o split).
    """
    common_ogs = disco_agg.index.intersection(ortho.index)

    n_ortho_total   = ortho.shape[0]
    n_common        = len(common_ogs)
    n_ortho_only    = n_ortho_total - n_common

    n_disco_forms   = disco_forms.shape[0]
    n_disco_bases   = disco_agg.shape[0]
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


def compute_summary_stats(per_og: pd.DataFrame, og_meta: pd.DataFrame,
                           ortho: pd.DataFrame,
                           n_og_raw: int = 0,
                           n_og_filtered: int = 0,
                           good_species_frac: float = 0.60,
                           mw_pval: float = np.nan,
                           n_disco_forms: int = 0,
                           n_disco_bases: int = 0,
                           n_extra_forms: int = 0) -> pd.DataFrame:
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
            "mw_pval_sc_vs_sp":            mw_pval,
            "median_ratio_sum":            sub["ratio_sum"].median(),
            "mean_ratio_sum":              sub["ratio_sum"].mean(),
            "frac_ratio_gt_0.5":           (sub["ratio_sum"] > 0.5).mean(),
            "median_disco_sum":            sub["disco_sum"].median(),
            "median_ortho_sum":            sub["ortho_sum"].median(),
            "median_delta_mean":           sub["delta_mean"].median(),
            "median_disco_frac_nz":        sub["disco_frac_nonzero"].median(),
            "median_ortho_frac_nz":        sub["ortho_frac_nonzero"].median(),
        })

    ortho_only_ogs = og_meta[og_meta["category"] == "ortho_only"].index
    ortho_only_ogs = [og for og in ortho_only_ogs if og in ortho.index]
    if ortho_only_ogs:
        sums = ortho.loc[ortho_only_ogs].sum(axis=1, skipna=True)
        frac_nz = (ortho.loc[ortho_only_ogs] > 0).mean(axis=1)
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


CAT_LABELS = {
    "single_complete": "Single-complete",
    "split":           "Split (≥2 forme)",
    "ortho_only":      "Solo ORTHO",
}


# ─────────────────────────────────────────────────────────────
# 3. Processing single OG
# ─────────────────────────────────────────────────────────────

def process_motif(disco_path: str, ortho_path: str, outdir: str):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    motif_name = Path(disco_path).stem.replace("totalscore_", "")

    print(f"[{motif_name}] Caricamento dati...")
    disco, ortho, n_og_raw, n_og_filtered = load_tables(disco_path, ortho_path)
    print(f"  DISCO: {disco.shape[0]} righe × {disco.shape[1]} specie")
    print(f"  ORTHO: {ortho.shape[0]} OG × {ortho.shape[1]} specie "
          f"(raw={n_og_raw}, filtrati={n_og_filtered})")

    print(f"[{motif_name}] Classificazione ortogruppi...")
    og_meta, disco_agg, disco_forms = classify_orthogroups(disco, ortho)
    cat_counts = og_meta["category"].value_counts()
    for cat, n in cat_counts.items():
        print(f"  {CAT_LABELS.get(cat, cat)}: {n} OG")

    print(f"[{motif_name}] Calcolo metriche comparative...")
    per_og = build_per_og_comparison(disco_agg, ortho, og_meta, disco_forms)

    ratio_sc = per_og[per_og["category"] == "single_complete"]["ratio_sum"].dropna()
    ratio_sp = per_og[per_og["category"] == "split"]["ratio_sum"].dropna()
    if len(ratio_sc) > 5 and len(ratio_sp) > 5:
        stat, pval = stats.mannwhitneyu(ratio_sc, ratio_sp, alternative="two-sided")
        print(f"  Mann-Whitney ratio (single vs split): U={stat:.0f}, p={pval:.4f}")
    else:
        pval = np.nan
        print("  Campione troppo piccolo per test statistico")

    n_disco_forms_cnt = disco_forms.shape[0]
    n_disco_bases_cnt = disco_agg.shape[0]
    n_extra_forms_cnt = n_disco_forms_cnt - n_disco_bases_cnt

    summary = compute_summary_stats(per_og, og_meta, ortho,
                                    n_og_raw=n_og_raw,
                                    n_og_filtered=n_og_filtered,
                                    mw_pval=pval,
                                    n_disco_forms=n_disco_forms_cnt,
                                    n_disco_bases=n_disco_bases_cnt,
                                    n_extra_forms=n_extra_forms_cnt)

    summary.to_csv(out / "summary_stats.tsv", sep="\t", float_format="%.4f")
    per_og.to_csv(out / "per_og_comparison.tsv", sep="\t", float_format="%.4f")
    og_meta.to_csv(out / "og_classification.tsv", sep="\t")
    print(f"[{motif_name}] Tabelle salvate in: {out}")

    return {"summary_stats": summary, "per_og_comparison": per_og, "og_meta": og_meta}


# ─────────────────────────────────────────────────────────────
# 4. Comparison between every single motif in one run
# ─────────────────────────────────────────────────────────────

def run_all(disco_dir: str, ortho_dir: str, outdir: str):
    disco_dir = Path(disco_dir)
    ortho_dir = Path(ortho_dir)
    out_root  = Path(outdir)
    out_root.mkdir(parents=True, exist_ok=True)

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
            process_motif(str(disco_path), str(ortho_path), str(motif_outdir))
            processed.append(motif_name)
        except Exception as e:
            print(f"  ERRORE durante il processing di {motif_name}: {e}")
            skipped.append(motif_name)

    print(f"\n{'='*60}")
    print(f"Processati: {len(processed)}  |  Saltati/Errori: {len(skipped)}")
    if skipped:
        print(f"Saltati: {skipped}")

    print("\nGenerazione batch_summary.tsv...")
    batch_tsv = out_root / "batch_summary.tsv"
    batch_summary(str(out_root), str(batch_tsv))

    print(f"\nDone. Output globale in: {out_root}")


def batch_summary(results_root: str, output_file: str = "batch_summary.tsv"):
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


# ─────────────────────────────────────────────────────────────
# 5. Setting command line for bash launch
# ─────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Confronto distribuzione punteggi dei motivi DISCO vs Orthofinder (Versione solo tabelle).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    p.add_argument("--outdir", help="Cartella di output")
    p.add_argument("--run-all", action="store_true",
                   help="Processa tutti i motivi trovati in --disco-dir e --ortho-dir")
    p.add_argument("--disco-dir", default="02_totalscore",
                   help="Cartella DISCO  (default: 02_totalscore)")
    p.add_argument("--ortho-dir", default="05_totalscore_orthofinder",
                   help="Cartella ORTHO  (default: 05_totalscore_orthofinder)")
    p.add_argument("--batch-summary", metavar="ROOT",
                   help="Aggrega i summary_stats.tsv in ROOT → batch_summary.tsv in --outdir")

    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()

    if args.run_all:
        if not args.outdir:
            print("ERRORE: --outdir e obbligatorio con --run-all")
            sys.exit(1)
        run_all(args.disco_dir, args.ortho_dir, args.outdir)

    elif args.batch_summary:
        out_tsv = Path(args.outdir) / "batch_summary.tsv" if args.outdir else "batch_summary.tsv"
        batch_summary(args.batch_summary, str(out_tsv))

    else:
        print(__doc__)