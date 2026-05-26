
#!/usr/bin/env python3
"""
cyanometdb_match.py

- Match MS1 summaries to CyanometDB entries by precursor m/z within a tolerance
- Build an "unknowns" sheet for rows with no library hit and >=2 diagnostic has_* flags
- Export CSV / Excel / optional unknowns heatmap PNG
"""

from __future__ import annotations

import os
from pathlib import Path
from datetime import datetime
from typing import Iterable, Optional, Sequence

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

# Optional for unknowns heatmap
try:
    import seaborn as sns  # type: ignore
    _HAS_SNS = True
except Exception:
    _HAS_SNS = False


# --------------------
# Defaults / constants
# --------------------
MS1_BASE_KEEP: Sequence[str] = (
    "cluster_id",
    "merged_precmz",
    "n_scans",
    "scan_nunique",
    "scan_ids",
    "ms1_scan_ids",
    "files",
    "source_file_<lambda>",
    "rt_min",
    "rt_median",
    "rt_max",
)

LIB_MZ_COL = "Monoisotopic mass [M+H]+"
MS1_MZ_COL = "merged_precmz"
CLASS_COLS = ("Class of compound", "Alternative class names")


# --------------------
# IO helpers
# --------------------
def read_any_table(path: str | os.PathLike) -> pd.DataFrame:
    p = Path(path)
    suf = p.suffix.lower()
    if suf == ".csv":
        return pd.read_csv(p)
    if suf in {".tsv", ".txt"}:
        return pd.read_csv(p, sep="\t")
    if suf == ".parquet":
        return pd.read_parquet(p)
    raise ValueError(f"Unsupported table format: {p.suffix}")


# --------------------
# Library loading
# --------------------
def load_library(
    xlsx_path: str | os.PathLike,
    *,
    class_filter: Optional[Iterable[str] | str] = None,
    class_cols: Sequence[str] = CLASS_COLS,
    mz_col: str = LIB_MZ_COL,
    sheet_index: int = 1,
) -> pd.DataFrame:
    """
    Load the CyanometDB Excel sheet.
    Optionally filter by class using substring matching across class columns.
    """
    try:
        lib = pd.read_excel(xlsx_path, sheet_name=sheet_index, engine="openpyxl")
    except Exception:
        lib = pd.read_excel(xlsx_path, sheet_name=sheet_index)

    lib.columns = [str(c).strip() for c in lib.columns]

    if mz_col not in lib.columns:
        raise KeyError(f"mz_col '{mz_col}' not found in library columns.")

    lib[mz_col] = pd.to_numeric(lib[mz_col], errors="coerce")
    lib = lib[lib[mz_col].notna()].copy()

    if class_filter:
        if isinstance(class_filter, str):
            targets = {class_filter.strip().lower()}
        else:
            targets = {str(s).strip().lower() for s in class_filter if s is not None}

        def _row_ok(row):
            hay = []
            for col in class_cols:
                if col in lib.columns and pd.notna(row.get(col)):
                    hay.append(str(row[col]).lower())
            return any(any(t in h for h in hay) for t in targets)

        lib = lib[lib.apply(_row_ok, axis=1)].copy()

    return lib.sort_values(mz_col).reset_index(drop=True)


# --------------------
# MS1 selection helper
# --------------------
def select_ms1_columns(
    ms1_df: pd.DataFrame,
    base_keep: Sequence[str] = MS1_BASE_KEEP,
    drop_empty_has: bool = True,
) -> pd.DataFrame:
    """
    Keep core metadata columns plus has_* columns present in this table.
    Optionally drop has_* columns that are entirely empty / false.
    """
    cols = [c for c in base_keep if c in ms1_df.columns]
    has_cols_all = [c for c in ms1_df.columns if str(c).strip().startswith("has_")]

    if drop_empty_has:
        has_cols = []
        for c in has_cols_all:
            col = ms1_df[c]
            if col.dtype == bool:
                keep = col.any()
            else:
                keep = (col.fillna(0) != 0).any()
            if keep:
                has_cols.append(c)
    else:
        has_cols = has_cols_all

    keep = list(dict.fromkeys(cols + has_cols))
    return ms1_df.loc[:, keep].copy() if keep else ms1_df.copy()


# --------------------
# Matching
# --------------------
def match_ms1_to_lib(
    ms1_df: pd.DataFrame,
    lib_df: pd.DataFrame,
    *,
    ms1_mz_col: str = MS1_MZ_COL,
    lib_mz_col: str = LIB_MZ_COL,
    tol_da: float = 0.05,
) -> pd.DataFrame:
    """
    Match MS1 features to library entries based on precursor m/z within ± tol_da.
    Returns a merged DataFrame with library annotations appended to each hit.
    """
    ms1_df = ms1_df.copy()
    ms1_df[ms1_mz_col] = pd.to_numeric(ms1_df[ms1_mz_col], errors="coerce")
    valid_ms1 = ms1_df[ms1_df[ms1_mz_col].notna()].copy()

    results = []

    for _, ms1_row in valid_ms1.iterrows():
        mz = ms1_row[ms1_mz_col]
        hits = lib_df[np.abs(lib_df[lib_mz_col] - mz) <= tol_da]

        if hits.empty:
            results.append({**ms1_row.to_dict(), "Compound identifier": np.nan})
        else:
            for _, hit in hits.iterrows():
                results.append({**ms1_row.to_dict(), **hit.to_dict()})

    return pd.DataFrame(results)


# --------------------
# Unknowns helpers
# --------------------
def _excel_engine() -> str:
    try:
        import xlsxwriter  # noqa: F401
        return "xlsxwriter"
    except Exception:
        print("xlsxwriter not installed — falling back to openpyxl.")
        return "openpyxl"


def build_unknowns_sheet(matches: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    """
    Build a table of unknown rows:
    - no Compound identifier
    - >=2 diagnostic has_* flags
    """
    universe = matches.copy()
    bool_cols = [c for c in universe.columns if str(c).strip().startswith("has_")]

    if "Compound identifier" not in universe.columns:
        return pd.DataFrame(), bool_cols

    unknown_mask = universe["Compound identifier"].isna()
    no_database_hits = universe.loc[unknown_mask].copy()

    if no_database_hits.empty or not bool_cols:
        return pd.DataFrame(), bool_cols

    no_database_hits[bool_cols] = no_database_hits[bool_cols].fillna(False).astype(bool)
    no_database_hits["n_diagnostic"] = no_database_hits[bool_cols].sum(axis=1)
    no_database_hits = no_database_hits.loc[no_database_hits["n_diagnostic"] >= 2].copy()

    if no_database_hits.empty:
        return pd.DataFrame(), bool_cols

    def make_pattern(row):
        return "+".join([col.replace("has_", "", 1) for col in bool_cols if row[col]])

    no_database_hits["fragment_pattern"] = no_database_hits.apply(make_pattern, axis=1)

    group_keys = []
    for maybe in ["files", "source_file_<lambda>"]:
        if maybe in no_database_hits.columns:
            group_keys.append(maybe)
            break
    group_keys += ["merged_precmz", "fragment_pattern"]

    agg_map = {c: "max" for c in bool_cols}

    present_scan_cols = [
        c for c in ["scan_ids", "ms1_scan_ids", "scan_number", "ms1_scan"]
        if c in no_database_hits.columns
    ]

    def uniq_join(x):
        vals = []
        for v in x.dropna():
            if isinstance(v, str):
                vals.extend([s.strip() for s in v.split(",") if s.strip()])
            else:
                vals.append(str(v))
        vals = sorted(set(vals))
        return ",".join(vals)

    for c in present_scan_cols:
        agg_map[c] = uniq_join

    grouped_unknowns = (
        no_database_hits
        .groupby(group_keys, dropna=False)
        .agg(agg_map)
        .reset_index()
        .sort_values(["merged_precmz"] + [g for g in group_keys if g != "merged_precmz"])
    )

    files_col = next((k for k in ["files", "source_file_<lambda>"] if k in grouped_unknowns.columns), None)

    def _files_to_str(v):
        if isinstance(v, (list, tuple)):
            return ",".join(map(str, v))
        return str(v)

    grouped_unknowns["row_label"] = grouped_unknowns.apply(
        lambda row: (
            f"Unknown | m/z={row['merged_precmz']:.4f}"
            f" | pattern={row['fragment_pattern']}"
            + (f" | file={_files_to_str(row[files_col])}" if files_col else "")
        ),
        axis=1,
    )

    return grouped_unknowns, bool_cols


# --------------------
# Matched tiles plot
# --------------------
def plot_matched_tiles(
    matched_only: pd.DataFrame,
    out_dir_ts: str,
    ts: str,
    *,
    bool_prefix: str = "has_",
    mz_col: str = "merged_precmz",
    file_col_candidates: Sequence[str] = ("files", "source_file_<lambda>"),
    compound_col: str = "Compound name",
    base_cmap: str = "tab20",
) -> Optional[str]:
    """
    Plot matched compounds as a has_* grid colored by contributing file(s).
    """
    if matched_only is None or matched_only.empty:
        print("No matched rows to visualize — skipping matched-compound grid.")
        return None

    if compound_col not in matched_only.columns:
        print(f"No '{compound_col}' column found — skipping matched-compound grid.")
        return None

    bool_cols = [c for c in matched_only.columns if str(c).strip().startswith(bool_prefix)]
    if not bool_cols:
        print(f"No columns starting with '{bool_prefix}' found — skipping matched-compound grid.")
        return None

    file_col = next((c for c in file_col_candidates if c in matched_only.columns), None)
    if file_col is None:
        print(f"No file provenance column found in {file_col_candidates} — skipping matched-compound grid.")
        return None

    df = matched_only.copy()

    def _files_to_list(v):
        if pd.isna(v):
            return []
        if isinstance(v, (list, tuple, set)):
            return list(v)
        if isinstance(v, str):
            return [x.strip() for x in v.split(",") if x.strip()]
        return [str(v)]

    df[file_col] = df[file_col].apply(_files_to_list)

    grouped = (
        df.groupby([compound_col, mz_col] + bool_cols, dropna=False)[file_col]
        .apply(lambda x: sorted(set(f for sub in x for f in sub)))
        .reset_index()
    )

    if grouped.empty:
        print("Nothing to plot after grouping matched compounds.")
        return None

    grouped = grouped.sort_values(mz_col, ascending=True).reset_index(drop=True)
    grouped["row_label"] = (
        grouped[compound_col].fillna("No match").astype(str)
        + " | m/z="
        + grouped[mz_col].map(lambda x: f"{float(x):.4f}" if pd.notna(x) else "nan")
    )

    all_files = sorted({f for L in grouped[file_col] for f in L})
    if not all_files:
        print("No file provenance found after grouping — nothing to color by.")
        return None

    cmap = cm.get_cmap(base_cmap, max(len(all_files), 2))
    file_to_color = {f: cmap(i) for i, f in enumerate(all_files)}

    fig, ax = plt.subplots(figsize=(12, 0.5 * len(grouped) + 2))

    for row_idx, row in grouped.iterrows():
        for col_idx, col in enumerate(bool_cols):
            if bool(row[col]):
                files_here = row[file_col]

                if not files_here:
                    ax.add_patch(
                        plt.Rectangle((col_idx, row_idx), 1, 1,
                                      facecolor=(0.85, 0.85, 0.85, 1.0),
                                      edgecolor="white")
                    )
                elif len(files_here) == 1:
                    ax.add_patch(
                        plt.Rectangle((col_idx, row_idx), 1, 1,
                                      facecolor=file_to_color[files_here[0]],
                                      edgecolor="white")
                    )
                else:
                    stripe_h = 1.0 / len(files_here)
                    for k, f in enumerate(files_here):
                        ax.add_patch(
                            plt.Rectangle(
                                (col_idx, row_idx + k * stripe_h),
                                1,
                                stripe_h,
                                facecolor=file_to_color.get(f, (0.85, 0.85, 0.85, 1.0)),
                                edgecolor="white",
                                linewidth=0.3,
                            )
                        )

    ax.set_xlim(0, len(bool_cols))
    ax.set_ylim(0, len(grouped))
    ax.set_xticks(np.arange(len(bool_cols)) + 0.5)
    ax.set_xticklabels(bool_cols, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(grouped)) + 0.5)
    ax.set_yticklabels(grouped["row_label"], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Diagnostic fragments / has_* flags")
    ax.set_ylabel("Matched compounds (sorted by m/z)")
    ax.set_title("Matched compounds: file-colored tiles per diagnostic fragment", fontsize=14)

    patches = [mpatches.Patch(color=file_to_color[f], label=f) for f in all_files]
    n = len(all_files)
    ncols = 1 if n <= 12 else 2 if n <= 30 else 3 if n <= 60 else 4

    ax.legend(
        handles=patches,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        fontsize=6,
        ncol=ncols,
        title="Files",
        borderaxespad=0.0,
    )

    fig.tight_layout()

    os.makedirs(out_dir_ts, exist_ok=True)
    out_png = os.path.join(out_dir_ts, f"matched_compound_tiles_{ts}.png")
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"Saved matched-compounds grid: {os.path.abspath(out_png)}")
    plt.close(fig)
    return out_png


# --------------------
# Output writer
# --------------------
def write_outputs(
    matches: pd.DataFrame,
    matched_only: pd.DataFrame,
    out_dir: str,
    ts: str | None = None,
    write_excel: bool = True,
    make_heatmap: bool = True,
):
    """
    Write outputs into a timestamped directory.
    Returns (out_dir_ts, paths_dict).
    """
    if ts is None:
        ts = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    out_dir_ts = f"{out_dir}_{ts}"
    os.makedirs(out_dir_ts, exist_ok=True)

    paths = {}

    csv_path = os.path.join(out_dir_ts, f"cyanometdb_matches_{ts}.csv")
    matched_only.to_csv(csv_path, index=False)
    print(f"Saved: {os.path.abspath(csv_path)} (rows: {len(matched_only)})")
    paths["matches_csv"] = csv_path

    unknowns, bool_cols = build_unknowns_sheet(matches)

    if write_excel:
        excel_path = os.path.join(out_dir_ts, f"cyanometdb_matches_{ts}.xlsx")
        engine = _excel_engine()
        with pd.ExcelWriter(excel_path, engine=engine) as xw:
            matched_only.to_excel(xw, index=False, sheet_name="matches")
            if not unknowns.empty:
                unknowns.to_excel(xw, index=False, sheet_name="unknowns_>=2diag")
            else:
                pd.DataFrame({"note": ["no unknowns with >=2 diagnostic ions"]}).to_excel(
                    xw, index=False, sheet_name="unknowns_>=2diag"
                )
        print(f"Excel written: {os.path.abspath(excel_path)}")
        paths["excel"] = excel_path

    if not unknowns.empty:
        unk_csv = os.path.join(out_dir_ts, f"unknown_features_with_scans_{ts}.csv")
        unknowns.to_csv(unk_csv, index=False)
        print(f"Exported unknown features with scans: {os.path.abspath(unk_csv)} (rows: {len(unknowns)})")
        paths["unknowns_csv"] = unk_csv

        if make_heatmap and _HAS_SNS and bool_cols:
            matrix = unknowns.set_index("row_label")[bool_cols]
            if not matrix.empty:
                fig, ax = plt.subplots(figsize=(12, 0.5 * len(matrix) + 2))
                sns.heatmap(
                    matrix.astype(int),
                    cmap=["lightgrey", "orange"],
                    cbar=False,
                    linewidths=0.5,
                    linecolor="white",
                    ax=ax,
                )
                ax.set_title("Unknown Features (≥2 diagnostic ions, per-file detail)", fontsize=14)
                ax.set_xlabel("Diagnostic Fragments")
                ax.set_ylabel("Unknown Precursor m/z (pattern, file)")
                plt.xticks(rotation=45, ha="right")
                plt.yticks(fontsize=7)
                fig.tight_layout()

                out_png = os.path.join(out_dir_ts, f"unknown_features_with_scans_{ts}.png")
                fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
                print(f"Saved figure: {os.path.abspath(out_png)}")
                plt.close(fig)
                paths["unknowns_png"] = out_png
        elif make_heatmap and not _HAS_SNS:
            print("seaborn not installed — skipping heatmap PNG for unknowns.")

    return out_dir_ts, paths
