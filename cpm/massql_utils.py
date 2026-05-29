# massql_utils.py
from typing import List, Tuple, Dict, Iterable, Optional
import os
import math
import re
import tempfile

import numpy as np
import pandas as pd
from massql import msql_engine, msql_fileloading


# -------------------------------------------------------------------------------------------
# Workaround for mzML files with scan IDs like: id="merged=4505 row=0"
# MassQL expects scan IDs like: id="scan=4505" so if it has row=0 it will break. this changes that so it's all formatted correctly
# -------------------------------------------------------------------------------------------
def _fix_mzml_scan_ids(input_file: str) -> str:
    if not input_file.lower().endswith(".mzml"):
        return input_file

    with open(input_file, "r", errors="ignore") as f:
        text = f.read()

    fixed_text = re.sub(
        r'id="merged=(\d+)\s+row=\d+"',
        r'id="scan=\1"',
        text
    )

    tmp = tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".mzML",
        delete=False
    )
    tmp.write(fixed_text)
    tmp.close()

    return tmp.name


# -------------------------------------------------------------------------------------------
# Keep only ms2 spectra that belong to features retained in ms1 feature table
# -------------------------------------------------------------------------------------------
def _filter_ms2_df_to_kept_features(
    ms2_df: pd.DataFrame,
    kept_feats: Optional[pd.DataFrame],
    *,
    mz_tol_da: float = 0.5,
    rt_tol_min: float = 0.3,
    ms2_precursor_col: str = "precmz",
    ms2_rt_col: str = "rt",
) -> pd.DataFrame:

    if ms2_df is None or ms2_df.empty:
        return ms2_df
    if kept_feats is None or len(kept_feats) == 0:
        return ms2_df

    if ms2_precursor_col not in ms2_df.columns or ms2_rt_col not in ms2_df.columns:
        return ms2_df

    if "merged_precmz" not in kept_feats.columns or "rt_median" not in kept_feats.columns:
        return ms2_df

    ms2_prec = ms2_df[ms2_precursor_col].astype(float).to_numpy()
    ms2_rt = ms2_df[ms2_rt_col].astype(float).to_numpy()

    kept_mz = kept_feats["merged_precmz"].astype(float).to_numpy()
    kept_rt = kept_feats["rt_median"].astype(float).to_numpy()

    keep_mask = np.zeros(len(ms2_df), dtype=bool)

    for mz0, rt0 in zip(kept_mz, kept_rt):
        keep_mask |= (
            (ms2_prec >= (mz0 - mz_tol_da)) &
            (ms2_prec <= (mz0 + mz_tol_da)) &
            (ms2_rt >= (rt0 - rt_tol_min)) &
            (ms2_rt <= (rt0 + rt_tol_min))
        )

    return ms2_df.loc[keep_mask].copy()



# -------------------------------------------------------------------------------------------
# Load a single file in
# -------------------------------------------------------------------------------------------
def load_file(
    input_file: str,
    *,
    kept_feature_table: Optional[pd.DataFrame] = None,
    mz_tol_da: float = 0.2,
    rt_tol_min: float = 0.3,
    blank_files: Optional[set] = None,
):

    base = os.path.basename(input_file)

    fixed_input_file = _fix_mzml_scan_ids(input_file)

    if blank_files and base in blank_files:
        ms1_df, ms2_df = msql_fileloading.load_data(fixed_input_file)
        return ms1_df.iloc[0:0].copy(), ms2_df.iloc[0:0].copy()

    ms1_df, ms2_df = msql_fileloading.load_data(fixed_input_file)

    ms1_df = ms1_df.copy()
    ms2_df = ms2_df.copy()

    if "source_file" not in ms1_df.columns:
        ms1_df["source_file"] = base
    if "source_file" not in ms2_df.columns:
        ms2_df["source_file"] = base

    if kept_feature_table is not None and len(kept_feature_table) > 0:
        ms2_df = _filter_ms2_df_to_kept_features(
            ms2_df,
            kept_feature_table,
            mz_tol_da=mz_tol_da,
            rt_tol_min=rt_tol_min,
            ms2_precursor_col="precmz",
            ms2_rt_col="rt",
        )

    return ms1_df, ms2_df


# -------------------------------------------------------------------------------------------
# Load multiple files
# -------------------------------------------------------------------------------------------
def load_files(
    input_files: List[str],
    *,
    kept_feature_table: Optional[pd.DataFrame] = None,
    mz_tol_da: float = 0.2,
    rt_tol_min: float = 0.3,
    blank_files: Optional[set] = None,
) -> Dict[str, object]:

    ms1_list, ms2_list = [], []
    per_file: Dict[str, Dict[str, pd.DataFrame]] = {}

    for f in input_files:
        ms1_df, ms2_df = load_file(
            f,
            kept_feature_table=kept_feature_table,
            mz_tol_da=mz_tol_da,
            rt_tol_min=rt_tol_min,
            blank_files=blank_files,
        )

        base = os.path.basename(f)
        per_file[base] = {"ms1": ms1_df, "ms2": ms2_df}
        ms1_list.append(ms1_df)
        ms2_list.append(ms2_df)

    ms1_all = pd.concat(ms1_list, ignore_index=True) if ms1_list else pd.DataFrame()
    ms2_all = pd.concat(ms2_list, ignore_index=True) if ms2_list else pd.DataFrame()

    return {"per_file": per_file, "merged": (ms1_all, ms2_all)}


# -------------------------------------------------------------------------------------------
# Build a MassQL query string
# -------------------------------------------------------------------------------------------
def build_massql_query(
    ions_mz: Iterable[float],
    tol_mz: float = 0.01,
    polarity: Optional[str] = None,
    rt_window: Optional[Tuple[float, float]] = None,
) -> str:

    clauses = []

    if polarity:
        clauses.append(f"POLARITY={polarity}")

    if rt_window:
        rtmin, rtmax = rt_window
        clauses.append(f"RTMIN={rtmin} AND RTMAX={rtmax}")

    for m in ions_mz:
        clauses.append(f"MS2PROD={m}:TOLERANCEMZ={tol_mz}")

    where = " AND ".join(clauses)

    return f"QUERY scaninfo(MS2DATA) WHERE {where}"


# -------------------------------------------------------------------------------------------
# Run across files
# -------------------------------------------------------------------------------------------
def run_across_files_individual(
    input_files: Iterable[str],
    ions: Iterable[float],
    tol_mz: float = 0.01,
    polarity: Optional[str] = None,
    rt_window: Optional[Tuple[float, float]] = None,
    *,
    data: Optional[Dict[str, object]] = None,
) -> pd.DataFrame:

    rows = []

    for f in input_files:
        base = os.path.basename(f)

        if data is not None:
            ms1_df = data["per_file"][base]["ms1"]
            ms2_df = data["per_file"][base]["ms2"]
        else:
            ms1_df, ms2_df = load_file(f)

        for ion in ions:
            q = build_massql_query(
                [ion],
                tol_mz=tol_mz,
                polarity=polarity,
                rt_window=rt_window,
            )

            df = msql_engine.process_query(
                q,
                f,
                ms1_df=ms1_df,
                ms2_df=ms2_df,
            )

            if df is not None and not df.empty:
                df = df.copy()
                df["ion"] = float(ion)
                df["source_file"] = base
                rows.append(df)

    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
# -------------------------------------------------------------------------------------------
# Add friendly labels for MP ions
def add_MP_labels(df: pd.DataFrame, ion_to_MP: Dict[float, str], label_col: str = "MP") -> pd.DataFrame:
    out = df.copy()
    if "ion" in out.columns:
        out[label_col] = out["ion"].map(lambda x: ion_to_MP.get(float(x), f"{x:.4f}"))
    elif "combo_label" in out.columns:
        def label_combo(s: str) -> str:
            parts = s.split("+")
            return "+".join([ion_to_MP.get(float(p), p) for p in parts])
        out[label_col] = out["combo_label"].map(label_combo)
    return out

#add  friendly labels for MC ions
def add_MC_labels(df: pd.DataFrame, ion_to_MC: Dict[float, str], label_col: str = "MC") -> pd.DataFrame:
    out = df.copy()
    if "ion" in out.columns:
        out[label_col] = out["ion"].map(lambda x: ion_to_MC.get(float(x), f"{x:.4f}"))
    elif "combo_label" in out.columns:
        def label_combo(s: str) -> str:
            parts = s.split("+")
            return "+".join([ion_to_MC.get(float(p), p) for p in parts])
        out[label_col] = out["combo_label"].map(label_combo)
    return out

#add friendly labels for AR
def add_AR_labels(df: pd.DataFrame, ion_to_AR: Dict[float, str], label_col: str = "AR") -> pd.DataFrame:
    out = df.copy()
    if "ion" in out.columns:
        out[label_col] = out["ion"].map(lambda x: ion_to_AR.get(float(x), f"{x:.4f}"))
    elif "combo_label" in out.columns:
        def label_combo(s: str) -> str:
            parts = s.split("+")
            return "+".join([ion_to_AR.get(float(p), p) for p in parts])
        out[label_col] = out["combo_label"].map(label_combo)
    return out

#add friendly labels for AB
def add_AB_labels(df: pd.DataFrame, ion_to_AB: Dict[float, str], label_col: str = "AB") -> pd.DataFrame:
    out = df.copy()
    if "ion" in out.columns:
        out[label_col] = out["ion"].map(lambda x: ion_to_AB.get(float(x), f"{x:.4f}"))
    elif "combo_label" in out.columns:
        def label_combo(s: str) -> str:
            parts = s.split("+")
            return "+".join([ion_to_AB.get(float(p), p) for p in parts])
        out[label_col] = out["combo_label"].map(label_combo)
    return out

#add friendly labels for MV
def add_MV_labels(df: pd.DataFrame, ion_to_MV: Dict[float, str], label_col: str = "MV") -> pd.DataFrame:
    out = df.copy()
    if "ion" in out.columns:
        out[label_col] = out["ion"].map(lambda x: ion_to_MV.get(float(x), f"{x:.4f}"))
    elif "combo_label" in out.columns:
        def label_combo(s: str) -> str:
            parts = s.split("+")
            return "+".join([ion_to_MV.get(float(p), p) for p in parts])
        out[label_col] = out["combo_label"].map(label_combo)
    return out

#add friendly labels for AG
def add_AG_labels(df: pd.DataFrame, ion_to_AG: Dict[float, str], label_col: str = "AG") -> pd.DataFrame:
    out = df.copy()
    if "ion" in out.columns:
        out[label_col] = out["ion"].map(lambda x: ion_to_AG.get(float(x), f"{x:.4f}"))
    elif "combo_label" in out.columns:
        def label_combo(s: str) -> str:
            parts = s.split("+")
            return "+".join([ion_to_AG.get(float(p), p) for p in parts])
        out[label_col] = out["combo_label"].map(label_combo)
    return out

#add friendly labels for MG
def add_MG_labels(df: pd.DataFrame, ion_to_MG: Dict[float, str], label_col: str = "MG") -> pd.DataFrame:
    out = df.copy()
    if "ion" in out.columns:
        out[label_col] = out["ion"].map(lambda x: ion_to_MG.get(float(x), f"{x:.4f}"))
    elif "combo_label" in out.columns:
        def label_combo(s: str) -> str:
            parts = s.split("+")
            return "+".join([ion_to_MG.get(float(p), p) for p in parts])
        out[label_col] = out["combo_label"].map(label_combo)
    return out
