# CPM_cli_04_14_2026_test.py

from __future__ import annotations

import argparse
from contextlib import ExitStack, nullcontext, redirect_stdout, redirect_stderr
from datetime import datetime
import importlib
from pathlib import Path
import shutil
import sys

import numpy as np
import pandas as pd
from pyteomics import mzml

from . import massql_utils as mu
from . import adduct_finder as af
from . import adduct_pipeline as ap

from .summary_builder import (
    make_summary_ind,
    explode_matches_to_per_file,
    add_ms1_auc_from_points,
    pool_auc_back_to_matches,
    make_scan_level_has_table,
)

from .qc_helpers import (
    QCConfig,
    load_metadata,
    load_ms1_points,
    compute_ref_auc_df,
    add_reference_normalization,
    apply_batch_correction_median_scaling,
    blank_filter_perfile_table,
    blank_filter_perfile_table_by_batch,
    save_qc_audit,
    flatten_subfolders,
)

from .rt_histograms import plot_rt_histograms
from .cyanopeptide_counts_plots import plot_indiv_counts
from .rt_mz_plot import plot_precursor_rt
from .plotting_ind_heatmap import plot_heatmaps
from .indiv_combo_dot_plot import plot_indiv_scatter
from .sum_intensity_from_scans import sum_intensities
from .ms2_tilemap_intensities import plot_has_tilemap
from .cyanometdb_match import (
    load_library,
    read_any_table,
    match_ms1_to_lib,
    select_ms1_columns,
    write_outputs,
    plot_matched_tiles,
)

BUNDLED_LIB_PATH = "data/CyanoMetDB_Version03.xlsx"

# ------------------------------------------------------------------
# Ion lists
# ------------------------------------------------------------------
MCIONS = [135.0804, 163.1113, 213.0870, 446.2286, 553.3093, 599.3552]
MPIONS = [184.06, 215.1192, 243.1127, 134.0961, 181.1331, 169.0967, 150.0912, 167.1178, 454.15]
ARIONS = [112.0964, 140.1066, 334.0838, 300.1232, 284.1268, 250.1440, 221.1646, 281.1914, 314.2199]
ABIONS = [114.0550, 164.1069, 192.1019, 233.1285, 263.1391, 362.2075]
MGIONS = [100.1122, 134.0727, 168.0341, 201.9955, 114.1278, 148.0888, 182.0494, 128.1423, 162.1039,
          196.0639, 142.1590, 176.1195, 210.0795]

# ------------------------------------------------------------------
# Ion labels dictionary
# ------------------------------------------------------------------
ION_TO_LABEL = {
    184.06:   "MP NMeTyrCl 184.06",
    215.1192: "MP AhpPhe-CO-H2O 215.12",
    243.1127: "MP AhpPhe-H2O 243.11",
    134.0961: "MP NMePhe 134.10",
    181.1331: "MP AhpLxx-CO-H2O 181.1331",
    169.0967: "MP AhpThr-CO-H2O 169.10",
    150.0912: "MP NMeTyr 150.0912",
    167.1178: "MP AhpVal-CO-H2O 167.12",
    454.15:   "MP Phe-Ahp-complete 454.15",

    135.0804: "MC adda fragment 135.0804",
    163.1113: "MC adda fragment 163.1113",
    213.0870: "MC Mdha 213.0870",
    446.2286: "MC adda (163) + Glu + Mdha + Ala",
    553.3093: "MC Adha + Ala + Leu + MeAsp + Arg",
    599.3552: "MC [MeAsp-Arg-Adda + H]+, [Asp-Har-Adda + H]+, or [Arg-Adda-Glu + H]+ 599.3552",

    112.0964: "AR Choi-h2o 112.0964",
    140.1066: "AR Choi 140.1066",
    334.0838: "AR Cl-Hpla-Tyr-CO 334.0838",
    300.1232: "AR Hpla-Tyr-CO 300.1232",
    284.1268: "AR Hlpa-Phe-CO 284.1268",
    250.1440: "AR Hlpa-Leu-CO 250.1440",
    221.1646: "AR Hlpa-Leu 221.1646",
    281.1914: "AR Agma 281.1914",
    314.2199: "AR OH-Choi-Agma 314.2199",

    114.0550: "AB CO + Arg- CH5H3 -CO 114.0550",
    164.1069: "AB NMe-HTyr 164.1069",
    192.1019: "AB CO-NMe-Hty 192.1019",
    233.1285: "AB Phe + MeAla +H 233.1285",
    263.1391: "AB HTyr +MeAla +H 263.1391",
    362.2075: "AB HTyr+ MeAla +Val +H 362.2075",

    100.1122: "MG Ahoa 100.1122",
    134.0727: "MG Ahoa (Cl) 134.0727",
    168.0341: "MG Ahoa (Di Cl) 168.0341",
    201.9955: "MG Ahoa (Tri Cl) 201.9955",
    114.1278: "MG NMe Ahoa 114.1278",
    148.0888: "MG NMe Ahoa (Cl) 148.0888",
    182.0494: "MG NMe Ahoa (Di Cl) 182.0494",
    128.1423: "MG Ahda 128.1423",
    162.1039: "MG Ahda (Cl) 162.1039",
    196.0639: "MG Ahda (Di Cl) 196.0639",
    142.1590: "MG NMe Ahda 142.1590",
    176.1195: "MG NMe Ahda (Cl) 176.1195",
    210.0795: "MG NMe Ahda (Di Cl) 210.0795",
}

# ------------------------------------------------------------------
# Class configs
# ------------------------------------------------------------------
CLASS_CONFIGS = {
    "AB": {
        "IONS": ABIONS,
        "ADD_LABELS_FN": mu.add_AB_labels,
        "LABEL_COL": "CyanopeptideClass_AB",
        "LIB_CLASS_FILTER": "Anabaenopeptin",
    },
    "MC": {
        "IONS": MCIONS,
        "ADD_LABELS_FN": mu.add_MC_labels,
        "LABEL_COL": "CyanopeptideClass_MC",
        "LIB_CLASS_FILTER": "Microcystin",
    },
    "AR": {
        "IONS": ARIONS,
        "ADD_LABELS_FN": mu.add_AR_labels,
        "LABEL_COL": "CyanopeptideClass_AR",
        "LIB_CLASS_FILTER": "Aeruginosin",
    },
    "MP": {
        "IONS": MPIONS,
        "ADD_LABELS_FN": mu.add_MP_labels,
        "LABEL_COL": "CyanopeptideClass_MP",
        "LIB_CLASS_FILTER": "Cyanopeptolin",
    },
    "MG": {
        "IONS": MGIONS,
        "ADD_LABELS_FN": mu.add_MG_labels,
        "LABEL_COL": "CyanopeptideClass_MG",
        "LIB_CLASS_FILTER": "Microginin",
    },
}


# ------------------------------------------------------------------
# Utilities
# ------------------------------------------------------------------
class Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for s in self.streams:
            s.write(data)
            s.flush()

    def flush(self):
        for s in self.streams:
            s.flush()


def _make_pipeline_log(output_root, class_tag=None):
    output_root = Path(output_root)
    pipelinelog_dir = output_root / "pipeline_log"
    pipelinelog_dir.mkdir(parents=True, exist_ok=True)

    ts = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    tag = f"_{class_tag}" if class_tag else ""
    log_path = pipelinelog_dir / f"pipeline_log{tag}_{ts}.txt"
    return pipelinelog_dir, log_path


def _parse_scan_number(spec: dict) -> int | None:
    sid = spec.get("id", "")
    if isinstance(sid, str) and "scan=" in sid:
        try:
            return int(sid.split("scan=")[-1].split()[0])
        except Exception:
            return None
    return None


def _get_rt_minutes(spec: dict, assume_time_unit: str = "min") -> float | None:
    rt = spec.get("scan start time", None)
    if rt is None:
        try:
            sc = spec.get("scanList", {}).get("scan", [])
            if isinstance(sc, dict):
                sc = [sc]
            if sc:
                rt = sc[0].get("scan start time", None)
        except Exception:
            rt = None

    if rt is None:
        return None

    rt = float(rt)
    if assume_time_unit.lower().startswith("sec"):
        rt = rt / 60.0
    return rt


def _guess_polarity(spec: dict) -> int | None:
    txt = str(spec).lower()
    if "positive scan" in txt:
        return 1
    if "negative scan" in txt:
        return -1
    return None


def extract_ms1_points_from_mzml(
    mzml_path: str | Path,
    *,
    mz_round: int | None = None,
    intensity_min: float = 0.0,
    assume_time_unit: str = "min",
    rt_window: tuple[float, float] | None = None,
) -> pd.DataFrame:
    mzml_path = str(mzml_path)
    src = Path(mzml_path).name

    out_source = []
    out_scan = []
    out_rt = []
    out_mz = []
    out_i = []
    out_ms1 = []
    out_pol = []

    with mzml.MzML(mzml_path) as reader:
        for spec in reader:
            mslevel = spec.get("ms level") or spec.get("msLevel") or spec.get("mslevel")
            if mslevel is None or int(mslevel) != 1:
                continue

            rt = _get_rt_minutes(spec, assume_time_unit=assume_time_unit)
            if rt is None:
                continue

            if rt_window is not None:
                rt_lo, rt_hi = rt_window
                if rt_lo is not None and rt < rt_lo:
                    continue
                if rt_hi is not None and rt > rt_hi:
                    continue

            scan_num = _parse_scan_number(spec)
            pol = _guess_polarity(spec)

            mzs = spec.get("m/z array", None)
            ints = spec.get("intensity array", None)
            if mzs is None or ints is None:
                continue

            mzs = np.asarray(mzs, dtype=float)
            ints = np.asarray(ints, dtype=float)

            if intensity_min > 0:
                keep = ints >= float(intensity_min)
                if not np.any(keep):
                    continue
                mzs = mzs[keep]
                ints = ints[keep]

            if mz_round is not None:
                mzs = np.round(mzs, int(mz_round))

            n = mzs.size
            if n == 0:
                continue

            out_source.extend([src] * n)
            out_scan.extend([scan_num] * n)
            out_rt.extend([rt] * n)
            out_mz.extend(mzs.tolist())
            out_i.extend(ints.tolist())
            out_ms1.extend([1] * n)
            out_pol.extend([pol] * n)

    return pd.DataFrame(
        {
            "source_file": out_source,
            "scan": out_scan,
            "rt": out_rt,
            "precmz": out_mz,
            "i": out_i,
            "mslevel": out_ms1,
            "polarity": out_pol,
        }
    )


def build_ms1_points(
    files: list[str],
    out_dir: Path,
    *,
    mz_round: int | None,
    intensity_min: float,
    assume_time_unit: str,
    rt_window: tuple[float, float] | None,
) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    ts = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    out_file = out_dir / f"ms1_points_{ts}.csv"

    all_dfs = []
    for p in files:
        print(f"[INFO] extracting MS1: {Path(p).name}")
        df = extract_ms1_points_from_mzml(
            p,
            mz_round=mz_round,
            intensity_min=intensity_min,
            assume_time_unit=assume_time_unit,
            rt_window=rt_window,
        )
        print(f"   rows: {len(df)}")
        all_dfs.append(df)

    ms1_points = pd.concat(all_dfs, ignore_index=True) if all_dfs else pd.DataFrame()
    ms1_points.to_csv(out_file, index=False)
    print(f"[DONE] wrote: {out_file}")
    return out_file


def discover_files(input_dir: Path, pattern: str = "*.mzML") -> list[str]:
    return sorted(str(p) for p in input_dir.glob(pattern))


def validate_class_labels():
    for class_tag, cfg in CLASS_CONFIGS.items():
        missing = [ion for ion in cfg["IONS"] if ion not in ION_TO_LABEL]
        if missing:
            print(f"[WARNING] {class_tag} missing labels for: {missing}")


# ------------------------------------------------------------------
# Main pipeline code
# ------------------------------------------------------------------
def run_class_pipeline(
    *,
    files: list[str],
    class_tag: str,
    metadata_path: Path | None,
    CyanoMetDBLibrary: Path | None,
    ms1_points_file: Path | None,
    output_root: Path,
    tol: float,
    polarity: str | None,
    rt_window: tuple[float, float] | None,
    ref_mz: float | None,
    ref_tol: float,
    ref_rt_window: tuple[float, float] | None,
    do_blank_filter: bool,
    do_batch_correct: bool,
    tol_da: float,
):
    cfg_class = CLASS_CONFIGS[class_tag]

    IONS = cfg_class["IONS"]
    ADD_LABELS_FN = cfg_class["ADD_LABELS_FN"]
    LABEL_COL = cfg_class["LABEL_COL"]
    LIB_CLASS_FILTER = cfg_class["LIB_CLASS_FILTER"]

    ts = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    run_dir = output_root / f"{class_tag}_run_{ts}"
    run_dir.mkdir(parents=True, exist_ok=True)

    raw_plots_dir = run_dir / "raw_ms2_plots"
    auc_dir = run_dir / "AUC_outputs"
    qc_dir = run_dir / "QC_audit"
    dotplot_dir = run_dir / "DotPlot_individual_ion_plots"
    adduct_dir = run_dir / "Adduct_and_summary_outputs"
    match_dir = run_dir / "CyanoMetDB_matches_out"

    for d in [
        raw_plots_dir,
        auc_dir,
        qc_dir,
        dotplot_dir,
        adduct_dir,
        match_dir,
        raw_plots_dir / "RT_histogram_plots",
        raw_plots_dir / "plots_diagnostic_ion_counts",
        raw_plots_dir / "Heatmap_Individual_ion_plots",
        raw_plots_dir / "Precursor_rt_plot",
    ]:
        d.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] Running class: {class_tag}")
    print(f"[INFO] Output directory: {run_dir.resolve()}")

    # ------------------------------------------------------------
    # 1) MassQL MS2 diagnostic ion search
    # ------------------------------------------------------------
    ind_hits = mu.run_across_files_individual(
        files,
        IONS,
        tol_mz=tol,
        polarity=polarity,
        rt_window=rt_window,
    )

    ind_hits_l = ADD_LABELS_FN(ind_hits, ION_TO_LABEL, label_col=LABEL_COL)

    ind_csv = run_dir / f"individual_hits_{class_tag}_{ts}.csv"
    ind_hits_l.to_csv(ind_csv, index=False)
    print(f"[DONE] wrote: {ind_csv}")

    expected_labels = [ION_TO_LABEL[ion] for ion in IONS if ion in ION_TO_LABEL]

    scan_has_df = make_scan_level_has_table(
        ind_hits_l,
        label_col=LABEL_COL,
        expected_labels=expected_labels,
    )

    scan_has_file = run_dir / f"scan_level_has_flags_{class_tag}_{ts}.csv"
    scan_has_df.to_csv(scan_has_file, index=False)
    print(f"[DONE] wrote: {scan_has_file}")

    # ------------------------------------------------------------
    # 2) Plots
    # ------------------------------------------------------------
    plot_rt_histograms(
        ind_hits_l,
        ION_TO_LABEL,
        out_dir_root=str(raw_plots_dir / "RT_histogram_plots"),
    )

    plot_indiv_counts(
        ind_hits_l,
        out_dir=str(raw_plots_dir / "plots_diagnostic_ion_counts"),
    )

    plot_precursor_rt(
        ind_hits_l,
        ion_to_label=ION_TO_LABEL,
        save=True,
        out_dir=str(raw_plots_dir / "Precursor_rt_plot"),
        show=False,
    )

    plot_heatmaps(
        ind_hits_l,
        ion_to_label=ION_TO_LABEL,
        save=True,
        out_dir=str(raw_plots_dir / "Heatmap_Individual_ion_plots"),
        fmt="png",
    )

    # ------------------------------------------------------------
    # 3) Build merged precursor summary
    # ------------------------------------------------------------
    indiv_summary = make_summary_ind(ind_hits_l, ion_to_label=ION_TO_LABEL)

    # ------------------------------------------------------------
    # 4) Load MS1 points
    # ------------------------------------------------------------
    if ms1_points_file is None:
        print("No MS1 points file provided. Automatically extracting MS1 from files...")

        ms1_points_file = build_ms1_points(
            files,
            output_root / "MS1_points",
            mz_round=None,
            intensity_min=0.0,
            assume_time_unit="min",
            rt_window=rt_window,
        )

        print(f"[INFO] built ms1_points_file: {ms1_points_file}")

    ms1_points = load_ms1_points(ms1_points_file)

    # ------------------------------------------------------------
    # 5) Metadata + QC config
    # ------------------------------------------------------------
    cfg = QCConfig(
        metadata_path=metadata_path,
        blank_col="sample_type",
        blank_value="blank",
        batch_col="batch",
        ref_mz=ref_mz,
        ref_tol=ref_tol,
        ref_rt_window=ref_rt_window if ref_mz is not None else None,
        apply_blank_filter=do_blank_filter,
        blank_ratio_thresh=3.0,
        blank_stat="median",
        apply_batch_correction=do_batch_correct,
        batch_use_nonblank_only=True,
    )

    metadata = load_metadata(cfg) if metadata_path is not None else None

    # ------------------------------------------------------------
    # 6) Compute MS1 AUC
    # ------------------------------------------------------------
    perfile_summary = explode_matches_to_per_file(indiv_summary)

    perfile_summary_auc = add_ms1_auc_from_points(
        perfile_summary,
        ms1_points,
        tol_mz=tol,
        rt_pad=0.25,
        polarity=polarity if isinstance(polarity, int) else None,
        intensity_col="i",
        mz_col="precmz",
    )

    # ------------------------------------------------------------
    # 6.1) Reference normalization
    # ------------------------------------------------------------
    if ref_mz is not None:
        ref_auc_df = compute_ref_auc_df(
            ms1_points,
            ref_mz=cfg.ref_mz,
            ref_tol=cfg.ref_tol,
            ref_rt_window=cfg.ref_rt_window,
            polarity=polarity if isinstance(polarity, int) else None,
        )
        perfile_raw = add_reference_normalization(
            perfile_summary_auc,
            ref_auc_df,
            ref_mz=cfg.ref_mz,
        )
    else:
        perfile_raw = perfile_summary_auc.copy()

    # ------------------------------------------------------------
    # 6.2) Blank filtering
    # ------------------------------------------------------------
    feature_cols = [c for c in ["merged_precmz", "rt_median", "charge"] if c in perfile_raw.columns]

    if do_blank_filter and metadata is not None:
        intensity_col = "ms1_auc_over_ref" if "ms1_auc_over_ref" in perfile_raw.columns else "ms1_auc"

        has_batch = cfg.batch_col in metadata.columns

        if has_batch:
            print(f"[INFO] Using batch-aware blank filter with batch column '{cfg.batch_col}'.")
            perfile_clean, keep_feats, removed_feats = blank_filter_perfile_table_by_batch(
                perfile_raw,
                metadata,
                cfg,
                intensity_col=intensity_col,
                feature_cols=feature_cols,
            )
        else:
            print(f"[INFO] No '{cfg.batch_col}' column in metadata — using non-batch blank filter.")
            perfile_clean, keep_feats, removed_feats = blank_filter_perfile_table(
                perfile_raw,
                metadata,
                cfg,
                intensity_col=intensity_col,
                feature_cols=feature_cols,
            )

        if "batch" in perfile_clean.columns:
            print("[DEBUG] dropping batch before batch correction")
            perfile_clean = perfile_clean.drop(columns=["batch"])

        print("[DEBUG] after drop_blank_rows")
        print("[DEBUG] perfile_clean columns:", perfile_clean.columns.tolist())
        print("[DEBUG] has batch?", "batch" in perfile_clean.columns)
        print("[DEBUG] batch-like columns:", [c for c in perfile_clean.columns if "batch" in c.lower()])
        print("[DEBUG] shape:", perfile_clean.shape)
    else:
        perfile_clean = perfile_raw.copy()
        keep_feats = None
        removed_feats = None

    # ------------------------------------------------------------
    # 6.3) Batch correction
    # ------------------------------------------------------------
    if do_batch_correct and metadata is not None:
        intensity_col = "ms1_auc_over_ref" if "ms1_auc_over_ref" in perfile_clean.columns else "ms1_auc"
        perfile_clean_bc = apply_batch_correction_median_scaling(
            perfile_clean,
            metadata,
            cfg,
            intensity_col=intensity_col,
        )
    else:
        perfile_clean_bc = perfile_clean.copy()

    # ------------------------------------------------------------
    # 6.4) Pooled tables
    # ------------------------------------------------------------
    pooled_raw = pool_auc_back_to_matches(perfile_raw)

    pooled_clean = (
        pool_auc_back_to_matches(perfile_clean_bc)
        if perfile_clean_bc is not None and not perfile_clean_bc.empty
        else pooled_raw.iloc[0:0].copy()
    )

    # remove feature-level MS2 flags from quant tables
    has_cols_raw = [c for c in perfile_raw.columns if c.startswith("has_")]
    perfile_raw = perfile_raw.drop(columns=has_cols_raw, errors="ignore")

    has_cols_clean = [c for c in perfile_clean_bc.columns if c.startswith("has_")]
    perfile_clean_bc = perfile_clean_bc.drop(columns=has_cols_clean, errors="ignore")

    # ------------------------------------------------------------
    # 7) Save AUC outputs
    # ------------------------------------------------------------
    raw_perfile_out = auc_dir / f"{class_tag}_perfile_RAW_{ts}.csv"
    raw_pooled_out = auc_dir / f"{class_tag}_pooled_RAW_{ts}.csv"
    clean_perfile_out = auc_dir / f"{class_tag}_perfile_CLEAN_{ts}.csv"
    clean_pooled_out = auc_dir / f"{class_tag}_pooled_CLEAN_{ts}.csv"

    perfile_raw.to_csv(raw_perfile_out, index=False)
    pooled_raw.to_csv(raw_pooled_out, index=False)
    perfile_clean_bc.to_csv(clean_perfile_out, index=False)
    pooled_clean.to_csv(clean_pooled_out, index=False)

    print(f"[DONE] Saved RAW per-file AUC table: {raw_perfile_out}")
    print(f"[DONE] Saved RAW pooled AUC table:  {raw_pooled_out}")
    print(f"[DONE] Saved CLEAN per-file AUC table: {clean_perfile_out}")
    print(f"[DONE] Saved CLEAN pooled AUC table:  {clean_pooled_out}")

    removed_feats_out = qc_dir / f"{class_tag}_blank_removed_features_{ts}.csv"

    # ------------------------------------------------------------
    # Build detailed removed-feature audit table
    # ------------------------------------------------------------
    if removed_feats is not None and not removed_feats.empty and metadata is not None:

        meta_key = "source_file" if "source_file" in metadata.columns else "filename"

        meta_for_merge = metadata[[meta_key, cfg.blank_col, cfg.batch_col]].copy()
        meta_for_merge[meta_key] = (
            meta_for_merge[meta_key]
            .astype(str)
            .str.strip()
            .apply(lambda x: Path(x).name.strip())
        )
        meta_for_merge = meta_for_merge.rename(columns={meta_key: "source_file"})

        tmp = perfile_raw.copy()
        tmp["source_file"] = (
            tmp["source_file"]
            .astype(str)
            .str.strip()
            .apply(lambda x: Path(x).name.strip())
        )

        tmp = tmp.merge(meta_for_merge, on="source_file", how="left")

        tmp["is_blank"] = (
            tmp[cfg.blank_col]
            .astype(str)
            .str.lower()
            .eq(str(cfg.blank_value).lower())
        )

        print(
            tmp.loc[
                tmp["source_file"].str.contains("blank", case=False, na=False),
                ["source_file", cfg.blank_col, cfg.batch_col]
            ]
            .drop_duplicates()
            .sort_values("source_file")
        )

        group_cols = feature_cols + ([cfg.batch_col] if cfg.batch_col in removed_feats.columns else [])

        def summarize_removed_feature(g):
            intensity_col_local = "ms1_auc_over_ref" if "ms1_auc_over_ref" in g.columns else "ms1_auc"

            blank = g.loc[g["is_blank"], intensity_col_local]
            samp = g.loc[~g["is_blank"], intensity_col_local]

            bmed = blank.median()
            smed = samp.median()

            ratio = np.nan
            if pd.notna(bmed) and bmed > 0:
                ratio = smed / bmed

            return pd.Series({
                "blank_files": ";".join(sorted(g.loc[g["is_blank"], "source_file"].dropna().unique())),
                "sample_files": ";".join(sorted(g.loc[~g["is_blank"], "source_file"].dropna().unique())),
                "blank_median_auc": bmed,
                "sample_median_auc": smed,
                "sample_blank_ratio": ratio,
            })

        print(
            tmp.loc[
                tmp["source_file"].str.contains("blank", case=False, na=False),
                ["source_file", cfg.blank_col, cfg.batch_col]
            ].drop_duplicates().sort_values("source_file")
        )

        removed_with_stats = (
            tmp.merge(removed_feats, on=group_cols, how="inner")
               .groupby(group_cols, dropna=False)
               .apply(summarize_removed_feature)
               .reset_index()
        )

        removed_with_stats.to_csv(removed_feats_out, index=False)
        print(f"[DONE] Saved detailed blank-REMOVED features table: {removed_feats_out}")

    else:
        print("[INFO] No removed-features table to save.")

    save_qc_audit(
        out_dir=str(qc_dir),
        perfile_raw=perfile_raw,
        perfile_clean=perfile_clean_bc,
        keep_table=keep_feats,
        removed_table=removed_feats,
        tag=class_tag,
    )

    # ------------------------------------------------------------
    # 8) Downstream scatter plot
    # ------------------------------------------------------------
    plot_indiv_scatter(
        indiv_summary,
        out_dir=str(dotplot_dir),
        fmt="png",
    )

    # ------------------------------------------------------------
    # 9) Adduct pipeline
    # ------------------------------------------------------------
    def make_summary_ind_with_labels(df, *args, **kwargs):
        return make_summary_ind(df, *args, ion_to_label=ION_TO_LABEL, **kwargs)

    merged_summary, merged_edges, G = ap.run_merged(
        ind_hits_l,
        make_summary_ind=make_summary_ind_with_labels,
        af_module=af,
        mz_col="merged_precmz",
        charge_col="charge",
        rt_col="rt_median",
        out_dir=str(adduct_dir),
        save_graph=True,
    )

    merged_summary_file = adduct_dir / f"indiv_merged_summary_{class_tag}_{ts}.csv"
    merged_summary.to_csv(merged_summary_file, index=False)
    print(f"[DONE] Saved merged summary: {merged_summary_file}")

    # ------------------------------------------------------------
    # 9.1) Sum MS2 scan intensities
    # ------------------------------------------------------------
    summary_with_i_df, summary_with_i_file = sum_intensities(
        summary_file=str(merged_summary_file),
        hits_file=str(ind_csv),
        output_dir=str(adduct_dir),
    )

    # ------------------------------------------------------------
    # 9.2) Plot MS2 intensity tilemap
    # ------------------------------------------------------------
    plot_has_tilemap(
        summary_file=summary_with_i_file,
        output_dir=str(adduct_dir),
        ion_to_label=ION_TO_LABEL,
    )

# ------------------------------------------------------------
# 10) CyanoMetDB matching
# ------------------------------------------------------------
    if CyanoMetDBLibrary is not None:
        lib_df = load_library(
            CyanoMetDBLibrary,
            class_filter=LIB_CLASS_FILTER,
            sheet_index=1,
        )

        # Match on merged summary (keeps has_* columns needed for unknowns)
        ms1_df_for_match = read_any_table(str(merged_summary_file))
        ms1_sel = select_ms1_columns(ms1_df_for_match)
        matches = match_ms1_to_lib(ms1_sel, lib_df, tol_da=tol_da)

        # Save RAW matches (before CLEAN filtering)
        ts_raw = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        out_dir_raw, paths_raw = write_outputs(
            matches,
            matches[matches["Compound identifier"].notna()],
            out_dir=str(match_dir),
            ts=f"{ts_raw}_RAW",
        )

        print("[INFO] Saved RAW CyanoMetDB matches.")
        # ------------------------------------------------------------
        # CLEAN filter: keep only CyanoMetDB rows whose features survived pooled_clean
        # ------------------------------------------------------------
        if pooled_clean is not None and (not pooled_clean.empty):
            key_cols = [c for c in ["merged_precmz", "rt_median"] if c in matches.columns and c in pooled_clean.columns]

            if key_cols:
                clean_keys = pooled_clean[key_cols].drop_duplicates()

                matches_clean = matches.merge(
                    clean_keys,
                    on=key_cols,
                    how="inner",
                )
                print(f"[INFO] Filtered CyanoMetDB results to CLEAN features: {len(matches)} -> {len(matches_clean)}")
            else:
                matches_clean = matches.copy()
                print("[INFO] No shared key columns between matches and pooled_clean — skipping CLEAN filter.")
        else:
            matches_clean = matches.copy()
            print("[INFO] pooled_clean not available — skipping CLEAN filter.")

        # Add clean quant columns back onto the cleaned match table
        if pooled_clean is not None and (not pooled_clean.empty) and ("merged_precmz" in pooled_clean.columns):
            quant_cols = [
                c for c in pooled_clean.columns
                if c in {"merged_precmz", "rt_median", "ms1_auc", "ms1_auc_over_ref", "ms1_auc_batchcorr", "ms1_auc_clean", "ms1_auc_sum_across_files"}
            ]
            quant_cols = list(dict.fromkeys(quant_cols))

            merge_keys = [c for c in ["merged_precmz", "rt_median"] if c in matches_clean.columns and c in pooled_clean.columns]

            if merge_keys:
                matches_clean = matches_clean.merge(
                    pooled_clean[quant_cols].drop_duplicates(subset=merge_keys),
                    on=merge_keys,
                    how="left",
                    suffixes=("", "_quant"),
                )
                print("[INFO] Merged CLEAN quant columns onto cleaned CyanoMetDB results from pooled_clean.")
            else:
                print("[INFO] No shared merge keys for CLEAN quant columns — skipping quant merge.")
        else:
            print("[INFO] pooled_clean not available — cleaned matches will not include clean quant columns.")

        matched_only = matches_clean[matches_clean["Compound identifier"].notna()].copy()
        putative_novel_only = matches_clean[matches_clean["Compound identifier"].isna()].copy()

        ts2 = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        out_dir_ts, paths = write_outputs(
            matches_clean,
            matched_only,
            out_dir=str(match_dir),
            ts=f"{ts2}_CLEAN",
        )

        putative_novel_out = Path(out_dir_ts) / f"putative_novel_only_{ts2}_CLEAN.csv"
        putative_novel_only.to_csv(putative_novel_out, index=False)
        print(f"[DONE] Saved CLEAN putative novel features: {putative_novel_out}")

        plot_matched_tiles(matched_only, out_dir_ts=out_dir_ts, ts=f"{ts2}_CLEAN")
    # flatten timestamped subfolders into base folders
    flatten_subfolders(run_dir)


# trying to merge these together so we don't have repeat folders created... not
# working. code runs but the folders are not merging.
def merge_timestamped_subfolders(run_dir: str | Path, remove_empty: bool = True):
    """
    Merge files from timestamped duplicate folders into their base folder.

    Example:
      Adduct_and_summary_outputs_26-03-22_18-23-02
      -> Adduct_and_summary_outputs

    It only acts when the base folder already exists.
    If a filename already exists, it appends '_dup1', '_dup2', etc.
    """
    run_dir = Path(run_dir)

    if not run_dir.exists():
        print(f"[INFO] merge step skipped; run_dir not found: {run_dir}")
        return

    folders = [p for p in run_dir.rglob("*") if p.is_dir()]
    folders_sorted = sorted(folders, key=lambda p: len(p.parts), reverse=True)

    for dup_dir in folders_sorted:
        name = dup_dir.name

        # Find base folder by stripping the last 2 underscore-separated chunks
        # Adduct_and_summary_outputs_26-03-22_18-23-02 -> Adduct_and_summary_outputs
        parts = name.split("_")
        if len(parts) < 3:
            continue

        candidate_base_name = "_".join(parts[:-2])
        if candidate_base_name == name:
            continue

        base_dir = dup_dir.parent / candidate_base_name

        if not base_dir.exists() or base_dir == dup_dir:
            continue

        print("[INFO] Merging:")
        print(f"       from: {dup_dir}")
        print(f"         to: {base_dir}")

        for item in dup_dir.iterdir():
            dest = base_dir / item.name

            if dest.exists():
                stem = dest.stem
                suffix = dest.suffix
                k = 1
                while True:
                    alt = base_dir / f"{stem}_dup{k}{suffix}"
                    if not alt.exists():
                        dest = alt
                        break
                    k += 1

            shutil.move(str(item), str(dest))
            print(f"       moved: {item.name} -> {dest.name}")

        if remove_empty:
            try:
                dup_dir.rmdir()
                print(f"       removed empty folder: {dup_dir}")
            except OSError:
                print(f"       folder not empty, kept: {dup_dir}")
    print(f"[DONE] Pipeline complete for {class_tag}")
    print(f"[DONE] All class-specific outputs are in: {run_dir.resolve()}")
    print(f"[DONE] MS1 points file used: {ms1_points_file}")
    return run_dir


def get_default_lib_path():
    """
    Get path to the bundled DB file.  It might be in two locations,
    depending if we're running out of the git repository (or similar) or if
    we're a properly installed package.

    The return value is a context manager.
    """
    fullpath = Path(BUNDLED_LIB_PATH)
    data_pkg = fullpath.parent.parts[0]  # min of two parts
    path = Path(*fullpath.parts[1:])
    for anchor in [data_pkg, 'cpm.' + data_pkg]:
        try:
            if importlib.resources.is_resource(anchor, path):
                return importlib.resources.path(anchor, path)
        except ModuleNotFoundError:
            pass
    print(f'[NOTICE] Could not find bundled {BUNDLED_LIB_PATH}')
    return nullcontext(None)


# ------------------------------------------------------------------
# Pipeline notebook
# ------------------------------------------------------------------
def run_pipeline_notebook(
    *,
    class_tag: str,
    files: list[str] | None = None,
    input_dir: str | Path | None = None,
    pattern: str = "*.mzML",
    metadata_path: str | Path | None = None,
    CyanoMetDBLibrary: str | Path | None = None,
    ms1_points_file: str | Path | None = None,
    output_root: str | Path = "pipeline_outputs",
    extract_ms1: bool = False,
    tol: float = 0.01,
    polarity: str | None = "POSITIVE",
    rt_window: tuple[float, float] | None = (2.0, 25.0),
    ref_rt_window: tuple[float, float] | None = (1.0, 25.0),
    ref_mz: float | None = None,
    ref_tol: float = 0.01,
    do_blank_filter: bool = True,
    do_batch_correct: bool = True,
    tol_da: float = 0.05,
    mz_round: int | None = None,
    intensity_min: float = 0.0,
    assume_time_unit: str = "min",
):
    output_root = Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    pipelinelog_dir, log_path = _make_pipeline_log(output_root, class_tag=class_tag)

    with open(log_path, "w") as log_f:
        tee_out = Tee(sys.stdout, log_f)
        tee_err = Tee(sys.stderr, log_f)

        with redirect_stdout(tee_out), redirect_stderr(tee_err):
            print("=" * 80)
            print("CPM PIPELINE AUDIT")
            print(f"class_tag: {class_tag}")
            print(f"output_root: {output_root}")
            print(f"log_path: {log_path}")
            print("=" * 80)

            validate_class_labels()

            print("rt_window:", rt_window, type(rt_window))
            print("ref_rt_window:", ref_rt_window, type(ref_rt_window))
            print("ms1_points_file:", ms1_points_file, type(ms1_points_file))

            if files is not None:
                files = [str(Path(f)) for f in files]
            elif input_dir is not None:
                files = discover_files(Path(input_dir), pattern=pattern)
            else:
                raise ValueError("Provide either files=[...] or input_dir=...")

            if not files:
                raise ValueError("No mzML files found.")

            print(f"[INFO] n_files: {len(files)}")
            for f in files:
                print(f"   - {f}")

            metadata_path = Path(metadata_path) if metadata_path is not None else None
            if CyanoMetDBLibrary is None:
                estack = ExitStack()  # left open until end-of-program
                CyanoMetDBLibrary = estack.enter_context(get_default_lib_path())
            elif isinstance(CyanoMetDBLibrary, str):
                CyanoMetDBLibrary = Path(CyanoMetDBLibrary)
            ms1_points_file = Path(ms1_points_file) if ms1_points_file is not None else None

            print(f"[INFO] metadata_path: {metadata_path}")
            print(f"[INFO] CyanoMetDBLibrary: {CyanoMetDBLibrary}")
            print(f"[INFO] extract_ms1: {extract_ms1}")

            if extract_ms1:
                ms1_points_file = build_ms1_points(
                    files,
                    output_root / "MS1_points",
                    mz_round=mz_round,
                    intensity_min=intensity_min,
                    assume_time_unit=assume_time_unit,
                    rt_window=rt_window,
                )
                print(f"[INFO] built ms1_points_file: {ms1_points_file}")

            if ms1_points_file is None:
                raise ValueError("ms1_points_file is required unless extract_ms1=True.")

            if class_tag == "ALL":
                class_tags_to_run = list(CLASS_CONFIGS.keys())
            else:
                if class_tag not in CLASS_CONFIGS:
                    raise ValueError(f"Unknown class_tag: {class_tag}")
                class_tags_to_run = [class_tag]

            print(f"[INFO] class_tags_to_run: {class_tags_to_run}")

            run_dirs = []
            for tag in class_tags_to_run:
                print(f"\n[INFO] starting class_tag={tag}")

                run_dir = run_class_pipeline(
                    files=files,
                    class_tag=tag,
                    metadata_path=metadata_path,
                    CyanoMetDBLibrary=CyanoMetDBLibrary,
                    ms1_points_file=ms1_points_file,
                    output_root=output_root,
                    tol=tol,
                    polarity=polarity,
                    rt_window=rt_window,
                    ref_rt_window=ref_rt_window,
                    ref_mz=ref_mz,
                    ref_tol=ref_tol,
                    do_blank_filter=do_blank_filter,
                    do_batch_correct=do_batch_correct,
                    tol_da=tol_da,
                )
                run_dirs.append(run_dir)
                print(f"[INFO] finished class_tag={tag} -> {run_dir}")

            result = {
                "class_tags_run": class_tags_to_run,
                "files": files,
                "ms1_points_file": ms1_points_file,
                "run_dirs": run_dirs,
                "pipeline_log_file": log_path,
                "pipeline_log_dir": pipelinelog_dir,
            }

            print("\n[DONE] Pipeline completed")
            print(f"[DONE] QC audit log saved to: {log_path}")

            return result


# ------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------
def main():
    estack = ExitStack()
    #  Obtain a pathlib.Path to the bundled DB/lib, the exit stack stays open
    #  until end of program.
    default_lib_path = estack.enter_context(get_default_lib_path())

    parser = argparse.ArgumentParser(description="Cyanopeptide pipeline CLI")

    parser.add_argument(
        "--class-tag",
        choices=list(CLASS_CONFIGS.keys()) + ["ALL"],
        required=True,
        help="Which class to run",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        help="Directory containing mzML files",
    )
    parser.add_argument(
        "--files",
        nargs="*",
        help="Explicit list of mzML files",
    )
    parser.add_argument(
        "--pattern",
        default="*.mzML",
        help="Glob pattern for mzML discovery inside --input-dir",
    )

    parser.add_argument("--metadata", type=Path, default=None)
    parser.add_argument(
        "--CyanoMetDBLibrary",
        type=Path,
        default=default_lib_path,
        required=(default_lib_path is None),
        help="CyanoMetDB library Excel file " + (
            "(required)" if default_lib_path is None
            else f"(default: bundled {default_lib_path.name})"
        ),
    )
    parser.add_argument("--ms1-points-file", type=Path, default=None)
    parser.add_argument("--output-root", type=Path, default=Path("pipeline_outputs"))

    parser.add_argument("--tol", type=float, default=0.01)
    parser.add_argument("--polarity", default="POSITIVE")
    parser.add_argument("--rt-min", type=float, default=2.0)
    parser.add_argument("--rt-max", type=float, default=25.0)

    parser.add_argument("--ref-rt-min", type=float, default=None)
    parser.add_argument("--ref-rt-max", type=float, default=None)

    parser.add_argument("--ref-mz", type=float, default=None)
    parser.add_argument("--ref-tol", type=float, default=0.01)

    parser.add_argument("--blank-filter", action="store_true")
    parser.add_argument("--batch-correct", action="store_true")

    parser.add_argument("--tol-da", type=float, default=0.05)

    parser.add_argument("--extract-ms1", action="store_true")
    parser.add_argument("--mz-round", type=int, default=None)
    parser.add_argument("--intensity-min", type=float, default=0.0)
    parser.add_argument("--assume-time-unit", default="min")

    args = parser.parse_args()

    validate_class_labels()

    if args.files:
        files = args.files
    elif args.input_dir:
        files = discover_files(args.input_dir, pattern=args.pattern)
    else:
        raise SystemExit("Provide either --files or --input-dir")

    if not files:
        raise SystemExit("No mzML files found.")

    args.output_root.mkdir(parents=True, exist_ok=True)

    rt_window = None if args.rt_min is None and args.rt_max is None else (args.rt_min, args.rt_max)

    ref_rt_window = (
        None
        if args.ref_rt_min is None and args.ref_rt_max is None
        else (args.ref_rt_min, args.ref_rt_max)
    )
    ms1_points_file = args.ms1_points_file
    if args.extract_ms1:
        ms1_points_file = build_ms1_points(
            files,
            args.output_root / "MS1_points",
            mz_round=args.mz_round,
            intensity_min=args.intensity_min,
            assume_time_unit=args.assume_time_unit,
            rt_window=rt_window
        )

    if args.class_tag == "ALL":
        class_tags_to_run = list(CLASS_CONFIGS.keys())
    else:
        class_tags_to_run = [args.class_tag]

    for class_tag in class_tags_to_run:
        run_class_pipeline(
            files=files,
            class_tag=class_tag,
            metadata_path=args.metadata,
            CyanoMetDBLibrary=args.CyanoMetDBLibrary,
            ms1_points_file=ms1_points_file,
            output_root=args.output_root,
            tol=args.tol,
            polarity=args.polarity,
            rt_window=rt_window,
            ref_mz=args.ref_mz,
            ref_tol=args.ref_tol,
            ref_rt_window=ref_rt_window,
            do_blank_filter=args.blank_filter,
            do_batch_correct=args.batch_correct,
            tol_da=args.tol_da,
        )


if __name__ == "__main__":
    main()
