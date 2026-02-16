"""Build donor/timepoint/vaccination cell-count summary for Sun et al. scRNA-seq.

This script uses the authors' demuxlet outputs and sample metadata to
recreate the mapping from individual cells to donors (S1–S21), timepoints
(Td0/Tm3), and vaccination status (BCG/CTL), then aggregates cell counts.

It purposely mirrors the logic in the R script
`scRNA_analyses/1_data_processing_clustering/1_create_seurat_objects`:

- demuxlet BEST calls are used to identify singlets vs doublets/ambiguous
- timepoint is inferred from capture ID with a special case for capture1
- vaccination is inferred from a hard-coded list of CTL donors
- donor IDs (S1–S21) come from `all_samples_meta_data.csv`

If a file `donor_sex.csv` is available with columns:

    donor,sex

it will be joined to the summary so that each donor carries a sex label.
If that file is not present, a template will be written that you can fill in
manually based on Table S1 / clinical metadata from the paper, and the
summary will contain `sex = NA` for all donors.

Outputs (written in the current working directory by default):

- scRNA_cell_metadata.tsv.gz
    One row per demuxlet *singlet* cell with columns:
      capture, barcode, best, timepoint, vaccination, unique_id,
      donor, sex (if available)

- scRNA_donor_timepoint_summary.csv
    Aggregated table with one row per donor × timepoint × vaccination:
      donor, sex, timepoint, vaccination, n_cells, n_captures,
      captures (semicolon-separated list)

Run from the repository root, e.g.:

    cd /Users/fsfatemi/sun_et_al
    python make_scRNA_metadata_summary.py

Dependencies:
    pip install pandas
"""

from __future__ import annotations

import os
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent


def load_all_samples_metadata() -> pd.DataFrame:
    """Load all_samples_meta_data.csv (sample → timepoint/donor/vaccination).

    File location (relative to repo root):
        BCG-humanBM-code/scRNA_analyses/1_data_processing_clustering/all_samples_meta_data.csv
    """

    meta_path = (
        ROOT
        / "BCG-humanBM-code"
        / "scRNA_analyses"
        / "1_data_processing_clustering"
        / "all_samples_meta_data.csv"
    )
    if not meta_path.exists():
        raise FileNotFoundError(f"Cannot find metadata file: {meta_path}")

    meta = pd.read_csv(meta_path)
    expected_cols = {"Sample", "timepoint", "donor", "vaccination"}
    missing = expected_cols.difference(meta.columns)
    if missing:
        raise ValueError(f"Metadata file {meta_path} is missing columns: {missing}")

    return meta


def iter_demuxlet_files():
    """Yield (capture, path) for each demuxlet CSV file.

    Files live in:
        BCG-humanBM-code/scRNA_analyses/1_data_processing_clustering/demuxlet_csv_outs/
    and are named like `capture1_demuxlet_out.csv`.
    """

    demux_dir = (
        ROOT
        / "BCG-humanBM-code"
        / "scRNA_analyses"
        / "1_data_processing_clustering"
        / "demuxlet_csv_outs"
    )
    if not demux_dir.exists():
        raise FileNotFoundError(f"Cannot find demuxlet output directory: {demux_dir}")

    for path in sorted(demux_dir.glob("capture*_demuxlet_out.csv")):
        name = path.name
        # e.g. capture1_demuxlet_out.csv → capture1
        capture = name.split("_")[0]
        yield capture, path


def classify_demuxlet_call(best: str) -> str:
    """Classify demuxlet BEST string as singlet/doublet/ambiguous.

    Mirrors grepl("SNG"/"DBL"/"AMB") logic in the R script.
    """

    if pd.isna(best):
        return "unknown"
    if best.startswith("SNG-"):
        return "singlet"
    if best.startswith("DBL-"):
        return "doublet"
    if best.startswith("AMB-"):
        return "ambiguous"
    # Fallback: unknown pattern
    return "unknown"


def assign_timepoint(capture: str, best: str) -> str:
    """Assign Td0/Tm3 following the R logic in 1_create_seurat_objects.

    In R:

        if batchID == "capture1":
            if BEST == "SNG-LB-SS-1S-RS-S8-CD34neg_S7_R1_001": timepoint <- "Td0"
            else:                                                  timepoint <- "Tm3"
        else if batchID in c("capture3", "capture5", ...):      timepoint <- "Td0"
             else:                                                 timepoint <- "Tm3"
    """

    if capture == "capture1":
        if best == "SNG-LB-SS-1S-RS-S8-CD34neg_S7_R1_001":
            return "Td0"
        return "Tm3"

    td0_captures = {
        "capture3",
        "capture5",
        "capture7",
        "capture9",
        "capture11",
        "capture13",
    }
    return "Td0" if capture in td0_captures else "Tm3"


def assign_vaccination(best: str) -> str:
    """Assign CTL/BCG using the CTL list from the R script.

    In R:

        CTL <- c(
          "SNG-LB-SS-1S-RS-S13-CD34neg_S12_R1_001",
          "SNG-LB-SS-1S-RS-S7-CD34neg_S6_R1_001",
          "SNG-LB-SS-1S-RS-S11-CD34neg_S10_R1_001",
          "SNG-LB-SS-1S-RS-S3-CD34neg_S3_R1_001",
          "SNG-LB-SS-1S-RS-S18-CD34neg_S17_R1_001")

        vaccination <- ifelse(BEST %in% CTL, "CTL", "BCG")
    """

    ctl_samples = {
        "SNG-LB-SS-1S-RS-S13-CD34neg_S12_R1_001",
        "SNG-LB-SS-1S-RS-S7-CD34neg_S6_R1_001",
        "SNG-LB-SS-1S-RS-S11-CD34neg_S10_R1_001",
        "SNG-LB-SS-1S-RS-S3-CD34neg_S3_R1_001",
        "SNG-LB-SS-1S-RS-S18-CD34neg_S17_R1_001",
    }
    return "CTL" if best in ctl_samples else "BCG"


def load_donor_sex(meta: pd.DataFrame) -> pd.DataFrame:
    """Load donor→sex mapping if provided; else create a template.

    The template lists all donors present in `meta['donor']` and writes
    donor_sex_template.csv in the current directory with `sex` set to NA.
    """

    donors = sorted(meta["donor"].unique())
    sex_path = ROOT / "donor_sex.csv"

    if sex_path.exists():
        donor_sex = pd.read_csv(sex_path)
        if "donor" not in donor_sex.columns or "sex" not in donor_sex.columns:
            raise ValueError(
                f"donor_sex.csv at {sex_path} must have columns 'donor' and 'sex'"
            )
        return donor_sex

    # No donor_sex.csv yet – create a template and return NA for all donors
    template_path = ROOT / "donor_sex_template.csv"
    if not template_path.exists():
        tmpl = pd.DataFrame({"donor": donors, "sex": ["NA"] * len(donors)})
        tmpl.to_csv(template_path, index=False)
        print(
            f"No donor_sex.csv found. Wrote template with donors to {template_path}. "
            "Fill this file with 'M'/'F' (or similar) and save as donor_sex.csv "
            "to include sex in downstream summaries."
        )

    return pd.DataFrame({"donor": donors, "sex": ["NA"] * len(donors)})


def build_cell_metadata() -> pd.DataFrame:
    """Construct per-cell metadata dataframe for all demuxlet singlets.

    Returns columns:
        capture, barcode, best, call_type, timepoint, vaccination,
        unique_id, donor, sex (if available)
    """

    meta = load_all_samples_metadata()
    donor_sex = load_donor_sex(meta)

    records = []

    for capture, path in iter_demuxlet_files():
        print(f"Reading demuxlet file for {capture}: {path}")
        df = pd.read_csv(path)
        if "BARCODE" not in df.columns or "BEST" not in df.columns:
            raise ValueError(f"Demuxlet file {path} missing BARCODE/BEST columns")

        df = df[["BARCODE", "BEST"]].copy()
        df["capture"] = capture
        df["call_type"] = df["BEST"].apply(classify_demuxlet_call)

        # Keep only singlets, mirroring the R code where only demuxlet_status
        # == "singlet" is retained for downstream analysis.
        df = df[df["call_type"] == "singlet"].copy()
        if df.empty:
            continue

        df["timepoint"] = df.apply(
            lambda row: assign_timepoint(row["capture"], row["BEST"]), axis=1
        )
        df["vaccination"] = df["BEST"].apply(assign_vaccination)
        df["unique_id"] = (
            df["BEST"] + "_" + df["timepoint"] + "_" + df["vaccination"]
        )

        records.append(df)

    if not records:
        raise RuntimeError("No singlet cells found in demuxlet outputs.")

    cells = pd.concat(records, ignore_index=True)

    # Join to sample-level metadata to get donor IDs
    cells = cells.merge(
        meta,
        how="left",
        left_on="unique_id",
        right_on="Sample",
    )

    # Drop cells without a matching sample row (e.g. missing S9 Td0 sample)
    before = len(cells)
    cells = cells[~cells["donor"].isna()].copy()
    after = len(cells)
    if after < before:
        print(f"Dropped {before - after} cells without donor metadata.")

    # Attach sex information by donor
    cells = cells.merge(donor_sex, how="left", on="donor")

    # Select and rename columns for a cleaner output
    out = cells[[
        "capture",
        "BARCODE",
        "BEST",
        "call_type",
        "timepoint_x",
        "vaccination_x",
        "unique_id",
        "donor",
        "sex",
    ]].copy()

    out.rename(
        columns={
            "BARCODE": "barcode",
            "BEST": "best",
            "timepoint_x": "timepoint",
            "vaccination_x": "vaccination",
        },
        inplace=True,
    )

    return out


def write_outputs(cell_meta: pd.DataFrame) -> None:
    """Write per-cell metadata and donor×timepoint summaries to disk."""

    cell_meta_path = ROOT / "scRNA_cell_metadata.tsv.gz"
    summary_path = ROOT / "scRNA_donor_timepoint_summary.csv"

    print(f"Writing per-cell metadata to {cell_meta_path}")
    cell_meta.to_csv(cell_meta_path, sep="\t", index=False, compression="gzip")

    # Summarize by donor × timepoint × vaccination
    grp_cols = ["donor", "sex", "timepoint", "vaccination"]
    grouped = cell_meta.groupby(grp_cols, dropna=False)

    summary = grouped.agg(
        n_cells=("barcode", "size"),
        n_captures=("capture", pd.Series.nunique),
        captures=("capture", lambda x: ";".join(sorted(set(x)))),
    ).reset_index()

    print(f"Writing donor/timepoint summary to {summary_path}")
    summary.to_csv(summary_path, index=False)


def main() -> None:
    cell_meta = build_cell_metadata()
    print(f"Constructed metadata for {len(cell_meta)} demuxlet singlet cells.")
    write_outputs(cell_meta)


if __name__ == "__main__":
    main()
