"""Infer donor-level sex from scRNA-seq expression.

This script uses the GEO-supplied 10X matrices for GSE248728 together with
the per-cell metadata in `scRNA_cell_metadata.tsv.gz` to infer each donor's
sex from gene expression.

Approach (expression-based, not using any external clinical metadata):

1. For each capture (capture1–capture14), load the 10X matrix
   (features, barcodes, matrix) from GSE248728_downloads/.
2. Restrict barcodes to demuxlet *singlet* cells present in
   `scRNA_cell_metadata.tsv.gz` and assigned to a donor.
3. For each donor, build a pseudobulk profile by summing counts across all
   their cells (across captures).
4. Compute, per donor:
     - Mean expression of key Y-chromosome genes (e.g. RPS4Y1, DDX3Y,
       EIF1AY, KDM5D, UTY).
     - Mean expression of XIST.
5. Classify donors using simple thresholds:

     - If Y-score is clearly > 0 and XIST is low → "M".
     - If Y-score is ~0 and XIST is high       → "F".
     - Otherwise                               → "undetermined".

The exact thresholds are heuristic but conservative; the result is written
to `donor_sex_inferred.csv`, which can then be reviewed and, if desired,
copied/edited into `donor_sex.csv` for use by `make_scRNA_metadata_summary.py`.

Run from the repository root, e.g.:

    cd /Users/fsfatemi/sun_et_al
    python infer_donor_sex_from_scRNA.py

Dependencies:
    pip install pandas scipy numpy

This script does *not* use any external clinical metadata; it is purely
expression-based, as requested.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy import io as spio


ROOT = Path(__file__).resolve().parent


def load_cell_metadata() -> pd.DataFrame:
    """Load per-cell metadata produced by make_scRNA_metadata_summary.

    Expects scRNA_cell_metadata.tsv.gz in the repo root.
    """

    meta_path = ROOT / "scRNA_cell_metadata.tsv.gz"
    if not meta_path.exists():
        raise FileNotFoundError(f"Cannot find {meta_path}; run make_scRNA_metadata_summary.py first.")

    with gzip.open(meta_path, "rt") as f:
        meta = pd.read_csv(f, sep="\t")

    # Keep only cells with an assigned donor (demuxlet singlets are already enforced)
    meta = meta[~meta["donor"].isna()].copy()
    return meta


def load_10x_capture(capture: str) -> Tuple[np.ndarray, List[str], List[str]]:
    """Load 10X matrix for a given capture from GSE248728_downloads.

    Returns (matrix, barcodes, gene_symbols).
    """

    cap_dir = ROOT / "GSE248728_downloads" / capture
    mat_path = cap_dir / f"GSE248728_{capture}_matrix.mtx.gz"
    bar_path = cap_dir / f"GSE248728_{capture}_barcodes.tsv.gz"
    feat_path = cap_dir / f"GSE248728_{capture}_features.tsv.gz"

    if not (mat_path.exists() and bar_path.exists() and feat_path.exists()):
        raise FileNotFoundError(f"Missing 10X files for {capture} in {cap_dir}")

    # MatrixMarket sparse matrix (genes × cells)
    mat = spio.mmread(str(mat_path)).tocsr()

    # Barcodes
    with gzip.open(bar_path, "rt") as f:
        barcodes = [line.strip() for line in f]

    # Features: gene_id, gene_name, feature_type
    gene_symbols: List[str] = []
    with gzip.open(feat_path, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            gene_symbols.append(parts[1])

    if mat.shape[0] != len(gene_symbols) or mat.shape[1] != len(barcodes):
        raise ValueError(
            f"Dimension mismatch for {capture}: matrix {mat.shape}, "
            f"{len(gene_symbols)} genes, {len(barcodes)} barcodes."
        )

    return mat, barcodes, gene_symbols


def build_donor_pseudobulk(
    cell_meta: pd.DataFrame,
    y_genes: List[str],
    xist_genes: List[str],
) -> Tuple[pd.DataFrame, Dict[str, float], Dict[str, float]]:
    """Construct donor-level pseudobulk and return expression scores.

    Returns:
        pseudobulk_counts: DataFrame (donor × gene)
        donor_y_score: mean expression of Y genes per donor
        donor_xist_score: mean expression of XIST per donor
    """

    # Store gene symbols from the first capture; gene order must match the
    # matrix rows. Duplicated gene symbols are allowed and will be handled
    # later by aggregating columns with the same name.
    gene_symbols_global: List[str] | None = None
    donors = sorted(cell_meta["donor"].unique())

    # We'll accumulate pseudobulk in a dict donor → vector of counts
    donor_counts: Dict[str, np.ndarray] = {}
    n_genes_global: int | None = None

    for capture in sorted(cell_meta["capture"].unique()):
        subset = cell_meta[cell_meta["capture"] == capture]
        if subset.empty:
            continue

        print(f"Processing {capture} with {len(subset)} cells...")
        mat, barcodes, gene_symbols = load_10x_capture(capture)

        # Build barcode index for this capture.
        # 10X barcodes are of the form "AAAC...-1", while the demuxlet-derived
        # metadata barcodes are the same strings *without* the "-1" suffix.
        # Mirror the R preprocessing (where "-1" is appended) by indexing
        # using the stripped version so metadata barcodes match.
        barcode_to_col = {}
        for i, bc in enumerate(barcodes):
            base = bc[:-2] if bc.endswith("-1") else bc
            # In case of any duplicates (should not happen), keep first
            if base not in barcode_to_col:
                barcode_to_col[base] = i

        # Initialize global gene dimension / symbols on first capture
        if n_genes_global is None:
            n_genes_global = mat.shape[0]
            gene_symbols_global = list(gene_symbols)
        else:
            if n_genes_global != mat.shape[0]:
                raise ValueError("Gene dimension differs across captures; this script assumes they match.")

        # For each donor in this capture, sum counts across their cells
        for donor in subset["donor"].unique():
            donor_cells = subset[subset["donor"] == donor]
            cols = []
            for bc in donor_cells["barcode"]:
                idx = barcode_to_col.get(bc)
                if idx is not None:
                    cols.append(idx)
            if not cols:
                continue

            # Sum selected columns (axis=1 => per gene)
            sub_mat = mat[:, cols]
            donor_vec = np.asarray(sub_mat.sum(axis=1)).ravel()

            if donor not in donor_counts:
                donor_counts[donor] = donor_vec
            else:
                donor_counts[donor] += donor_vec

    if not donor_counts:
        raise RuntimeError("No donor counts accumulated; check metadata and barcodes.")

    # Build pseudobulk DataFrame
    assert n_genes_global is not None
    if gene_symbols_global is None:
        raise RuntimeError("No gene symbols recorded; this should not happen.")

    donor_vecs = [donor_counts[d] for d in donors]
    lengths = {v.shape[0] for v in donor_vecs}
    if len(lengths) != 1:
        raise ValueError(f"Inconsistent gene vector lengths across donors: {lengths}")
    n_genes = lengths.pop()
    if n_genes != len(gene_symbols_global):
        raise ValueError(
            f"Numeric gene dimension ({n_genes}) does not match number of gene symbols ({len(gene_symbols_global)})."
        )

    pseudobulk_raw = pd.DataFrame(
        np.vstack(donor_vecs),
        index=donors,
        columns=gene_symbols_global,
    )

    # Aggregate any duplicated gene symbols by summing across columns with
    # the same name so that downstream indexing by gene symbol is well-defined.
    pseudobulk = pseudobulk_raw.groupby(pseudobulk_raw.columns, axis=1).sum()

    # Compute expression scores
    donor_y_score: Dict[str, float] = {}
    donor_xist_score: Dict[str, float] = {}

    # Intersect requested genes with available genes
    available_genes = set(pseudobulk.columns)
    y_used = [g for g in y_genes if g in available_genes]
    xist_used = [g for g in xist_genes if g in available_genes]

    if not y_used:
        print("Warning: none of the specified Y-chromosome genes found in features; Y score will be 0.")
    if not xist_used:
        print("Warning: XIST not found in features; XIST score will be 0.")

    for donor in donors:
        row = pseudobulk.loc[donor]
        y_val = float(row[y_used].mean()) if y_used else 0.0
        xist_val = float(row[xist_used].mean()) if xist_used else 0.0
        donor_y_score[donor] = y_val
        donor_xist_score[donor] = xist_val

    return pseudobulk, donor_y_score, donor_xist_score


def classify_sex(donor_y: Dict[str, float], donor_xist: Dict[str, float]) -> pd.DataFrame:
    """Classify donors as M/F/undetermined based on Y and XIST scores.

    Uses simple heuristics on log1p-normalized scores.
    """

    donors = sorted(donor_y.keys())

    records = []
    for d in donors:
        y_raw = donor_y[d]
        x_raw = donor_xist[d]

        # log1p to reduce scale effects
        y = np.log1p(y_raw)
        x = np.log1p(x_raw)

        # Heuristic thresholds; these can be tuned after inspection
        if y > 1.0 and x < 0.5:
            inferred = "M"
        elif y < 0.2 and x > 1.0:
            inferred = "F"
        else:
            inferred = "undetermined"

        records.append(
            {
                "donor": d,
                "y_score_raw": y_raw,
                "xist_score_raw": x_raw,
                "y_score_log1p": y,
                "xist_score_log1p": x,
                "sex_inferred": inferred,
            }
        )

    return pd.DataFrame.from_records(records)


def main() -> None:
    cell_meta = load_cell_metadata()

    # Define marker genes: robust chrY genes and XIST
    y_genes = ["RPS4Y1", "DDX3Y", "EIF1AY", "KDM5D", "UTY", "ZFY"]
    xist_genes = ["XIST"]

    pseudobulk, donor_y, donor_xist = build_donor_pseudobulk(
        cell_meta, y_genes=y_genes, xist_genes=xist_genes
    )

    inferred = classify_sex(donor_y, donor_xist)

    out_path = ROOT / "donor_sex_inferred.csv"
    print(f"Writing inferred donor sex to {out_path}")
    inferred.to_csv(out_path, index=False)

    # Also print a small summary to stdout
    print("\nInferred donor sex:")
    print(inferred[["donor", "sex_inferred"]].to_string(index=False))


if __name__ == "__main__":
    main()
