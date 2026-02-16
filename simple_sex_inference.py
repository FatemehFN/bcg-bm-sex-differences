"""Simple donor sex inference from scRNA-seq expression using established markers.

This script infers donor sex by aggregating expression of sex-linked genes across
all cells for each donor, then applying simple classification rules.

Strategy:
1. Load cell metadata to map cells to donors
2. For each capture, load 10X data and aggregate by donor
3. Use Y-chromosome genes (male markers) and XIST (female marker) for classification
4. Apply conservative thresholds to classify as M/F/undetermined

Dependencies: pip install pandas numpy scipy
Run from repo root: python simple_sex_inference.py
"""

import gzip
import pandas as pd
import numpy as np
from pathlib import Path
from scipy.io import mmread
from typing import Dict, List, Tuple


def load_cell_metadata() -> pd.DataFrame:
    """Load per-cell metadata linking barcodes to donors."""
    with gzip.open('scRNA_cell_metadata.tsv.gz', 'rt') as f:
        meta = pd.read_csv(f, sep='\t')
    return meta[~meta['donor'].isna()].copy()


def load_10x_single_capture(capture: str) -> Tuple[np.ndarray, List[str], List[str]]:
    """Load 10X data for one capture."""
    cap_dir = Path('GSE248728_downloads') / capture
    
    # Load matrix (genes x cells)
    mat_path = cap_dir / f'GSE248728_{capture}_matrix.mtx.gz'
    mat = mmread(str(mat_path)).tocsr()
    
    # Load barcodes
    bar_path = cap_dir / f'GSE248728_{capture}_barcodes.tsv.gz'
    with gzip.open(bar_path, 'rt') as f:
        barcodes = [line.strip() for line in f]
    
    # Load gene symbols
    feat_path = cap_dir / f'GSE248728_{capture}_features.tsv.gz'
    gene_symbols = []
    with gzip.open(feat_path, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            gene_symbols.append(parts[1])  # gene symbol is 2nd column
    
    return mat, barcodes, gene_symbols


def aggregate_donor_expression(meta: pd.DataFrame) -> Dict[str, pd.Series]:
    """Build donor-level expression profiles by summing across all their cells."""
    
    # Sex marker genes (confirmed available from previous check)
    y_genes = ['SRY', 'RPS4Y1', 'ZFY', 'USP9Y', 'DDX3Y', 'UTY', 'KDM5D', 'EIF1AY']
    sex_markers = ['XIST'] + y_genes
    
    donor_profiles = {}
    
    for capture in meta['capture'].unique():
        print(f"Processing {capture}...")
        
        cap_meta = meta[meta['capture'] == capture]
        if cap_meta.empty:
            continue
            
        try:
            mat, barcodes, gene_symbols = load_10x_single_capture(capture)
        except Exception as e:
            print(f"  Error loading {capture}: {e}")
            continue
        
        # Create gene symbol -> row index map
        gene_to_idx = {gene: i for i, gene in enumerate(gene_symbols)}
        
        # Get indices for sex marker genes that exist
        marker_indices = []
        marker_names = []
        for gene in sex_markers:
            if gene in gene_to_idx:
                marker_indices.append(gene_to_idx[gene])
                marker_names.append(gene)
        
        if not marker_indices:
            print(f"  No sex markers found in {capture}")
            continue
        
        # Extract sex marker expression matrix
        sex_mat = mat[marker_indices, :]  # shape: (n_sex_markers, n_cells)
        
        # Map 10X barcodes to metadata barcodes (strip "-1" suffix)
        bc_map = {}
        for i, bc in enumerate(barcodes):
            base_bc = bc[:-2] if bc.endswith('-1') else bc
            bc_map[base_bc] = i
        
        # Aggregate by donor
        for donor in cap_meta['donor'].unique():
            donor_cells = cap_meta[cap_meta['donor'] == donor]
            
            # Find column indices for this donor's cells
            col_indices = []
            for bc in donor_cells['barcode']:
                if bc in bc_map:
                    col_indices.append(bc_map[bc])
            
            if not col_indices:
                continue
            
            # Sum expression across this donor's cells
            donor_counts = np.array(sex_mat[:, col_indices].sum(axis=1)).flatten()
            donor_series = pd.Series(donor_counts, index=marker_names)
            
            # Add to accumulated profile
            if donor in donor_profiles:
                donor_profiles[donor] += donor_series
            else:
                donor_profiles[donor] = donor_series.copy()
    
    return donor_profiles


def classify_sex_simple(donor_profiles: Dict[str, pd.Series]) -> pd.DataFrame:
    """Classify donor sex using simple marker-based rules."""
    
    results = []
    
    for donor, profile in donor_profiles.items():
        
        # Get XIST expression
        xist_expr = profile.get('XIST', 0)
        
        # Get Y-chromosome gene expression (sum of all Y genes)
        y_genes = ['SRY', 'RPS4Y1', 'ZFY', 'USP9Y', 'DDX3Y', 'UTY', 'KDM5D', 'EIF1AY']
        y_expr = sum(profile.get(gene, 0) for gene in y_genes if gene in profile.index)
        
        # Log-transform for better thresholding
        xist_log = np.log1p(xist_expr)
        y_log = np.log1p(y_expr)
        
        # Classification rules based on observed data patterns
        # High Y genes + low XIST = Male
        # High XIST + low Y genes = Female
        if y_log > 9.0 and xist_log < 6.0:
            sex_pred = 'M'
            confidence = 'high' if y_log > 10.0 and xist_log < 5.0 else 'medium'
        elif xist_log > 9.0 and y_log < 7.0:
            sex_pred = 'F'
            confidence = 'high' if xist_log > 9.5 and y_log < 6.5 else 'medium'
        else:
            sex_pred = 'undetermined'
            confidence = 'low'
        
        results.append({
            'donor': donor,
            'xist_raw': xist_expr,
            'y_genes_raw': y_expr,
            'xist_log1p': xist_log,
            'y_genes_log1p': y_log,
            'sex_inferred': sex_pred,
            'confidence': confidence
        })
    
    return pd.DataFrame(results).sort_values('donor')


def main():
    """Main execution function."""
    print("Loading cell metadata...")
    meta = load_cell_metadata()
    print(f"Found {len(meta)} cells from {meta['donor'].nunique()} donors")
    
    print("Aggregating donor expression profiles...")
    donor_profiles = aggregate_donor_expression(meta)
    print(f"Built profiles for {len(donor_profiles)} donors")
    
    print("Classifying sex...")
    results = classify_sex_simple(donor_profiles)
    
    # Save full results
    results.to_csv('donor_sex_inference_detailed.csv', index=False)
    print("Detailed results saved to: donor_sex_inference_detailed.csv")
    
    # Save simple donor/sex mapping
    simple = results[['donor', 'sex_inferred']].rename(columns={'sex_inferred': 'sex'})
    simple.to_csv('donor_sex.csv', index=False)
    print("Simple mapping saved to: donor_sex.csv")
    
    # Print summary
    print("\nSex inference summary:")
    print(results.groupby(['sex_inferred', 'confidence']).size())
    
    print("\nDetailed results:")
    pd.set_option('display.max_columns', None)
    print(results)


if __name__ == '__main__':
    main()