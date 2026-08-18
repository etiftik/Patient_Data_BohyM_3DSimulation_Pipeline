"""
BoHyM 3D Gene Enrichment Analysis
=================================

This script performs three related RNA-seq analyses focused on cytoskeleton-
associated genes:

1. Compare two sample groups in a bulk tumor RNA-seq dataset.


Main steps
----------
- Load raw count data
- Convert counts to CPM
- Perform gene-wise Mann-Whitney U tests
- Correct p-values for multiple testing
- Restrict results to cytoskeleton-related Gene Ontology terms
- Annotate significant genes using MyGene.info
- Export selected results as CSV files


"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

import gseapy as gp
import mygene


# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------

DATA_DIR = Path(".")
OUTPUT_DIR = Path("results")
OUTPUT_DIR.mkdir(exist_ok=True)

BULK_TUMOR_FILE = DATA_DIR / "BulkTumor_RNAseq_Collaboration.xlsx"
GSE166044_FILE = DATA_DIR / "GSE166044_raw_counts_GRCh38.p13_NCBI.tsv"

ORGANISM = "Human"
GO_LIBRARY = "GO_Cellular_Component_2023"

# Count filtering.
# Setting both to 0 effectively keeps all genes.
MIN_COUNTS = 0
MIN_SAMPLE_FRACTION = 0

# Statistical thresholds.
ALPHA_ANALYSIS_1 = 0.05
ALPHA_ANALYSIS_2 = 0.01
FDR_ANALYSIS_3 = 0.001
LOG2FC_ANALYSIS_3 = -1


# ---------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------

def calculate_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    """Convert a raw count matrix to counts per million (CPM)."""
    total_reads = counts.sum(axis=0)

    # Avoid division by zero for empty libraries.
    total_reads = total_reads.replace(0, np.nan)

    return counts.divide(total_reads, axis=1) * 1e6


def filter_genes_by_count(
    counts: pd.DataFrame,
    min_counts: float = MIN_COUNTS,
    min_sample_fraction: float = MIN_SAMPLE_FRACTION,
) -> pd.DataFrame:
    """
    Keep genes that reach `min_counts` in at least `min_sample_fraction`
    of samples.
    """
    genes_above_threshold = (counts >= min_counts).sum(axis=1)
    required_samples = int(np.ceil(min_sample_fraction * counts.shape[1]))

    return counts.loc[genes_above_threshold >= required_samples].copy()


def mann_whitney_differential_expression(
    cpm: pd.DataFrame,
    group1_columns,
    group2_columns,
    log2fc_direction: str = "group2_over_group1",
    correction_method: str = "fdr_bh",
) -> pd.DataFrame:
    """
    Perform a gene-wise Mann-Whitney U test between two sample groups.

    Parameters
    ----------
    cpm
        CPM-normalized gene-expression matrix. Rows are genes.
    group1_columns, group2_columns
        Column names belonging to the two groups.
    log2fc_direction
        Either "group2_over_group1" or "group1_over_group2".
    correction_method
        Multiple-testing correction method accepted by statsmodels.

    Returns
    -------
    pd.DataFrame
        Original CPM matrix plus group means, log2 fold change, raw p-values,
        adjusted p-values, and a significance flag.
    """
    group1 = cpm.loc[:, group1_columns]
    group2 = cpm.loc[:, group2_columns]

    group1_mean = group1.mean(axis=1)
    group2_mean = group2.mean(axis=1)

    if log2fc_direction == "group2_over_group1":
        numerator = group2_mean
        denominator = group1_mean
    elif log2fc_direction == "group1_over_group2":
        numerator = group1_mean
        denominator = group2_mean
    else:
        raise ValueError(
            "log2fc_direction must be 'group2_over_group1' "
            "or 'group1_over_group2'."
        )

    # A small pseudocount prevents divide-by-zero and log2(0).
    pseudocount = 1e-12
    log2fc = np.log2(
        (numerator + pseudocount) / (denominator + pseudocount)
    )

    p_values = []
    for gene in cpm.index:
        _, p_value = stats.mannwhitneyu(
            group1.loc[gene].values,
            group2.loc[gene].values,
            alternative="two-sided",
        )
        p_values.append(p_value)

    rejected, adjusted_p_values, _, _ = multipletests(
        p_values,
        alpha=0.05,
        method=correction_method,
    )

    results = cpm.copy()
    results["group1_mean_cpm"] = group1_mean
    results["group2_mean_cpm"] = group2_mean
    results["log2FC"] = log2fc
    results["p_value"] = p_values
    results["adjusted_p_value"] = adjusted_p_values
    results["significant_after_correction"] = rejected

    return results


def get_cytoskeleton_gene_set() -> set[str]:
    """Build the cytoskeleton-related gene set from selected GO CC terms."""
    go_cc = gp.get_library(
        name=GO_LIBRARY,
        organism=ORGANISM,
    )

    go_terms = [
        "Actin Cytoskeleton (GO:0015629)",
        "Actin Filament (GO:0005884)",
        "Cell-Substrate Junction (GO:0030055)",
        "Cortical Actin Cytoskeleton (GO:0030864)",
        "Cortical Cytoskeleton (GO:0030863)",
        "Cytoskeleton (GO:0005856)",
        "Focal Adhesion (GO:0005925)",
        "Myosin Filament (GO:0032982)",
    ]

    genes = []
    for term in go_terms:
        genes.extend(go_cc[term])

    return set(genes)


def annotate_genes_with_mygene(gene_symbols) -> pd.DataFrame:
    """Retrieve gene names and summaries from MyGene.info."""
    mg = mygene.MyGeneInfo()

    gene_info = mg.querymany(
        list(gene_symbols),
        scopes="symbol",
        fields="name,summary",
        species="human",
    )

    records = []
    for entry in gene_info:
        symbol = entry.get("query")
        records.append(
            {
                "gene_symbol": symbol,
                "full_name": entry.get("name", "NA"),
                "function": entry.get("summary", "NA"),
            }
        )

    return pd.DataFrame(records).set_index("gene_symbol")


def keep_cytoskeleton_genes(
    results: pd.DataFrame,
    cytoskeleton_genes: set[str],
) -> pd.DataFrame:
    """Restrict results to genes in the cytoskeleton-related GO gene set."""
    subset = results.loc[results.index.isin(cytoskeleton_genes)].copy()
    return subset.sort_values("log2FC")


def annotate_and_merge(results: pd.DataFrame) -> pd.DataFrame:
    """Add MyGene.info annotations to a differential-expression table."""
    annotations = annotate_genes_with_mygene(results.index)
    return annotations.join(results, how="right")


# ---------------------------------------------------------------------
# Dataset loading
# ---------------------------------------------------------------------

def load_bulk_tumor_dataset() -> pd.DataFrame:
    """
    Load the bulk tumor RNA-seq workbook and return the gene-count matrix.

    The original notebook used the first data column as the gene identifier
    and excluded it from the sample matrix.
    """
    data = pd.read_excel(BULK_TUMOR_FILE, index_col=0)

    counts = data.iloc[:, 1:].copy()
    counts.index = data.iloc[:, 0]
    counts.index.name = "GeneID"

    return counts


# ---------------------------------------------------------------------
# Analysis 1: bulk tumor dataset
# ---------------------------------------------------------------------

def run_analysis_1(
    counts: pd.DataFrame,
    cytoskeleton_genes: set[str],
) -> pd.DataFrame:
    """Compare the two predefined groups in the bulk tumor dataset."""
    filtered_counts = filter_genes_by_count(counts)
    cpm = calculate_cpm(filtered_counts)

  
    group1_columns = list(cpm.columns[0:20])
    group2_columns = list(cpm.columns[21:])

    results = mann_whitney_differential_expression(
        cpm,
        group1_columns,
        group2_columns,
        log2fc_direction="group2_over_group1",
        correction_method="fdr_by",
    )

    significant = results.loc[
        results["p_value"] < ALPHA_ANALYSIS_1
    ].copy()

    significant.to_csv(
        OUTPUT_DIR / "analysis1_significant_genes.csv"
    )

    cytoskeleton_results = keep_cytoskeleton_genes(
        significant,
        cytoskeleton_genes,
    )

    cytoskeleton_results.to_csv(
        OUTPUT_DIR / "analysis1_cytoskeleton_genes.csv"
    )

    annotated = annotate_and_merge(cytoskeleton_results)
    annotated.to_csv(
        OUTPUT_DIR / "analysis1_bohym_input.csv"
    )

    return annotated




# ---------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------

def main():
    """Run all three analyses."""
    print("Loading cytoskeleton-related GO gene set...")
    cytoskeleton_genes = get_cytoskeleton_gene_set()
    print(f"Loaded {len(cytoskeleton_genes):,} unique cytoskeleton-related genes.")

    print("\nLoading bulk tumor dataset...")
    bulk_counts = load_bulk_tumor_dataset()
    print(f"Bulk tumor matrix: {bulk_counts.shape}")

    print("\nRunning analysis 1...")
    analysis1_results = run_analysis_1(
        bulk_counts,
        cytoskeleton_genes,
    )
    print(f"Analysis 1 retained {len(analysis1_results):,} cytoskeleton genes.")


if __name__ == "__main__":
    main()
