"""Simple promiscuity filtering for ProteomeScan pipeline."""

import json
import logging
from typing import Set

logger = logging.getLogger(__name__)


def load_non_promiscuous_targets(promis_file: str, threshold: str = "25%_22") -> Set[str]:
    """Load non-promiscuous target genes from JSON file.

    Parameters
    ----------
    promis_file : str
        Path to JSON file containing promiscuity threshold data.
    threshold : str, optional
        Key name for threshold to use from JSON data, default "25%_22".

    Returns
    -------
    set of str
        Set of gene names classified as non-promiscuous. Returns empty
        set if loading fails.

    Notes
    -----
    JSON file should contain threshold keys mapping to lists of gene names.
    Returns the gene list for the specified threshold key.

    Examples
    --------
    >>> targets = load_non_promiscuous_targets('promis_thresholds_may19.json')
    >>> print(f"Found {len(targets)} non-promiscuous targets")
    """
    try:
        with open(promis_file, 'r') as f:
            promis_data = json.load(f)

        non_promiscuous = set(promis_data.get(threshold, []))
        logger.info(f"Loaded {len(non_promiscuous)} non-promiscuous targets using {threshold} threshold")
        return non_promiscuous

    except Exception as e:
        logger.error(f"Failed to load promiscuity data: {e}")
        return set()


def is_gene_promiscuous(gene_name: str, non_promiscuous_targets: Set[str]) -> bool:
    """Check if gene is promiscuous (not in non-promiscuous set).

    Parameters
    ----------
    gene_name : str
        Gene name to check.
    non_promiscuous_targets : set of str
        Set of non-promiscuous gene names.

    Returns
    -------
    bool
        True if gene is promiscuous (not in non-promiscuous set),
        False if gene is non-promiscuous (in the set).

    Notes
    -----
    Conservative approach: genes not explicitly listed as non-promiscuous
    are considered promiscuous.

    Examples
    --------
    >>> non_promiscuous = load_non_promiscuous_targets('promis_thresholds_may19.json')
    >>> is_promiscuous = is_gene_promiscuous('BRAF', non_promiscuous)
    >>> print(f"BRAF is promiscuous: {is_promiscuous}")
    """
    return gene_name not in non_promiscuous_targets