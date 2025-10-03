"""Simple promiscuity filtering for ProteomeScan pipeline."""

import json
import logging
from typing import Set

logger = logging.getLogger(__name__)


def load_non_promiscuous_targets(promis_file: str, threshold: str = "25%_22") -> Set[str]:
    """Load non-promiscuous targets from JSON file.

    Args:
        promis_file: Path to promis_thresholds_may19.json
        threshold: Threshold key (default: "25%_22")

    Returns:
        Set of non-promiscuous gene names
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
    """Check if a gene is promiscuous (not in non-promiscuous list).

    Args:
        gene_name: Gene name to check
        non_promiscuous_targets: Set of non-promiscuous genes

    Returns:
        True if gene is promiscuous, False otherwise
    """
    return gene_name not in non_promiscuous_targets