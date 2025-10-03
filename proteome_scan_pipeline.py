#!/usr/bin/env python3
"""ProteomeScan pipeline for molecular docking and binding analysis.

This module provides a complete pipeline for proteome scanning that:
    - Reads inputs from separate input folder (SDF files, JSON with ligand and PDB information)
    - Sets up gene directories and downloads PDB files
    - Performs molecular docking using AutoDock Vina
    - Analyzes binding poses using fpocket
    - Filters results based on pocket coverage

The pipeline is designed to be modular and uses existing utilities from the proteome_scan package.
"""

import os
import sys
import json
import pandas as pd
import shutil
import argparse
from pathlib import Path
import logging
from typing import Dict, List, Any, Optional, Tuple

from proteome_scan import (
    download_pdb_file,
    setup_gene_from_config,
    vina_docking,
    run_pose_analysis
)
from proteome_scan.gene_guided_docking_utils.promiscuity_filter import (
    load_non_promiscuous_targets,
    is_gene_promiscuous
)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger: logging.Logger = logging.getLogger(__name__)

def run_docking_for_gene(
    gene_name: str,
    ligand_name: str,
    ligand_address: str,
    scan_dir: str
) -> bool:
    """Run molecular docking for a specific gene-ligand pair.

    This function performs AutoDock Vina docking between a ligand and all PDB structures
    associated with a given gene. It processes results and saves the top scoring
    complexes for further analysis.

    Args:
        gene_name: Name of the target gene (e.g., 'BRAF', 'MEK1').
        ligand_name: Name of the ligand compound (e.g., 'Dabrafenib').
        ligand_address: File path to the ligand SDF file.
        scan_dir: Root directory for pipeline results.

    Returns:
        True if docking completed successfully, False otherwise.

    Raises:
        Exception: If docking process encounters critical errors.

    Note:
        Requires a PDB metadata CSV file at {scan_dir}/{gene_name}/{gene_name}_pdbs.csv.
        Creates complex PDB files in {scan_dir}/{gene_name}/complexes/.
        Saves results to {scan_dir}/{ligand_name}/top_score_{gene_name}_{ligand_name}.csv.
    """
    from rdkit import Chem

    try:
        logger.info(f"Docking {gene_name} with {ligand_name}")

        # Check if results already exist
        results_file: str = f"{scan_dir}/{ligand_name}/top_score_{gene_name}_{ligand_name}.csv"
        if os.path.exists(results_file):
            logger.info(f"Skipping {gene_name} {ligand_name} - results already exist")
            return True

        # Create directories
        os.makedirs(f"{scan_dir}/{gene_name}/complexes", exist_ok=True)

        # Check for PDBs CSV
        pdbs_file_path: str = f"{scan_dir}/{gene_name}/{gene_name}_pdbs.csv"
        if not os.path.exists(pdbs_file_path):
            logger.error(f"PDBs metadata file not found: {pdbs_file_path}")
            return False

        # Load PDB metadata
        try:
            df: pd.DataFrame = pd.read_csv(pdbs_file_path, index_col='id')
        except (pd.errors.EmptyDataError, pd.errors.ParserError) as e:
            logger.error(f"Failed to parse PDBs metadata file {pdbs_file_path}: {e}")
            return False
        docking_scores: Dict[str, List[Optional[float]]] = {}

        # Run docking for each PDB
        for pdb_id in list(df.index):
            pdb_address: str = os.path.join(scan_dir, df.loc[pdb_id]['path'])

            try:
                if not os.path.exists(pdb_address):
                    raise FileNotFoundError(f"PDB file not found: {pdb_address}")
                if not os.path.exists(ligand_address):
                    raise FileNotFoundError(f"Ligand file not found: {ligand_address}")

                logger.info(f"  Docking PDB {pdb_id}...")
                docking_result = vina_docking(pdb_address, ligand_address, scan_dir)
                complex: Optional[Tuple[Any, ...]] = docking_result[0]
                score: List[Optional[float]] = docking_result[1]

            except (FileNotFoundError, ValueError, RuntimeError) as e:
                logger.warning(f"  Docking failed for {pdb_id}: {e}")
                score: List[Optional[float]] = [None]
                complex: Optional[Tuple[Any, ...]] = None
            except Exception as e:
                logger.error(f"  Unexpected error during docking for {pdb_id}: {e}")
                score: List[Optional[float]] = [None]
                complex: Optional[Tuple[Any, ...]] = None

            # Save complex if successful
            if complex is not None:
                combined_mol = Chem.CombineMols(complex[0][0], complex[0][1])
                complex_path: str = f"{scan_dir}/{gene_name}/complexes/complex_{gene_name}_{pdb_id}_{ligand_name}.pdb"
                Chem.rdmolfiles.MolToPDBFile(combined_mol, complex_path)
                logger.info(f"  Saved complex: complex_{gene_name}_{pdb_id}_{ligand_name}.pdb")

            docking_scores[pdb_id] = score

        # Process docking results
        docking_results_df: pd.DataFrame = pd.DataFrame({
            'id': list(docking_scores.keys()),
            'scores': list(docking_scores.values())
        })
        docking_results_df.dropna(inplace=True)

        if len(docking_results_df) == 0:
            logger.error(f"No successful docking results for {gene_name}")
            return False

        docking_results_df['top_score'] = docking_results_df['scores'].apply(lambda x: x[0])
        docking_results_df.set_index('id', inplace=True)

        merged_results_df: pd.DataFrame = pd.merge(
            df, docking_results_df, how='left', left_index=True, right_index=True
        )
        sorted_results_df: pd.DataFrame = merged_results_df[
            ['path', 'top_score', 'scores']
        ].sort_values(by='top_score')

        # Take top 2 results or all if fewer than 2
        top_results_df: pd.DataFrame = (
            sorted_results_df.iloc[[0, 1]] if len(sorted_results_df) > 2
            else sorted_results_df
        )
        top_results_df['gene_name'] = gene_name

        # Save results
        output_dir: str = f"{scan_dir}/{ligand_name}"
        os.makedirs(output_dir, exist_ok=True)

        output_file: str = f"{output_dir}/top_score_{gene_name}_{ligand_name}.csv"
        top_results_df.to_csv(output_file)

        logger.info(f"Docking completed for {gene_name}: {len(top_results_df)} results")
        return True

    except Exception as e:
        logger.error(f"Error in docking {gene_name} {ligand_name}: {e}")
        return False


def run_pipeline(config_file: str, sdf_dir: str) -> None:
    """Execute the complete ProteomeScan pipeline from configuration file.

    This is the main pipeline orchestrator that coordinates all pipeline steps:
        1. Gene setup and PDB download
        2. Molecular docking execution
        3. Pose analysis and pocket coverage calculation
        4. Results filtering and file organization

    Args:
        config_file: Path to JSON configuration file containing work_dir, genes, and ligands.
        sdf_dir: Directory containing ligand SDF files.

    Raises:
        FileNotFoundError: If config file or SDF directory doesn't exist.
        json.JSONDecodeError: If config file is malformed.

    Example:
        >>> run_pipeline('config.json', 'data/ligands/processed/')

    Note:
        Config file format:
        {
            "work_dir": "pipeline_results",
            "genes": [
                {"gene_name": "BRAF", "pdb_ids": ["4XV2", "5VAL"]},
                {"gene_name": "MEK1", "pdb_ids": ["3EQH"]}
            ],
            "ligands": [
                {"ligand_name": "Dabrafenib", "sdf_file": "Processed_Dabrafenib.sdf"}
            ]
        }
    """
    logger.info("=== PROTEOME SCAN PIPELINE ===")

    # Load configuration
    try:
        with open(config_file, 'r', encoding='utf-8') as f:
            config: Dict[str, Any] = json.load(f)
    except FileNotFoundError:
        logger.error(f"Configuration file not found: {config_file}")
        raise
    except json.JSONDecodeError as e:
        logger.error(f"Invalid JSON in configuration file {config_file}: {e}")
        raise

    work_dir: str = config.get('work_dir', 'pipeline_results')

    try:
        genes: List[Dict[str, Any]] = config['genes']
        ligands: List[Dict[str, Any]] = config['ligands']
    except KeyError as e:
        logger.error(f"Missing required key in configuration file: {e}")
        raise ValueError(f"Configuration file must contain '{e.args[0]}' key") from e

    # Create work directory
    Path(work_dir).mkdir(exist_ok=True)

    # Step 1: Setup genes
    logger.info("Step 1: Setting up genes and downloading PDBs")
    valid_genes: List[Dict[str, Any]] = []

    for gene_config in genes:
        if setup_gene_from_config(gene_config, work_dir):
            valid_genes.append(gene_config)
            logger.info(f"✓ {gene_config['gene_name']}: Setup successful")
        else:
            logger.error(f"✗ {gene_config['gene_name']}: Setup failed")

    if not valid_genes:
        logger.error("No valid genes. Stopping.")
        return

    # Step 2: Run docking
    logger.info("Step 2: Running docking")

    for ligand_config in ligands:
        ligand_name: str = ligand_config['ligand_name']
        sdf_filename: str = ligand_config['sdf_file']
        ligand_path: str = os.path.join(sdf_dir, sdf_filename)

        if not os.path.exists(ligand_path):
            logger.error(f"Ligand SDF file not found: {ligand_path}")
            continue

        logger.info(f"Processing {ligand_name}")

        for gene_config in valid_genes:
            gene_name: str = gene_config['gene_name']
            run_docking_for_gene(gene_name, ligand_name, ligand_path, work_dir)

    # Step 2.5: Load promiscuity filter
    logger.info("Step 2.5: Loading promiscuity filter")
    promis_file: str = os.path.join(work_dir, "../proteome_scan/promis_thresholds_may19.json")
    if not os.path.exists(promis_file):
        # Try alternative path
        promis_file = "proteome_scan/promis_thresholds_may19.json"

    non_promiscuous_targets: set = load_non_promiscuous_targets(promis_file, threshold="25%_22")
    if not non_promiscuous_targets:
        logger.warning("Could not load promiscuity data - proceeding without filtering")
        non_promiscuous_targets = set()

    # Step 3: Run pose analysis (with promiscuity filtering)
    logger.info("Step 3: Running pose analysis with promiscuity filtering")

    all_results: List[Dict[str, Any]] = []

    for ligand_config in ligands:
        ligand_name: str = ligand_config['ligand_name']

        for gene_config in valid_genes:
            gene_name: str = gene_config['gene_name']

            # Check promiscuity before pose analysis
            if non_promiscuous_targets and is_gene_promiscuous(gene_name, non_promiscuous_targets):
                logger.info(f"Skipping pose analysis for promiscuous gene: {gene_name}")
                continue

            complexes_dir: Path = Path(work_dir) / gene_name / "complexes"

            if complexes_dir.exists():
                for complex_file in complexes_dir.glob(f"complex_{gene_name}_*_{ligand_name}.pdb"):
                    # Extract PDB ID from filename
                    parts: List[str] = complex_file.stem.split('_')
                    if len(parts) >= 3:
                        pdb_id: str = parts[2]

                        results_dir: Path = Path(work_dir) / "pose_analysis"
                        results_dir.mkdir(exist_ok=True)

                        logger.info(f"Running pose analysis for non-promiscuous gene: {gene_name}")
                        result: Optional[Dict[str, Any]] = run_pose_analysis(
                            str(complex_file), str(results_dir), gene_name, pdb_id, ligand_name
                        )
                        if result:
                            all_results.append(result)

    # Save final results
    if all_results:
        results_df: pd.DataFrame = pd.DataFrame(all_results)
        results_file: Path = Path(work_dir) / "pipeline_results.csv"
        results_df.to_csv(results_file, index=False)

        # Filter complexes with >50% pocket coverage
        high_coverage_df: pd.DataFrame = results_df[results_df['total_coverage'] > 50.0].copy()

        if len(high_coverage_df) > 0:
            # Save high coverage complexes
            high_coverage_dir: Path = Path(work_dir) / "results" / "high_coverage_complexes"
            high_coverage_dir.mkdir(parents=True, exist_ok=True)

            high_coverage_file: Path = high_coverage_dir / "high_coverage_complexes.csv"
            high_coverage_df.to_csv(high_coverage_file, index=False)

            # Copy actual complex PDB files to high coverage directory
            for _, row in high_coverage_df.iterrows():
                gene_name: str = row['gene_name']
                pdb_id: str = row['pdb_id']
                ligand_name: str = row['ligand_name']

                # Find original complex file
                complex_file: Path = Path(work_dir) / gene_name / "complexes" / f"complex_{gene_name}_{pdb_id}_{ligand_name}.pdb"
                if complex_file.exists():
                    # Copy to high coverage directory
                    target_file: Path = high_coverage_dir / f"complex_{gene_name}_{pdb_id}_{ligand_name}_coverage_{row['total_coverage']:.1f}pct.pdb"
                    shutil.copy2(complex_file, target_file)
                    logger.info(f"Saved high coverage complex: {target_file.name}")

            logger.info("=== HIGH COVERAGE COMPLEXES (>50%) ===")
            logger.info(f"Found {len(high_coverage_df)} complexes with >50% pocket coverage:")
            for _, row in high_coverage_df.iterrows():
                logger.info(f"{row['gene_name']} {row['pdb_id']} {row['ligand_name']}: {row['total_coverage']:.1f}% coverage")

            logger.info(f"High coverage results saved to: {high_coverage_file}")
            logger.info(f"Complex PDB files copied to: {high_coverage_dir}")
        else:
            logger.warning("No complexes found with >50% pocket coverage")

        logger.info("=== ALL RESULTS (NON-PROMISCUOUS ONLY) ===")
        for _, row in results_df.iterrows():
            logger.info(f"{row['gene_name']} {row['pdb_id']} {row['ligand_name']}: {row['total_coverage']:.1f}% total, {row['top1_coverage']:.1f}% top1")

        # Show promiscuity filtering summary
        if non_promiscuous_targets:
            total_genes = len(valid_genes)
            analyzed_genes = len(set(results_df['gene_name']))
            skipped_genes = total_genes - analyzed_genes
            logger.info(f"=== PROMISCUITY FILTERING SUMMARY ===")
            logger.info(f"Total genes: {total_genes}")
            logger.info(f"Non-promiscuous genes analyzed: {analyzed_genes}")
            logger.info(f"Promiscuous genes skipped: {skipped_genes}")

        logger.info(f"All results saved to: {results_file}")
    else:
        logger.warning("No final results")

    logger.info("=== PIPELINE COMPLETED ===")


def main() -> None:
    """Command-line interface for the ProteomeScan pipeline.

    Parses command-line arguments and executes the pipeline with the specified
    configuration file and SDF directory.

    Command-line Arguments:
        --config: Path to JSON configuration file (required).
        --sdf-dir: Directory containing ligand SDF files (required).

    Example:
        python proteome_scan_pipeline.py --config config.json --sdf-dir data/ligands/processed/

    Raises:
        SystemExit: If required files/directories are not found or pipeline fails.
    """
    parser: argparse.ArgumentParser = argparse.ArgumentParser(
        description='ProteomeScan Pipeline - Molecular docking and pose analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    python proteome_scan_pipeline.py --config config.json --sdf-dir data/ligands/processed/

Config file should contain gene and ligand information in JSON format.
See documentation for detailed configuration format.
        """
    )
    parser.add_argument(
        '--config',
        required=True,
        type=str,
        help='JSON configuration file with genes and ligands information'
    )
    parser.add_argument(
        '--sdf-dir',
        required=True,
        type=str,
        help='Directory containing ligand SDF files'
    )

    args: argparse.Namespace = parser.parse_args()

    # Validate input arguments
    if not os.path.exists(args.config):
        logger.error(f"Configuration file not found: {args.config}")
        sys.exit(1)

    if not os.path.isdir(args.sdf_dir):
        logger.error(f"SDF directory not found or not a directory: {args.sdf_dir}")
        sys.exit(1)

    # Execute pipeline
    try:
        run_pipeline(args.config, args.sdf_dir)
        logger.info("Pipeline completed successfully!")
    except (FileNotFoundError, ValueError, json.JSONDecodeError) as e:
        logger.error(f"Pipeline configuration error: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Pipeline failed with unexpected error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
