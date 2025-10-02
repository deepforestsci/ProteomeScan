#!/usr/bin/env python3
"""
Simple general ProteomeScan pipeline - Refactored to use existing utilities

This module provides a complete pipeline for proteome scanning that:
- Reads inputs from separate input folder (SDF files, JSON with ligand and PDB information)
- Sets up gene directories and downloads PDB files
- Performs molecular docking using Vina
- Analyzes binding poses using fpocket
- Filters results based on pocket coverage

"""

import os
import sys
import json
import pandas as pd
import shutil
import argparse
from pathlib import Path
import logging
from typing import Dict, List, Any, Optional, Tuple, Union

# Import utilities from the proteome_scan package
from proteome_scan import (
    download_pdb_file,
    setup_gene_from_config,
    vina_docking,
    run_pose_analysis
)

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger: logging.Logger = logging.getLogger(__name__)

# download_pdb_file is now imported from proteome_scan.gene_pdb_utils

# vina_docking is now imported from proteome_scan.gene_guided_docking_utils

def run_docking_for_gene(
    gene_name: str,
    ligand_name: str,
    ligand_address: str,
    scan_dir: str
) -> bool:
    """
    Run molecular docking for a specific gene-ligand pair.

    This function performs Vina docking between a ligand and all PDB structures
    associated with a given gene. It processes results and saves the top scoring
    complexes for further analysis.

    Args:
        gene_name (str): Name of the target gene (e.g., 'BRAF', 'MEK1')
        ligand_name (str): Name of the ligand compound (e.g., 'Dabrafenib')
        ligand_address (str): File path to the ligand SDF file
        scan_dir (str): Root directory for pipeline results

    Returns:
        bool: True if docking completed successfully, False otherwise

    Raises:
        Exception: If docking process encounters critical errors

    Note:
        - Requires a PDB metadata CSV file at {scan_dir}/{gene_name}/{gene_name}_pdbs.csv
        - Creates complex PDB files in {scan_dir}/{gene_name}/complexes/
        - Saves results to {scan_dir}/{ligand_name}/top_score_{gene_name}_{ligand_name}.csv
    """
    from rdkit import Chem

    try:
        logger.info(f"Docking {gene_name} with {ligand_name}")

        # Check if results already exist
        if os.path.exists(f"{scan_dir}/{ligand_name}/top_score_{gene_name}_{ligand_name}.csv"):
            logger.info(f"Skipping {gene_name} {ligand_name} - results already exist")
            return True

        # Create directories
        os.makedirs(f"{scan_dir}/{gene_name}/complexes", exist_ok=True)

        # Check for PDbs CSV
        if not os.path.exists(f"{scan_dir}/{gene_name}/{gene_name}_pdbs.csv"):
            logger.error(f"PDbs file not found for {gene_name}")
            return False

        # Load PDB info
        df = pd.read_csv(f"{scan_dir}/{gene_name}/{gene_name}_pdbs.csv", index_col='id')
        docking_scores = {}

        # Run docking for each PDB
        for pdb_id in list(df.index):
            pdb_address = os.path.join(scan_dir, df.loc[pdb_id]['path'])

            try:
                assert os.path.exists(pdb_address), f"PDB file not found for {pdb_id}"
                assert os.path.exists(ligand_address), f"Ligand file not found"

                logger.info(f"  Docking PDB {pdb_id}...")
                complex, score = vina_docking(pdb_address, ligand_address, scan_dir)

            except Exception as e:
                logger.warning(f"  Docking failed for {pdb_id}: {e}")
                score = [None]
                complex = None

            # Save complex if successful
            if complex is not None:
                complex_mol = Chem.CombineMols(complex[0][0], complex[0][1])
                complex_path = f"{scan_dir}/{gene_name}/complexes/complex_{gene_name}_{pdb_id}_{ligand_name}.pdb"
                Chem.rdmolfiles.MolToPDBFile(complex_mol, complex_path)
                logger.info(f"  Saved complex: complex_{gene_name}_{pdb_id}_{ligand_name}.pdb")

            docking_scores[pdb_id] = score

        # Process results
        df_docked = pd.DataFrame({'id': docking_scores.keys(), 'scores': docking_scores.values()})
        df_docked.dropna(inplace=True)

        if len(df_docked) == 0:
            logger.error(f"No successful docking results for {gene_name}")
            return False

        df_docked['top_score'] = df_docked['scores'].apply(lambda x: x[0])
        df_docked = df_docked.set_index('id')

        merged_df = pd.merge(df, df_docked, how='left', left_index=True, right_index=True)
        top_df = merged_df[['path', 'top_score', 'scores']].sort_values(by='top_score')

        # Take top results
        if len(top_df) > 2:
            top2_df = top_df.iloc[[0,1]]
        else:
            top2_df = top_df

        top2_df['gene_name'] = [gene_name]*len(top2_df)

        # Save results
        os.makedirs(f"{scan_dir}/{ligand_name}", exist_ok=True)
        top2_df.to_csv(f"{scan_dir}/{ligand_name}/top_score_{gene_name}_{ligand_name}.csv")

        logger.info(f"Docking completed for {gene_name}: {len(top2_df)} results")
        return True

    except Exception as e:
        logger.error(f"Error in docking {gene_name} {ligand_name}: {e}")
        return False

# run_pose_analysis is now imported from proteome_scan.pose_binding_analysis

# setup_gene_from_config is now imported from proteome_scan.gene_pdb_utils

def run_pipeline(config_file: str, sdf_dir: str) -> None:
    """
    Execute the complete ProteomeScan pipeline from configuration file.

    This is the main pipeline orchestrator that coordinates all pipeline steps:
    1. Gene setup and PDB download
    2. Molecular docking execution
    3. Pose analysis and pocket coverage calculation
    4. Results filtering and file organization

    Args:
        config_file (str): Path to JSON configuration file containing:
            - work_dir: Output directory path
            - genes: List of gene configurations with names and PDB IDs
            - ligands: List of ligand configurations with names and SDF files
        sdf_dir (str): Directory containing ligand SDF files

    Returns:
        None

    Raises:
        FileNotFoundError: If config file or SDF directory doesn't exist
        json.JSONDecodeError: If config file is malformed

    Example:
        >>> run_pipeline('config.json', 'data/ligands/processed/')

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
    logger.info("=== SIMPLE PROTEOME SCAN PIPELINE ===")

    # Load configuration
    with open(config_file, 'r') as f:
        config: Dict[str, Any] = json.load(f)

    work_dir: str = config.get('work_dir', 'pipeline_results')
    genes: List[Dict[str, Any]] = config['genes']
    ligands: List[Dict[str, Any]] = config['ligands']

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
        ligand_name = ligand_config['ligand_name']
        sdf_filename = ligand_config['sdf_file']
        ligand_path = os.path.join(sdf_dir, sdf_filename)

        if not os.path.exists(ligand_path):
            logger.error(f"Ligand SDF not found: {ligand_path}")
            continue

        logger.info(f"Processing {ligand_name}")

        for gene_config in valid_genes:
            gene_name = gene_config['gene_name']
            run_docking_for_gene(gene_name, ligand_name, ligand_path, work_dir)

    # Step 3: Run pose analysis
    logger.info("Step 3: Running pose analysis")

    all_results: List[Dict[str, Any]] = []

    for ligand_config in ligands:
        ligand_name = ligand_config['ligand_name']

        for gene_config in valid_genes:
            gene_name = gene_config['gene_name']
            complexes_dir = Path(work_dir) / gene_name / "complexes"

            if complexes_dir.exists():
                for complex_file in complexes_dir.glob(f"complex_{gene_name}_*_{ligand_name}.pdb"):
                    # Extract PDB ID from filename
                    parts = complex_file.stem.split('_')
                    if len(parts) >= 3:
                        pdb_id = parts[2]

                        results_dir = Path(work_dir) / "pose_analysis"
                        results_dir.mkdir(exist_ok=True)

                        result: Optional[Dict[str, Any]] = run_pose_analysis(
                            str(complex_file), str(results_dir), gene_name, pdb_id, ligand_name
                        )
                        if result:
                            all_results.append(result)

    # Save final results
    if all_results:
        results_df: pd.DataFrame = pd.DataFrame(all_results)
        results_file = Path(work_dir) / "pipeline_results.csv"
        results_df.to_csv(results_file, index=False)

        # Filter complexes with >50% pocket coverage
        high_coverage_df: pd.DataFrame = results_df[results_df['total_coverage'] > 50.0].copy()

        if len(high_coverage_df) > 0:
            # Save high coverage complexes
            high_coverage_dir = Path(work_dir) / "results" / "high_coverage_complexes"
            high_coverage_dir.mkdir(parents=True, exist_ok=True)

            high_coverage_file = high_coverage_dir / "high_coverage_complexes.csv"
            high_coverage_df.to_csv(high_coverage_file, index=False)

            # Copy actual complex PDB files to high coverage directory
            for _, row in high_coverage_df.iterrows():
                gene_name = row['gene_name']
                pdb_id = row['pdb_id']
                ligand_name = row['ligand_name']

                # Find original complex file
                complex_file = Path(work_dir) / gene_name / "complexes" / f"complex_{gene_name}_{pdb_id}_{ligand_name}.pdb"
                if complex_file.exists():
                    # Copy to high coverage directory
                    target_file = high_coverage_dir / f"complex_{gene_name}_{pdb_id}_{ligand_name}_coverage_{row['total_coverage']:.1f}pct.pdb"
                    shutil.copy2(complex_file, target_file)
                    logger.info(f"Saved high coverage complex: {target_file.name}")

            logger.info(f"=== HIGH COVERAGE COMPLEXES (>50%) ===")
            logger.info(f"Found {len(high_coverage_df)} complexes with >50% pocket coverage:")
            for _, row in high_coverage_df.iterrows():
                print(f"{row['gene_name']} {row['pdb_id']} {row['ligand_name']}: {row['total_coverage']:.1f}% coverage")

            logger.info(f"High coverage results saved to: {high_coverage_file}")
            logger.info(f"Complex PDB files copied to: {high_coverage_dir}")
        else:
            logger.warning("No complexes found with >50% pocket coverage")

        logger.info("=== ALL RESULTS ===")
        for _, row in results_df.iterrows():
            print(f"{row['gene_name']} {row['pdb_id']} {row['ligand_name']}: {row['total_coverage']:.1f}% total, {row['top1_coverage']:.1f}% top1")

        logger.info(f"All results saved to: {results_file}")
    else:
        logger.warning("No final results")

    logger.info("=== PIPELINE COMPLETED ===")


def main() -> None:
    """
    Command-line interface for the ProteomeScan pipeline.

    Parses command-line arguments and executes the pipeline with the specified
    configuration file and SDF directory.

    Command-line Arguments:
        --config: Path to JSON configuration file (required)
        --sdf-dir: Directory containing ligand SDF files (required)

    Example:
        python proteome_scan_pipeline.py --config config.json --sdf-dir data/ligands/processed/

    Raises:
        SystemExit: If required files/directories are not found or pipeline fails
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
        logger.error(f"Config file not found: {args.config}")
        sys.exit(1)

    if not os.path.exists(args.sdf_dir):
        logger.error(f"SDF directory not found: {args.sdf_dir}")
        sys.exit(1)

    # Execute pipeline
    try:
        run_pipeline(args.config, args.sdf_dir)
        logger.info("Pipeline completed successfully!")
    except Exception as e:
        logger.error(f"Pipeline failed with error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
