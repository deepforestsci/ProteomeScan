#!/usr/bin/env python3
"""ProteomeScan pipeline for molecular docking and binding analysis.

This module provides a complete pipeline for proteome scanning that performs
molecular docking between ligands and protein targets, followed by binding
pose analysis and filtering based on pocket coverage.

The pipeline consists of four main stages:
1. Gene setup and PDB file download from the Protein Data Bank
2. Molecular docking using AutoDock Vina
3. Binding pose analysis using fpocket for cavity detection
4. Result filtering based on pocket coverage and promiscuity

The pipeline is designed to be modular and uses existing utilities from the
proteome_scan package for maximum flexibility and reusability.

Examples
--------
Run the pipeline from command line:

    $ python proteome_scan_pipeline.py --config config.json --sdf-dir data/ligands/

Run programmatically:

    >>> from proteome_scan_pipeline import run_pipeline
    >>> run_pipeline('config.json', 'data/ligands/processed/')

Notes
-----
The pipeline requires:
- A JSON configuration file specifying genes and ligands
- SDF files for ligands in the specified directory
- Network connectivity for PDB file downloads
- AutoDock Vina and fpocket must be installed and accessible

See Also
--------
proteome_scan.vina_docking : Core docking functionality
proteome_scan.analyse_pose_main : Pose analysis with fpocket using multi_pose approach
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
    vina_docking
)
from proteome_scan.pose_binding_analysis.multi_pose_run import get_overall_ligand_interactions, get_total_top_n_bucket_percentages
from proteome_scan.pose_binding_analysis.analyse_pose_script import main as analyse_pose_main
from proteome_scan.gene_guided_docking_utils.promiscuity_filter import (
    load_promiscuous_targets,
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

    This function performs AutoDock Vina docking between a ligand and all PDB
    structures associated with a given gene. It processes the docking results,
    saves the top scoring complexes as PDB files, and generates a CSV summary
    of the docking scores.

    Parameters
    ----------
    gene_name : str
        Name of the target gene (e.g., 'BRAF', 'MEK1'). Must correspond to
        a directory with PDB metadata.
    ligand_name : str
        Name of the ligand compound (e.g., 'Dabrafenib'). Used for naming
        output files and directories.
    ligand_address : str
        Full file path to the ligand SDF file containing the 3D structure
        of the small molecule.
    scan_dir : str
        Root directory for pipeline results. All gene subdirectories and
        output files will be created here.

    Returns
    -------
    bool
        True if docking completed successfully for at least one PDB structure,
        False if all docking attempts failed or critical errors occurred.

    Raises
    ------
    Exception
        If docking process encounters critical errors that prevent execution.
        Individual PDB docking failures are logged but do not raise exceptions.

    Notes
    -----
    Requires a PDB metadata CSV file at {scan_dir}/{gene_name}/{gene_name}_pdbs.csv,
    valid SDF file at ligand_address, and AutoDock Vina installation.

    Creates complex PDB files in {scan_dir}/{gene_name}/complexes/ and results CSV
    at {scan_dir}/{ligand_name}/top_score_{gene_name}_{ligand_name}.csv.

    The results CSV contains the top 2 docking scores (or all if fewer than 2 PDBs)
    sorted by binding affinity, with more negative scores indicating stronger binding.

    Examples
    --------
    >>> success = run_docking_for_gene(
    ...     'BRAF', 'Dabrafenib', 'ligands/Dabrafenib.sdf', 'results/'
    ... )
    >>> if success:
    ...     print("Docking completed successfully")
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
            pdb_address: str = os.path.join(scan_dir, str(df.loc[pdb_id, 'path']))

            # Initialize result variables
            complex_result: Optional[Tuple[Any, ...]] = None
            score_result: List[Optional[float]] = [None]

            try:
                if not os.path.exists(pdb_address):
                    raise FileNotFoundError(f"PDB file not found: {pdb_address}")
                if not os.path.exists(ligand_address):
                    raise FileNotFoundError(f"Ligand file not found: {ligand_address}")

                logger.info(f"  Docking PDB {pdb_id}...")
                docking_result = vina_docking(pdb_address, ligand_address, scan_dir)
                complex_result = docking_result[0]
                score_result = docking_result[1]

            except (FileNotFoundError, ValueError, RuntimeError) as e:
                logger.warning(f"  Docking failed for {pdb_id}: {e}")
                score_result = [None]
                complex_result = None
            except Exception as e:
                logger.error(f"  Unexpected error during docking for {pdb_id}: {e}")
                score_result = [None]
                complex_result = None

            # Save complex if successful
            if complex_result is not None:
                combined_mol = Chem.CombineMols(complex_result[0][0], complex_result[0][1])
                complex_path: str = f"{scan_dir}/{gene_name}/complexes/complex_{gene_name}_{pdb_id}_{ligand_name}.pdb"
                Chem.rdmolfiles.MolToPDBFile(combined_mol, complex_path)
                logger.info(f"  Saved complex: complex_{gene_name}_{pdb_id}_{ligand_name}.pdb")

            docking_scores[pdb_id] = score_result

        # Process docking results
        docking_results_df: pd.DataFrame = pd.DataFrame({
            'id': list(docking_scores.keys()),
            'scores': list(docking_scores.values())
        })
        docking_results_df.dropna(inplace=True)

        if len(docking_results_df) == 0:
            logger.error(f"No successful docking results for {gene_name}")
            return False

        docking_results_df['top_score'] = docking_results_df['scores'].apply(lambda x: x[0] if x and len(x) > 0 else None)
        docking_results_df.set_index('id', inplace=True)

        merged_results_df: pd.DataFrame = pd.merge(
            df, docking_results_df, how='left', left_index=True, right_index=True
        )
        selected_columns_df = merged_results_df[['path', 'top_score', 'scores']]
        sorted_results_df: pd.DataFrame = selected_columns_df.sort_values('top_score')  # type: ignore[call-overload]

        # Take top 2 results or all if fewer than 2
        if len(sorted_results_df) >= 2:
            top_results_df = sorted_results_df.iloc[[0, 1]].copy()
        else:
            top_results_df = sorted_results_df.copy()

        top_results_df.loc[:, 'gene_name'] = gene_name

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

    This is the main pipeline orchestrator that coordinates all pipeline steps
    in sequence. The pipeline performs comprehensive molecular docking analysis
    with automated filtering for promiscuous targets and pocket coverage analysis.

    The pipeline executes the following stages:
    1. Gene setup and PDB file download from the Protein Data Bank
    2. Molecular docking execution using AutoDock Vina
    3. Promiscuity filtering to exclude non-selective targets
    4. Pose analysis and pocket coverage calculation using fpocket
    5. Results filtering and high-coverage complex identification
    6. File organization and summary report generation

    Parameters
    ----------
    config_file : str
        Path to JSON configuration file containing work_dir, genes, and ligands
        specifications. Must be a valid JSON file with required keys.
    sdf_dir : str
        Directory containing ligand SDF files. All SDF files referenced in the
        configuration must exist in this directory.

    Raises
    ------
    FileNotFoundError
        If config file or SDF directory doesn't exist, or if required SDF
        files are missing from the specified directory.
    json.JSONDecodeError
        If config file contains malformed JSON that cannot be parsed.
    ValueError
        If config file is missing required keys ('genes' or 'ligands').

    Notes
    -----
    The pipeline generates individual gene directories with PDB metadata and complex
    files, pose_analysis/ directory with fpocket cavity analysis results,
    results/high_coverage_complexes/ with complexes having >50% pocket coverage,
    and pipeline_results.csv with comprehensive analysis results.

    Only non-promiscuous genes (based on internal promiscuity thresholds) are
    subjected to pose analysis to focus on selective drug targets.

    Examples
    --------
    Basic pipeline execution:

    >>> run_pipeline('config.json', 'data/ligands/processed/')

    The configuration file should follow this format:

    >>> config = {
    ...     "work_dir": "pipeline_results",
    ...     "genes": [
    ...         {"gene_name": "BRAF", "pdb_ids": ["4XV2", "5VAL"]},
    ...         {"gene_name": "MEK1", "pdb_ids": ["3EQH"]}
    ...     ],
    ...     "ligands": [
    ...         {"ligand_name": "Dabrafenib", "sdf_file": "Processed_Dabrafenib.sdf"}
    ...     ]
    ... }

    See Also
    --------
    run_docking_for_gene : Individual gene-ligand docking execution
    proteome_scan.analyse_pose_main : Pocket coverage analysis using multi_pose approach
    proteome_scan.setup_gene_from_config : Gene setup and PDB download
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
    # Use the correct path to the promis_thresholds_may19.json file
    import proteome_scan
    proteome_scan_dir = os.path.dirname(proteome_scan.__file__)
    promis_file: str = os.path.join(proteome_scan_dir, "promis_thresholds_may19.json")

    if not os.path.exists(promis_file):
        logger.error(f"Promiscuity file not found at: {promis_file}")
        logger.warning("Proceeding without promiscuity filtering - all genes will be processed")
        promiscuous_targets = set[str]()
    else:
        logger.info(f"Loading promiscuity data from: {promis_file}")
        promiscuous_targets: set[str] = load_promiscuous_targets(promis_file, threshold="25%_22")
        if not promiscuous_targets:
            logger.warning("Could not load promiscuity data - proceeding without filtering")
            promiscuous_targets = set[str]()

    # Step 3: Run pose analysis (with promiscuity filtering)
    logger.info("Step 3: Running pose analysis with promiscuity filtering")

    all_results: List[Dict[str, Any]] = []

    for ligand_config in ligands:
        ligand_name_analysis: str = ligand_config['ligand_name']

        for gene_config in valid_genes:
            gene_name_analysis: str = gene_config['gene_name']

            # Check promiscuity before pose analysis
            if promiscuous_targets and is_gene_promiscuous(gene_name_analysis, promiscuous_targets):
                logger.info(f"Skipping pose analysis for promiscuous gene: {gene_name_analysis}")
                continue

            complexes_dir: Path = Path(work_dir) / gene_name_analysis / "complexes"

            if complexes_dir.exists():
                for complex_file in complexes_dir.glob(f"complex_{gene_name_analysis}_*_{ligand_name_analysis}.pdb"):
                    # Extract PDB ID from filename
                    parts: List[str] = complex_file.stem.split('_')
                    if len(parts) >= 3:
                        pdb_id_analysis: str = parts[2]

                        results_dir: Path = Path(work_dir) / "pose_analysis"
                        results_dir.mkdir(exist_ok=True)

                        logger.info(f"Running pose analysis for non-promiscuous gene: {gene_name_analysis}")

                        # Use multi_pose_run approach directly
                        try:
                            # Run analysis using analyse_pose_main (same as multi_pose_run does)
                            analyse_pose_main(str(complex_file), str(results_dir), is_clean_up=True)

                            # Extract metrics using exact same logic as multi_pose_run.py
                            complex_name = f"complex_{gene_name_analysis}_{pdb_id_analysis}_{ligand_name_analysis}"
                            result_file = os.path.join(results_dir, complex_name, "full_analysis.csv")

                            if os.path.exists(result_file):
                                df = pd.read_csv(result_file)
                                if len(df) > 0:
                                    # Use multi_pose_run.py functions directly
                                    total_coverage = get_overall_ligand_interactions(df)
                                    top1_coverage = get_total_top_n_bucket_percentages(df, 1)  # Top 1 pocket by druggability

                                    result = {
                                        'gene_name': gene_name_analysis,
                                        'pdb_id': pdb_id_analysis,
                                        'ligand_name': ligand_name_analysis,
                                        'total_coverage': total_coverage,
                                        'top1_coverage': top1_coverage,
                                        'num_pockets': len(df),
                                        'status': 'SUCCESS'
                                    }
                                    all_results.append(result)
                                    logger.info(f"  {gene_name_analysis}: {total_coverage:.1f}% total, {top1_coverage:.1f}% top1")
                        except Exception as e:
                            logger.error(f"Pose analysis failed for {gene_name_analysis}: {e}")

    # Save final results
    if all_results:
        results_df: pd.DataFrame = pd.DataFrame(all_results)
        results_file: Path = Path(work_dir) / "pipeline_results.csv"
        results_df.to_csv(results_file, index=False)

        # Filter complexes with >50% pocket coverage
        high_coverage_df: pd.DataFrame = results_df[results_df['total_coverage'] > 50.0].copy()  # type: ignore[assignment]

        if len(high_coverage_df) > 0:
            # Save high coverage complexes
            high_coverage_dir: Path = Path(work_dir) / "results" / "high_coverage_complexes"
            high_coverage_dir.mkdir(parents=True, exist_ok=True)

            high_coverage_file: Path = high_coverage_dir / "high_coverage_complexes.csv"
            high_coverage_df.to_csv(high_coverage_file, index=False)

            # Copy actual complex PDB files to high coverage directory
            for _, row in high_coverage_df.iterrows():
                row_gene_name: str = str(row['gene_name'])
                row_pdb_id: str = str(row['pdb_id'])
                row_ligand_name: str = str(row['ligand_name'])

                # Find original complex file
                source_complex_file: Path = Path(work_dir) / row_gene_name / "complexes" / f"complex_{row_gene_name}_{row_pdb_id}_{row_ligand_name}.pdb"
                if source_complex_file.exists():
                    # Copy to high coverage directory
                    target_file: Path = high_coverage_dir / f"complex_{row_gene_name}_{row_pdb_id}_{row_ligand_name}_coverage_{row['total_coverage']:.1f}pct.pdb"
                    shutil.copy2(source_complex_file, target_file)
                    logger.info(f"Saved high coverage complex: {target_file.name}")

            logger.info("=== HIGH COVERAGE COMPLEXES (>50%) ===")
            logger.info(f"Found {len(high_coverage_df)} complexes with >50% pocket coverage:")
            for _, row in high_coverage_df.iterrows():  # type: ignore[misc]
                logger.info(f"{row['gene_name']} {row['pdb_id']} {row['ligand_name']}: {row['total_coverage']:.1f}% coverage")

            logger.info(f"High coverage results saved to: {high_coverage_file}")
            logger.info(f"Complex PDB files copied to: {high_coverage_dir}")
        else:
            logger.warning("No complexes found with >50% pocket coverage")

        logger.info("=== ALL RESULTS (NON-PROMISCUOUS ONLY) ===")
        for _, row in results_df.iterrows():  # type: ignore[misc]
            logger.info(f"{row['gene_name']} {row['pdb_id']} {row['ligand_name']}: {row['total_coverage']:.1f}% total, {row['top1_coverage']:.1f}% top1")

        # Show promiscuity filtering summary
        if promiscuous_targets:
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

    Provides a command-line entry point for executing the complete ProteomeScan
    pipeline. Parses command-line arguments, validates input files and directories,
    and executes the pipeline with comprehensive error handling and logging.

    This function serves as the main entry point when the script is executed
    directly from the command line. It handles argument parsing, input validation,
    and provides informative error messages for common failure scenarios.

    Parameters
    ----------
    None
        Command-line arguments are parsed from sys.argv using argparse.

    Command-line Arguments
    ----------------------
    --config : str, required
        Path to JSON configuration file containing gene and ligand specifications.
        Must be a valid, readable JSON file with required structure.
    --sdf-dir : str, required
        Directory containing ligand SDF files referenced in the configuration.
        Must be an existing directory with appropriate read permissions.

    Raises
    ------
    SystemExit
        Exit code 1 if required files/directories are not found, configuration
        is invalid, or pipeline execution fails with unrecoverable errors.

    Notes
    -----
    The function provides detailed error messages and logging throughout execution.
    Pipeline progress and results are logged to stdout with timestamps and
    appropriate log levels (INFO, WARNING, ERROR).

    All pipeline outputs are saved to the work directory specified in the
    configuration file, allowing for easy result inspection and analysis.

    Examples
    --------
    Execute pipeline from command line:

    .. code-block:: bash

        python proteome_scan_pipeline.py --config config.json --sdf-dir data/ligands/processed/

    View help and usage information:

    .. code-block:: bash

        python proteome_scan_pipeline.py --help

    See Also
    --------
    run_pipeline : Core pipeline execution function
    argparse.ArgumentParser : Command-line argument parsing
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
