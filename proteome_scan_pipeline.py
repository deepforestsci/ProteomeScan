#!/usr/bin/env python3
"""
Simple general ProteomeScan pipeline
Reads inputs from separate input folder:
- SDF files
- JSON with ligand and PDB information
"""

import os
import sys
import json
import pandas as pd
import subprocess
import shutil
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def download_pdb_file(pdb_id, target_dir):
    """Download PDB file from RCSB"""
    pdb_file = target_dir / f"{pdb_id}.pdb"

    if pdb_file.exists():
        logger.info(f"PDB file {pdb_id}.pdb already exists")
        return True

    logger.info(f"Downloading PDB {pdb_id}...")

    # Try wget first, then curl
    for cmd_template in [
        "wget -O {output} https://files.rcsb.org/download/{pdb_id}.pdb",
        "curl -o {output} https://files.rcsb.org/download/{pdb_id}.pdb"
    ]:
        try:
            cmd = cmd_template.format(output=pdb_file, pdb_id=pdb_id)
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

            if result.returncode == 0 and pdb_file.exists():
                logger.info(f"Downloaded {pdb_id}.pdb")
                return True
        except Exception as e:
            continue

    logger.error(f"Failed to download {pdb_id}.pdb")
    return False

def vina_docking(raw_pdb_address, raw_ligand_address, scan_dir, exhaustiveness=32, num_modes=8):
    """Vina docking function - copied from working code"""
    import datetime as dt
    from deepchem.dock.pose_generation import VinaPoseGenerator
    from rdkit import Chem

    dir_name = "vina_docking_"+str(dt.datetime.now().isoformat())
    tmp = os.path.join(scan_dir, dir_name)
    os.mkdir(tmp)

    pg = VinaPoseGenerator()
    complex, scores = pg.generate_poses(
                        molecular_complex=(raw_pdb_address, raw_ligand_address),
                        exhaustiveness=exhaustiveness,
                        num_modes=num_modes,
                        out_dir=tmp,
                        generate_scores=True)
    return complex, scores

def run_docking_for_gene(gene_name, ligand_name, ligand_address, scan_dir):
    """Run docking for a gene-ligand pair"""
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

def run_pose_analysis(complex_file, results_dir, gene_name, pdb_id, ligand_name):
    """Run pose analysis using fpocket"""
    try:
        logger.info(f"Running pose analysis for {gene_name}_{pdb_id}_{ligand_name}")

        if not os.path.exists(complex_file):
            logger.error(f"Complex file not found: {complex_file}")
            return None

        # Import pose analysis functions
        sys.path.append('proteome_scan/pose_binding_analysis')
        from analyse_pose_script import main as analyze_pose

        # Run analysis
        analyze_pose(complex_file, results_dir, is_clean_up=True)

        # Check results
        complex_name = f"complex_{gene_name}_{pdb_id}_{ligand_name}"
        result_file = os.path.join(results_dir, complex_name, "result.csv")

        if os.path.exists(result_file):
            result_df = pd.read_csv(result_file)

            if len(result_df) > 0:
                total_coverage = result_df['% Ligand inside pocket'].sum()
                top1_coverage = result_df['% Ligand inside pocket'].max()

                logger.info(f"  {gene_name}: {total_coverage:.1f}% total, {top1_coverage:.1f}% top1")

                return {
                    'gene_name': gene_name,
                    'pdb_id': pdb_id,
                    'ligand_name': ligand_name,
                    'total_coverage': min(100.0, total_coverage),
                    'top1_coverage': min(100.0, top1_coverage),
                    'num_pockets': len(result_df),
                    'status': 'SUCCESS'
                }

        logger.warning(f"No valid results for {gene_name}")
        return None

    except Exception as e:
        logger.error(f"Pose analysis failed for {gene_name}: {e}")
        return None

def setup_gene_from_config(gene_config, work_dir):
    """Setup gene directory from config"""
    gene_name = gene_config['gene_name']
    pdb_ids = gene_config['pdb_ids']

    logger.info(f"Setting up {gene_name} with PDBs: {pdb_ids}")

    # Create directories
    gene_dir = Path(work_dir) / gene_name
    gene_dir.mkdir(parents=True, exist_ok=True)

    complexes_dir = gene_dir / "complexes"
    complexes_dir.mkdir(exist_ok=True)

    # Download PDBs
    successful_pdbs = []
    for pdb_id in pdb_ids:
        if download_pdb_file(pdb_id, gene_dir):
            successful_pdbs.append(pdb_id)

    if not successful_pdbs:
        logger.error(f"No PDBs downloaded for {gene_name}")
        return False

    # Create CSV with PDB info - keep it simple like other files
    pdb_data = []
    for pdb_id in successful_pdbs:
        pdb_data.append({
            'id': pdb_id,
            'path': f"{gene_name}/{pdb_id}.pdb"
        })

    pdb_df = pd.DataFrame(pdb_data)
    pdb_df.set_index('id', inplace=True)
    pdb_df.to_csv(gene_dir / f"{gene_name}_pdbs.csv")

    logger.info(f"Setup complete for {gene_name}: {len(successful_pdbs)} PDBs")
    return True

def run_pipeline(config_file, sdf_dir):
    """Run the complete pipeline from config file and SDF directory"""

    logger.info("=== SIMPLE PROTEOME SCAN PIPELINE ===")

    # Load config
    with open(config_file, 'r') as f:
        config = json.load(f)

    work_dir = config.get('work_dir', 'pipeline_results')
    genes = config['genes']
    ligands = config['ligands']

    # Create work directory
    Path(work_dir).mkdir(exist_ok=True)

    # Step 1: Setup genes
    logger.info("Step 1: Setting up genes and downloading PDBs")
    valid_genes = []

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

    all_results = []

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

                        result = run_pose_analysis(str(complex_file), str(results_dir), gene_name, pdb_id, ligand_name)
                        if result:
                            all_results.append(result)

    # Save final results
    if all_results:
        results_df = pd.DataFrame(all_results)
        results_file = Path(work_dir) / "pipeline_results.csv"
        results_df.to_csv(results_file, index=False)

        # Filter complexes with >50% pocket coverage
        high_coverage_df = results_df[results_df['total_coverage'] > 50.0].copy()

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

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Simple ProteomeScan Pipeline')
    parser.add_argument('--config', required=True, help='JSON config file with genes and ligands info')
    parser.add_argument('--sdf-dir', required=True, help='Directory containing SDF files')

    args = parser.parse_args()

    if not os.path.exists(args.config):
        logger.error(f"Config file not found: {args.config}")
        sys.exit(1)

    if not os.path.exists(args.sdf_dir):
        logger.error(f"SDF directory not found: {args.sdf_dir}")
        sys.exit(1)

    run_pipeline(args.config, args.sdf_dir)