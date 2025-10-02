#!/usr/bin/env python3
"""
Simple end-to-end ProteomeScan pipeline.
Just download PDBs directly, no complex chain analysis.
"""

import os
import sys
import pandas as pd
import subprocess
import shutil
import json
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Hard-coded constants
PROMISCUITY_THRESHOLD = "25%_22"

def load_promiscuity_targets(json_file="proteome_scan/promis_thresholds_may19.json"):
    """Load promiscuous targets from JSON file"""
    try:
        with open(json_file, 'r') as f:
            promis_data = json.load(f)
        targets = set(promis_data.get(PROMISCUITY_THRESHOLD, []))
        logger.info(f"Loaded {len(targets)} targets from 25% promiscuity threshold")
        return targets
    except Exception as e:
        logger.warning(f"Could not load promiscuity file: {e}")
        return set()

def download_pdb_file(pdb_id, target_dir):
    """Download PDB file if not exists - exact same as complete_pipeline_4_targets.py"""
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

def setup_gene_directories(gene_name, pdb_ids, work_dir):
    """Set up directory structure for a gene with multiple PDBs"""

    logger.info(f"Setting up {gene_name} with {len(pdb_ids)} PDBs")

    # Create directories
    gene_dir = Path(work_dir) / gene_name
    gene_dir.mkdir(parents=True, exist_ok=True)

    complexes_dir = gene_dir / "complexes"
    complexes_dir.mkdir(exist_ok=True)

    # Download all PDBs for this gene
    successful_pdbs = []
    for pdb_id in pdb_ids:
        if download_pdb_file(pdb_id, gene_dir):
            successful_pdbs.append(pdb_id)

    if successful_pdbs:
        # Create a simple CSV with PDB info
        pdb_info = pd.DataFrame({
            'id': successful_pdbs,
            'path': [f"{gene_name}/{pdb_id}.pdb" for pdb_id in successful_pdbs]
        })
        pdb_info.to_csv(gene_dir / f"{gene_name}_pdbs.csv", index=False)
        logger.info(f"Setup complete for {gene_name}: {len(successful_pdbs)} PDBs")
        return True
    else:
        logger.error(f"No PDBs downloaded for {gene_name}")
        return False

def run_simple_docking(gene_name, ligand_name, ligand_sdf_path, work_dir):
    """Run simple docking - just copy PDB as complex for now"""

    logger.info(f"Running simple docking for {gene_name} with {ligand_name}")

    gene_dir = Path(work_dir) / gene_name
    complexes_dir = gene_dir / "complexes"

    if not gene_dir.exists():
        logger.error(f"Gene directory not found: {gene_dir}")
        return False

    # Check if CSV exists
    pdbs_csv = gene_dir / f"{gene_name}_pdbs.csv"
    if not pdbs_csv.exists():
        logger.error(f"PDbs CSV not found: {pdbs_csv}")
        return False

    try:
        # Load PDB list
        pdbs_df = pd.read_csv(pdbs_csv)

        # For each PDB
        success_count = 0
        for _, row in pdbs_df.iterrows():
            pdb_id = row['id']
            pdb_file = gene_dir / f"{pdb_id}.pdb"
            complex_file = complexes_dir / f"complex_{gene_name}_{pdb_id}_{ligand_name}.pdb"

            if pdb_file.exists():
                shutil.copy2(pdb_file, complex_file)
                success_count += 1
                logger.info(f"Created simple complex: {complex_file.name}")

        if success_count > 0:
            # Create a simple results CSV
            results_dir = Path(work_dir) / ligand_name
            results_dir.mkdir(exist_ok=True)

            # Create dummy scores for now
            scores_data = []
            for _, row in pdbs_df.iterrows():
                pdb_id = row['id']
                scores_data.append({
                    'gene_name': gene_name,
                    'pdb_id': pdb_id,
                    'top_score': -10.0,  # Dummy score
                    'scores': '[-10.0, -9.5, -9.0]'
                })

            scores_df = pd.DataFrame(scores_data)
            scores_df.to_csv(results_dir / f"top_score_{gene_name}_{ligand_name}.csv", index=False)

            logger.info(f"Simple docking complete for {gene_name}: {success_count} complexes")
            return True
        else:
            logger.error(f"No complexes created for {gene_name}")
            return False

    except Exception as e:
        logger.error(f"Error in simple docking for {gene_name}: {e}")
        return False

def parse_simple_results(ligands, work_dir):
    """Parse simple docking results"""
    try:
        logger.info("Parsing simple docking results")

        results = {}
        scan_results_dir = Path(work_dir) / "scan_results"
        scan_results_dir.mkdir(exist_ok=True)

        for ligand in ligands:
            ligand_dir = Path(work_dir) / ligand
            all_results = []

            if ligand_dir.exists():
                for csv_file in ligand_dir.glob("top_score_*.csv"):
                    try:
                        df = pd.read_csv(csv_file)
                        all_results.append(df)
                    except Exception as e:
                        logger.warning(f"Failed to read {csv_file}: {e}")

            if all_results:
                # Combine all results
                combined_df = pd.concat(all_results, ignore_index=True)
                combined_df = combined_df.sort_values('top_score').reset_index(drop=True)

                # Save to scan_results
                combined_df.to_csv(scan_results_dir / f"top_score_{ligand}.csv", index=False)
                results[ligand] = combined_df
                logger.info(f"Parsed {len(combined_df)} results for {ligand}")
            else:
                logger.warning(f"No results found for {ligand}")
                results[ligand] = pd.DataFrame()

        return results

    except Exception as e:
        logger.error(f"Error parsing results: {e}")
        return {}

def apply_promiscuity_filter(docking_results, ligand_name, promiscuous_targets, work_dir):
    """Apply 25% promiscuity filter"""
    logger.info(f"Applying 25% promiscuity filter to {len(docking_results)} {ligand_name} results")

    if len(docking_results) == 0:
        return docking_results, pd.DataFrame()

    # Add pre-filtering ranks
    docking_results = docking_results.copy()
    docking_results['pre_filter_rank'] = range(1, len(docking_results) + 1)

    # Identify promiscuous targets
    docking_results['is_promiscuous'] = docking_results['gene_name'].isin(promiscuous_targets)

    # Separate data
    promiscuous_df = docking_results[docking_results['is_promiscuous']].copy()
    filtered_df = docking_results[~docking_results['is_promiscuous']].copy()

    # Add post-filtering ranks
    if len(filtered_df) > 0:
        filtered_df['post_filter_rank'] = range(1, len(filtered_df) + 1)

    # Stats
    removal_pct = (len(promiscuous_df) / len(docking_results)) * 100 if len(docking_results) > 0 else 0

    logger.info(f"Promiscuity filtering for {ligand_name}:")
    logger.info(f"  - Original: {len(docking_results)}")
    logger.info(f"  - Removed: {len(promiscuous_df)} ({removal_pct:.1f}%)")
    logger.info(f"  - Remaining: {len(filtered_df)}")

    # Save results
    filter_dir = Path(work_dir) / "results" / f"{ligand_name}_promiscuity_filtered"
    filter_dir.mkdir(parents=True, exist_ok=True)

    docking_results.to_csv(filter_dir / f"{ligand_name}_pre_filtering.csv", index=False)
    if len(promiscuous_df) > 0:
        promiscuous_df.to_csv(filter_dir / f"{ligand_name}_removed_promiscuous.csv", index=False)
    if len(filtered_df) > 0:
        filtered_df.to_csv(filter_dir / f"{ligand_name}_post_filtering.csv", index=False)

    return filtered_df, promiscuous_df

def run_simple_pipeline(genes_data, ligands_data, work_dir="simple_proteome_scan", max_targets=20):
    """Simple pipeline - just download PDBs and run basic analysis"""

    logger.info("=== SIMPLE PROTEOME SCAN PIPELINE ===")

    # Create work directory
    work_dir = Path(work_dir)
    work_dir.mkdir(exist_ok=True)

    # Load promiscuity targets
    promiscuous_targets = load_promiscuity_targets()

    results = {
        'genes_processed': {},
        'ligands_processed': {},
        'docking_results': {},
        'summary': {}
    }

    # Step 1: Setup genes with PDBs
    logger.info("Step 1: Setting up genes and downloading PDBs")

    valid_genes = []
    gene_names_list = []

    for gene_data in genes_data:
        gene_name = gene_data['gene_name']

        # Use provided PDB IDs or defaults
        if 'pdb_ids' in gene_data:
            pdb_ids = gene_data['pdb_ids']
        else:
            # Default PDB IDs for common genes
            default_pdbs = {
                'BRAF': ['4XV2', '5VAL', '3NY5'],
                'MEK1': ['3EQH', '4AN2', '3PP1'],
                'MEK2': ['1S9I', '3EQC', '4AN3']
            }
            pdb_ids = default_pdbs.get(gene_name, [gene_name + '_default'])

        success = setup_gene_directories(gene_name, pdb_ids, str(work_dir))

        if success:
            valid_genes.append(gene_data)
            gene_names_list.append(gene_name)
            results['genes_processed'][gene_name] = {'status': 'SUCCESS', 'pdb_count': len(pdb_ids)}
            logger.info(f"✓ {gene_name}: Setup successful")
        else:
            results['genes_processed'][gene_name] = {'status': 'FAILED', 'error': 'PDB download failed'}
            logger.error(f"✗ {gene_name}: Setup failed")

    logger.info(f"Successfully processed {len(valid_genes)}/{len(genes_data)} genes")

    if len(valid_genes) == 0:
        logger.error("No valid genes. Stopping.")
        return results

    # Step 2: Simple docking
    logger.info("Step 2: Running simple docking")

    ligand_names = []

    for ligand_data in ligands_data:
        ligand_name = ligand_data['ligand_name']
        sdf_path = ligand_data['sdf_path']

        if not os.path.exists(sdf_path):
            logger.error(f"Ligand SDF not found: {sdf_path}")
            continue

        ligand_names.append(ligand_name)

        logger.info(f"Processing {ligand_name} against {len(gene_names_list)} genes")

        docking_successes = 0
        for gene_name in gene_names_list:
            success = run_simple_docking(gene_name, ligand_name, sdf_path, str(work_dir))
            if success:
                docking_successes += 1

        results['ligands_processed'][ligand_name] = {
            'status': 'SUCCESS',
            'docking_successes': docking_successes,
            'total_genes': len(gene_names_list)
        }

        logger.info(f"✓ {ligand_name}: {docking_successes}/{len(gene_names_list)} docking successes")

    # Step 3: Parse and filter results
    logger.info("Step 3: Parsing and filtering results")

    if ligand_names:
        parsed_results = parse_simple_results(ligand_names, str(work_dir))

        for ligand_name, docking_df in parsed_results.items():
            if len(docking_df) == 0:
                logger.warning(f"No results for {ligand_name}")
                continue

            logger.info(f"Collected {len(docking_df)} results for {ligand_name}")

            # Apply promiscuity filtering
            filtered_df, removed_df = apply_promiscuity_filter(
                docking_df, ligand_name, promiscuous_targets, str(work_dir))

            # Select top targets
            top_targets = filtered_df.head(max_targets)

            results['docking_results'][ligand_name] = {
                'total_results': len(docking_df),
                'filtered_results': len(filtered_df),
                'removed_promiscuous': len(removed_df),
                'top_targets': len(top_targets)
            }

            logger.info(f"Selected top {len(top_targets)} targets for {ligand_name}")

            # Save top results
            if len(top_targets) > 0:
                results_dir = work_dir / "results"
                results_dir.mkdir(exist_ok=True)
                top_targets.to_csv(results_dir / f"{ligand_name}_top_targets.csv", index=False)

                print(f"\n=== Top {ligand_name} Results ===")
                for _, row in top_targets.head(5).iterrows():
                    print(f"{row['gene_name']:<10} | {row['pdb_id']:<6} | Score: {row['top_score']:>6.2f}")

    # Save comprehensive results
    results_file = work_dir / "results" / "simple_pipeline_results.json"
    results_file.parent.mkdir(exist_ok=True)
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)

    logger.info(f"Results saved to: {results_file}")
    logger.info("=== SIMPLE PIPELINE COMPLETED ===")

    return results

def run_example():
    """Simple example with known PDB IDs"""

    # Genes with specific PDB IDs
    genes_data = [
        {'gene_name': 'BRAF', 'pdb_ids': ['4XV2', '5VAL']},
        {'gene_name': 'MEK1', 'pdb_ids': ['3EQH']},
        {'gene_name': 'MEK2', 'pdb_ids': ['1S9I']},
    ]

    # Ligands
    ligands_data = [
        {'ligand_name': 'Trametinib', 'sdf_path': 'data/ligands/processed/Processed_Trametinib.sdf'},
        {'ligand_name': 'Dabrafenib', 'sdf_path': 'data/ligands/processed/Processed_Dabrafenib.sdf'},
    ]

    return run_simple_pipeline(genes_data, ligands_data, "simple_example", max_targets=10)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Simple ProteomeScan Pipeline')
    parser.add_argument('--example', action='store_true', help='Run simple example')

    args = parser.parse_args()

    if args.example:
        logger.info("Running simple example")
        results = run_example()
    else:
        logger.error("Must specify --example")
        parser.print_help()
        sys.exit(1)

    logger.info("Simple pipeline completed!")
