import os
import shutil
import json
import pandas as pd
from proteome_scan import get_cleaned_pdbs, run_docking, parse_results
from proteome_scan.post_scan_analysis import filter_promiscuous_targets, run_multi_pose_analysis


with open('data/ligands/ligand_SMILES.json', 'r') as f:
    smiles_data = json.load(f)['ligand_smiles']
    f.close()

with open('data/ligands/ligands_known_targets.json', 'r') as f:
    ligands_known_targets = json.load(f)
    f.close()

proteome_scan_gene_df = pd.read_csv("data/ProteomeScan_7657_genes.csv")


def run_proteome_scan(ligands, gene_names, scan_dir):
    """
    Runs the proteome scan for the given ligands and gene names
    and stores top results in `{scan_dir}/scan_results`.
    """
    os.makedirs(scan_dir, exist_ok=True)
    # Add known target genes to the gene_names list
    ligand_smiles = []
    for ligand in ligands:
        gene_names.extend(ligands_known_targets[ligand]['target_genes'])
        ligand_smiles.append(smiles_data[ligand])

    # get ligand sdf addresses
    ligand_sdf_addresses = {ligand: f"data/ligands/processed/Processed_{ligand}.sdf" for ligand in ligands}

    # get pdbs and their metadata for each gene
    gene_metadata = {}
    for gene_name in gene_names:
        gene_metadata[gene_name] = {}
        s = proteome_scan_gene_df.loc[
        proteome_scan_gene_df['Gene Names (primary)'].eq(gene_name), 'Entry']
        if not s.empty:
            entry_id = s.iloc[0]
        else:
            raise ValueError(f"Gene {gene_name} not found in ProteomeScan gene dataframe")
        gene_metadata[gene_name]['entry_id'] = entry_id
        if os.path.exists(f"./{scan_dir}/{gene_name}"):
            print(f"skipping {gene_name} because it already exists in {scan_dir}")
            df_metadata = pd.read_csv(f"./{scan_dir}/{gene_name}/{gene_name}_pdbs.csv", index_col='id')
        else:
            print(f"fetching pdbs for {gene_name}")
            df_metadata = get_cleaned_pdbs(gene_name, entry_id)
            # move gene data to scan_dir
            print(f"moving {gene_name} data to {scan_dir}")
            shutil.move(f"./{gene_name}", f"./{scan_dir}/{gene_name}")
        gene_metadata[gene_name]['pdbs_metadata'] = df_metadata.to_dict('records')


    for gene_name in gene_names:
        for ligand in ligands:
            print(f"Docking selected {gene_name} PDBs with {ligand}")
            if run_docking(gene_name, ligand, ligand_sdf_addresses[ligand], scan_dir):
                print(f"Docking {gene_name} with {ligand} completed")
            else:
                print(f"Docking {gene_name} with {ligand} failed")

    parse_results(ligands, scan_dir) # stores results in scan_dir/scan_results

    # filter promiscuous targets
    thresholds = [
                    (15, 1), # in top 15% of all targets and common in 1 target
                ]
    filtered_results_dict = filter_promiscuous_targets(thresholds, scan_dir)
    print(filtered_results_dict)

    # run pose binding analysis (example for one gene)
    pose_analysis_gene = "GBA3"
    complexes_dir = os.path.join(scan_dir, pose_analysis_gene, "complexes")
    complexes = [os.path.join(complexes_dir, f) for f in os.listdir(complexes_dir)]
    df_results = run_multi_pose_analysis(complexes, np=8)
    os.makedirs("pose_analysis", exist_ok=True)
    df_results.to_csv(os.path.join("pose_analysis", f"{pose_analysis_gene}_pose_analysis.csv"), index=False)
    
    print("Proteome scan completed")


if __name__ == "__main__":
    ligands = ['Trametinib', 'Tucatinib']
    gene_names = ["GBA3", "SLC7A11", "FABP2", "CYP1A1"]
    scan_dir = "proteome_scan_test1"
    run_proteome_scan(ligands, gene_names, scan_dir)
