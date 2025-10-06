import os
import re
import pandas as pd
from pymol import cmd, stored
import shutil
import logging

# Set up logging
logger = logging.getLogger(__name__)


def separate_protein_ligand(pose_path, run_dir):
    # Load the complex PDB file
    cmd.load(pose_path, "complex")

    # Select the ligand(s) with residue names LIG or UNL
    cmd.select("ligand", "resn LIG+UNL")
    # cmd.select("ligand", "resn UNL")

    # Select the protein by excluding the ligand and water molecules
    cmd.select("protein_1", "polymer and not resn LIG+UNL and not resn HOH")

    # Save the ligand and protein selections to separate PDB files
    ligand_path = os.path.join(run_dir, "ligand.pdb")
    protein_path = os.path.join(run_dir, "protein.pdb")
    cmd.save(ligand_path, "ligand")
    cmd.save(protein_path, "protein_1")

    # Optional: Clean up selections
    # cmd.delete("ligand")
    # cmd.delete("protein_1")
    # cmd.delete("complex")
    cmd.reinitialize()
    return ligand_path, protein_path

def run_fpocket(protein_path):
    os.system(f"fpocket -f {protein_path}")

def parse_pocket_data(file_path):
    with open(file_path, 'r') as file:
        content = file.read()

    # Split the content into blocks for each pocket
    pocket_blocks = re.split(r'Pocket\s+\d+\s*:', content)
    pocket_headers = re.findall(r'Pocket\s+\d+\s*:', content)

    pockets = []

    for i, (header, block) in enumerate(zip(pocket_headers, pocket_blocks[1:])):  # Skip the first split as it will be empty
        pocket_data = {}
        pocket_data['pocket_id'] = i+1
        lines = block.strip().split('\n')
        for line in lines:
            if ':' in line:
                key, value = line.split(':', 1)
                key = key.strip()
                value = value.strip()
                # Attempt to convert numerical values to float or int
                try:
                    if '.' in value:
                        value = float(value)
                    else:
                        value = int(value)
                except ValueError:
                    pass  # Keep as string if conversion fails
                pocket_data[key] = value
        pockets.append(pocket_data)

    return pockets

def analyse_overlaps(ligand_path, pockets_out_path, pockets_df):
    # Load structures
    cmd.load(ligand_path)
    cmd.load(pockets_out_path)
    # Identify STP residues
    stored.list = []
    cmd.iterate("(resn STP)", "stored.list.append(resi)")
    lastSTP = stored.list[-1]

    # Hide lines for STP residues
    cmd.hide("lines", "resn STP")

    # Create and visualize pocket spheres
    for my_index in range(1, int(lastSTP) + 1):
        pocket_name = f"pocket{my_index}"
        cmd.select(pocket_name, f"resn STP and resi {my_index}")
        cmd.show("spheres", pocket_name)
        cmd.set("sphere_scale", 0.3, pocket_name)
        cmd.set("sphere_transparency", 0.1, pocket_name)
        cmd.color(my_index, pocket_name)

    overlap_within_alpha_spheres = []
    ratio_pocket_filled = []
    for my_index in range(1,int(lastSTP)+1): 
        cmd.do(f"select overlapping{my_index}, ligand within 3 of pocket{my_index}")
        cmd.do(f"select p_overlapping{my_index}, pocket{my_index} within 3 of ligand")
        overlap_within_alpha_spheres.append(cmd.count_atoms(f"overlapping{my_index}"))
        pocket_alphasphere_count = cmd.count_atoms(f"pocket{my_index}")
        ratio_p_filled = cmd.count_atoms(f"p_overlapping{my_index}")/pocket_alphasphere_count
        ratio_pocket_filled.append(ratio_p_filled)

    # Calculate VDW overlap between each pocket and the ligand
    VDWoverlap_data = []
    for my_index in range(1, int(lastSTP) + 1):
        pocket_name = f"pocket{my_index}"
        overlap_value = cmd.overlap("ligand", pocket_name)
        VDWoverlap_data.append(overlap_value)
        print(f"VDW overlap between {pocket_name} and ligand: {overlap_value:.3f}")
    ligand_atoms_count = cmd.count_atoms("ligand")
    pockets_df['Ligand atoms overlap count'] = pd.Series(overlap_within_alpha_spheres)
    pockets_df['% Ligand inside pocket'] = pockets_df['Ligand atoms overlap count']/ligand_atoms_count*100
    pockets_df['% Pocket Filled'] = pd.Series(ratio_pocket_filled)*100
    pockets_df['Ligand_VDWoverlap'] = pd.Series(VDWoverlap_data)
    analysis_df = pockets_df[['pocket_id', 'Ligand_VDWoverlap', '% Ligand inside pocket', '% Pocket Filled', 'Ligand atoms overlap count', 'Number of Alpha Spheres', 'Druggability Score', 'Pocket_Druggability_Rank', 'Volume', 'Mean local hydrophobic density']]
    cmd.reinitialize()
    return pockets_df, analysis_df

def main(pose_path, results_dir, is_clean_up=False):
    run_name = os.path.basename(pose_path).split('.')[0]
    run_dir = os.path.join(os.getcwd(), run_name)
    os.makedirs(run_dir, exist_ok=True)
    ligand_path, protein_path = separate_protein_ligand(pose_path, run_dir)
    run_fpocket(protein_path)
    pockets_out_path = os.path.join(run_dir,"protein_out", "protein_out.pdb")
    pockets_info_path = os.path.join(run_dir,"protein_out", "protein_info.txt")
    parsed_pockets = parse_pocket_data(pockets_info_path)
    if not parsed_pockets:
        if is_clean_up:
            shutil.rmtree(run_dir)
        raise ValueError("No pockets found!")
    pockets_df = pd.DataFrame(parsed_pockets)
    pockets_df['Pocket_Druggability_Rank'] = pockets_df['Druggability Score'].rank(ascending=False, method='min')
    pockets_df, analysis_df = analyse_overlaps(ligand_path, pockets_out_path, pockets_df)
    final_df = analysis_df[analysis_df['Ligand_VDWoverlap']>0].sort_values('Ligand_VDWoverlap')  # type: ignore[call-overload]
    output_dir = os.path.join(results_dir, run_name)
    os.makedirs(output_dir)
    pockets_df.to_csv(os.path.join(output_dir, "full_analysis.csv"), index=False)
    final_df.to_csv(os.path.join(output_dir, "result.csv"), index=False)

    # clean-up
    if is_clean_up:
        shutil.rmtree(run_dir)

def run_pose_analysis(complex_file, results_dir, gene_name, pdb_id, ligand_name):
    """Analyze protein-ligand complex using fpocket cavity detection.

    Parameters
    ----------
    complex_file : str
        Path to protein-ligand complex PDB file. Ligand must have residue
        name 'LIG' or 'UNL'.
    results_dir : str
        Directory where analysis results will be saved.
    gene_name : str
        Gene name for result identification.
    pdb_id : str
        PDB identifier for result identification.
    ligand_name : str
        Ligand name for result identification.

    Returns
    -------
    dict or None
        Result dictionary with keys: gene_name, pdb_id, ligand_name,
        total_coverage, top1_coverage, num_pockets, status.
        Returns None if analysis fails.

    Notes
    -----
    Separates protein and ligand using PyMOL, runs fpocket on protein to detect
    cavities, then calculates ligand overlap with detected pockets.

    total_coverage is sum of '% Ligand inside pocket' (capped at 100%).
    top1_coverage is max of '% Ligand inside pocket' (capped at 100%).

    Creates subdirectory in results_dir named complex_{gene_name}_{pdb_id}_{ligand_name}.

    Examples
    --------
    >>> result = run_pose_analysis(
    ...     'complex_BRAF_4XV2_Dabrafenib.pdb',
    ...     './pose_analysis/',
    ...     'BRAF',
    ...     '4XV2',
    ...     'Dabrafenib'
    ... )
    >>> if result:
    ...     print(f"Coverage: {result['total_coverage']:.1f}%")
    """
    try:
        logger.info(f"Running pose analysis for {gene_name}_{pdb_id}_{ligand_name}")

        if not os.path.exists(complex_file):
            logger.error(f"Complex file not found: {complex_file}")
            return None

        # Run analysis using the main function
        main(complex_file, results_dir, is_clean_up=True)

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


if __name__ == "__main__":
    pose_path = "" # complex pdb file
    results_dir = "./results"
    main(pose_path, results_dir, is_clean_up=True)
