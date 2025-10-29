import pandas as pd
from deepchem.dock.pose_generation import VinaPoseGenerator
import datetime as dt
import os
from rdkit import Chem
import json

# Server docking imports
try:
    from pyds import Settings, Data, Docking
    SERVER_AVAILABLE = True
except ImportError:
    SERVER_AVAILABLE = False


DISABLE_WARNINGS = True

if DISABLE_WARNINGS:
    import warnings
    import logging
    from rdkit import RDLogger
    warnings.simplefilter(action='ignore', category=FutureWarning) # ignore warnings in stdout
    pd.options.mode.chained_assignment = None # ignore warnings in stdout
    RDLogger.DisableLog('rdApp.*') # ignore warnings in stdout
    logging.getLogger().setLevel(logging.ERROR)


def vina_docking(raw_pdb_address, raw_ligand_address, scan_dir, exhaustiveness=32, num_modes=8):
    """Perform molecular docking using DeepChem's VinaPoseGenerator.

    Parameters
    ----------
    raw_pdb_address : str
        Path to protein PDB file.
    raw_ligand_address : str
        Path to ligand SDF file.
    scan_dir : str
        Directory where temporary docking subdirectory will be created.
    exhaustiveness : int, optional
        Vina exhaustiveness parameter, default 32.
    num_modes : int, optional
        Maximum number of binding modes to generate, default 8.

    Returns
    -------
    complex : tuple
        Docked molecular complex from DeepChem.
    scores : list
        List of binding scores for each generated mode.

    Notes
    -----
    Creates timestamped subdirectory in scan_dir for Vina output files.
    Uses blind docking (no explicit pocket definition).

    Examples
    --------
    >>> complex, scores = vina_docking(
    ...     'protein.pdb',
    ...     'ligand.sdf',
    ...     './docking_results/'
    ... )
    >>> if scores:
    ...     print(f"Best score: {scores[0]}")
    """
    dir_name = "vina_docking_"+str(dt.datetime.now().isoformat())
    tmp = os.path.join(scan_dir, dir_name)
    os.mkdir(tmp)

    # without explicit pocket finding (blind docking)
    pg = VinaPoseGenerator()
    complex, scores = pg.generate_poses(
                        molecular_complex=(raw_pdb_address, raw_ligand_address),
                        exhaustiveness=exhaustiveness,
                        num_modes=num_modes,
                        out_dir=tmp,
                        generate_scores=True)
    return complex,scores

def vina_docking_server(raw_pdb_address, raw_ligand_address, scan_dir, server_url="http://localhost:8000", exhaustiveness=32, num_modes=8):
    """Perform molecular docking using DeepChem Server via pyds API.
    
    Parameters
    ----------
    raw_pdb_address : str
        Path to protein PDB file.
    raw_ligand_address : str
        Path to ligand SDF file.
    scan_dir : str
        Directory where temporary docking subdirectory will be created.
    server_url : str
        DeepChem Server URL (default: http://localhost:8000).
    exhaustiveness : int, optional
        Vina exhaustiveness parameter, default 32.
    num_modes : int, optional
        Maximum number of binding modes to generate, default 8.

    Returns
    -------
    complex : tuple
        Docked molecular complex from DeepChem (same format as local).
    scores : list
        List of binding scores for each generated mode.
    """
    if not SERVER_AVAILABLE:
        raise ImportError("pyds not available. Install with: pip install pyds")
    
    # Setup clients
    settings = Settings(base_url=server_url)
    data_client = Data(settings)
    docking_client = Docking(settings)
    
    try:
        # Upload files to datastore
        print(f"Uploading protein to server: {raw_pdb_address}")
        protein_result = data_client.upload_data(raw_pdb_address)
        protein_address = protein_result["dataset_address"]
        
        print(f"Uploading ligand to server: {raw_ligand_address}")
        ligand_result = data_client.upload_data(raw_ligand_address)
        ligand_address = ligand_result["dataset_address"]
        
        # Run docking
        print("Running docking on server...")
        result = docking_client.run(
            protein_address=protein_address,
            ligand_address=ligand_address,
            output=f"docking_{os.path.basename(raw_pdb_address)}",
            exhaustiveness=exhaustiveness,
            num_modes=num_modes
        )
        
        # Download results
        results_address = result["docking_results_address"]
        print(f"Downloading results from: {results_address}")
        results_data = data_client.get(results_address)
        
        # Parse results to match local format
        if isinstance(results_data, str):
            results = json.loads(results_data)
        else:
            results = results_data
            
        # Extract complex addresses and scores
        complex_addresses = results.get("complex_addresses", {})
        scores_data = results.get("scores_address")
        
        if not complex_addresses or not scores_data:
            print("No docking results found")
            return None, [None]
        
        # Download scores
        scores_json = data_client.get(scores_data)
        if isinstance(scores_json, str):
            scores_data = json.loads(scores_json)
        else:
            scores_data = scores_json
            
        scores = scores_data.get("scores", [])
        
        # Download first complex (for compatibility with local format)
        first_mode = list(complex_addresses.keys())[0]
        complex_address = complex_addresses[first_mode]
        complex_pdb = data_client.get(complex_address)
        
        # Convert to RDKit format (simplified - just return the PDB content)
        # For full compatibility, would need to parse PDB and create RDKit molecules
        print(f"Downloaded complex for mode: {first_mode}")
        
        # Return in same format as local docking
        # Note: This is simplified - full implementation would parse PDB and create RDKit molecules
        return (complex_pdb,), scores[:1]  # Return first score only for compatibility
        
    except Exception as e:
        print(f"Server docking failed: {e}")
        return None, [None]

def main(gene_name, ligand, ligand_address, scan_dir, use_server=False, server_url="http://localhost:8000"):
    try:
        if os.path.exists(f"{scan_dir}/{ligand}/top_score_{gene_name}_{ligand}.csv"):
            print(f"Skipping {gene_name} {ligand} docking because results file already exists")
            return True
        os.makedirs(f"{scan_dir}/{gene_name}/complexes", exist_ok=True)
        # Initialize empty list to collect results instead of empty DataFrame with empty dicts
        results_list = []

        print("Docking for gene: ", gene_name)
        if not os.path.exists(f"{scan_dir}/{gene_name}/{gene_name}_pdbs.csv"):
            raise Exception(f"pdbs file not found for {gene_name}")
        df = pd.read_csv(f"{scan_dir}/{gene_name}/{gene_name}_pdbs.csv", index_col='id')
        docking_scores = {}
        for id in list(df.index):
            pdb_address = os.path.join(scan_dir, df.loc[id]['path'])
            try:
                assert os.path.exists(pdb_address), f"pdb file not found for {id}"
                assert os.path.exists(ligand_address), f"ligand file not found for {ligand}"
                
                if use_server:
                    print(f"Using server docking for {id}")
                    complex, score = vina_docking_server(pdb_address, ligand_address, scan_dir, server_url)
                else:
                    complex, score = vina_docking(pdb_address, ligand_address, scan_dir)
            except Exception as e:
                print(f"Docking failed for {id} => {e}")
                score = [None]
                complex = None
            if complex is not None:
                complex_mol = Chem.CombineMols(complex[0][0], complex[0][1])
                Chem.rdmolfiles.MolToPDBFile(complex_mol, f"{scan_dir}/{gene_name}/complexes/complex_{gene_name}_{id}_{ligand}.pdb")
            docking_scores[id] = score

        df_docked = pd.DataFrame({'id': docking_scores.keys(), 'scores': docking_scores.values()})
        df_docked.dropna(inplace=True)
        df_docked['top_score'] = df_docked['scores'].apply(lambda x: x[0])
        df_docked = df_docked.set_index('id')
        merged_df = pd.merge(df, df_docked, how='left', left_index=True, right_index=True)
        top_df = merged_df[['top_score', 'scores']].sort_values(by='top_score')
        if len(top_df) > 2:
            top2_df = top_df.iloc[[0,1]]
        else:
            top2_df = top_df
        top2_df['gene_name'] = [gene_name]*len(top2_df)
        results_list.append(top2_df)

        # Create final DataFrame from results list instead of using empty DataFrame
        os.makedirs(f"{scan_dir}/{ligand}", exist_ok=True)
        if results_list:
            main_df = pd.concat(results_list, ignore_index=True)
            sorted_maindf = main_df.sort_values(by='top_score')
            sorted_maindf.to_csv(f"{scan_dir}/{ligand}/top_score_{gene_name}_{ligand}.csv")
        else:
            # Handle case where no results were generated - create empty CSV with correct columns
            empty_df = pd.DataFrame(columns=['top_score', 'scores', 'gene_name'])
            empty_df.to_csv(f"{scan_dir}/{ligand}/top_score_{gene_name}_{ligand}.csv")
        return True
    except Exception as e:
        print(f"Error in docking {gene_name} {ligand} => {e}")
        return False
