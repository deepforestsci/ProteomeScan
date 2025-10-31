import pandas as pd
from deepchem.dock.pose_generation import VinaPoseGenerator
import datetime as dt
import os
from rdkit import Chem
import json
import requests

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


def vina_docking(raw_pdb_address: str,
                raw_ligand_address: str,
                scan_dir: str,
                exhaustiveness: int = 32,
                num_modes: int = 8) -> tuple[list, list]:
    """
    Run VINA docking for a given protein pdb file and ligand sdf file.

    Parameters
    ----------
    raw_pdb_address: str
        Path to the protein pdb file.
    raw_ligand_address: str
        Path to the ligand sdf file.
    scan_dir: str
        Path to the scan directory.
    exhaustiveness: int
        Exhaustiveness for VINA docking.
    num_modes: int
        Number of modes for VINA docking.

    Returns
    -------
    complex: list
        List of complex molecules.
    scores: list
        List of scores.

    Examples
    --------
    >>> vina_docking("1abc.pdb", "ligand.sdf", "./scan_dir", exhaustiveness=32, num_modes=8)
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
    return complex, scores

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
    
    # Configure settings properly (profile and project required by pyds API)
    settings = Settings()
    settings.set_profile("default")  # Default profile
    settings.set_project("proteomescan")  # Default project for ProteomeScan
    settings.set_base_url(server_url)
    
    # Setup clients with properly configured settings
    data_client = Data(settings)
    docking_client = Docking(settings)
    
    # Get profile and project for HTTP requests (pyds doesn't have get_data method)
    profile = settings.get_profile()
    project = settings.get_project()
    
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
        
        # Retrieve results: Since pyds doesn't have get_data() and server may not have endpoint,
        # we try to access datastore directly if available (for local server)
        results_address = result["docking_results_address"]
        print(f"Retrieving results from: {results_address}")
        
        # Try to access server datastore directly (works if server is local and we can import)
        try:
            from deepchem_server.core import config
            datastore = config.get_datastore()
            if datastore is not None:
                # Use server's internal datastore API
                results = datastore.get(results_address, kind='data', fetch_sample=False)
                if isinstance(results, str):
                    results = json.loads(results)
            else:
                raise ImportError("Datastore not available")
        except (ImportError, Exception):
            # Fallback: Return addresses only, let user handle retrieval
            print("WARNING: Cannot retrieve data automatically. Results addresses returned.")
            print(f"Results address: {results_address}")
            # Return minimal info - user will need server-side access to retrieve
            return None, [None]
            
        # Extract complex addresses and scores
        complex_addresses = results.get("complex_addresses", {})
        scores_address = results.get("scores_address")
        
        if not complex_addresses or not scores_address:
            print("No docking results found")
            return None, [None]
        
        # Retrieve scores
        print(f"Retrieving scores from: {scores_address}")
        try:
            scores_data = datastore.get(scores_address, kind='data', fetch_sample=False)
            if isinstance(scores_data, str):
                scores_data = json.loads(scores_data)
            scores = scores_data.get("scores", [])
        except Exception as e:
            print(f"Failed to retrieve scores: {e}")
            scores = []
        
        # Retrieve first complex PDB and convert to RDKit format to match local function
        first_mode = list(complex_addresses.keys())[0]
        complex_address = complex_addresses[first_mode]
        print(f"Retrieving complex from: {complex_address}")
        try:
            complex_data = datastore.get(complex_address, kind='data', fetch_sample=False)
            
            # Handle mdtraj.Trajectory - convert to PDB string first
            import tempfile
            if hasattr(complex_data, 'save_pdb'):
                with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp:
                    tmp_path = tmp.name
                complex_data.save_pdb(tmp_path)
                complex_pdb_path = tmp_path
            elif isinstance(complex_data, str):
                # If string, write to temp file for RDKit parsing
                with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp:
                    tmp.write(complex_data)
                    complex_pdb_path = tmp.name
            else:
                raise ValueError(f"Unexpected complex data type: {type(complex_data)}")
            
            # Parse PDB to extract protein and ligand separately
            # The complex PDB contains both protein (ATOM) and ligand (HETATM) records
            # We need to separate them to match the local function's return format
            
            # Read PDB file line by line to separate protein and ligand
            with open(complex_pdb_path, 'r') as f:
                pdb_lines = f.readlines()
            
            # Separate ATOM (protein) and HETATM (ligand) records
            protein_lines = []
            ligand_lines = []
            
            for line in pdb_lines:
                if line.startswith('ATOM'):
                    protein_lines.append(line)
                elif line.startswith('HETATM'):
                    ligand_lines.append(line)
                elif line.startswith('END') or line.startswith('TER'):
                    protein_lines.append(line)
            
            # Write separate temp files for protein and ligand
            with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp_protein:
                tmp_protein.write(''.join(protein_lines))
                protein_path = tmp_protein.name
            
            with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp_ligand:
                tmp_ligand.write(''.join(ligand_lines))
                ligand_path = tmp_ligand.name
            
            # Parse protein and ligand molecules separately
            protein_mol = Chem.rdmolfiles.MolFromPDBFile(protein_path)
            ligand_mol = Chem.rdmolfiles.MolFromPDBFile(ligand_path)
            
            # Clean up temp files
            os.unlink(complex_pdb_path)
            os.unlink(protein_path)
            os.unlink(ligand_path)
            
            if protein_mol is None or ligand_mol is None:
                print("Failed to parse protein or ligand from complex PDB")
                return None, [None]
            
        except Exception as e:
            print(f"Failed to retrieve/parse complex: {e}")
            import traceback
            traceback.print_exc()
            return None, [None]
        
        print(f"Downloaded and parsed complex for mode: {first_mode}")
        
        # Return in same format as local docking: ((protein_mol, ligand_mol),)
        return ((protein_mol, ligand_mol),), scores[:1]
        
    except Exception as e:
        print(f"Server docking failed: {e}")
        return None, [None]

def main(gene_name, ligand, ligand_address, scan_dir, use_server=False, server_url="http://localhost:8000"):
    try:
        if os.path.exists(f"{scan_dir}/{ligand}/top_score_{gene_name}_{ligand}.csv"):
            print(f"Skipping {gene_name} {ligand} docking because results file already exists")
            return True
        os.makedirs(f"{scan_dir}/{gene_name}/complexes", exist_ok=True)
        main_df = pd.DataFrame({'chains': {},
                                'resolution': {},
                                'coverage': {},
                                'top_score': {},
                                'scores': {},
                                'gene_name': {}})

        print("Docking for gene: ", gene_name)
        if not os.path.exists(f"{scan_dir}/{gene_name}/{gene_name}_pdbs.csv"):
            raise Exception(f"pdbs file not found for {gene_name}")
        df = pd.read_csv(f"{scan_dir}/{gene_name}/{gene_name}_pdbs.csv", index_col='id')
        docking_scores = {}
        for id in list(df.index):
            pdb_address = os.path.join(scan_dir, df.loc[id]['path'])
            try:
                assert os.path.exists(pdb_address), f"pdb file not found for {id}"
                print(ligand_address)
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
        top_df = merged_df[['chains', 'resolution','coverage', 'top_score', 'scores']].sort_values(by='top_score')
        if len(top_df) > 2:
            top2_df = top_df.iloc[[0,1]]
        else:
            top2_df = top_df
        top2_df['gene_name'] = [gene_name]*len(top2_df)
        main_df = pd.concat([main_df, top2_df])

        sorted_maindf = main_df.sort_values(by='top_score')
        os.makedirs(f"{scan_dir}/{ligand}", exist_ok=True)
        sorted_maindf.to_csv(f"{scan_dir}/{ligand}/top_score_{gene_name}_{ligand}.csv")
        return True
    except Exception as e:
        print(f"Error in docking {gene_name} {ligand} => {e}")
        return False
