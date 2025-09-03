import pandas as pd
from deepchem.dock.pose_generation import VinaPoseGenerator
import datetime as dt
import os
from rdkit import Chem


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

def main(gene_name, ligand, ligand_address, scan_dir):
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
