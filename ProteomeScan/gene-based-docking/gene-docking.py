import warnings
warnings.simplefilter(action='ignore', category=FutureWarning) # ignore warnings in stdout

import pandas as pd
pd.options.mode.chained_assignment = None # ignore warnings in stdout
from deepchem.dock.pose_generation import VinaPoseGenerator
import datetime as dt
import os
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*') # ignore warnings in stdout
import zipfile

import logging
logging.getLogger().setLevel(logging.ERROR)

def unzip_file(zip_path, extract_to):
    # Ensure the extract_to directory exists
    if not os.path.exists(extract_to):
        os.makedirs(extract_to)

    # Unzip the file
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(extract_to)

def vina_docking(raw_pdb_address, raw_ligand_address):
    dir_name = "vina_docking_"+str(dt.datetime.now().isoformat())
    tmp = os.path.join(root, dir_name)
    os.mkdir(tmp)

    # without explicit pocket finding (blind docking)
    pg = VinaPoseGenerator()
    _, scores = pg.generate_poses(
                        molecular_complex=(raw_pdb_address, raw_ligand_address),
                        exhaustiveness=32,
                        num_modes=8,
                        out_dir=tmp,
                        generate_scores=True)
    return scores

if __name__ == "__main__":
    zip_path = './dependencies/genes-PDB.zip'
    extract_to = './dependencies'

    unzip_file(zip_path, extract_to)

    root = "dependencies"

    ligand = os.environ["LIGAND"] # input
    gene_name = os.environ["GENE_NAME"] # input
    raw_ligand_address = f"./dependencies/Processed_{ligand}.sdf"

    main_df = pd.DataFrame({'chains': {},
                            'resolution': {},
                            'coverage': {},
                            'top_score': {},
                            'scores': {},
                            'gene_name': {}})

    print(gene_name)
    df = pd.read_csv(f"./dependencies/genes-PDB/{gene_name}/{gene_name}_pdbs.csv", index_col='id')
    docking_scores = {}
    for id in list(df.index):
        raw_pdb_address = os.path.join("./dependencies/genes-PDB", df.loc[id]['path'])
        try:
            score = vina_docking(raw_pdb_address, raw_ligand_address)
        except Exception as e:
            print(f"Docking failed for {id} => {e}")
            score = [None]
        docking_scores[id] = score

    df_docked = pd.DataFrame({'id': docking_scores.keys(), 'scores': docking_scores.values()})
    df_docked.dropna(inplace=True)
    df_docked['top_score'] = df_docked['scores'].apply(lambda x: x[0])
    df_docked = df_docked.set_index('id')
    merged_df = pd.merge(df, df_docked, how='left', left_index=True, right_index=True)
    top_df = merged_df[['chains', 'resolution','coverage', 'top_score', 'scores']].sort_values(by='top_score')
    if len(top_df) >2:
        top2_df = top_df.iloc[[0,1]]
    else:
        top2_df = top_df
    top2_df['gene_name'] = [gene_name]*len(top2_df)
    main_df = pd.concat([main_df, top2_df])

    sorted_maindf = main_df.sort_values(by='top_score')
    sorted_maindf.to_csv(f"top_score_{gene_name}_{ligand}.csv")
