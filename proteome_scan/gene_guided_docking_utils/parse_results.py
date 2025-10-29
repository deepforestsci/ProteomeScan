import pandas as pd
import os

def concat_csv_from_folder(folder_path: str) -> pd.DataFrame:
    """
    Concatenate all csv files in a folder into a single dataframe.

    Parameters
    ----------
    folder_path: str
        Path to the folder containing the csv files.

    Returns
    -------
    main_df: pd.DataFrame
        Concatenated dataframe of all csv files in the folder.
    """
    main_df = pd.DataFrame()
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        if filename.endswith('.csv'):
            temp_df = pd.read_csv(file_path)
            main_df = pd.concat([main_df, temp_df], ignore_index=True)

    return main_df

def parse_results(ligands: list[str], scan_dir: str) -> None:
    """
    Parse the results of the docking runs for each ligand and save the top score csv files to the scan_results folder.

    Parameters
    ----------
    ligands: list[str]
        List of ligands.
    scan_dir: str
        Path to the scan directory.
    """
    os.makedirs(os.path.join(scan_dir, "scan_results"), exist_ok=True)
    for ligand in ligands:
        raw_results_folder_path = f'{scan_dir}/{ligand}'
        df = concat_csv_from_folder(raw_results_folder_path)
        df = df.sort_values(by="top_score")
        df = df.drop_duplicates("gene_name", keep='first')
        df = df.rename(columns={'Unnamed: 0': 'pdb_id'})
        df.to_csv(os.path.join(scan_dir, "scan_results", f"top_score_{ligand}.csv"), index=False)
