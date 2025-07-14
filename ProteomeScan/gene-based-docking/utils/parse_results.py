import pandas as pd
import os

def concat_csv_from_folder(folder_path):
    main_df = pd.DataFrame()
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        if filename.endswith('.csv'):
            temp_df = pd.read_csv(file_path)
            main_df = pd.concat([main_df, temp_df], ignore_index=True)

    return main_df


all_ligands = [
                'Trametinib',
                'Binimetinib',
                'Dabrafenib',
                'Vemurafenib',
                'TAK-632',
                'Tucatinib',
                'Regorafenib',
                'SOS1-IN-11',
                'SOS1-IN-15',
                'Palbociclib',
                'BAY-293',
                'BI-3406',
                'Alpelisib',
                'Erlotinib',
                'Olaparib',
                'RMC6236',
                'Sapanisertib',
                'Selinexor',
                'SN-38',
                'Vociprotafib'
            ]

for ligand in all_ligands:
    raw_results_folder_path = f'{ligand}_docking_scores'
    df = concat_csv_from_folder(raw_results_folder_path)
    df = df.sort_values(by="top_score")
    df = df.drop_duplicates("gene_name", keep='first')
    df.to_csv(os.path.join("top_scores", f"top_score_{ligand}.csv"), index=False)
