import os
import pandas as pd
import json
from typing import Dict, List


def get_promiscuous_targets(thresholds: List[tuple], scan_dir: str) -> dict:
    """
    Get the promiscuous targets for the given thresholds and scan results directory
    
    Parameters
    ----------
    thresholds: List[tuple]
        list of tuples of (m, n)
    scan_results_dir: str
        directory containing the scan results

    Returns
    -------
    data_promis: Dict
        dictionary of (target_threshold, n) and the promiscuous targets
    """
    scan_results_dir: str = os.path.join(scan_dir, "scan_results")
    data_promis: Dict = {}
    for target_threshold, n in thresholds:
        target_threshold/=100
        df_list: List[pd.DataFrame] = []
        for filename in os.listdir(scan_results_dir):
            file_path: str = os.path.join(scan_results_dir, filename)

            # Check if the file is a CSV
            if filename.endswith('.csv'):
                # Read the CSV file into a DataFrame
                temp_df: pd.DataFrame = pd.read_csv(file_path)
                df_list.append(temp_df[:int(len(temp_df)*target_threshold)+1])

        f_df: pd.DataFrame = pd.concat(df_list)['gene_name'].value_counts().reset_index()
        f_df.rename(columns={'index': 'gene name', 'gene_name': 'occurrence in top'}, inplace=True)
        common_occurence: List[str] = f_df[f_df['occurrence in top']>=n]['gene name'].to_list()
        raw_key: tuple = (str(int(target_threshold*100)), str(n))
        key: str = "%_".join(raw_key)
        data_promis[key] = common_occurence

    os.makedirs(os.path.join(scan_dir, "promiscuity_analysis"), exist_ok=True)
    save_promis_path: str = os.path.join(scan_dir, "promiscuity_analysis", "promiscuous_targets.json")
    with open(save_promis_path, "w") as f:
        json.dump(data_promis, f, indent=2)
        f.close()

    return data_promis


def filter_promiscuous_targets(thresholds: list[tuple], scan_dir: str) -> dict:
    """
    Filter the promiscuous targets for the given thresholds and scan results directory

    Parameters
    ----------
    thresholds: list[tuple]
        list of tuples of (m, n)
    scan_dir: str
        directory containing the scan results

    Returns
    -------
    filtered_results_dict: dict
        dictionary of (target_threshold, n) and the filtered results directory
    """
    scan_results_dir: str = os.path.join(scan_dir, "scan_results")
    data_promis: Dict = get_promiscuous_targets(thresholds, scan_dir)

    filtered_results_dict: Dict = {}
    for threshold, promiscuous_targets in data_promis.items():
        filtered_results_dir: str = os.path.join(scan_dir, f"scan_results_promiscuity_filtered_{threshold}")
        os.makedirs(filtered_results_dir, exist_ok=True)
        for filename in os.listdir(scan_results_dir):
            file_path: str = os.path.join(scan_results_dir, filename)

            # Check if the file is a CSV
            if filename.endswith('.csv'):
                # Read the CSV file into a DataFrame
                temp_df: pd.DataFrame = pd.read_csv(file_path)
                temp_df = temp_df[~temp_df['gene_name'].isin(promiscuous_targets)]
                temp_df.to_csv(os.path.join(filtered_results_dir, filename), index=False)
        filtered_results_dict[threshold] = filtered_results_dir
    return filtered_results_dict
