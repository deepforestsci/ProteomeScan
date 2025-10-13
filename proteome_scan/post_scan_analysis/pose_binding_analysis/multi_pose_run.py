import os
import pandas as pd
import multiprocessing
from proteome_scan.post_scan_analysis.pose_binding_analysis.analyse_pose_script import main as analyse_pose

temp_results_dir = "/temp_results"

def run_analyse_pose(complex_path):
    try:
        complex_path_ = os.path.join(root, complex_path)
        analyse_pose(pose_path=complex_path_, results_dir=temp_results_dir, is_clean_up=True)
        return True
    except Exception as e:
        print(f"error {complex_path}: {e}")
        return False

def get_overall_ligand_interactions(df):
    percents = df['% Ligand inside pocket']
    total_percent = sum(percents) if sum(percents) <= 100 else 100
    return total_percent

def get_total_top_n_bucket_percentages(df, n):
    """
    Extract the percentages of liquid in the top-n ranked buckets based on scores.
    """
    sum_ = df[df['Pocket_Druggability_Rank']<=n]['% Ligand inside pocket'].to_numpy().sum()
    sum_ = min(sum_, 100)
    return sum_

def run_multi_pose_analysis(complexes, np=4):
    global temp_results_dir

    with multiprocessing.Pool(processes=np) as pool:
        results = pool.map(run_analyse_pose, complexes)

    failed_poses = []
    for c, r in zip(complexes, results):
        if not r:
            failed_poses.append(c)

    data = []
    for complex_path in complexes:
        datapoint = {}
        run_name = os.path.basename(complex_path).split('.')[0]
        datapoint['complex'] = run_name
        if complex_path in failed_poses:
            datapoint['total % Ligand inside pocket'] = None
            for n in [1, 5, 10]:
                datapoint[f'total % Ligand inside pockets (top{n} pockets)'] = None
            data.append(datapoint)
            continue
        results_path = os.path.join(temp_results_dir, run_name, "full_analysis.csv")
        df = pd.read_csv(results_path)
        total_percent = get_overall_ligand_interactions(df)
        datapoint['total % Ligand inside pockets'] = total_percent

        for n in [1, 5, 10]:
            datapoint[f'total % Ligand inside pockets (top{n} pockets)'] = get_total_top_n_bucket_percentages(df, n)
        data.append(datapoint)

    df_results = pd.DataFrame(data)
    # df_results['gene_name'] = df_results['complex'].apply(lambda x: x.split("_")[1])
    # df_results['ligand'] = df_results['complex'].apply(lambda x: x.split("_")[2])
    return df_results


if __name__ == '__main__':
    root = "complexes"
    complexes = os.listdir(root)
    df_results = run_multi_pose_analysis(complexes, np=8)
    df_results.to_csv("pose_analysis.csv", index=False)
