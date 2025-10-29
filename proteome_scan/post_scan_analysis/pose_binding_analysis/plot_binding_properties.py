import matplotlib.pyplot as plt
import pandas as pd
import os

def plot_grouped_bars(df: pd.DataFrame, columns: list, title: str, root: str) -> None:
    """
    Plot the grouped bars for a given dataframe for each pocket id in complex.

    Parameters
    ----------
    df: pd.DataFrame
        Dataframe to plot.
    columns: list
        List of columns to plot.
    title: str
        Title of the plot.
    root: str
        Root directory to save the plot.

    Returns
    -------
    None
    """
    if len(columns) != 3:
        raise ValueError("This function requires exactly 3 columns.")

    x = range(len(df))
    width = 0.25

    plt.figure(figsize=(10, 6))
    plt.bar([i - width for i in x], df[columns[0]], width=width, label=columns[0])
    plt.bar(x, df[columns[1]], width=width, label=columns[1])

    plt.scatter(x, df[columns[2]], color='black', label=columns[2], zorder=5)

    plt.xlabel("Pocket Id")
    plt.ylabel("Value")
    plt.title(f"{title} - Ligand Pocket Analysis")
    plt.xticks(ticks=x, labels=x)
    plt.legend()
    plt.grid(True, axis='y')
    plt.tight_layout()
    plt.savefig(os.path.join(root,f"{title}.png"))
    plt.close()


if __name__ == "__main__":
    root = "complexes"
    complexes = os.listdir(root)
    results_dir = "/results" # from multi-pose run

    for complex_path in complexes:
        datapoint = {}
        run_name = os.path.basename(complex_path).split('.')[0]
        results_path = os.path.join(results_dir, run_name, "full_analysis.csv")

        df = pd.read_csv(results_path)
        df1 = df[['% Ligand inside pocket', 'Pocket_Druggability_Rank', 'Druggability Score']].sort_values(by='Pocket_Druggability_Rank')
        max_druggability_score = df1['Druggability Score'].max()
        df1['Relative Druggability Score'] = df1['Druggability Score'].apply(lambda x: x/max_druggability_score)
        df1['Rank Score'] = df1['Pocket_Druggability_Rank'].apply(lambda x: (len(df1)-x)/(len(df1)-1))
        df1['% Ligand inside pocket'] = df1['% Ligand inside pocket']/100
        df2 = df1[df1['% Ligand inside pocket']>=0]
        plot_grouped_bars(df2, ['% Ligand inside pocket', 'Relative Druggability Score', 'Rank Score'], title=run_name, root="binding_image_analysis")
