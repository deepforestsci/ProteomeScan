# ProteomeScan
ProteomeScan is a Python toolkit for large-scale human proteome scan screening against a given set of drugs/ligands.

## Requirements

- Python 3.10 and above
- Create a conda environment and install dependencies from `requirements.txt`
- Internet connection (for UniProt/PDBe/PDBe-KB access)
- Your ligand files (see `/data` folder below)

## Data Folder

- `data/gene_selection_raw`: Contains files used to filter Human Proteome and shortlist gene names for ProteomeScan.
- `data/ligands/processed/`: **Place your ligand SDF files here** before running the pipeline. Each ligand should be provided as its own file named as `Processed_<LigandName>.sdf`.
- If you have raw Ligand SMILES, use `data/ligands/prepare_ligands.py`.

## Quick Start

**To run the full pipeline, simply execute:**

```bash
python example.py
```

`example.py` is provided as an executable, end-to-end pipeline with example gene and ligand lists. You can customize these lists or provide your own, and place the corresponding ligand `.sdf` files in `data/ligands/processed/`.

What happens:
- For each gene, the tool fetches optimal (curated, best) PDBs by scraping UniProt, PDBe, and PDBe-KB.
- Each ligand is docked against optimal PDBs for each gene target.
- All results and logs are saved in the provided output directory.
- Post-scan run includes promiscuity filtering and pose analysis.

## Pseudo-Code

For workflow logic at a high level, see the pseudo-code below:

```
Require: ligands list, gene names list

N_ligands ← length(ligands list)
N_gene_names ← length(gene names list)
Optimal PDBs ← ∅
for each gene name in gene names list do
    Optimal PDBs ← get_optimal_cleaned pdbs(gene name)
end for
Scores ← ∅
for each ligand in ligands list do
    for each gene name in gene names list do
        Scores[ligand][gene name] ← run_gene_based_docking(ligand, gene name)
    end for
end for
Ranked_Scores ← ∅
for each ligand in ligands list do
    Ranked_Scores[ligand] ← Sort(Scores[ligand])
end for
Promiscuous targets ← Analyse_Target_Promiscuity(N_ligands, N_gene_names, Ranked_Scores)
Updated_Promiscuous_targets ← [ ]
for each target in Promiscuous targets do
    if Check_binding(target, ligands list) then
        Append target to Updated_Promiscuous_targets
    end if
end for
for each ligand in ligands list do
    Remove Updated_Promiscuous_targets from Ranked_Scores[ligand]
end for
return Ranked_Scores per ligand
```