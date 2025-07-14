# ProteomeScan
ProteomeScan

## Pseudo-Code

Algorithm: Proteome Scan Algorithm

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