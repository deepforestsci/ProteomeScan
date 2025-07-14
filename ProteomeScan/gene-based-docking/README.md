# Gene-based Docking

## Pseudo-Code
### Algorithm: Run Gene-Based Docking

```
Input: ligand, gene name

Ligand_pdb ← get_pdbs(ligand) {Convert SMILES SDF then PDB}
(PDB_ids, Cleaned_pdb_path_list) ← get_pdbs(gene name)
Result ← {}
for (pdb_id, protein_pdb) in (PDB_ids, Cleaned_pdb_path_list) do
    Scores ← VINA(protein_pdb, Ligand_pdb, exhaustiveness=32, mode=8)
    Top_score ← min(Scores)
    Result[pdb_id] ← Top_score
end for
(Best_PDB_id, Best_Top_score) ← entry in Result with minimum score
return (Best_PDB_id, Best_Top_score)
```
