# Get Optimal Cleaned PDBs (per Gene)

## Pseudo-Code
### Algorithm: Get Optimal Cleaned PDBs

```
Required: gene name

query url ← BuildQuery(UniProt, gene name)
response ← GET(query url)
if response.status = 200 then
    data ← ParseData(response.body)
    canonical protein ID ← Extract(data)
else
    canonical protein ID ← null
end if
if canonical protein ID̸ = null then
    PDB_metadata ← [ ]
    query url ← BuildQuery(PDBe-KB, canonical protein ID)
    response ← GET(query url)
    if response.status = 200 then
        data ← ParseData(response.body)
        PDB_metadata ← Extract(data) {[pdb_id, seq start, seq end, resolution]}
        for all metadata in PDB_metadata do
            metadata.coverage ← metadata.end_position − metadata.start_position
        end for
    end if
end if
cleaned_PDB_list ← [ ]
while cleaned_PDB_list is empty and PDB_metadata is not empty do
    selected_PDBs ← Select_Optimal_PDBs(PDB_metadata)
    downloaded_PDBs ← [ ]
    for all pdb_id in selected_PDBs do
        download_url ← BuildQuery(PDBe, pdb_id)
        Path ← "∼ /{gene name}/{pdb_id}"
        if WGET(download_url, Path) then
            append (pdb_id, Path) to downloaded_PDBs
        else
            remove pdb_id from PDB_metadata
        end if
    end for
    if downloaded_PDBs != [ ] then
        for all (pdb_id, Path) in downloaded_PDBs do
            Cleaned_Path ← "∼ /{gene name}/clean {pdb_id}"
            if Clean_PDB(Path) then
                append (pdb_id, Cleaned_Path) to cleaned_PDB_list
            else
                remove pdb_id from PDB_metadata
            end if
        end for
    end if
end while
return cleaned_PDB_list
```

## Pseudo-Code Helper Functions

### Algorithm: Select_Optimal_PDBs

```
Input: A list of PDB metadata, each containing:
    • pdb_id
    • start_position (on protein sequence)
    • end_position (on protein sequence)
    • resolution
    • coverage
sorted_list ← Sort the list by start_position, then by end_position
selected_ranges ← [ ]
Set last_end_position ← −∞
for all range in sorted_list do
    if range.start_position > last_end_position then
        Add range to selected_ranges
        last_end_position ← range.end_position
    else
        Let last ← last element in selected_ranges
        if range.resolution < last.resolution then
            Replace last with range in selected_ranges
        else if abs(range.resolution − last.resolution) < threshold and range.coverage > last.coverage then
            Replace last with range in selected_ranges
        end if
        last_end_position ← max(last_end_position, range.end_position)
    end if
end for
selected_PDBs ← Map selected_ranges to their corresponding pdb_ids
return selected_PDBs
```

### Algorithm: Clean_PDB

```
Input: original PDB file path and output cleaned path
Result ← PDBFixer(Path, Cleaned_Path) 
    {
        Procedure using PDBFixer:
        • Remove chains not produced by the gene
        • Remove heterogens
        • Remove water molecules
        • Add hydrogens at pH 7
    }
if Result then
    return True
else
    return False
end if
```
