# ProteomeScan Pipeline

Complete end-to-end molecular docking pipeline for proteome-wide compound screening.

## Overview

ProteomeScan performs automated molecular docking of small molecules against protein targets, with integrated promiscuity filtering and binding pose analysis. The pipeline downloads PDB structures, runs AutoDock Vina docking, filters promiscuous targets, and analyzes binding poses using fpocket.

## Requirements

- Python 3.7+
- AutoDock Vina
- fpocket
- PyMOL
- Required Python packages: pandas, rdkit, deepchem, subprocess

## Quick Start

1. Prepare input configuration file (`pipeline_input/config.json`):
```json
{
  "work_dir": "pipeline_results",
  "genes": [
    {"gene_name": "BRAF", "pdb_ids": ["4XV2", "5VAL"]},
    {"gene_name": "MEK1", "pdb_ids": ["3EQH"]}
  ],
  "ligands": [
    {"ligand_name": "Trametinib", "sdf_file": "Processed_Trametinib.sdf"},
    {"ligand_name": "Dabrafenib", "sdf_file": "Processed_Dabrafenib.sdf"}
  ]
}
```

2. Place ligand SDF files in `pipeline_input/` directory

3. Run the pipeline:
```bash
python proteome_scan_simple_pipeline.py --config pipeline_input/config.json
```

## Pipeline Steps

1. **PDB Download**: Downloads protein structures from RCSB PDB database
2. **Molecular Docking**: Runs AutoDock Vina docking for each ligand-protein pair
3. **Promiscuity Filtering**: Removes targets with >25% promiscuity threshold
4. **Pose Analysis**: Analyzes binding poses using fpocket for pocket coverage
5. **Results Filtering**: Outputs complexes with >50% pocket coverage

## Input Files

- **Config JSON**: Specifies genes (with PDB IDs) and ligands (with SDF files)
- **SDF Files**: Small molecule structures in SDF format
- **Promiscuity Data**: `proteome_scan/promis_thresholds_may19.json`

## Output Structure

```
pipeline_results/
├── pdbs/                    # Downloaded PDB structures
├── docking/                 # Docking results and complexes
│   └── [gene_name]/
│       ├── complexes/       # Generated protein-ligand complexes
│       └── results/         # Docking scores and poses
├── pose_analysis/           # fpocket analysis results
└── results/                 # Final filtered results
    ├── promiscuity_filtered/ # After promiscuity filtering
    └── high_coverage_complexes/ # Complexes with >50% pocket coverage
```

## Key Output Files

- `[ligand]_top_targets.csv`: Top-ranked targets per ligand
- `[ligand]_promiscuity_filtered.csv`: Results after promiscuity filtering
- `high_coverage_complexes.csv`: Final complexes with >50% pocket coverage
- `pipeline_results.json`: Complete pipeline execution summary

## Example Usage

```bash
# Run with example configuration
python proteome_scan_simple_pipeline.py --config pipeline_input/config.json

# Check results
ls pipeline_results/results/high_coverage_complexes/
```

## Notes

- Pipeline performs real molecular docking calculations (not simulation)
- Uses 25% promiscuity threshold for target filtering
- Requires sufficient computational resources for large-scale screening
- All complexes with >50% pocket coverage are automatically saved to results