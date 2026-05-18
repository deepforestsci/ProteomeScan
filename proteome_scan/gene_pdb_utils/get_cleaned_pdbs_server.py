"""
Server-backed pdb_clean for ProteomeScan.

Calls the deepchem-server pdb_clean primitive via pyds, then downloads all
artifacts (cleaned PDB files + metadata CSV) to a local gene directory.
Returns a DataFrame identical in structure to the local get_cleaned_pdbs().
"""
from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path
from typing import Optional

import pandas as pd

try:
    from pyds import Settings, Data
    from pyds.primitives.proteome_scan import PdbClean
    PYDS_AVAILABLE = True
except ImportError:
    PYDS_AVAILABLE = False


def get_cleaned_pdbs_server(
    gene_name: str,
    entry_id: str,
    scan_id: str = "default_scan",
    output: Optional[str] = None,
    server_url: str = "http://localhost:8000",
    profile: str = "default",
    project: str = "proteomescan",
    min_res_val: float = 2.5,
) -> pd.DataFrame:
    """Download and clean PDB structures for a gene via the deepchem-server.

    Submits a pdb_clean job to the server, waits for completion, then
    downloads all artefacts into the same local directory layout that the
    local get_cleaned_pdbs() produces.  The rest of the ProteomeScan pipeline
    (docking, parse_results) is therefore unaffected.

    Parameters
    ----------
    gene_name : str
        Gene symbol, e.g. 'GBA3'.
    entry_id : str
        UniProt accession, e.g. 'Q9Y4K5'.
    scan_id : str
        Unique run identifier used for server-side caching (pass the scan
        directory name to namespace artefacts per run).
    output : str, optional
        Datastore output prefix.  Defaults to gene_name.
    server_url : str
        Base URL of the running deepchem-server instance.
    profile : str
        pyds profile name (must exist in the server's profile registry).
    project : str
        pyds project name.
    min_res_val : float
        Maximum acceptable resolution in Angstroms.  Structures coarser than
        this are excluded from the selection.  Default 2.5.

    Returns
    -------
    pd.DataFrame
        Subset of PDB metadata for the selected structures, with an added
        ``path`` column containing the local path to each cleaned PDB file.
        Also saved to ``./{gene_name}/{gene_name}_pdbs.csv``.

    Raises
    ------
    ImportError
        If pyds is not installed.
    RuntimeError
        If the server call fails or no cleaned PDBs are returned.
    """
    if not PYDS_AVAILABLE:
        raise ImportError(
            "pyds is not installed; cannot use server-backed pdb_clean. "
            "Install it or set use_server=False to use the local implementation."
        )

    if output is None:
        output = gene_name

    settings = Settings(profile=profile, project=project, base_url=server_url)
    pdb_clean_client = PdbClean(settings)
    data_client = Data(settings)

    result = pdb_clean_client.run(
        gene_name=gene_name,
        entry_id=entry_id,
        scan_id=scan_id,
        output=output,
        min_res_val=min_res_val,
        profile_name=profile,
        project_name=project,
    )
    summary_address = result["pdb_clean_results_address"]

    gene_dir = Path(f"./{gene_name}")
    gene_dir.mkdir(exist_ok=True)

    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        data_client.download_data(
            summary_address, tmp_path,
            profile_name=profile, project_name=project,
        )
        with open(tmp_path) as f:
            summary = json.load(f)
    finally:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)

    cleaned_pdb_addresses: dict = summary["cleaned_pdb_addresses"]
    pdbs_metadata_csv_address: str = summary["pdbs_metadata_csv_address"]

    if not cleaned_pdb_addresses:
        raise RuntimeError(
            f"Server returned no cleaned PDB files for gene {gene_name} (entry {entry_id})"
        )

    csv_path = str(gene_dir / f"{gene_name}_pdbs.csv")
    data_client.download_data(
        pdbs_metadata_csv_address, csv_path,
        profile_name=profile, project_name=project,
    )
    df = pd.read_csv(csv_path, index_col="id")

    local_paths: dict[str, str] = {}
    for pdb_id, ds_address in cleaned_pdb_addresses.items():
        local_pdb = str(gene_dir / f"cleaned_g_{gene_name}_p_{pdb_id}.pdb")
        data_client.download_data(
            ds_address, local_pdb,
            profile_name=profile, project_name=project,
        )
        local_paths[pdb_id] = local_pdb

    df = df.loc[df.index.isin(local_paths)].copy()
    df["path"] = [local_paths[pid] for pid in df.index]
    df.to_csv(csv_path)

    return df
