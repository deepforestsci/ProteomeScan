import requests
import os
import gzip
import json
from tqdm import tqdm
import re
import pandas as pd

import os
from concurrent.futures import ThreadPoolExecutor


def run_on_multiple_threads(fn, values, max_workers):
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        results = list(tqdm(executor.map(fn, values), total=len(values), desc="Processing"))
    return list(results)
# %%
def get_canon_pdb(gene_name):
    print("---------")
    print(gene_name)

    p_id = None
    search_results = requests.get(f"https://rest.uniprot.org/uniprotkb/search?compressed=true&format=json&query=%28%28gene%3A{gene_name}%29%29+AND+%28model_organism%3A9606%29+AND+%28reviewed%3Atrue%29&size=500", timeout=10)
    # return search_results
    if search_results.status_code == 200:
        data = gzip.decompress(search_results.content)
        text = data.decode('utf-8')
        site_meta_data = json.loads(text)
        
        results = []
        for result in site_meta_data['results']:
            if 'UniProtKB reviewed' in result['entryType']:
                print(result['entryType'])
                p_id = result['primaryAccession']
                protein_name = result['proteinDescription']['recommendedName']['fullName']['value']
                results.append((protein_name, p_id))
        print(results)
        return results

# %%
def get_pdbs_df(protein_id):
    # pdb_list = []
    search_results = requests.get(f"https://www.ebi.ac.uk/pdbe/graph-api/uniprot/protvista/unipdb/{protein_id}")
    data = search_results.json()
    # for pdb_struct in data[protein_id]['tracks'][0]['data']:
    def get_pdb(pdb_struct):
        output = {'type': 'PDB'}
        pdb_id = pdb_struct['label']['id']
        output['id'] = pdb_id
        output['resolution'] = pdb_struct['label']['resolution']
        if not output['resolution']:
            output['resolution'] = None
        start = pdb_struct['locations'][0]['fragments'][0]['start']
        end = pdb_struct['locations'][0]['fragments'][-1]['end']
        pattern = r'href="([^"]+)"'
        # Search for the pattern in the raw data
        match = re.search(pattern, pdb_struct['locations'][0]['fragments'][0]['tooltipContent'])

        # Extract the link if found
        if match:
            def get_chains(match):
                href_link = match.group(1)
                # print("Extracted href link:", href_link)
                search_results = requests.get(f"https://www.ebi.ac.uk/pdbe/graph-api/pdbe_pages/protvista/chains/{pdb_id}/{href_link.split('/')[-1]}")
                data = search_results.json()
                # print(data[pdb_id]['tracks'][0]['data'][0]['chainId'])
                chains = set()
                for data_ in data[pdb_id]['tracks'][0]['data']:
                    # print(data_['label'])
                    text = data_['label']
                    pattern = r"\(auth (.*)\)"
                    match = re.search(pattern, text)

                    if match:
                        auth_chain = match.group(1)
                        chains.add(auth_chain)
                chains_string = '/'.join(chains)+"="+str(start)+"-"+str(end)
                return chains_string
            fetched = False
            count = 0
            while not fetched and count<20:
                print(f"trying get chains for {pdb_id}: {count}")
                try:
                    output['chains'] = get_chains(match)
                    fetched = True
                except json.decoder.JSONDecodeError:
                    print("error occured in get chains")
                    count+=1
            if not fetched:
                raise Exception("get chains failed")
        else:
            output['chains'] = None
            print(f"No href link found for chains ({pdb_id})")
        # print(output)
        # pdb_list.append(output)
        return output

    pdb_list = run_on_multiple_threads(get_pdb, data[protein_id]['tracks'][0]['data'], max_workers=5)

    pdbs_df = pd.DataFrame(pdb_list)
    pdbs_df['id'] = pdbs_df['id'].str.upper()
    pdbs_df[['chain_type', 'chain_start', 'chain_end']] = pdbs_df['chains'].str.extract(r'(.+)?=(\d+)-(\d+)')
    pdbs_df['coverage'] = pdbs_df.apply(lambda x: int(x.chain_end)-int(x.chain_start) ,axis=1)
    pdbs_df['chain_type'] = pdbs_df['chain_type'].apply(lambda x: x.split('/'))
    pdbs_df['chain_start'] = pdbs_df['chain_start'].astype(int)
    pdbs_df['chain_end'] = pdbs_df['chain_end'].astype(int)
    pdbs_df = pdbs_df.dropna()
    pdbs_df = pdbs_df.set_index('id', drop=False)
    return pdbs_df

# %%
def get_protein_details(protein_id):
    search_results = requests.get(f"https://www.ebi.ac.uk/proteins/api/proteins/{protein_id}")
    data = search_results.json()
    length, seq = data['sequence']['length'], data['sequence']['sequence']
    return length, seq


# %%
import matplotlib.pyplot as plt
import numpy as np

def visualize_ranges(ranges, selected_ranges, MAX):
    plt.figure(figsize=(15, 5))
    for r in ranges:
        plt.plot([r.start/MAX*10, r.end/MAX*10], [1, 1], lw=6, alpha=0.2, label=f'Range({r.start}, {r.end}, {r.error_score})' if r not in selected_ranges else None)

    # Highlight selected ranges
    for r in selected_ranges:
        plt.plot([r.start/MAX*10, r.end/MAX*10], [1.05, 1.05], lw=6, color='red', label=f'Selected: Range({r.start}, {r.end}, {r.error_score})')
        # pass

    plt.yticks([])
    plt.title('Range Selection Visualization')
    plt.xlabel('Range')
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.xlim(0, 15)  # Adjust based on your range data
    plt.ylim(0, 2)
    plt.axhline(y=1, color='gray', lw=0.5)
    plt.savefig(f"visualize_ranges.png")


# %%
class Range:
    def __init__(self, id, start, end, error_score):
        self.id = id
        self.start = start
        self.end = end
        self.error_score = error_score
        self.coverage = end-start

    def __repr__(self):
        return f"Range({self.id}, {self.start}, {self.end}, {self.error_score}, coverage = {self.coverage})"

def select_ranges(ranges):
    # Sort ranges by starting point, then by end point
    sorted_ranges = sorted(ranges, key=lambda r: (r.start, r.end))
    print(sorted_ranges)
    selected_ranges = []
    last_end = float('-inf')

    for current_range in sorted_ranges:
        if current_range.start > last_end:
            # No overlap
            selected_ranges.append(current_range)
            last_end = current_range.end
        else:
            # Overlapping ranges
            if current_range.error_score < selected_ranges[-1].error_score or abs(current_range.error_score - selected_ranges[-1].error_score) < 0.5:
                if current_range.coverage > selected_ranges[-1].coverage:
                    selected_ranges.append(current_range)

                last_end = max(last_end, current_range.end)
    return selected_ranges

def get_optimal_pdbs_df(df, seq_length, min_res_val=2.5):
    df['range_obj'] = df.apply(lambda x: Range(x.id, x.chain_start, x.chain_end, x.resolution), axis=1)
    min_res_for_max_cov = min(df[df['coverage']==max(df['coverage'])]['resolution'].values)

    df_1 = df[df['resolution']<=min_res_val]
    df_2 = df[df['resolution']<=min_res_for_max_cov]

    all_data = df['range_obj'].to_list()
    data1 = df_1['range_obj'].to_list()
    data2 = df_2['range_obj'].to_list()
    print("set1 and set2 length for finding optimal res/coverage pdbs: ", len(data1) , len(data2))
    
    selected_ranges1 = select_ranges(data1)
    print("Selected high res Ranges:", selected_ranges1)

    selected_ranges2 = select_ranges(data2)
    print("Selected overall Ranges:", selected_ranges2)

    for i in selected_ranges1:
        print(i)

    for i in selected_ranges1:
        print(i.id,end=",")

    for i in selected_ranges2:
        print(i)

    for i in selected_ranges2:
        print(i.id,end=",")
    
    selected_ranges = list(set(selected_ranges1+selected_ranges2))

    # visualize_ranges(all_data, selected_ranges, seq_length)

    ids = [i.id for i in selected_ranges]
    selected_pdbs = df.loc[ids]
    return selected_pdbs

# %%
def download_pdbs(gene_name, pdb_id_list):
    failed_pdbs = []
    for pdb_id in pdb_id_list:
        if not os.path.isdir(f'{gene_name}'):
            os.mkdir(f'{gene_name}')
        if os.path.isfile(os.path.join(gene_name, f"g_{gene_name}_p_{pdb_id}.pdb")):
            print(os.path.join(gene_name, f"g_{gene_name}_p_{pdb_id}.pdb"), " exists!")
            continue
        res = os.system(f'cd {gene_name}; wget -O g_{gene_name}_p_{pdb_id}.pdb https://www.ebi.ac.uk/pdbe/entry-files/download/pdb{pdb_id.lower()}.ent')
        if res !=0:
            if os.path.isfile(os.path.join(gene_name, f"g_{gene_name}_p_{pdb_id}.pdb")):
                os.remove(os.path.join(gene_name, f"g_{gene_name}_p_{pdb_id}.pdb"))
            failed_pdbs.append(pdb_id)
    return failed_pdbs

# %%
from Bio import PDB

def get_chain_ids(pdb_file):
    # Create a PDB parser
    parser = PDB.PDBParser(QUIET=True)

    # Parse the PDB file
    structure = parser.get_structure('PDB_structure', pdb_file)

    # Get chain IDs
    chain_ids = set()  # Using a set to avoid duplicates
    for model in structure:
        for chain in model:
            chain_ids.add(chain.id)

    return list(chain_ids)

# %%
def pdb_cleaner(gene_name, id, remove_chains=[]):
    # remove_chains = []
    replace_nonstandard_residues = None
    # fix_add_missing_atoms = True
    remove_heterogens = True
    remove_water = True
    add_hydrogens = True
    pH = 7.0
    pdb_path = f"./{gene_name}/g_{gene_name}_p_{id}.pdb"

    try:
        from rdkit import Chem  # type: ignore
        from rdkit.Chem import AllChem  # type: ignore
    except ModuleNotFoundError:
        raise ImportError("This function requires RDKit to be installed.")
    try:
        from pdbfixer.pdbfixer import PDBFixer, proteinResidues, dnaResidues  # type: ignore
    except:  # noqa
        raise ImportError("This function requires pdbfixer to be installed.")
    try:
        from openmm.app import PDBFile  # type: ignore
    except:  # noqa
        raise ImportError("This function requires OpenMM to be installed.")

    if remove_chains == 'None':
        remove_chains = None

    if True:
        fixer = PDBFixer(pdb_path)
        output_path = f"./{gene_name}/cleaned_g_{gene_name}_p_{id}.pdb"

        try:
            if remove_chains:
                fixer.removeChains(chainIds=remove_chains)
            if replace_nonstandard_residues:
                fixer.findMissingResidues()
                fixer.findNonstandardResidues()
                fixer.replaceNonstandardResidues()
                fixer.findMissingAtoms()
                fixer.addMissingAtoms()
            if remove_heterogens and not remove_water:
                fixer.removeHeterogens(True)
            if remove_heterogens and remove_water:
                fixer.removeHeterogens(False)
            # if not replace_nonstandard_residues and fix_add_missing_atoms:
            #     fixer.findMissingResidues()
            #     fixer.findMissingAtoms()
            #     fixer.addMissingAtoms()
            if add_hydrogens:
                try:
                    fixer.addMissingHydrogens(pH)
                except:
                    fixer.findMissingResidues()
                    fixer.findMissingAtoms()
                    fixer.addMissingAtoms()
                    fixer.addMissingHydrogens(pH)

            PDBFile.writeFile(fixer.topology, fixer.positions,
                                open(output_path, 'w'))
            return output_path
        except Exception as e:  # noqa
            print(f'failed PDB fixing: {pdb_path} => {e}')
            return None


def get_cleaned_pdbs(gene_name, entry_id):
    protein_id = entry_id

    MAX_COVERAGE, prot_seq = get_protein_details(protein_id)
    print(f"Fetching pdbs for canonical protein {protein_id}")
    pdbs_df = get_pdbs_df(protein_id)

    # selected_pdbs = get_optimal_pdbs_df(pdbs_df, MAX_COVERAGE, min_res_val=2.5)
    # id_chain_map = selected_pdbs['chain_type'].to_dict()

    # failed_pdbs = download_pdbs(gene_name, list(id_chain_map.keys()))
    failed_pdbs = None
    while failed_pdbs != []:
        print("retrying get_optimal_pdbs")
        if failed_pdbs is not None:
            pdbs_df = pdbs_df.drop(index=failed_pdbs)
        if len(pdbs_df) == 0:
            raise OverflowError("no pdbs left to check, likely over sized pdbs")
        selected_pdbs = get_optimal_pdbs_df(pdbs_df, MAX_COVERAGE, min_res_val=2.5)
        id_chain_map = selected_pdbs['chain_type'].to_dict()
        failed_pdbs = download_pdbs(gene_name, list(id_chain_map.keys()))
        if failed_pdbs != []:
            print(f"failed_pdbs: {failed_pdbs}")
            continue

        for id, chains in id_chain_map.items():
            all_chains = get_chain_ids(f"./{gene_name}/g_{gene_name}_p_{id}.pdb")
            chains_to_remove = list(set(all_chains).difference(set(chains)))
            id_chain_map[id] = id_chain_map[id], chains_to_remove
        
        cleaned_pdbs_path = []
        uncleanable_pdbs = []
        for id, chains in id_chain_map.items():
            out = pdb_cleaner(gene_name, id=id, remove_chains=chains[-1])
            cleaned_pdbs_path.append(out)
            if out is None:
                uncleanable_pdbs.append(id)
        print(uncleanable_pdbs)
        failed_pdbs = uncleanable_pdbs.copy()
        print(f"failed_pdbs: {failed_pdbs}")
    
    print("pdbs download successful")

    if uncleanable_pdbs:
        print(f"uncleanable_pdbs: {uncleanable_pdbs}")

    if not cleaned_pdbs_path:
        raise Exception("No clean pdbs")

    selected_pdbs['path'] = cleaned_pdbs_path
    selected_pdbs.to_csv(f"./{gene_name}/{gene_name}_pdbs.csv")
    return selected_pdbs


def run(item):
    entry_id = item['Entry']
    primary_genes = item['Gene Names (primary; single)']
    gene_name = primary_genes.split('; ')[0].strip()
    print(gene_name)
    if os.path.exists(f"/home/ubuntu/download_pdbs/{gene_name}"):
        print(f"skipping {gene_name}")
        return 'skipped'
    try:
        output = get_cleaned_pdbs(gene_name, entry_id)
    except OverflowError as e:
        print(e)
        return None
    return output

if __name__ == "__main__":
    df = pd.read_csv("ProteomeScan/gene_selection/experimental/valid_for_pipeline_7681_human_prot.csv")
    data = df[['Entry', 'Gene Names (primary; single)']].to_dict('index')

    to_run = list(data.values())

    print(len(to_run))

    from concurrent.futures import ProcessPoolExecutor
    from tqdm import tqdm

    with ProcessPoolExecutor(max_workers=5) as executor:
        results = list(tqdm(executor.map(run, to_run), desc="Parallel Processing", total=len(to_run)))
