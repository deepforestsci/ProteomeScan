import requests
import gzip
import json
from multiprocessing import Pool
from tqdm import tqdm


def check_kinase(gene_name):
    url = (
        "https://rest.uniprot.org/uniprotkb/search"
        f"?compressed=true&format=json"
        f"&query=((gene:{gene_name}) AND (organism_id:9606)"
        f" AND (keyword:KW-0418)) AND (reviewed:true)&size=500"
    )
    r = requests.get(url)
    r.raise_for_status()
    decompressed = gzip.decompress(r.content)
    data = json.loads(decompressed.decode("utf-8"))
    return gene_name, bool(data.get("results"))

def worker(genes):
    return [check_kinase(g) for g in genes]

def run_parallel(gene_list, processes=None):
    # Split genes into roughly equal chunks per process
    n = processes or None  # None → default to cpu_count()
    with Pool(processes=n) as pool:
        results = list(tqdm(pool.imap_unordered(check_kinase, gene_list),
                            total=len(gene_list),
                            desc="Checking kinases"))
    return dict(results)

if __name__ == "__main__":
    gene_list = ['MAP2K1', 'TOP1']  # example
    results = run_parallel(gene_list)