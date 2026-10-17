"""Fetch CHS-family sequences from UniProt plus the reference (alfalfa CHS2, P30074)."""
import requests

URL = "https://rest.uniprot.org/uniprotkb/search"
QUERY = '(protein_name:"chalcone synthase") AND (length:[350 TO 450]) AND (reviewed:true)'

def get_fasta(query, size=200):
    r = requests.get(URL, params={"query": query, "format": "fasta", "size": size}, timeout=60)
    r.raise_for_status()
    return r.text

if __name__ == "__main__":
    seqs = get_fasta(QUERY)
    ref = get_fasta("accession:P30074", size=1)
    with open("data/chs_all.fasta", "w") as f:
        # reference first so it is easy to find in the alignment
        f.write(ref)
        f.write(seqs)
    print("Wrote data/chs_all.fasta:", (ref + seqs).count(">"), "sequences")
