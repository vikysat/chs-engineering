"""Download 1CGK, list residues within CUTOFF angstrom of Cys164 SG, and write chain A sequence."""
import pandas as pd
from Bio.PDB import PDBList, PDBParser, NeighborSearch, is_aa
from Bio.SeqUtils import seq1

CUTOFF = 8.0
CATALYTIC_RESSEQ = 164

PDBList().retrieve_pdb_file("1CGK", pdir="data", file_format="pdb")
structure = PDBParser(QUIET=True).get_structure("1CGK", "data/pdb1cgk.ent")
chain = structure[0]["A"]

residues = [r for r in chain if is_aa(r, standard=True)]
seq = "".join(seq1(r.get_resname()) for r in residues)
seq_idx = {r.get_id()[1]: i + 1 for i, r in enumerate(residues)}  # PDB resseq -> 1-based seq index

cys = chain[CATALYTIC_RESSEQ]
assert cys.get_resname() == "CYS", f"Residue {CATALYTIC_RESSEQ} is {cys.get_resname()}, not CYS"
sg = cys["SG"]

atoms = [a for r in residues for a in r]
near = {a.get_parent() for a in NeighborSearch(atoms).search(sg.coord, CUTOFF)}
rows = [{"resseq": r.get_id()[1], "seq_idx": seq_idx[r.get_id()[1]], "aa": seq1(r.get_resname()),
         "dist_SG": round(min(a - sg for a in r), 2)} for r in near]
df = pd.DataFrame(rows).sort_values("dist_SG")
df.to_csv("data/pocket.csv", index=False)

with open("data/1CGK_A.fasta", "w") as f:
    f.write(f">1CGK_A\n{seq}\n")
print(f"{len(df)} residues within {CUTOFF} A of Cys{CATALYTIC_RESSEQ}; wrote data/pocket.csv")
