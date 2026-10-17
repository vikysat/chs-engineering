"""Superimpose a predicted model on 1CGK and compare pocket residues and confidence.
Usage: python src/07_compare_models.py data/models/Q8H4L3_model.pdb Q8H4L3
"""
import sys
import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import PDBParser, PDBIO, Superimposer, is_aa
from Bio.SeqUtils import seq1

path, name = sys.argv[1], sys.argv[2]

def residues(p):
    return [r for r in PDBParser(QUIET=True).get_structure("x", p)[0]["A"] if is_aa(r, standard=True)]

ref, mod = residues("data/pdb1cgk.ent"), residues(path)
sr = "".join(seq1(r.get_resname()) for r in ref)
sm = "".join(seq1(r.get_resname()) for r in mod)

al = Align.PairwiseAligner()
al.substitution_matrix = substitution_matrices.load("BLOSUM62")
al.open_gap_score, al.extend_gap_score = -10, -0.5
a = al.align(sr, sm)[0]
m = {t0 + k: q0 + k for (t0, t1), (q0, q1) in zip(*a.aligned) for k in range(t1 - t0)}

sup = Superimposer()
sup.set_atoms([ref[i]["CA"] for i in m], [mod[j]["CA"] for j in m.values()])
sup.apply([at for r in mod for at in r])
ident = sum(sr[i] == sm[j] for i, j in m.items()) / len(m)

pocket = pd.read_csv("data/pocket.csv", keep_default_na=False)
rows, pl = [], []
for _, p in pocket.iterrows():
    i = int(p.seq_idx) - 1
    if i in m:
        j = m[i]
        pl.append(mod[j]["CA"].get_bfactor())
        rows.append((p.resseq, sr[i], mod[j].get_id()[1], sm[j], round(ref[i]["CA"] - mod[j]["CA"], 2),
                     "" if sr[i] == sm[j] else "DIFFERENT"))
df = pd.DataFrame(rows, columns=["1CGK_res", "1CGK_aa", f"{name}_res", f"{name}_aa", "CA_dist", "note"])
df.to_csv(f"data/pocket_compare_{name}.csv", index=False)

io = PDBIO(); io.set_structure(mod[0].get_parent().get_parent()); io.save(f"data/models/{name}_aligned.pdb")

print(f"{name}: {len(m)} residues aligned, identity {ident:.0%}, CA RMSD {sup.rms:.2f} A")
print(f"pLDDT whole chain {sum(r['CA'].get_bfactor() for r in mod)/len(mod):.1f}, pocket {sum(pl)/len(pl):.1f}")
print(f"Max pocket CA distance {df.CA_dist.max()} A")
print(df[df.note != ""].to_string(index=False))
