"""Keep only sequences that have the catalytic triad (C164, H303, N336 in alfalfa CHS2 numbering).

Run first:  mafft --auto data/chs_all.fasta > data/chs_aln.fasta
"""
from Bio import AlignIO

TRIAD = {164: "C", 303: "H", 336: "N"}  # 1-based ungapped positions in the reference
REF_ID = "P30074"

aln = AlignIO.read("data/chs_aln.fasta", "fasta")
ref = next(r for r in aln if REF_ID in r.id)
ref_seq = str(ref.seq)

# map ungapped reference position -> alignment column
col_of = {}
pos = 0
for col, aa in enumerate(ref_seq):
    if aa != "-":
        pos += 1
        col_of[pos] = col

for p, aa in TRIAD.items():
    got = ref_seq[col_of[p]]
    if got != aa:
        print(f"WARNING: reference has {got} at {p}, expected {aa}. Check numbering (e.g. Met1 offset).")

kept = [r for r in aln if all(str(r.seq)[col_of[p]] == aa for p, aa in TRIAD.items())]
print(f"{len(kept)}/{len(aln)} sequences keep the triad")

seen = set()
with open("data/chs_triad.fasta", "w") as f:
    for r in kept:
        acc = r.id.split("|")[1] if "|" in r.id else r.id   # e.g. P30074
        if acc in seen:
            continue
        seen.add(acc)
        f.write(f">{acc}\n{str(r.seq)}\n")   # keep gaps: IQ-TREE needs an alignment
