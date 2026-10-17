"""Score every single substitution at pocket positions with ESM-2 (masked marginal: log p(mut) - log p(wt)).

Positive score = ESM finds the mutant more plausible than wild type. It is a prior, not a kinetics prediction.
"""
import pandas as pd
import torch
import esm

AAS = "ACDEFGHIKLMNPQRSTVWY"
# catalytic triad is excluded from mutagenesis
PROTECTED = {164, 303, 336}

seq = "".join(l.strip() for l in open("data/1CGK_A.fasta") if not l.startswith(">"))
pocket = pd.read_csv("data/pocket.csv")
pocket = pocket[~pocket.resseq.isin(PROTECTED)]

model, alphabet = esm.pretrained.esm2_t12_35M_UR50D()
model.eval()
conv = alphabet.get_batch_converter()
_, _, tokens = conv([("wt", seq)])

rows = []
for _, p in pocket.iterrows():
    i = int(p.seq_idx)                      # 1-based; token index is i because of BOS at 0
    masked = tokens.clone()
    masked[0, i] = alphabet.mask_idx
    with torch.no_grad():
        logp = torch.log_softmax(model(masked)["logits"][0, i], dim=-1)
    wt = seq[i - 1]
    for aa in AAS:
        if aa != wt:
            rows.append({"resseq": int(p.resseq), "wt": wt, "mut": aa,
                         "variant": f"{wt}{int(p.resseq)}{aa}",
                         "esm_score": float(logp[alphabet.get_idx(aa)] - logp[alphabet.get_idx(wt)])})

df = pd.DataFrame(rows).sort_values("esm_score", ascending=False)
df.to_csv("data/esm_scores.csv", index=False)
print(df.head(15).to_string(index=False))
