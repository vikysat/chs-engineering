"""Make figures: library size vs screening burden, ESM variant heatmap, phylogenetic tree."""
import importlib.util, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

spec = importlib.util.spec_from_file_location("nnk", "src/05_nnk.py")
nnk = importlib.util.module_from_spec(spec); spec.loader.exec_module(nnk)
os.makedirs("figures", exist_ok=True)

# 1. library size vs clones to screen
sites = [1, 2, 3, 4]
clones = [nnk.clones_for_coverage(nnk.library_size(n)) for n in sites]
fig, ax = plt.subplots(figsize=(5, 3.5))
ax.bar([str(n) for n in sites], clones, color="#3b6ea5")
ax.set_yscale("log"); ax.set_xlabel("NNK sites randomised"); ax.set_ylabel("Clones for 95% coverage")
ax.axhline(10_000, ls="--", color="grey"); ax.text(-0.4, 12_000, "~10k clones: typical plate-based screen", fontsize=7)
ax.set_title("Screening burden vs library size"); fig.tight_layout()
fig.savefig("figures/nnk_library.png", dpi=200)

# 2. ESM heatmap (needs data/esm_scores.csv from step 04)
if os.path.exists("data/esm_scores.csv"):
    import pandas as pd
    df = pd.read_csv("data/esm_scores.csv", keep_default_na=False)
    df["pos"] = df.wt + df.resseq.astype(str)
    mat = df.pivot(index="pos", columns="mut", values="esm_score")
    mat = mat.loc[sorted(mat.index, key=lambda s: int(s[1:]))]
    fig, ax = plt.subplots(figsize=(7, max(4, 0.22 * len(mat))))
    im = ax.imshow(mat.values, cmap="RdBu", vmin=-6, vmax=6, aspect="auto")
    ax.set_xticks(range(len(mat.columns)), mat.columns); ax.set_yticks(range(len(mat.index)), mat.index, fontsize=6)
    fig.colorbar(im, label="ESM-2 log-odds (mut vs wt)"); ax.set_title("Pocket single-mutant scan")
    fig.tight_layout(); fig.savefig("figures/esm_heatmap.png", dpi=200)
else:
    print("skip heatmap: run step 04 first")

# 3. tree (needs IQ-TREE output)
tree_file = "data/chs_triad.fasta.treefile"
if os.path.exists(tree_file):
    from Bio import Phylo
    fig = plt.figure(figsize=(8, 10)); ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(Phylo.read(tree_file, "newick"), axes=ax, do_show=False, label_func=lambda c: None)
    fig.savefig("figures/chs_tree.png", dpi=200)
else:
    print("skip tree: run IQ-TREE first")
