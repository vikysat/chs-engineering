"""NNK library size and clones needed for a target coverage."""
import math

NNK_CODONS = 32  # 32 codons encode all 20 amino acids + 1 stop (TAG)
NNK_AA = 20

def library_size(n_sites, per_site=NNK_CODONS):
    return per_site ** n_sites

def clones_for_coverage(n_variants, coverage=0.95):
    """Clones to screen so each variant is seen with probability >= coverage (uniform sampling)."""
    return math.ceil(math.log(1 - coverage) / math.log(1 - 1 / n_variants))

if __name__ == "__main__":
    print("sites | codon combos | clones @95% | amino-acid combos | clones @95%")
    for n in range(1, 5):
        c, a = library_size(n), library_size(n, NNK_AA)
        print(f"{n:5d} | {c:12,d} | {clones_for_coverage(c):11,d} | {a:17,d} | {clones_for_coverage(a):11,d}")
