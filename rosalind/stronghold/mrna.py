from collections import Counter
from toolkit.structures import RNA_Codons
from toolkit.utilities import unpackFASTAToStr

protein = unpackFASTAToStr("../compute_data/rosalind_mrna.txt")

permutations = 1

for base in protein:
    permutations *= Counter(RNA_Codons.values())[base]

permutations *= 3  # To account for stop codon

print(permutations % 1000000)
