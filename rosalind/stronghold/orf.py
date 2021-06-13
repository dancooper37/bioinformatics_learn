from toolkit.dna import proteinScan
from toolkit.utilities import unpackFASTAToStr

DNAStr = unpackFASTAToStr("../compute_data/rosalind_orf.txt")

for protein in proteinScan(DNAStr, 0, 0, True):
    print(f"{protein}")