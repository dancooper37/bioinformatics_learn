from toolkit.utilities import unpackFASTAToDict
from toolkit.dna import consensusNuc

dnaMatrix = []

dataDict = unpackFASTAToDict("../compute_data/rosalind_cons.txt")

for value in dataDict.values():
    dnaMatrix.append(value)  # Unpacks the sequences from fasta into a list of strings

dnaVector = []  # Stores the nth base from each sequence, cleared after use
consensus = []  # List of most frequent base in each position over all sequences

for base in range(len(dnaMatrix[0])):
    for seq in range(len(dnaMatrix)):
        dnaVector.append(dnaMatrix[seq][base])
        if seq == len(dnaMatrix) - 1:  # Once the vector has been fully generated for each of the sequences
            consensus.append(consensusNuc(dnaVector))  # Find the most common base from the vector
            dnaVector.clear()  # Empty the list

print("".join(consensus))













