from collections import Counter
from toolkit.structures import *


def validate(dna_seq):
    """Checks if string is valid DNA sequence"""
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotides:
            return False
    return tmpseq


def countNucFrequency(seq):
    """Counts nucleotide occurence in a sequence"""
    tmpFreqDict = {"A": 0, "C": 0, "G": 0, "T": 0}
    for nuc in seq:
        tmpFreqDict[nuc] += 1
    return tmpFreqDict
    # More pythonic approach
    # return dict(collections.Counter(seq))


def consensusNuc(vector):
    """Finds the most common nucleotide from a vector/list"""
    tmpFreqDict = {"A": 0, "C": 0, "G": 0, "T": 0}
    for nuc in vector:
        tmpFreqDict[nuc] += 1
    return max(tmpFreqDict, key=tmpFreqDict.get)


def transcribeToRNA(seq):
    """DNA to RNA translation. Replaces Thymine with Uracil"""
    return seq.replace("T", "U")


def reverseCompliment(seq):
    """Finds complimentary base pairs and reverses list"""
    return "".join([DNA_ReverseCompliment[nuc] for nuc in seq])[::-1]  # Reverses the list
    # More pythonic approach
    # mapping = str.maketrans("ATCG", "TAGC")
    # return seq.translate(mapping)[::-1]


def percentGC(seq):
    """GC content in DNA/RNA sequence"""
    return round((seq.count("C") + seq.count("G")) / len(seq) * 100)


def percentATU(seq):
    """AT/U content in DNA/RNA sequence"""
    return round((seq.count("A") + seq.count("T") + seq.count("U")) / len(seq) * 100)


def subseqGCPercent(seq, k=20):
    """GC Content in a DNA/RNA sub-sequence length k. k = 20 by default"""
    res = []
    for i in range(0, len(seq) - k + 1, k):
        subseq = seq[i:i + k]
        res.append(percentGC(subseq))
    return res


def subseqATUPercent(seq, k=20):
    """ATU Content in a DNA/RNA sub-sequence length k. k = 20 by default"""
    res = []
    for i in range(0, len(seq) - k + 1, k):
        subseq = seq[i:i + k]
        res.append(percentATU(subseq))
    return res


def translateDNASeq(seq, init_pos=0):
    """Translates DNA sequence into amin acid sequence"""
    return [DNA_Codons[seq[pos:pos + 3]] for pos in range(init_pos, len(seq) - 2, 3)]


def translateRNASeq(seq, init_pos=0):
    """Translates DNA sequence into amin acid sequence"""
    return [RNA_Codons[seq[pos:pos + 3]] for pos in range(init_pos, len(seq) - 2, 3)]


def codonFrequency(seq, aminoacid):
    """Provides frequency of each codon encoding a given amino acid in DNA sequence"""
    tmpList = []
    for i in range(0, len(seq) - 2, 3):
        if DNA_Codons[seq[i:i + 3]] == aminoacid:
            tmpList.append(seq[i:i + 3])

    freqDict = dict(Counter(tmpList))
    totalWeight = sum(freqDict.values())
    for seq in freqDict:
        freqDict[seq] = round(freqDict[seq] / totalWeight, 2)
    return freqDict


def genReadingFrames(seq):
    """Generate the six reading frames of a DNA sequence, including the reverse compliment"""
    frames = [
        translateDNASeq(seq, 0),
        translateDNASeq(seq, 1),
        translateDNASeq(seq, 2),
        translateDNASeq(reverseCompliment(seq), 0),
        translateDNASeq(reverseCompliment(seq), 1),
        translateDNASeq(reverseCompliment(seq), 2)
    ]
    return frames


def proteinListFromRF(aaSeq):
    """Finds all possible proteins from an AA sequence in reading frame. Returns list of possible proteins"""
    currentProtein = []
    proteins = []
    for aa in aaSeq:
        # Stop protein list if _ marker found
        if aa[:1] == "_":  # index used to align with structures.py format"
            if currentProtein:
                for p in currentProtein:
                    proteins.append(p)
                currentProtein = []
        else:
            # Start protein list if methionine found
            if aa[:1] == "M":
                currentProtein.append("")
            for i in range(len(currentProtein)):
                currentProtein[i] += aa
    return proteins


def proteinScan(seq, startReadPos = 0, endReadPos = 0, ordered = False):
    """Find all possible unique proteins for all open reading frames"""
    if endReadPos > startReadPos:
        rfs = genReadingFrames(seq[startReadPos:endReadPos])
    else:
        rfs = genReadingFrames(seq)

    res = []
    for rf in rfs:
        proteins = proteinListFromRF(rf)
        for p in proteins:
            res.append(p)

    res = list(dict.fromkeys(res))  # Removes duplicates

    if ordered:
        return sorted(res, key = len, reverse = True)

    return res



