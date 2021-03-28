def readFile(filePath):
    """Reads file and returns list of lines"""
    with open(filePath, "r") as f:
        return [l.strip() for l in f.readlines()]


def gc_content(seq):
    """GC content in DNA/RNA sequence"""
    return (seq.count("C") + seq.count("G")) / len(seq) * 100


# Unpacks file contents into list
# FASTAFile = readFile("../test_data/gc.txt.txt")
FASTAFile = readFile("../compute_data/rosalind_gc.txt")
# Dictionary for labels + data
FASTADict = {}
# String for current label
FASTALabel = ""

# Converts FASTA list file data into dictionary
for line in FASTAFile:
    if ">" in line:
        FASTALabel = line
        FASTADict[FASTALabel] = ""
    else:
        FASTADict[FASTALabel] += line

# Uses dictionary comprehension to generate new dictionary with GC content values
RESULTDict = {key: gc_content(value) for (key, value) in FASTADict.items()}

# Finds max value in values() of results dictionary
MaxGCKey = max(RESULTDict, key=RESULTDict.get)

print(f"{MaxGCKey[1:]}\n{RESULTDict[MaxGCKey]}")
