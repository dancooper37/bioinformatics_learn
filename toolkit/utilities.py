def colored(seq):
    """Colors nucleotides in terminal for easier reading"""
    bcolors = {
        "A": "\033[92m",
        "C": "\033[94m",
        "G": "\033[93m",
        "T": "\033[91m",
        "U": "\033[91m",
        "reset": "\033[0;0m"
    }

    tmpStr = ""

    for nuc in seq:
        if nuc in bcolors:
            tmpStr += bcolors[nuc] + nuc
        else:
            tmpStr += bcolors["reset"] + nuc

    return tmpStr + "\033[0;0m"


def readFile(filePath):
    """Reads file and returns list of lines"""
    with open(filePath, "r") as f:
        return [l.strip() for l in f.readlines()]


def unpackFASTAToDict(url):
    """Unpacks the contents of a FASTA file into a dictionary"""
    # Unpacks file contents into list
    FASTAFile = readFile(url)
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

    return FASTADict


def unpackFASTAToStr(url):
    """Unpacks the contents of a FASTA file into a single string"""
    # Unpacks file contents into list
    FASTAFile = readFile(url)
    FASTAStr = ""

    # Converts FASTA list file data into dictionary
    for line in FASTAFile:
        if ">" in line:
            FASTALabel = line
        else:
            FASTAStr += line

    return FASTAStr
