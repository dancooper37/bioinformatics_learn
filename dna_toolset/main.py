from DNAToolkit import *
import random

# Creates random DNA sequence for testing
randDNAStr = "".join([random.choice(Nucleotides) for nuc in range(50)])

DNAStr = validateSeq(randDNAStr)
print(countNucFrequency(DNAStr))