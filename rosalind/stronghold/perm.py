from math import factorial
from itertools import permutations


n = 7
nList = []

for i in range(n):
    nList.append(i + 1)

print(factorial(n))

permGen = permutations(nList)

for permList in list(permGen):
    print(" ".join(map(str, permList)))

