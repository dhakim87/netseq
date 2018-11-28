from candidateAssignment import *


#Test Calc Mass
testStart = [DEFAULT_ALPHABET.sortedAminos[0], DEFAULT_ALPHABET.sortedAminos[1]]
print("Calc Mass: " + str(calcMass(testStart)))

#Test Generate Candidates Single
testStart = [DEFAULT_ALPHABET.sortedAminos[0], DEFAULT_ALPHABET.sortedAminos[1]]
testNeighbors = generateCandidatesSingle(testStart, 2);
print("||L|| = " + str(len(DEFAULT_ALPHABET.sortedAminos)))
print(testNeighbors)
print(len(testNeighbors))

#Test Generate Candidates
testStart = [
[DEFAULT_ALPHABET.sortedAminos[0], DEFAULT_ALPHABET.sortedAminos[1]],
[DEFAULT_ALPHABET.sortedAminos[2], DEFAULT_ALPHABET.sortedAminos[3]]
]
testNeighbors = generateCandidates(0, 0, testStart, 2, 10000);
print("||L|| = " + str(len(DEFAULT_ALPHABET.sortedAminos)))
print(testNeighbors)
print(len(testNeighbors))
