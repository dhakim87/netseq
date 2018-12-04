from candidateAssignment import *
from Composition import *

#Test Generate Candidates Single V3
testStart = [DEFAULT_ALPHABET.sortedAminos[0], DEFAULT_ALPHABET.sortedAminos[1]]
transitionPaths = Composition.generatePaths(Composition.fromAminoArr(testStart), 114.0428, 2, 0.02)
testNeighbors = generateCandidatesSingle3(testStart, 2, transitionPaths)
print("||L|| = " + str(len(DEFAULT_ALPHABET.sortedAminos)))
print(testNeighbors)
print(len(testNeighbors))
print(transitionPaths)


#Test Calc Mass
testStart = [DEFAULT_ALPHABET.sortedAminos[0], DEFAULT_ALPHABET.sortedAminos[1]]
print("Calc Mass: " + str(calcMass(testStart)))

#Test Generate Candidates Single
testStart = [DEFAULT_ALPHABET.sortedAminos[0], DEFAULT_ALPHABET.sortedAminos[1]]
testNeighbors = generateCandidatesSingle2(testStart, 2);
print("||L|| = " + str(len(DEFAULT_ALPHABET.sortedAminos)))
print(testNeighbors)
print(len(testNeighbors))

#Test Generate Candidates
testStart = [
[DEFAULT_ALPHABET.sortedAminos[0], DEFAULT_ALPHABET.sortedAminos[1]],
[DEFAULT_ALPHABET.sortedAminos[2], DEFAULT_ALPHABET.sortedAminos[3]]
]
testNeighbors = generateCandidates(0, 0, testStart, 2, .02);
print("||L|| = " + str(len(DEFAULT_ALPHABET.sortedAminos)))
print(testNeighbors)
print(len(testNeighbors))
