import Queue
from eq import fuzzyEqual
from globals import DEFAULT_ALPHABET
import scoring
import parsers
from SpectralGraph import *
from Composition import *
from collections import defaultdict

#There are several places where we will branch and bound and keep only the M highest candidates by score.  You can control M here.
M = 1000

#General Idea:  We will go with option 1 as it makes more sense as a module within the system.
#Option 1.  Create module that assigns candidates to all nodes in the network:
#    assignCandidates(SG, initialCandidates)
#    {
#        annotatedNodes = initialCandidates
#        while "more nodes to annotate"
#            generate all possible candidates at unannotated nodes
#            score all candidates individually by combination of pairwise score (MST score?) with known candidates and score on associated spectra
#            Sort all candidates by score, take top m,
#            Take best candidate by score and set node to annotate as node containing this candidate
#            Assign all candidates from top m associated with this node as candidates of this node
#
#        return annotatedNodes
#    }
#
#Option 2.  Create module that creates a single candidate network:
#    assignCandidates(SG, initialCandidates)
#    {
#        annotatedNodes = initialCandidates
#        while "more nodes to annotate"
#            generate all possible candidates at unannotated nodes
#            score all candidates individually by combination of pairwise score (MST score?) with known annotations and score on associated spectra
#            Take best candidate by score, add (candidate, nodeContainingCandidate) to annotatedNodes
#
#        Print network chosen and see how it lines up with our existing scoring mechanisms.
#
#        //How does this handle multiple candidates in the initial input?
#    }

"""Calculate the mass of a list of Aminos by summing up their individual masses"""
def calcMass(candidate):
    totalMass = 0
    for amino in candidate:
        totalMass += amino.mass
    return totalMass

"""Generate candidates, where each candidate is a list of aminos, from a starting list of aminos.
   Candidates are generated as all lists within a cyclic edit distance of maximumLDist"""
def generateCandidatesSingle2(neighbor, maximumLDist, alphabet = DEFAULT_ALPHABET):
    #Ugh!
    #Can formulate as a graph search in an absurdly high dimensional space.  Ugh!
    #Use breadth first search.  Store everything in spin invariant canonical.
    nodeToDist = {}
    nodeToVal = {}
    nodeToVisited = {}
    
    val = alphabet.toSpinInvariant(neighbor)
    strVal = str(val)
    nodeToDist[strVal] = 0
    nodeToVal[strVal] = val
    
    """helper function to add neighbors to the visit queue and update necessary mappings"""
    def evaluateNeighbor(newVal, newDist, nodeToDist, nodeToVal, toVisit):
        newVal = alphabet.toSpinInvariant(newVal)
        newStr = str(newVal)
        if newStr in nodeToDist:
            nodeToDist[newStr] = min(newDist, nodeToDist[newStr])
        else:
            nodeToDist[newStr] = newDist
        nodeToVal[newStr] = newVal
        toVisit.put(newStr)
    
    toVisit = Queue.Queue()
    toVisit.put(strVal)
    #TODO FIXME HACK:  Rather than organizing this as a breadth first search across individual transitions, if we can generate a list of multi-step paths of ADD/SWP/DEL operations
    #that add to the necessary change in mass, we can rewrite to loop over the paths rather than doing an exhaustive search of paths, many of which will end in invalid masses.
    while not toVisit.empty():
        curStr = toVisit.get()
        if curStr in nodeToVisited:
            continue
        nodeToVisited[curStr] = True
        
        curVal = nodeToVal[curStr]
        curDist = nodeToDist[curStr]
        if curDist == maximumLDist:
            continue #Neighbors will be out of bounds of this search, skip

        #Apply all possible transitions on this node.
        #DEL
        for position in range(0, len(curVal)):
            newVal = list(curVal)
            del newVal[position]
            evaluateNeighbor(newVal, curDist + 1, nodeToDist, nodeToVal, toVisit)
        #SWP
        for position in range(0, len(curVal)):
            for swapTo in alphabet.sortedAminos:
                newVal = list(curVal)
                newVal[position] = swapTo
                evaluateNeighbor(newVal, curDist + 1, nodeToDist, nodeToVal, toVisit)
        #ADD
        for position in range(0, max(len(curVal), 1)):
            #Note that inserting at position 0 and adding to end of list are the same because these are cyclopeptides.  Note also that both length 0 and length 1 lists have a single insertion point.
            for toAdd in alphabet.sortedAminos:
                newVal = list(curVal)
                newVal.insert(position, toAdd)
                evaluateNeighbor(newVal, curDist + 1, nodeToDist, nodeToVal, toVisit)

    resultsList = []
    for strVal in nodeToVal:
        resultsList.append(nodeToVal[strVal])

    print(" Neighbors of: " + str(neighbor) + " " + str(len(resultsList)))
    return resultsList


def generateCandidatesSingle3(neighbor, maximumLDist, transitionPaths, alphabet = DEFAULT_ALPHABET):
    val = alphabet.toSpinInvariant(neighbor)
    
    outputList = []
    for transitionPath in transitionPaths:
        activeList = [val]
        for transition in transitionPath:
            newActiveList = []
            for active in activeList:
                if transition.action == 0: #ADD
                    for position in range(0, max(len(active), 1)):
                        newActive = list(active)
                        newActive.insert(position, transition.amino1)
                        newActiveList.append(newActive)
                if transition.action == 1: #DEL
                    for position in range(0, len(active)):
                        if active[position] != transition.amino1:
                            continue
                        newActive = list(active)
                        del newActive[position]
                        newActiveList.append(newActive)
                if transition.action == 2: #SWP
                    for position in range(0, len(active)):
                        if active[position] != transition.amino1:
                            continue
                        newActive = list(active)
                        newActive[position] = transition.amino2
                        newActiveList.append(newActive)
            #Update the active list for each transition in the path
            activeList = alphabet.removeStrDups(newActiveList)
        #After running all transitions in the transition path, save the candidates to the output
        for candidate in activeList:
            outputList.append(alphabet.toSpinInvariant(candidate))

    #Finally, after running all transition paths, clean up the output list and return it
    outputList = alphabet.removeStrDups(outputList)
    return outputList;


"""Generate candidates, where each candidate is a list of aminos, from a starting set of lists of aminos with the same total approximate mass.
   Candidates are generated as all lists within a cyclic edit distance of maximumLDist, then are filtered down to be within epsilon of the target mass.  """
def generateCandidates(targetMass, neighborMass, neighborCandidates, maximumLDist, massEpsilon, alphabet = DEFAULT_ALPHABET):
    print("Generate Candidates: " + str(neighborMass) + " -> " + str(targetMass))
    #Ugh, there are a lot of ways we could generate candidates depending on the assumptions we're willing to make.  The problem seems tractable with a fixed alphabet, working at the level of masses directly we could generate a heck of a lot of potential candidates and it isn't clear to me which of those make biological sense.  So we'll start with a fixed alphabet for now and keep this function separated so we can try different implementations.

    #Option 1:  Absolute worst possible correct algorithm:  Enumerate all strings between length min neighborCandidate length - maximumLDist and length max neighborCandidate length + maximumLDist.  Check all strings for Ldist to all neighbor candidates and proximity to target mass.
    #Ugh, this is basically a password cracker.  Depending on length of alphabet, this could take quite a while.  Could potentially speed up the Ldist check based on some tree.

    #Option 2:  Second worst possible correct algorithm:  For each neighbor candidate, generate all strings within L dist, filter out all strings with the wrong target mass.
    #The generation step is tough, and will potentially generate the same string multiple times.  This should still be significantly faster than checking a string for its L dist to each neighbor candidate.

    #Option 3:  Potentially better algorithm:  Enumerate all possible changes (ADD/SWAP/DEL) with L dist 1 and associate them with their change in mass.  Find all sets of changes with L dist <= maximum L dist and change in mass equal to targetMass - neighborMass.  Permute all possible change sets to a canonical form to first apply all ADDs, then all SWAPs, then all DELs.  For each neighborCandidate, apply all possible change sets (with each individual step potentially occurring at each position of the candidate).
    #Ugh, this is possibly better, but quite hard to implement correctly and loses a lot of its potential due to the multiple ways to get to each candidate.

    #Option 4:  Potential compromise:  Project all candidates into a "composition space" where only their compositions, not their order, is maintained.  This space is an ||A|| dimensional space where each axis is the number of amino acids of each label in the alphabet.  L dist in this space is a taxicab distance, so we can very quickly generate the set of possible compositions that are within a maximum L dist and narrow in on compositions that allow for the target mass.  We can then check all permutations of these compositions for L dist to our existing candidates.  If our candidate set contains all permutations of every composition, then this is the correct answer, as there is no reason to check for L dist.

    #I think we'll go with option 2 for now and attempt option 3 later.
    #For each neighbor, generate all candidates, union them together, filter out anything with the wrong mass.
#    allCandidates = []
#    totalComplete = 0
#    for neighbor in neighborCandidates:
#        candidatesFromNeighbor = generateCandidatesSingle2(neighbor, maximumLDist)
#        print("Unfiltered: " + str(len(candidatesFromNeighbor)))
#        filtered = 0
#        for candidate in candidatesFromNeighbor:
#            if fuzzyEqual(calcMass(candidate), targetMass, massEpsilon):
#                filtered = filtered + 1
#                allCandidates.append(candidate)
#        print ("Filtered: " + str(filtered))
#        totalComplete += 1
#        print(str(totalComplete) + "/" + str(len(neighborCandidates)))
#
#    allCandidates = alphabet.removeStrDups(allCandidates)
#
#    return allCandidates

    #Now an attempt at option 3!
    allCandidates = []
    compositionToNeighbor = defaultdict(list)
    for neighbor in neighborCandidates:
        comp = Composition.fromAminoArr(neighbor)
        compositionToNeighbor[comp].append(neighbor)
    
    allPossibleTransitionPaths = None
    for comp in compositionToNeighbor:
        mass = comp.calcMass()
        allPossibleTransitionPaths = Composition.generatePossiblePaths(mass, targetMass, maximumLDist, massEpsilon * 10, alphabet) #TODO FIXME HACK:  Why do people keep passing candidates that are nowhere near the mass epsilon tolerances???
        print("Possible transition paths: " + str(len(allPossibleTransitionPaths)) + " within tolerance: " + str(massEpsilon * 10))
        break
    
    totalComplete = 0
    for composition in compositionToNeighbor:
        neighborList = compositionToNeighbor[composition]
#        (transitionPaths, endingComps) = Composition.generatePaths(composition, targetMass, maximumLDist, massEpsilon, alphabet)
        transitionPaths = composition.filterPaths(allPossibleTransitionPaths, targetMass, massEpsilon)
        
        for neighbor in neighborList:
            candidatesFromNeighbor = generateCandidatesSingle3(neighbor, maximumLDist, transitionPaths, alphabet)
            for candidate in candidatesFromNeighbor:
                allCandidates.append(candidate)
        totalComplete += len(neighborList)

        print("Starting Composition Mass: " + str(composition.calcMass()))
        print("Paths: " + str(transitionPaths))
#        print("Filts: " + str(transitionPaths2))
#        print("Comps: " + str(endingComps))
        print("Num transition paths: " + str(len(transitionPaths)))
#        print("Num transition paths2: " + str(len(transitionPaths2)))
#        if len(transitionPaths) != len(transitionPaths2):
#            print("-----")
#            print(transitionPaths)
#            print("VS")
#            print(transitionPaths2)
#            for comp in endingComps:
#                print(comp.calcMass())
#            print("Starting Composition Mass: " + str(composition.calcMass()))
#            print("Claimed node mass: " + str(neighborMass))
#            print("-----")
#            raise Exception("Shouldn't Happen")
#        print("Num final comps: " + str(len(endingComps)))
        print(str(totalComplete) + "/" + str(len(neighborCandidates)))

    allCandidates = alphabet.removeStrDups(allCandidates)
    return allCandidates

def combineCandidateSets(listOfCandidateSets):
    #Potentially intersect sets from all neighbors if we believe this will give us consistency, potentially union all the sets if we just want as many candidates as we can get, potentially some weird partial operation where a candidate makes the cut if it is suggested by at least 50% of neighbors...
    #We'll start simple and do a union.
    x = []
    for candidateSet in listOfCandidateSets:
        for candidate in candidateSet:
            x.append(candidate)
    x = DEFAULT_ALPHABET.removeStrDups(x)

    return x

def generateAllCandidatesByNode(SG, unannotated, annotated, lDist = 2, epsilon = .02):
    allCandidatesByNode = {}
    for key in unannotated:
        neighboringNodes = SG.nodes[key].neighbors
        #Find all the neighboring nodes
        #For each annotated neighboring node, generate a candidate set for this node as the set of candidates that are within L dist of x and a matching mass
        listOfCandidateSets = []
        for neighbor in neighboringNodes:
            neighborID = neighbor.ID
            if neighborID not in annotated:
                continue #We haven't annotated this node yet, skip.
            listOfCandidateSets.append(generateCandidates(SG.nodes[key].mass-1, SG.nodes[neighborID].mass-1, annotated[neighborID], lDist, epsilon))
        #Potentially intersect sets from all neighbors if we believe this will give us consistency, potentially union all the sets if we just want as many candidates as we can get, potentially some weird partial operation where a candidate makes the cut if it is suggested by at least 50% of neighbors...
        finalCandidateSet = combineCandidateSets(listOfCandidateSets)
        allCandidatesByNode[key] = finalCandidateSet

    return allCandidatesByNode

#TODO FIXME HACK:  argh, need to check for connectedness or we'll report nodes we can't possibly annotate
def getUnannotated(all, annotated):
    unannotated = []
    for key in all:
        if key not in annotated and hasattr(all[key], "allSpectra"):
            unannotated.append(key)
    return unannotated

def assignCandidates(SG, initialCandidates):
    annotatedNodes = initialCandidates
    
    while len(annotatedNodes) < len(SG.nodes): #more nodes to annotate
        unannotatedList = getUnannotated(SG.nodes, annotatedNodes)

        #Generate all possible candidates
        allCandidatesByNode = generateAllCandidatesByNode(SG, unannotatedList, annotatedNodes)

        #Evaluate all candidates by their score and maximum possible pairwise score
        candidateScores = []
        
        print("Scoring candidate sequences...")
        for nodeID in allCandidatesByNode:
            node = SG.nodes[nodeID]
            print(nodeID + " " + str(node.mass))
            nodeExperimentalSpectra = list(parseMGF.findMasses(x[1]) for x in node.allSpectra)[0]
            i = 0
            for candidate in allCandidatesByNode[nodeID]:
                candidateTheoreticalSpectra = parseMGF.genTheoreticalMasses(DEFAULT_ALPHABET.toMassList(candidate))
                score = scoring.scoreSingleAnnotation(nodeExperimentalSpectra, candidateTheoreticalSpectra, .02)
                candidateScores.append((node, candidate, len(score)))
                i += 1
                if (i % 1000) == 0:
                    print(str(i) + "/" + len(allCandidatesByNode[nodeID]))

        candidateScores.sort(key = lambda x:x[2], reverse = True)
        
        if len(candidateScores) <= 0:
            break
    
        chosen = candidateScores[0]
        associatedTuples = []
        for tuple in candidateScores:
            if tuple[0] == chosen[0]:
                associatedTuples.append(tuple)
        annotatedNodes[chosen[0].ID] = []

        #Take Top M - This is branch and bound - The larger this M, the less chance of dropping the correct answer, but the longer we will take to run the algorithm.
        i = 0
        for tuple in associatedTuples:
            annotatedNodes[chosen[0].ID].append(tuple[1])
            i += 1
            if i == M:
                break

        expectedAnswer = DEFAULT_ALPHABET.toSpinInvariant(DEFAULT_ALPHABET.toCanonicalAminoArr(SG.nodes[chosen[0].ID].CORRECT_ANSWER.massList))

        print("Annotating node: " + str(chosen[0]))
        print("Candidates at node: " + str(len(annotatedNodes[chosen[0].ID])))
        print("Expected Answer at node: " + str(expectedAnswer))
        print("Expected Answer Weight: " + str(sum(SG.nodes[chosen[0].ID].CORRECT_ANSWER.massList)))
        
        foundIndex = None
        i = 0
        for annotation in annotatedNodes[chosen[0].ID]:
            if str(annotation) == str(expectedAnswer):
                foundIndex = i
                break
            i += 1
        
        print("Expected Answer Position: " + str(foundIndex))
        print(str(annotatedNodes[chosen[0].ID]))

        raw_input("Okay")


#MAIN PROGRAM START
print("\n\n---START---\n\n");

spectralGraphNodes = parsers.readSpectralGraphNodes("../data/small_surugamide_clusters_send2Daniel.txt")
spectralGraphEdges = parsers.readPairsInfo("../data/edges_e014513b89b34710a4909dadaee466bb.pairsinfo")

#Debug the data read from the files
for node in spectralGraphNodes:
    print(node["clusterIndex"] + ": " + node["parentMass"])

for edge in spectralGraphEdges:
    print(edge["A"] + " -> " + edge["B"])

SG = SpectralGraph(spectralGraphNodes, spectralGraphEdges, DEFAULT_ALPHABET)
SG.retrieveCandidates("../data/surugamide/", DEFAULT_ALPHABET)

raw_input("Starting Reading In Spectra")
import parseMGF
SG.nodes["69"].allSpectra = parseMGF.readMGF("../data/surugamide/spectra/cyclonovo_785.5337524.mgf")
SG.nodes["73"].allSpectra = parseMGF.readMGF("../data/surugamide/spectra/cyclonovo_799.5462036.mgf")
SG.nodes["86"].allSpectra = parseMGF.readMGF("../data/surugamide/spectra/cyclonovo_856.5719605.mgf")
SG.nodes["90"].allSpectra = parseMGF.readMGF("../data/surugamide/spectra/cyclonovo_870.5863648.mgf")
SG.nodes["95"].allSpectra = parseMGF.readMGF("../data/surugamide/spectra/cyclonovo_884.6019287.mgf")
SG.nodes["105"].allSpectra = parseMGF.readMGF("../data/surugamide/spectra/cyclonovo_898.6170044.mgf")
SG.nodes["110"].allSpectra = parseMGF.readMGF("../data/surugamide/spectra/cyclonovo_912.6369019.mgf")
SG.nodes["113"].allSpectra = parseMGF.readMGF("../data/surugamide/spectra/cyclonovo_914.61.mgf")
raw_input("Setting expected answers")
SG.nodes["69"].setCorrectAnswer(SpectralGraphNodeCandidate(SG.nodes["69"], 0, None, [113.084, 71.0371, 99.0684, 113.084, 128.095, 147.068, 113.084], "FAKE", 7000, -1, DEFAULT_ALPHABET))
SG.nodes["73"].setCorrectAnswer(SpectralGraphNodeCandidate(SG.nodes["73"], 1, None, [71.0371, 113.084, 113.084, 128.095, 113.084, 147.068, 113.084], "FAKE", 7000, -1, DEFAULT_ALPHABET))
SG.nodes["86"].setCorrectAnswer(SpectralGraphNodeCandidate(SG.nodes["86"], 2, None, [99.0684, 71.0371, 99.0684, 99.0684, 128.095, 99.0684, 147.068, 113.084], "FAKE", 7000, -1,DEFAULT_ALPHABET))
SG.nodes["90"].setCorrectAnswer(SpectralGraphNodeCandidate(SG.nodes["90"], 3, None, [99.0684, 71.0371, 99.0684, 99.0684, 128.095, 113.084, 147.068, 113.084], "FAKE", 7000, -1,DEFAULT_ALPHABET))
SG.nodes["95"].setCorrectAnswer(SpectralGraphNodeCandidate(SG.nodes["95"], 4, None, [99.0684, 71.0371, 99.0684, 113.084, 128.095, 113.084, 147.068, 113.084], "FAKE", 7000, -1,DEFAULT_ALPHABET))
SG.nodes["105"].setCorrectAnswer(SpectralGraphNodeCandidate(SG.nodes["105"], 5, None, [99.0684, 71.0371, 113.084, 113.084, 128.095, 113.084, 147.068, 113.084], "FAKE", 7000, -1,DEFAULT_ALPHABET))
SG.nodes["110"].setCorrectAnswer(SpectralGraphNodeCandidate(SG.nodes["110"], 6, None, [113.084, 71.0371, 113.084, 113.084, 128.095, 113.084, 147.068, 113.084], "FAKE", 7000, -1,DEFAULT_ALPHABET))
SG.nodes["113"].setCorrectAnswer(SpectralGraphNodeCandidate(SG.nodes["113"], 7, None, [99.0684, 71.0371, 113.084, 113.084, 128.095, 113.084, 163.063, 113.084], "FAKE", 7000, -1,DEFAULT_ALPHABET))



START_NODE_ID = "110" #TODO FIXME HACK:  Input this

initialCandidates = {}
initialCandidates[START_NODE_ID] = []

temp = []
for c in SG.nodes[START_NODE_ID].candidates:
    aminoArr = DEFAULT_ALPHABET.toCanonicalAminoArr(c.massList)
    spinInvariantAminoArr = DEFAULT_ALPHABET.toSpinInvariant(aminoArr)
    temp.append(spinInvariantAminoArr)
    
    #TEST SCORING
#    candidateTheoreticalSpectra = parseMGF.genTheoreticalMasses(c.massList)
#    nodeExperimentalSpectra = list(parseMGF.findMasses(x[1]) for x in SG.nodes[START_NODE_ID].allSpectra)[0]
#    score = scoring.scoreSingleAnnotation(nodeExperimentalSpectra, candidateTheoreticalSpectra, .02)
#    print(score)
#    print("Expected Score: " + str(c.score))
#    print("Computed Score: " + str(len(score)))
#    raw_input("Ok")

temp = DEFAULT_ALPHABET.removeStrDups(temp)
temp2 = []

nodeExperimentalSpectra = list(parseMGF.findMasses(x[1]) for x in SG.nodes[START_NODE_ID].allSpectra)[0]
for aminoArr in temp:
    candidateTheoreticalSpectra = parseMGF.genTheoreticalMasses(DEFAULT_ALPHABET.toMassList(aminoArr))
    score = scoring.scoreSingleAnnotation(nodeExperimentalSpectra, candidateTheoreticalSpectra, .02)
    temp2.append((aminoArr, score))

temp2.sort(key = lambda x:x[1], reverse = True) #sort candidates by score descending
#temp2 = temp2[0:M] #Take top M elements  #Ugh... this trashes the necessary candidate even at M=100... what's wrong here?

initialCandidates[START_NODE_ID] = []
for tuple in temp2:
    initialCandidates[START_NODE_ID].append(tuple[0])

print("Num initial candidates: " + str(len(initialCandidates[START_NODE_ID])))

raw_input("Done checking scoring module")

assignCandidates(SG, initialCandidates)
