import queue

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
#            score all candidates individually by combination of pairwise score (MST score?) with known candidates and score on associated spectra
#            Take best candidate by score, add (candidate, nodeContainingCandidate) to annotatedNodes
#
#        Print network chosen and see how it lines up with our existing scoring mechanisms.
#
#        //How does this handle multiple candidates in the initial input?
#    }

def calcMass(candidate):
    totalMass = 0
    for amino in candidate:
        totalMass += amino.mass
    return totalMass


def generateCandidates(neighbor, maximumLDist, alphabet = DEFAULT_ALPHABET):
    #Ugh!
    #Can formulate as a graph search in an absurdly high dimensional space.  Ugh!
    #Use breadth first search.  Store everything in spin invariant canonical.
    nodeToDist = {}
    nodeToVal = {}
    nodeToVisited = {}
    
    val = AminoLanguage.toSpinInvariant(AminoLanguage.toCanonicalAminoArr(neighbor))
    strVal = str(val)
    nodeToDist[strVal] = 0
    nodeToVal[strVal] = val
    
    def evaluateNeighbor(newVal, newDist, nodeToDist, nodeToVal, toVisit):
        newVal = AminoLanguage.toSpinInvariantCanonicalNameArr(newVal)
        newStr = str(newVal)
        if newStr in nodeToDist:
            nodeToDist[newStr] = min(newDist, nodeToDist[newStr])
        else:
            nodeToDist[newStr] = newDist
        nodeToVal[newStr] = newVal
        toVisit.put(newStr)

    
    toVisit = queue.Queue()
    toVisit.put(strVal)
    while len(toVisit) > 0:
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
        for position in range(0, len(curVal)):
            #Note that inserting at position 0 and adding to end of list are the same because these are cyclopeptides.
            for toAdd in alphabet.sortedAminos:
                newVal = list(curVal)
                newVal.insert(position, toAdd)
                evaluateNeighbor(newVal, curDist + 1, nodeToDist, nodeToVal, toVisit)

    resultsList = []
    for strVal in nodeToVal:
        resultsList.append(nodeToVal[strVal])
    return resultsList


def generateCandidates(targetMass, neighborMass, neighborCandidates, maximumLDist, massEpsilon):
    #Ugh, there are a lot of ways we could generate candidates depending on the assumptions we're willing to make.  The problem seems tractable with a fixed alphabet, working at the level of masses directly we could generate a heck of a lot of potential candidates and it isn't clear to me which of those make biological sense.  So we'll start with a fixed alphabet for now and keep this function separated so we can try different implementations.

    #Option 1:  Absolute worst possible correct algorithm:  Enumerate all strings between length min neighborCandidate length - maximumLDist and length max neighborCandidate length + maximumLDist.  Check all strings for Ldist to all neighbor candidates and proximity to target mass.
    #Ugh, this is basically a password cracker.  Depending on length of alphabet, this could take quite a while.  Could potentially speed up the Ldist check based on some tree.

    #Option 2:  Second worst possible correct algorithm:  For each neighbor candidate, generate all strings within L dist, filter out all strings with the wrong target mass.
    #The generation step is tough, and will potentially generate the same string multiple times.  This should still be significantly faster than checking a string for its L dist to each neighbor candidate.

    #Option 3:  Potentially better algorithm:  Enumerate all possible changes (ADD/SWAP/DEL) with L dist 1 and associate them with their change in mass.  Find all sets of changes with L dist <= maximum L dist and change in mass equal to targetMass - neighborMass.  Permute all possible change sets to first apply all ADDs, then all SWAPs, then all DELs.  For each neighborCandidate, apply all possible change sets (with each individual step potentially occurring at each position of the candidate).
    #Ugh, this is possibly better, but quite hard to implement correctly and loses a lot of its potential due to the multiple ways to get to each candiate.

    #Option 4:  Potential compromise:  Project all candidates into a "composition space" where only their compositions, not their order, is maintained.  This space is an ||A|| dimensional space where each axis is the number of amino acids of each label in the alphabet.  L dist in this space is a taxicab distance, so we can very quickly generate the set of possible compositions that are within a maximum L dist and narrow in on compositions that allow for the target mass.  We can then check all permutations of these compositions for L dist to our existing candidates.  If our candidate set contains all permutations of every composition, then this is the correct answer, as there is no reason to check for L dist.

    #I think we'll go with option 2 for now.
    #For each neighbor, generate all candidates, union them together, filter out anything with the wrong mass.
    allCandidates = []
    for neighbor in neighborCandidates:
        candidatesFromNeighbor = generateCandidates(neighbor, maximumLDist)
        for candidate in candidatesFromNeighbor:
            allCandidates.append(candidate)

    allCandidates = AminoLanguage.removeStrDups(allCandidates)

    filteredCandidates = []
    for candidate in allCandidates:
        if fuzzyEqual(calcMass(candidate), targetMass, massEpsilon):
            filteredCandidates.add(candidate)

    return filteredCandidates

def combineCandidateSets(listOfCandidateSets):
    #Potentially intersect sets from all neighbors if we believe this will give us consistency, potentially union all the sets if we just want as many candidates as we can get, potentially some weird partial operation where a candidate makes the cut if it is suggested by at least 50% of neighbors...
    #We'll start simple and do a union.
    x = set([])
    for candidateSet in listOfCandidateSets:
        for candidate in candidateSet:
            x.add(candidate)
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
            listOfCandidateSets.append(generateCandidates(SG.nodes[key].mass, annotated[neighborID], lDist, epsilon))
        #Potentially intersect sets from all neighbors if we believe this will give us consistency, potentially union all the sets if we just want as many candidates as we can get, potentially some weird partial operation where a candidate makes the cut if it is suggested by at least 50% of neighbors...
        finalCandidateSet = combineCandidateSets(listOfCandidateSets)
        
        allCandidatesByNode[key] = []


def getUnannotated(all, annotated):
    unannotated = []
    for key in all:
        if key not in annotated:
            unannotated.append(key)
    return unannotated

def assignCandidates(SG, initialCandidates):
    annotatedNodes = initialCandidates
    while len(annotatedNodes) < len(SG): #more nodes to annotate
        unannotatedList = getUnannotated(SG, annotatedNodes)

        #Generate all possible candidates
        allCandidatesByNode = generateAllCandidates(SG, unannotatedList, annotatedNodes)

