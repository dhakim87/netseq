"""
    Edges must have fields A and B that point to their nodes
    Nodes must have an ID field with suitable equality and hashing implementations
"""
def spanningTree(edgeList, weightFunc, maximal=False):
    #Kruskal's Algorithm:
    #Sort edges by their weight.  (If maximal is true, sort desc rather than asc)
    #Loop through edges, add edge to result if endpoints of edge are in different sets.
    def findSetRepresentative(lookupTable, candidate):
        id = candidate.ID
        while id in lookupTable and id != lookupTable[id]:
            id = lookupTable[id]
        return id

    def areInSameSet(lookupTable, candidateA, candidateB):
        return findSetRepresentative(lookupTable, candidateA) == findSetRepresentative(lookupTable, candidateB)

    def mergeSets(lookupTable, candidateA, candidateB):
        aRepr = findSetRepresentative(lookupTable, candidateA)
        bRepr = findSetRepresentative(lookupTable, candidateB)
        newRepr = min(aRepr, bRepr)
        id = candidateA.ID
        while id in lookupTable and id != lookupTable[id]:
            temp = id
            id = lookupTable[id]
            lookupTable[temp] = newRepr
        id = candidateB.ID
        while id in lookupTable and id != lookupTable[id]:
            temp = id
            id = lookupTable[id]
            lookupTable[temp] = newRepr
        lookupTable[aRepr] = newRepr
        lookupTable[bRepr] = newRepr

    sortedEdgeList = sorted(edgeList, key=weightFunc, reverse=maximal)
    lookupTable = {}
    MST = []

    for potentialEdge in sortedEdgeList:
        if not areInSameSet(lookupTable, potentialEdge.A, potentialEdge.B):
            MST.append(potentialEdge)
            mergeSets(lookupTable, potentialEdge.A, potentialEdge.B)

    return MST
