import os
import parsers

"""A node in the spectral network graph.  These nodes contain an adjacency list and a set of candidate cyclopeptides"""
class SpectralGraphNode:
    def __init__(self, tsvRow):
        self.tsvRow = tsvRow
        self.ID = tsvRow["clusterIndex"]
        self.mass = float(tsvRow["parentMass"])
        self.neighbors = []
        self.candidates = []

    def addNeighbor(self, neighbor):
        self.neighbors.append(neighbor)

    def addCandidate(self, c):
        self.candidates.append(c)

    def __str__(self):
        return "(" + self.ID + ", " + str(self.mass) + ")"

    def __repr__(self):
        return "(" + self.ID + ", " + str(self.mass) + ")"


    def setCorrectAnswer(self, c):
        self.CORRECT_ANSWER = c

"""An edge in the candidate graph"""
class CandidateEdge:
    def __init__(self, candA, candB, cyclicEditDist):
        self.A = candA
        self.B = candB
        self.editDist = cyclicEditDist

    def follow(self, fromCandidate):
        if fromCandidate is self.A:
            return self.B
        if fromCandidate is self.B:
            return self.A

"""A candidate cyclopeptide sequence explaining the mass measured at a SpectralGraphNode."""
class SpectralGraphNodeCandidate:
    def __init__(self, parentNode, ID, tsvRow, massList, graphFileName, score, rank, L):
        self.tsvRow = tsvRow
        self.massList = massList
        self.graphFile = graphFileName
        canonNameArr = L.toCanonicalNameArr(massList)
        self.name = L.toSpinInvariant(canonNameArr)
        self.counts = L.toCanonicalCounts(self.name)
        self.ID = ID
        self.parentNode = parentNode
        self.neighbors = []
        self.score = float(score)
        self.rank = rank

    def genSimpleName(self):
        s = ""
        for i in range(len(self.name)-1):
            amino = self.name[i]
            s += amino[0:3] + "-"
        s += self.name[-1][0:3]
        return s


"""The spectral network graph itself.  This is an undirected graph made of SpectralGraphNode instances."""
class SpectralGraph:
    def __init__(self, spectralGraphNodes, spectralGraphEdges, L):
        self.nodes = {}
        for nodeRow in spectralGraphNodes:
            node = SpectralGraphNode(nodeRow)
            self.nodes[node.ID] = node
        for edgeRow in spectralGraphEdges:
            aID = edgeRow["A"]
            bID = edgeRow["B"]
            
            if aID in self.nodes and bID in self.nodes:
                self.nodes[aID].addNeighbor(self.nodes[bID])
                self.nodes[bID].addNeighbor(self.nodes[aID])
            else:
                print("WARN: Could not add edge: " + aID + " -> " + bID)
        self.L = L

    def __str__(self):
        s = ""
        s += "Nodes: \n"
        for nodeKey in self.nodes:
            s += str(self.nodes[nodeKey])
            s += "\n"
        s += "Edges: \n"
        for nodeKey in self.nodes:
            s += str(self.nodes[nodeKey]) + "\n"
            for neighbor in self.nodes[nodeKey].neighbors:
                s += "\t-> " + str(neighbor) + "\n"
#        s += "Candidates: \n"
#        for nodeKey in self.nodes:
#            s += str(self.nodes[nodeKey]) + "\n"
#            for c in self.nodes[nodeKey].candidates:
#                s += "\t" + str(c.massList) + " (" + str(L.toCanonicalNameArr(c.massList)) + ")\n"
        return s
        
    """Retrieve candidate peptides for each mass from the filesystem"""
    def retrieveCandidates(self, folder, L):
        idGen = 0
        if not folder.endswith("/"):
            folder += "/"
        networkCandidates = []
        files = os.listdir(folder)
        for s in files:
            if s.endswith("_all_psms.concise.tsv"):
                candidatesFile = s
                peptideMass = s[0:len(s)-len("_all_psms.concise.tsv")]
                graphsFolder = "graphs_" + peptideMass + "_dir"
                networkCandidates.append((peptideMass, folder + candidatesFile, folder + graphsFolder))

        #print(networkCandidates)

        for candidate in networkCandidates:
            SGNode = self.findNodeByMass(float(candidate[0]), 0.01)
            candidateTSV = parsers.readCandidates(candidate[1])
            #columns = ["graphFile", "mgfFile", "score", "pValue"]
            scoreSet = set([])
            for row in candidateTSV:
                scoreSet.add(float(row["score"]))
            sortedUniques = sorted(list(scoreSet), reverse=True)
            scoreToRank = {}
            for i in range(len(sortedUniques)):
                scoreToRank[sortedUniques[i]] = i
            for row in candidateTSV:
                graphFile = row["graphFile"]
                graphFileName = os.path.basename(graphFile)
                graphFileNameRelativePath = candidate[2] + "/" + graphFileName
                massList = parsers.parseGraphFile(graphFileNameRelativePath)
                SGNode.addCandidate(SpectralGraphNodeCandidate(SGNode, idGen, row, massList, graphFileName, row["score"], scoreToRank[float(row["score"])], L))
                idGen += 1

    """A O(N) lookup for nodes in the graph by their mass (within a small tolerance)"""
    def findNodeByMass(self, mass, epsilon):
        result = None
        for key in self.nodes:
            n = self.nodes[key]
            if n.mass + epsilon >= mass and n.mass - epsilon <= mass:
                assert result is None, "Found 2+ matching nodes for mass: " + str(mass) + " with epsilon: " + str(epsilon)
                result = n
        return result
