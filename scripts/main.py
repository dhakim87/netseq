from collections import defaultdict
#import networkx as nx
import matplotlib.pyplot as plt
import levenshtein as Lev
import spectralPlotting
import sys
import os

SHOW_GRAPHS = False  #Whether or not to show graphs with networkx and matplotlib (VERY SLOW)
PRINT_CANDIDATES = False #Whether or not to print out the candidates at various points in the algorithm (Useful for debugging but a ton of text)
MAX_MUTATIONS_PER_EDGE = 2 #How many mutations we allow per edge in the network graph, this determines where we add edges in the candidate graph.  (>1 VERY SLOW)
DEFAULT_U_CONSISTENCY_THRESH = 0 #Override: U=xxx (>0 VERY SLOW)
DEFAULT_CANDIDATE_SCORE_THRESH = 0 #Override: MIN_SCORE=xxx

OUTPUT = "short" #Supports "mass", "amino", "filename", or "short"

"""An Amino Acid"""
class Amino:
    def __init__(self, tsvRow):
        self.tsvRow = tsvRow
        self.name = tsvRow["uc"]
        self.mass = float(tsvRow["mass"])

    def __str__(self):
        return str(self.mass) + " (" + self.name + ")";

"""A grab bag of amino acids.  ID'd to allow fast hash/eq lookups.  You must manually intern these objects to use their faster behavior."""
class AminoComposition:
    idGen = 0
    def __init__(self, aminoDict):
        self.ID = AminoComposition.idGen
        AminoComposition.idGen += 1
        self.aminoDict = aminoDict

    def __hash__(self):
        return hash(self.ID)

    def __eq__(self, other):
        return self.__class__ == other.__class__ and self.ID == other.ID

"""The language of amino acids used in construction of cyclopeptides"""
class AminoLanguage:
    MASS_EPSILON = 0.02 #Da -- This is the maximum difference between a mass and the amino acid we will label that mass as.  This should be defined based on the accuracy/precision of the mass spectrometer in use.
    fastLookup = {} #Caches and returns results of the toCanonicalName function
    compositionSet = {} #Allows interning of AminoCompositions for faster lookup and equality

    def __init__(self, sortedAminos):
        self.sortedAminos = sortedAminos

    """The canonical name of a mass is a unique string defining the set of all amino acids within MASS_EPSILON of that mass"""
    def toCanonicalName(self, mass):
        s = ""
        if mass in self.fastLookup:
            return self.fastLookup[mass]
        
        #todo - use binary search
        for amino in self.sortedAminos:
            if mass > amino.mass + AminoLanguage.MASS_EPSILON:
                continue
            if mass < amino.mass - AminoLanguage.MASS_EPSILON:
                break
            if s == "":
                s += amino.name
            else:
                s += "/" + amino.name
        if s == "":
            s = "?!?" #If you want to keep running with this unknown mass and not corrupt your results, add it to the amino list with a unique 3 letter code.
            raise Exception("Unexpected Amino Weight")
            
        AminoLanguage.fastLookup[mass] = s
        return s

    """The canonical name of an array of masses is an array of strings that each uniquely represent amino acids within MASS_EPSILON of each mass"""
    def toCanonicalNameArr(self, massList):
        result = []
        for m in massList:
            result.append(self.toCanonicalName(m))
        return result

    """The set of canonical name array spins is what you get by rotating the strings in a canonical name array.  Because we are working with cyclopeptides, all of these spins are considered identical"""
    def generateAllCanonicalNameSpins(self, canonicalNameArr):
        allSpins = []
        for i in range(len(canonicalNameArr)):
            allSpins.append(canonicalNameArr[i:] + canonicalNameArr[:i])
        return allSpins

    """To quickly identify spins of the same cyclopeptide, all spin sets can be converted to a spinInvariantCanonicalNameArr.  This is the rotation of the cyclopeptide that results in the least value string, sorted lexicographically.  For maintaining this interface, all that matters is that all rotations of a canonical name arr are assigned the same string value, but cyclopeptides that are not considered equal must generate different spin invariant canonical names"""
    def toSpinInvariantCanonicalNameArr(self, canonicalNameArr):
        allSpins = self.generateAllCanonicalNameSpins(canonicalNameArr)
        allSpins.sort(key=lambda x: str(x))
        return allSpins[0]

    """A composition is a dictionary from amino acids to counts.  By interning one of these dictionaries, you can later quickly compare them for equality or use them as keys in hashmaps"""
    def internComposition(self, aminoDict):
        #Interns compositions like they are strings
        list = []
        for key in sorted(aminoDict.keys()):
            list.append((key, aminoDict[key]))
        strVal = str(list)
        if strVal not in AminoLanguage.compositionSet:
            AminoLanguage.compositionSet[strVal] = AminoComposition(aminoDict)
        return AminoLanguage.compositionSet[strVal]

    """Counts the amino acids in a canonical name array and returns an interned AminoComposition object containing these counts"""
    def toCanonicalCounts(self, canonicalNameArr):
        aminoCounter = defaultdict(int)
        for amino in canonicalNameArr:
            aminoCounter[amino] += 1
        return self.internComposition(aminoCounter)
    
    def couldBeRelated(self, aCounts, bCounts, maxIndels):
        delta = {}
        allKeys = set([])
        for key in aCounts.aminoDict:
            allKeys.add(key)
        for key in bCounts.aminoDict:
            allKeys.add(key)

        for key in allKeys:
            diff = bCounts.aminoDict[key] - aCounts.aminoDict[key]
            if diff != 0:
                delta[key] = diff
        
        #Now that we have a set of adds and deletes, let's ensure that it is possible to do in maxIndels.  If adds or deletes of any amino acid requires more than maxIndels,
        #we can't possibly do it, so we can fail out early.
        
        for key in delta:
            if abs(delta[key]) > maxIndels:
                return False
        
        #Imagine these compositions of A and B existing in an ||L|| dimensional space with axes named by each amino acid.
        #We have computed the difference between these two compositions stored it into this delta dictionary.
        #Since we are counting indels, ADDs and DELs shift a point by 1 along some axis.
        #SWAPs are diagonal lines in this space, a combined DEL+ADD.
        #The shortest path between two compositions is the path with the most SWAPs, as that can be as low as half the length of path created with ADDs and DELs.
        #This means that the minimum possible path length between two compositions is the sum of the absolute values of all values in the delta dictionary, divided by 2.
        
        addDelLength = 0
        for key in delta:
            addDelLength += abs(delta[key])
        minPathLength = addDelLength / 2.0
        if minPathLength > maxIndels:
            return False
        return True
    
    """ Removes duplicate entities in a list.  Duplicate is determined by hash and equality of string value of each item """
    def removeStrDups(self, l):
        dupRemover = set([])
        resultList = []
        for item in l:
            itemStr = str(item)
            if itemStr in dupRemover:
                continue
            dupRemover.add(itemStr)
            resultList.append(item)
        return resultList
    
    """A Levenshtein distance based method for adding edges pairwise between candidates in two nodes."""
    def findRelations2(self, SGNodeA, SGNodeB):
        i = 0
        for candA in SGNodeA.candidates:
            spins = self.generateAllCanonicalNameSpins(candA.name)
            i += 1
            print(str(i) + " / " + str(len(SGNodeA.candidates)))
            for candB in SGNodeB.candidates:
                if not self.couldBeRelated(candA.counts, candB.counts, MAX_MUTATIONS_PER_EDGE):
                    continue
                
                closestDist = max(len(candA.name), len(candB.name), MAX_MUTATIONS_PER_EDGE + 1)
                for candASpin in spins:
                    dist = Lev.compute2(candASpin, candB.name, MAX_MUTATIONS_PER_EDGE)
                    if dist is not None:
                        closestDist = min(closestDist, dist)
                if closestDist <= MAX_MUTATIONS_PER_EDGE:
                    edge = CandidateEdge(candA, candB, closestDist)
                    candA.neighbors.append(edge)
                    candB.neighbors.append(edge)

    """Computes cyclic levenshtein distance ignoring any maximum delta"""
    def computeCyclicEditDist(self, candA, candB):
        spins = self.generateAllCanonicalNameSpins(candA.name)
        maxDist = max(len(candA.name), len(candB.name))
        closestDist = maxDist
        for candASpin in spins:
            dist = Lev.compute2(candASpin, candB.name, maxDist)
            if dist is not None:
                closestDist = min(closestDist, dist)
        return closestDist

"""Parse a tab separated value file into an array of rows.  Any rows in the file that do not have enough entries to fill all columns will be skipped and warnings will be printed to console.  You may optionally skip the first line of the file (to remove header info) by passing True as the third argument."""
def parseTSV(file, columns, skipFirstLine):
    firstLine = skipFirstLine
    
    arr = []
    arrIndex = 0
    with open(file, 'r') as f:
        for line in f:
            row = []
            if (firstLine):
                firstLine = False
                continue
            
            ss = line.split("\t")
            
            if (len(ss) < len(columns)):
                print("WARNING: COULD NOT PARSE ROW: " + (str)(len(ss)))
                continue
            for i in range(0, len(columns)):
                row.append(ss[i])

            arr.append(row)
            arrIndex = arrIndex + 1
    return arr;

"""Converts an array of objects and an array of column headers into a dictionary"""
def rowToDict(row, columns):
    x = {}
    for index, val in enumerate(row):
        x[columns[index]] = row[index];
    return x;

"""Converts an array of rows and an array of column headers into an array of dictionaries"""
def arrToDicts(arr, columns):
    allRows = []
    for row in arr:
        allRows.append(rowToDict(row, columns))
    return allRows;

"""Reads a spectral graph nodes file (no file ending atm)"""
def readSpectralGraphNodes(file):
    columns = ["clusterIndex", "numSpectra", "parentMass", "precursorCharge", "precursorMass", "sumPrecursorIntensity", "G1", "G2", "G3", "G4", "G5", "G6", "allFiles", "allGroups", "defaultGroups", "rtMean", "rtStdErr", "proteoSAFeClusterLink", "uniqueFileSourcesCount", "evenOdd", "libraryID", "numberOrganismIDs", "allOrganisms", "componentindex"]
    allEntries = parseTSV(file, columns, True)
    allEntriesAsDictionaries = arrToDicts(allEntries, columns)
    return allEntriesAsDictionaries

"""Reads a pairs info file (.pairsinfo)"""
def readPairsInfo(file):
    columns = ["A", "B", "deltaMass", "ignore1", "cosine", "ignore3", "ignore4"]
    allEntries = parseTSV(file, columns, False)
    allEntriesAsDictionaries = arrToDicts(allEntries, columns)
    return allEntriesAsDictionaries

"""Reads a candidate file (.concise.tsv)"""
def readCandidates(file):
    columns = ["graphFile", "mgfFile", "score", "pValue"]
    allEntries = parseTSV(file, columns, False)
    allEntriesAsDictionaries = arrToDicts(allEntries, columns)
    return allEntriesAsDictionaries

"""Reads an amino acid language file (aa_mass_list.txt)"""
def readAminoLanguage(file):
    columns = ["lc", "uc", "composition", "mass"]
    allEntries = parseTSV(file, columns, False)
    allEntriesAsDictionaries = arrToDicts(allEntries, columns)
    
    sortedAminos = []
    for row in allEntriesAsDictionaries:
        sortedAminos.append(Amino(row))
    
    sortedAminos.sort(key = lambda x: x.mass)
    L = AminoLanguage(sortedAminos)
    return L

"""Reads a graph file containing a cyclopeptide.  Ensures the file is actually circular, then returns in the form of a mass list.  (.graph)"""
def parseGraphFile(file):
    hasReadComponents = False
    componentsToRead = -1
    hasReadBonds = False
    bondsToRead = -1
    massList = []
    temp = 0
    with open(file, 'r') as f:
        for line in f:
            ss = line.split()
            if not hasReadComponents:
                componentsToRead = int(ss[4])
                hasReadComponents = True
                continue
            if componentsToRead > 0:
                massList.append(float(ss[2]))
                componentsToRead = componentsToRead - 1
                continue
            if not hasReadBonds:
                bondsToRead = int(ss[4])
                hasReadBonds = True
                continue
            if bondsToRead > 0:
                assert int(ss[0]) == temp, "Unsupported: Peptide is not circular"
                temp += 1
                bondsToRead -= 1
                assert int(ss[2]) == temp or (int(ss[2]) == 0 and bondsToRead == 0), "Unsupported: Peptide is not circular"
    return massList

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
        self.name = L.toSpinInvariantCanonicalNameArr(canonNameArr)
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
        if PRINT_CANDIDATES:
            s += "Candidates: \n"
            for nodeKey in self.nodes:
                s += str(self.nodes[nodeKey]) + "\n"
                for c in self.nodes[nodeKey].candidates:
                    s += "\t" + str(c.massList) + " (" + str(L.toCanonicalNameArr(c.massList)) + ")\n"
        return s

    """A O(N) lookup for nodes in the graph by their mass (within a small tolerance)"""
    def findNodeByMass(self, mass, epsilon):
        result = None
        for key in self.nodes:
            n = self.nodes[key]
            if n.mass + epsilon >= mass and n.mass - epsilon <= mass:
                assert result is None, "Found 2+ matching nodes for mass: " + str(mass) + " with epsilon: " + str(epsilon)
                result = n
        return result

"""Retrieve candidate peptides for each mass from the filesystem"""
def retrieveCandidates(SG, folder, L):
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
        SGNode = SG.findNodeByMass(float(candidate[0]), 0.01)
        candidateTSV = readCandidates(candidate[1])
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
            massList = parseGraphFile(graphFileNameRelativePath)
            SGNode.addCandidate(SpectralGraphNodeCandidate(SGNode, idGen, row, massList, graphFileName, row["score"], scoreToRank[float(row["score"])], L))
            idGen += 1

"""Use networkx and matplotlib to display the spectral graph"""
def plotSpectralGraph(SG):
    if not SHOW_GRAPHS:
        return
    G = nx.Graph()
    for nodeKey in SG.nodes:
        node = SG.nodes[nodeKey]
        G.add_node(str(node))
        for neighbor in node.neighbors:
            G.add_edge(str(node), str(neighbor))

    plt.subplot(111)
    nx.draw_networkx(G)
    plt.show()

"""Use networkx and matplotlib to display the candidate graph"""
def plotCandidateGraph(SG):
    if not SHOW_GRAPHS:
        return
    G = nx.Graph()
    for nodeKey in SG.nodes:
        node = SG.nodes[nodeKey]
        for candidate in node.candidates:
            G.add_node(node.ID + ":" + str(candidate.ID))
            for edge in candidate.neighbors:
                neighborCandidate = edge.follow(candidate)
                G.add_edge(node.ID + ":" + str(candidate.ID), neighborCandidate.parentNode.ID + ":" + str(neighborCandidate.ID))
    plt.subplot(111)
    nx.draw_networkx(G)
    plt.show()

#MAIN PROGRAM START
print("\n\n---START---\n\n");

#Read command line args
args = {}
for arg in sys.argv[1:]:
    ss = arg.split("=")
    if len(ss) == 2:
        args[ss[0]] = ss[1]


L = readAminoLanguage("../data/aa_mass_list.txt");

#Override Thresholds based on user input
uConsistency = DEFAULT_U_CONSISTENCY_THRESH
if "U" in args:
    uConsistency = int(args["U"])

print("G-u Consistency: u = " + str(uConsistency))

candidateScoreThresh = DEFAULT_CANDIDATE_SCORE_THRESH
if "MIN_SCORE" in args:
    candidateScoreThresh = int(args["MIN_SCORE"])

print("Candidate Min Score: " + str(candidateScoreThresh))


raw_input("Hit Enter To Start")

spectralGraphNodes = readSpectralGraphNodes("../data/small_surugamide_clusters_send2Daniel.txt")
spectralGraphEdges = readPairsInfo("../data/edges_e014513b89b34710a4909dadaee466bb.pairsinfo")

#Debug the data read from the files
for node in spectralGraphNodes:
    print(node["clusterIndex"] + ": " + node["parentMass"])

for edge in spectralGraphEdges:
    print(edge["A"] + " -> " + edge["B"])

SG = SpectralGraph(spectralGraphNodes, spectralGraphEdges, L)

retrieveCandidates(SG, "../data/surugamide/", L)
#TODO FIXME HACK:  THESE ARE THE ANSWERS FOR SURUGAMIDE.  NEED TO BUILD A FILE FORMAT TO STORE THESE AND READ IT IN FOR AN ARBITRARY NETWORK
SG.nodes["69"].setCorrectAnswer(SpectralGraphNodeCandidate(SG.nodes["69"], 0, None, [113.084, 71.0371, 99.0684, 113.084, 128.095, 147.068, 113.084], "FAKE", 7000, -1, L))
SG.nodes["73"].setCorrectAnswer(SpectralGraphNodeCandidate(SG.nodes["73"], 1, None, [71.0371, 113.084, 113.084, 128.095, 113.084, 147.068, 113.084], "FAKE", 7000, -1, L))
SG.nodes["86"].setCorrectAnswer(SpectralGraphNodeCandidate(SG.nodes["86"], 2, None, [99.0684, 71.0371, 99.0684, 99.0684, 128.095, 99.0684, 147.068, 113.084], "FAKE", 7000, -1,L))
SG.nodes["90"].setCorrectAnswer(SpectralGraphNodeCandidate(SG.nodes["90"], 3, None, [99.0684, 71.0371, 99.0684, 99.0684, 128.095, 113.084, 147.068, 113.084], "FAKE", 7000, -1,L))
SG.nodes["95"].setCorrectAnswer(SpectralGraphNodeCandidate(SG.nodes["95"], 4, None, [99.0684, 71.0371, 99.0684, 113.084, 128.095, 113.084, 147.068, 113.084], "FAKE", 7000, -1,L))
SG.nodes["105"].setCorrectAnswer(SpectralGraphNodeCandidate(SG.nodes["105"], 5, None, [99.0684, 71.0371, 113.084, 113.084, 128.095, 113.084, 147.068, 113.084], "FAKE", 7000, -1,L))
SG.nodes["110"].setCorrectAnswer(SpectralGraphNodeCandidate(SG.nodes["110"], 6, None, [113.084, 71.0371, 113.084, 113.084, 128.095, 113.084, 147.068, 113.084], "FAKE", 7000, -1,L))
SG.nodes["113"].setCorrectAnswer(SpectralGraphNodeCandidate(SG.nodes["113"], 7, None, [99.0684, 71.0371, 113.084, 113.084, 128.095, 113.084, 163.063, 113.084], "FAKE", 7000, -1,L))

#TODO FIXME HACK:  THESE ARE THE FILTERED SPECTRA FOR SURUGAMIDE.  NEED TO BUILD A FILE FORMAT TO STORE THESE AND READ IT IN FOR AN ARBITRARY NETWORK
#63: 771.517 - I don't have this one
#69: 785.534
#73: 799.54
#86: 856.571
#90: 870.586
#95: 884.599
#105: 898.616
#110: 912.632
#113: 914.61

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
raw_input("Finished Reading In Spectra")

def wipeEdge(SG, aID, bID):
    print("WIPING EDGE: " + str(SG.nodes[aID]) + " <-> " + str(SG.nodes[bID]) + " FOR TEST PURPOSES")
    SG.nodes[aID].neighbors.remove(SG.nodes[bID])
    SG.nodes[bID].neighbors.remove(SG.nodes[aID])

print("TEST TEST TEST!!!")
#wipeEdge(SG, "86", "95") #L2
#wipeEdge(SG, "113", "95") #L2
#wipeEdge(SG, "69", "95") #L2
#wipeEdge(SG, "105", "90") #L2
wipeEdge(SG, "69", "73") #L3
wipeEdge(SG, "105", "69") #L3
raw_input("HELLO!  Finished wiping edges for test purposes!")

#END TEST TEST TEST


#for nodeID in ["69", "73", "86", "90", "95", "105", "110", "113"]:
#    theoreticalSpectra = parseMGF.genTheoreticalMasses(SG.nodes[nodeID].CORRECT_ANSWER.massList)
#    print("\nNODE: " + str(SG.nodes[nodeID]))
#    peaksAccountedFor = 0
#    peaksNotAccountedFor = 0
#    for mass in SG.nodes[nodeID].allSpectra[0][1]:
#        accountedFor = False
#        for otherMass in theoreticalSpectra:
#            if abs(mass - otherMass) < .02:
#                accountedFor = True
#                break
#        if accountedFor:
#            peaksAccountedFor += 1
#        else:
#            peaksNotAccountedFor += 1
#
#    print("Peaks Accounted For: " + str(peaksAccountedFor))
#    print("Peaks Not Accounted For: " + str(peaksNotAccountedFor))
#    print("Total Peaks: " + str(peaksAccountedFor + peaksNotAccountedFor))
#
#
#    spectralPlotting.lineSpectrumMatch(
#        SG.nodes[nodeID].allSpectra[0][0]["PEPMASS"],
#        SG.nodes[nodeID].allSpectra[0][1],
#        theoreticalSpectra
#    )

raw_input("Finished Plotting Theoretical Vs Experimental Spectra")

#maximumExperimentalPairwiseScore = 0
#maximumTheoreticalPairwiseScore = 0
#for nodeID in ["69", "73", "86", "90", "95", "105", "110", "113"]:
#    node = SG.nodes[nodeID]
#    for neighbor in SG.nodes[nodeID].neighbors:
#        if node.mass > neighbor.mass:
#            continue
#        nodeExperimentalSpectraList = list(parseMGF.findMasses(x[1]) for x in node.allSpectra)
#        neighborExperimentalSpectraList = list(parseMGF.findMasses(x[1]) for x in neighbor.allSpectra)
#        massDeltas = spectralPlotting.computeDeltaMassDicts(nodeExperimentalSpectraList, neighborExperimentalSpectraList)
#
#        candidateA = node.CORRECT_ANSWER
#        candidateB = neighbor.CORRECT_ANSWER
#        candidateA.theoreticalSpectra = parseMGF.genTheoreticalMasses(candidateA.massList)
#        candidateB.theoreticalSpectra = parseMGF.genTheoreticalMasses(candidateB.massList)
#
#        #Remember that we are creating B-A, must be consistent for node and candidate.
#        candidateDeltas = spectralPlotting.computeDeltaMassDict(candidateA.theoreticalSpectra, candidateB.theoreticalSpectra)
#
#        experimentalTotal = 0
#        for key in massDeltas[0]:
#            count = massDeltas[0][key]
#            if abs(key - (neighbor.mass - node.mass)) < .02:
#                print("Experimental: " + str(key) + " Count: " + str(count))
#                experimentalTotal += count
#        maximumExperimentalPairwiseScore += experimentalTotal
#
#        theoreticalTotal = 0
#        for key in candidateDeltas:
#            count = candidateDeltas[key]
#            if abs(key - (neighbor.mass - node.mass)) < .02:
#                print("Theoretical:  " + str(key) + " Count: " + str(count))
#                theoreticalTotal += count
#        maximumTheoreticalPairwiseScore += theoreticalTotal
#
#        print("Node: " + str(node) + " Neighbor: " + str(neighbor))
#        print("Exp Total: " + str(experimentalTotal) + " Theoretical Total: " + str(theoreticalTotal))
#
#        cyclicLDist = str(L.computeCyclicEditDist(node.CORRECT_ANSWER, neighbor.CORRECT_ANSWER))
#        spectralPlotting.scatterSpectras(node.CORRECT_ANSWER.genSimpleName() + ":" + str(node.mass) + " - " + neighbor.CORRECT_ANSWER.genSimpleName() + ":" + str(neighbor.mass) + "(" + cyclicLDist + ")" + " n*m=" + str(len(massDeltas)), massDeltas, str(neighbor.mass - node.mass), candidateDeltas)
#
#print("Maximum experimental pairwise score: " + str(maximumExperimentalPairwiseScore))
#print("Maximum theoretical pairwise score: " + str(maximumTheoreticalPairwiseScore))

raw_input("Finished plotting pairwise diffs")

#print("Hacky print out of third file")
#eLookup = {}
#for e in spectralGraphEdges:
#    if e["A"] not in eLookup:
#        eLookup[e["A"]] = {}
#    if e["B"] not in eLookup:
#        eLookup[e["B"]] = {}
#    eLookup[e["A"]][e["B"]] = e
#    eLookup[e["B"]][e["A"]] = e
#
#filename = "edges.tsv"
#with open(filename, 'w') as f:
#    f.write("edge\tdeltaMass\tcosine\tcyclicLDist\tsequences\n")
#    for key in SG.nodes:
#        node = SG.nodes[key]
#        for neighbor in node.neighbors:
#            if hasattr(node, "CORRECT_ANSWER") and hasattr(neighbor, "CORRECT_ANSWER") and node.mass >= neighbor.mass:
#                f.write(str(node.mass)+"-"+str(neighbor.mass))
#                f.write("\t")
#                f.write(eLookup[node.ID][neighbor.ID]["deltaMass"])
#                f.write("\t")
#                f.write(eLookup[node.ID][neighbor.ID]["cosine"])
#                f.write("\t")
#                f.write(str(L.computeCyclicEditDist(node.CORRECT_ANSWER, neighbor.CORRECT_ANSWER)))
#                f.write("\t")
#                cA = node.CORRECT_ANSWER
#                for i in range(len(cA.name)):
#                    amino = cA.name[i]
#                    f.write(amino[0:3])
#                f.write("-")
#                cB = neighbor.CORRECT_ANSWER
#                for i in range(len(cB.name)):
#                    amino = cB.name[i]
#                    f.write(amino[0:3])
#                f.write("\n")
#


raw_input("Finished!")

#for key in ["69", "73", "86", "90", "95", "105", "110", "113"]:
#    SG.nodes[key].candidates = [SG.nodes[key].CORRECT_ANSWER]
#Ugh, half the edges in this spectral graph are wrong.  Let's see what happens if we remove them...

#Debug the graph created
print(SG)
plotSpectralGraph(SG)

#STEP 1: Remove Duplicate Candidates -- These are candidates within a node that have an identical spin invariant canonical form.
raw_input("Ready Step 1: Remove Duplicate Candidates\n")
#Remove duplicates from each node
for nodeKey in SG.nodes:
    node = SG.nodes[nodeKey]
    duplicates = []
    uniques = []
    spinInvariantNames = set([])
    for candidate in node.candidates:
        if str(candidate.name) in spinInvariantNames:
            duplicates.append(candidate)
        else:
            uniques.append(candidate)
            spinInvariantNames.add(str(candidate.name))
    print(str(node) + ": Duplicates: " + str(len(duplicates)) + " Remaining: " + str(len(uniques)))
    node.candidates = uniques

#STEP 2: Remove candidates with score < candidateScoreThresh
raw_input("Ready Step 2: Filter Low Scoring Candidates\n")
if candidateScoreThresh > 0:
    for nodeKey in SG.nodes:
        node = SG.nodes[nodeKey]
        highScorers = []
        for candidate in node.candidates:
            if candidate.score >= candidateScoreThresh:
                highScorers.append(candidate)
        print(str(node) + ": Low Scorers: " + str(len(node.candidates) - len(highScorers))  + " Remaining: " + str(len(highScorers)))
        node.candidates = highScorers

#STEP 2: Remove Empty Nodes -- Nodes with no candidates supplied immediately prevent us from constructing a G Consistent Graph.  So we remove these nodes and solve the pieces of the graph that we can solve.
raw_input("Ready Step 2: Filter nodes without candidates\n")
#Remove nodes with no candidates from the graph
filteredNodes = {}
for nodeKey in SG.nodes:
    node = SG.nodes[nodeKey]
    if len(node.candidates) == 0:
        print("WARNING: UNLINKING NODE WITHOUT CANDIDATES: " + str(node))
        #Erk, empty node.
        for neighbor in node.neighbors:
            print("\t REMOVING EDGE: " + str(node) + " <-> " + str(neighbor))
            neighbor.neighbors.remove(node)
    else:
        filteredNodes[nodeKey] = node
SG.nodes = filteredNodes
plotSpectralGraph(SG)

#STEP 3: Build the Candidate Graph -- Here we identify indels between candidates of adjacent nodes.  These indels form the edges of a new graph existing between candidates.
raw_input("Ready Step 3: Add candidate edges\n")
relationCount = {}
relationCount2 = {}
for nodeKey in SG.nodes:
    node = SG.nodes[nodeKey]
    for neighbor in node.neighbors:
        if node.mass < neighbor.mass: #Assume A has higher mass
            continue
        L.findRelations2(node, neighbor)

print("Relation Counts: ")
for key in relationCount2:
    print(str(key) + " -> " + str(relationCount2[key]))

#STEP 4: As we are searching for G Consistent subgraphs in the candidate graph, we can iteratively remove any candidates which we can prove cannot be part of such a subgraph.  This is an optimization pre-pass and consists of repeated quick tests to filter out candidates.
raw_input("Ready Step 4: Filter candidate graph for nodes that cannot be in G Consistent Set (PREVENTS G-u CONSISTENT SEARCH LATER!)\n")

#Remove a candidate from its neighbors' adjacency lists.  Does not remove the candidate from the parent node's candidate list.
def unlinkCandidateEdges(c):
    for edge in c.neighbors:
        edge.follow(c).neighbors.remove(edge)

#This filter strips out candidates that lack necessary edges to candidates in surrounding nodes.  By repeatedly filtering until no more candidates can be filtered, we are left with a new guarantee on viability of the remaining candidates:
#   Any candidate that survives convergence of this filter is part of a connected component in the candidate graph that spans that candidate's parent node's connected component in the spectral graph.  In this way, running this filter to convergence is similar to running a BFS that makes the same guarantee.  Unlike a pure BFS for connected components, this filter also removes candidates with non-viable edges that are connected to viable connected components in the candidate graph.
def filterCandidatesBySurroundingEdges():
    candidatesToUnlink = []
    filteredCandidateDict = {}
    for nodeKey in SG.nodes:
        node = SG.nodes[nodeKey]
        filteredCandidates = []
        for candidate in node.candidates:
            if len(candidate.neighbors) < len(node.neighbors):
                candidatesToUnlink.append(candidate)
                continue
            nodeLevelNeighbors = set([])
            for candidateNeighborEdge in candidate.neighbors:
                candidateNeighbor = candidateNeighborEdge.follow(candidate)
                nodeLevelNeighbors.add(candidateNeighbor.parentNode.ID)
            if len(nodeLevelNeighbors) < len(node.neighbors):
                candidatesToUnlink.append(candidate)
                continue
            filteredCandidates.append(candidate)
        print(str(node) + ": " + str(len(node.candidates)) + " -> " + str(len(filteredCandidates)))
        filteredCandidateDict[nodeKey] = filteredCandidates

    for nodeKey in SG.nodes:
        node = SG.nodes[nodeKey]
        node.candidates = filteredCandidateDict[nodeKey]
    for candidate in candidatesToUnlink:
        unlinkCandidateEdges(candidate)

    return len(candidatesToUnlink)

if uConsistency == 0:
    removeMoreCandidates = True
    while removeMoreCandidates:
        removeMoreCandidates = (filterCandidatesBySurroundingEdges() > 0)
else:
    print("SKIPPING PRE FILTERING OPTIMIZATION: u =" + str(uConsistency))

#STEP 5:  Enumerate G Consistent Subgraphs in the candidate graph!  This is the tough one.  This is done as a DFS on the spectral graph SG while maintaining and updating an exponentially growing list (wavefront) of viable candidate sequences.  Each time an SGNode is visited, every candidate of the new node is checked to see if it can be added to the ends of any entries in the wavefront.  For every candidate that can be added to a wavefront, the wavefront forks: A new wavefront is generated with entries for each viable candidate.
raw_input("Ready Step 5: Enumerate G Consistent Candidate Subgraphs\n")

def DFSHelper(stack, visited):
    connectedComponent = []
    while len(stack) > 0:
        c = stack.pop()
        if c.ID in visited:
            continue
        visited.add(c.ID)
        connectedComponent.append(c)
        for neighbor in c.neighbors:
            if neighbor.ID not in visited:
                stack.append(neighbor)
    return connectedComponent

#Returns a lookup table from nodeID to set of nodeIDs representing each connected component.
def DFSNode(SG):
    connectedComponents = []
    stack = []
    visited = set([])
    for nodeKey in SG.nodes:
        node = SG.nodes[nodeKey]
        if node.ID in visited:
            continue
        stack.append(node)
        connectedComponents.append(DFSHelper(stack, visited))
    connectionDict = {}
    for cc in connectedComponents:
        componentSet = set([])
        for node in cc:
            componentSet.add(node.ID)
        for node in cc:
            connectionDict[node.ID] = componentSet
    return connectionDict

class CandidateNetwork:
    def __init__(self, score, candidates, u):
        self.score = score
        self.candidates = candidates
        self.remainingBrokenEdges = u

    def canAfford(self, brokenEdges):
        return self.remainingBrokenEdges >= brokenEdges

    def copyAndAppend(self, visitedNode, newCandidate, brokenEdges):
        newDictionary = dict(self.candidates)
        newDictionary[visitedNode.ID] = newCandidate
        return CandidateNetwork(self.score + newCandidate.score, newDictionary, self.remainingBrokenEdges - brokenEdges)

def fuzzyEqual(a, b, eps):
    return abs(b-a) < eps

def pairwiseScoreMultiAnnotation(S1, S2, T1, T2, deltaMass, eps):
    if deltaMass <= 0:
        raise Exception("Delta Mass Must Be > 0")

    #Define fuzzyEqual(a,b) as abs(b-a) < eps
    #Return count of all tuples (s1, s2, t1, t2) such that:
    #   s1 is an element of S1,
    #   s2 is an element of S2,
    #   t1 is an element of T1,
    #   t2 is an element of T2,
    #   fuzzyEqual(s2 - s1, deltaMass),
    #   fuzzyEqual(t2 - t1, deltaMass),
    #   fuzzyEqual(t1, s1) and
    #   fuzzyEqual(t2, s2)
    results = []
    sPairs = defaultdict(list)
    tPairs = defaultdict(list)
    for s1 in S1:
        for s2 in S2:
            if s2 <= s1:
                continue
            if fuzzyEqual(s2 - s1, deltaMass, eps):
                sPairs[s1].append(s2)
    for t1 in T1:
        for t2 in T2:
            if t2 <= t1:
                continue
            if fuzzyEqual(t2 - t1, deltaMass, eps):
                tPairs[t1].append(t2)
    for s1 in sPairs:
        for s2 in sPairs[s1]:
            for t1 in tPairs:
                for t2 in tPairs[t1]:
                    if fuzzyEqual(t1, s1, eps) and fuzzyEqual(t2, s2, eps):
                        results.append((s1, s2, t1, t2))
    return results

def pairwiseScoreSingleAnnotation(S1, S2, T1, T2, deltaMass, eps):
    #Define fuzzyEqual(a,b) as abs(b-a) < eps
    #Return count of all tuples (t1, t2) such that:
    #   t1 is an element of T1,
    #   t2 is an element of T2,
    #   fuzzyEqual(t2-t1, deltaMass)
    #   There exists s1 is an element of S1 such that fuzzyEqual(s1,t1)
    #   There exists s2 is an element of S2 such that fuzzyEqual(s2,t2)
    results = []
    tPairs = defaultdict(list)
    for t1 in T1:
        for t2 in T2:
            if t2 <= t1:
                continue
            if fuzzyEqual(t2 - t1, deltaMass, eps):
                tPairs[t1].append(t2)

    for t1 in tPairs:
        found = False
        temp = None
        for s1 in S1:
            if fuzzyEqual(t1, s1, eps):
                found = True
                temp = s1
                break
        if found:
            for t2 in tPairs[t1]:
                for s2 in S2:
                    if fuzzyEqual(t2, s2, eps):
                        results.append((temp, s2, t1, t2))
                        break

    return results

#By filtering candidates by their surrounding edges until we stopped removing candidates
#we are now guaranteed that all remaining candidates exist in connected components which span the containing connected components of their parent nodes in the graph.
#Now the awful exponential part: And remember that everything has to be done separately per connected component of the SG node graph.
def GConsistentSearch(rootNode, u):
    #Implement as a search on the SpectralGraph, not the candidate graph.
    #Maintain a wavefront (list) of all sequences of candidates up to the current node.
    #When visiting a node, for each sequence identify which candidates in this node are possible   Generate a new sequence for each of these new candidates and update.  Storing as a tree would save space if you want.

    wavefront = []
    for c in rootNode.candidates:
        wavefront.append(CandidateNetwork(c.score, {c.parentNode.ID: c}, u))
    stack = []
    for neighbor in rootNode.neighbors:
        stack.append(neighbor)
    visited = set([rootNode.ID])
    visitOrder = [rootNode]

    print("ROOTED AT: " + str(rootNode))
    print("INITIAL FRONT: " + str(len(wavefront)))

    while len(stack) > 0:
        n = stack.pop()
        if n.ID in visited:
            continue
        print("VISITING: " + str(n))
        visitOrder.append(n)
        visited.add(n.ID)
        nextWavefront = []
        
        neighboringCandidatesLookup = {}
        neighboringCandidatesList = []
        for front in wavefront:
            neighboringCandidates = []
            for visitedNeighbor in n.neighbors:
                if visitedNeighbor.ID not in visited:
                    continue
                neighboringCandidates.append(front.candidates[visitedNeighbor.ID])

            #Optimization: Due to the sparseness of the spectral graph (and the redundancy of the wavefront), many fronts will have the same neighboring candidate set,
            #pull these sets out and produce a lookup back to the fronts so we can
            #compute against neighboring candidate sets instead of against fronts.
            neighboringCandidatesStr = str(neighboringCandidates)
            if neighboringCandidatesStr not in neighboringCandidatesLookup:
                neighboringCandidatesLookup[neighboringCandidatesStr] = []
                neighboringCandidatesList.append(neighboringCandidates)
            neighboringCandidatesLookup[neighboringCandidatesStr].append(front)
        
        if u == 0:
            #Version 1:  Full G Consistency.  Check each set of neighboring candidates against the edges in each candidate of node.  If the candidate has edges to every
            #neighbor in the set, it is valid, and all fronts associated with this candidate set can fork and add this new candidate.
            for neighboringCandidates in neighboringCandidatesList:
                for candidate in n.candidates:
                    isValidCandidate = True
                    for nc in neighboringCandidates:
                        if nc.ID not in candidate.neighborSet: #TODO FIXME HACK:  Should I keep this scratch space neighborSet in a separate dictionary?
                            isValidCandidate = False
                    if isValidCandidate:
                        for front in neighboringCandidatesLookup[str(neighboringCandidates)]:
                            nextWavefront.append(front.copyAndAppend(n, candidate, 0))
        else:
            #Version 2:  G-u Consistency.  Check each set of neighboring candidates against the edges in each candidate of node.  Count how many edges the candidate has
            #to neighbors in the set.  For each front associated with this candidate set, check if it has enough slack leftover to fork to this new candidate, and if so, do so.
            for neighboringCandidates in neighboringCandidatesList:
                for candidate in n.candidates:
                    brokenEdges = 0
                    for nc in neighboringCandidates:
                        if nc.ID not in candidate.neighborSet:
                            brokenEdges += 1
                    for front in neighboringCandidatesLookup[str(neighboringCandidates)]:
                        if front.canAfford(brokenEdges):
                            nextWavefront.append(front.copyAndAppend(n, candidate, brokenEdges))

        for neighborNode in n.neighbors:
            if neighborNode.ID not in visited:
                stack.append(neighborNode)
        wavefront = nextWavefront
        print("NEW FRONT: " + str(len(wavefront)))

    return (visitOrder, wavefront)

handled = set([])

#set up scratch (Yuck, python lets me use fields in arbitrary objects for scratch space...  How to reorganize this?)
for nodeKey in SG.nodes:
    node = SG.nodes[nodeKey]
    for c in node.candidates:
        c.neighborSet = set([])
        for candidateEdge in c.neighbors:
            cn = candidateEdge.follow(c)
            c.neighborSet.add(cn.ID)

for nodeKey in SG.nodes:
    node = SG.nodes[nodeKey]
    for c in node.candidates:
        c.theoreticalSpectra = None
        for candidateEdge in c.neighbors:
            candidateEdge.pairwiseScore = None

#Run a separate GConsistent search for each connected component of the spectral graph, possibilities for each of these are independent and should be treated as such.
SGGroups = DFSNode(SG)
answers = []
for nodeKey in SG.nodes:
    node = SG.nodes[nodeKey]
    if node.ID in handled:
        continue
    group = SGGroups[node.ID]
    for nodeID in group:
        handled.add(nodeID)
    print("SEARCHING GROUP (REPRESENTATIVE #" + node.ID + ") COUNT: " + str(len(group)))
    allSequences = GConsistentSearch(node, uConsistency)
    answers.append(allSequences)

    #Calculate metadata about each answer

    #Compute total edit distance
    for candidateNetwork in allSequences[1]:
        cyclicLevenshteinSum = 0
        for nodeID in group:
            for neighborNode in SG.nodes[nodeID].neighbors:
                if nodeID < neighborNode.ID:
                    continue
                candidateA = candidateNetwork.candidates[nodeID]
                candidateB = candidateNetwork.candidates[neighborNode.ID]
                for edge in candidateA.neighbors:
                    if edge.follow(candidateA) is candidateB:
                        cyclicLevenshteinSum += edge.editDist
        candidateNetwork.totalEditDistance = cyclicLevenshteinSum

    #Compute theoretical spectra and pairwise score
    for candidateNetwork in allSequences[1]:
        totalPairwiseScore = 0
        totalPairwiseScore2 = 0
        totalPairwiseScore3 = 0
        
        for nodeID in group:
            node = SG.nodes[nodeID]
            for neighborNode in SG.nodes[nodeID].neighbors:
                if node.mass >= neighborNode.mass:
                    continue
                candidateA = candidateNetwork.candidates[nodeID]
                candidateB = candidateNetwork.candidates[neighborNode.ID]
                edgeAB = None
                for edge in candidateA.neighbors:
                    if edge.follow(candidateA) is candidateB:
                        edgeAB = edge
                        break
                        
                nodeExperimentalSpectra = list(parseMGF.findMasses(x[1]) for x in node.allSpectra)[0]
                neighborExperimentalSpectra = list(parseMGF.findMasses(x[1]) for x in neighborNode.allSpectra)[0]

                if edgeAB.pairwiseScore is None:
                    if candidateA.theoreticalSpectra is None:
                        candidateA.theoreticalSpectra = parseMGF.genTheoreticalMasses(candidateA.massList)
                    if candidateB.theoreticalSpectra is None:
                        candidateB.theoreticalSpectra = parseMGF.genTheoreticalMasses(candidateB.massList)
                    
                    #Remember that we are creating B-A, must be consistent for node and candidate.
                    massDeltaValue = neighborNode.mass - node.mass

                    #Version 1:  Score by overlap in deltas between experimental and theoretical spectra
                    s1s2t1t2 = pairwiseScoreMultiAnnotation(
                        nodeExperimentalSpectra,
                        neighborExperimentalSpectra,
                        candidateA.theoreticalSpectra,
                        candidateB.theoreticalSpectra,
                        massDeltaValue,
                        .02)
                    edgeAB.pairwiseScore = len(s1s2t1t2)

                    #Version 2:  Score by intersection of
                    t1t2 = pairwiseScoreSingleAnnotation(
                        nodeExperimentalSpectra,
                        neighborExperimentalSpectra,
                        candidateA.theoreticalSpectra,
                        candidateB.theoreticalSpectra,
                        massDeltaValue,
                        .02)
                    edgeAB.pairwiseScore2 = len(t1t2)
                    
                    if len(s1s2t1t2) != len(t1t2):
                        print("Oops! \nQuads: " + str(s1s2t1t2))
                        print("Pairs: " + str(t1t2))
                    
                        
                    #Version 3:  Score only by theoretical spectra
                    candidateDeltas = spectralPlotting.computeDeltaMassDict(candidateA.theoreticalSpectra, candidateB.theoreticalSpectra)
                    edgeAB.pairwiseScore3 = spectralPlotting.associateMassDeltasEps(candidateDeltas, massDeltaValue, .02)
                
                totalPairwiseScore += edgeAB.pairwiseScore
                totalPairwiseScore2 += edgeAB.pairwiseScore2
                totalPairwiseScore3 += edgeAB.pairwiseScore3

        candidateNetwork.pairwiseScore = totalPairwiseScore
        candidateNetwork.pairwiseScore2 = totalPairwiseScore2
        candidateNetwork.pairwiseScore3 = totalPairwiseScore3

    #Double check results before sending to bahar :D
    totalEdgesSG = 0
    for nodeID in group:
        totalEdgesSG += len(SG.nodes[nodeID].neighbors)

    for candidateNetwork in allSequences[1]:
        candidateSet = set([])
        totalEdgesCounted = 0
        for nodeID in group:
            candidateSet.add(candidateNetwork.candidates[nodeID].ID)
        for nodeID in group:
            for candidateEdge in candidateNetwork.candidates[nodeID].neighbors:
                neighbor = candidateEdge.follow(candidateNetwork.candidates[nodeID])
                if neighbor.ID in candidateSet:
                    totalEdgesCounted += 1
        if totalEdgesSG != totalEdgesCounted + 2 * (uConsistency - candidateNetwork.remainingBrokenEdges):
            print Exception("OH NOES! Expected: " + str(totalEdgesSG) + " But Got: " + str(totalEdgesCounted) + " + " + str(candidateNetwork.remainingBrokenEdges))
    print("OH YEAH!")

#Write search results to a file, sort by total candidate scores.
for answer in answers:
    scoreX = []
    pairwiseScoreY = []
    pairwiseScore2Y = []
    pairwiseScore3Y = []
    correctRow = None
    
    representativeNodeID = answer[0][0].ID
    for node in answer[0]:
        if node.ID < representativeNodeID:
            representativeNodeID = node.ID

    foundAnswer = False
    filename = "group_" + representativeNodeID + ".cdt"
    with open(filename, 'w') as f:
        #Write header masses
        for node in answer[0]:
            f.write(str(node.mass))
            f.write("\t")
        f.write("score\t")
        f.write("pairwiseScore\t")
        f.write("totalEditDist\t")
        f.write("correct\n")

        #Write rows
        rowNumber = 0
        answer[1].sort(key=lambda x: -x.score)
        for row in answer[1]:
            rowNumber += 1
            correctCount = 0
            for node in answer[0]:
                candidate = row.candidates[node.ID]
                if candidate.name == node.CORRECT_ANSWER.name:
                    correctCount += 1
                if OUTPUT == "mass":
                    for i in range(len(candidate.massList)-1):
                        amino = str(candidate.massList[i])
                        f.write(amino + "-")
                    f.write(str(candidate.massList[-1]))
                if OUTPUT == "amino":
                    f.write(str(candidate.name))
                if OUTPUT == "filename":
                    f.write(str(candidate.graphFile))
                if OUTPUT == "short":
                    for i in range(len(candidate.name)-1):
                        amino = candidate.name[i]
                        f.write(amino[0:3] + "-")
                    f.write(candidate.name[-1][0:3])
                    f.write(":S"+str(int(candidate.score)) + ":R" + str(int(candidate.rank)))
                f.write("\t")
            f.write(str(row.score))
            f.write("\t")
            f.write(str(row.pairwiseScore))
            f.write("\t")
            f.write(str(row.totalEditDistance))
            f.write("\t")
            if correctCount == len(answer[0]):
                print("FILE: " + filename + " CORRECT ANSWER: ROW: " + str(rowNumber) + " SCORE: " + str(row.score) + " EDIT DIST: " + str(row.totalEditDistance))
                correctRow = row
                f.write("*")
                foundAnswer = True
            else:
                f.write(str(correctCount))

            f.write("\n")
            
            scoreX.append(row.score)
            pairwiseScoreY.append(row.pairwiseScore)
            pairwiseScore2Y.append(row.pairwiseScore2)
            pairwiseScore3Y.append(row.pairwiseScore3)

    if not foundAnswer:
        print("ANSWER NOT FOUND!!!")

    plt.scatter(scoreX, pairwiseScoreY)
    plt.title(filename)
    plt.xlabel("Score")
    plt.ylabel("Pairwise Score")
    if correctRow is not None:
        plt.axvline(x=correctRow.score, c='k', linestyle='dashed')
        plt.axhline(y=correctRow.pairwiseScore, c='k', linestyle='dashed')
        plt.axhline(y=0, c='k')
    plt.show()

    plt.scatter(scoreX, pairwiseScore2Y)
    plt.title(filename)
    plt.xlabel("Score")
    plt.ylabel("Pairwise Score 2")
    if correctRow is not None:
        plt.axvline(x=correctRow.score, c='k', linestyle='dashed')
        plt.axhline(y=correctRow.pairwiseScore2, c='k', linestyle='dashed')
        plt.axhline(y=0, c='k')
    plt.show()

    plt.scatter(scoreX, pairwiseScore3Y)
    plt.title(filename)
    plt.xlabel("Score")
    plt.ylabel("Pairwise Score 3")
    if correctRow is not None:
        plt.axvline(x=correctRow.score, c='k', linestyle='dashed')
        plt.axhline(y=correctRow.pairwiseScore3, c='k', linestyle='dashed')
        plt.axhline(y=0, c='k')
    plt.show()
