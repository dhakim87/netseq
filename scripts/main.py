__author__ = "Daniel Hakim, PevznerLab"

from collections import defaultdict
from parsers import *
from eq import fuzzyEqual
#import networkx as nx
import matplotlib.pyplot as plt
import levenshtein as Lev
import spectralPlotting
import sys
import globals
import scoring
import MST
from AminoLanguage import *
from SpectralGraph import *

SHOW_GRAPHS = False  #Whether or not to show graphs with networkx and matplotlib (VERY SLOW)
PRINT_CANDIDATES = False #Whether or not to print out the candidates at various points in the algorithm (Useful for debugging but a ton of text)
MAX_MUTATIONS_PER_EDGE = 2 #How many mutations we allow per edge in the network graph, this determines where we add edges in the candidate graph.  (>1 VERY SLOW)
DEFAULT_U_CONSISTENCY_THRESH = 0 #Override: U=xxx (>0 VERY SLOW)
DEFAULT_CANDIDATE_SCORE_THRESH = 0 #Override: MIN_SCORE=xxx

OUTPUT = "short" #Supports "mass", "amino", "filename", or "short"

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

"""A Levenshtein distance based method for adding edges pairwise between candidates in two nodes."""
def findRelations2(SGNodeA, SGNodeB, maxMutationsPerEdge):
    i = 0
    for candA in SGNodeA.candidates:
        spins = L.generateAllCanonicalNameSpins(candA.name)
        i += 1
        print(str(i) + " / " + str(len(SGNodeA.candidates)))
        for candB in SGNodeB.candidates:
            if not L.couldBeRelated(candA.counts, candB.counts, maxMutationsPerEdge):
                continue
            
            closestDist = max(len(candA.name), len(candB.name), maxMutationsPerEdge + 1)
            for candASpin in spins:
                dist = Lev.compute2(candASpin, candB.name, maxMutationsPerEdge)
                if dist is not None:
                    closestDist = min(closestDist, dist)
            if closestDist <= maxMutationsPerEdge:
                edge = CandidateEdge(candA, candB, closestDist)
                candA.neighbors.append(edge)
                candB.neighbors.append(edge)

#MAIN PROGRAM START
print("\n\n---START---\n\n");

#Read command line args
args = {}
for arg in sys.argv[1:]:
    ss = arg.split("=")
    if len(ss) == 2:
        args[ss[0]] = ss[1]

L = globals.DEFAULT_ALPHABET

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
SG.retrieveCandidates("../data/surugamide/", L)

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
        findRelations2(node, neighbor, MAX_MUTATIONS_PER_EDGE)

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
        candidateEdgeList = []
        
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
                    s1s2t1t2 = scoring.pairwiseScoreMultiAnnotation(
                        nodeExperimentalSpectra,
                        neighborExperimentalSpectra,
                        candidateA.theoreticalSpectra,
                        candidateB.theoreticalSpectra,
                        massDeltaValue,
                        .02)
                    edgeAB.pairwiseScore = len(s1s2t1t2)

                    #Version 2:  Score by intersection of
                    t1t2 = scoring.pairwiseScoreSingleAnnotation(
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
                candidateEdgeList.append(edgeAB)

        candidateNetwork.pairwiseScore = totalPairwiseScore
        candidateNetwork.pairwiseScore2 = totalPairwiseScore2
        candidateNetwork.pairwiseScore3 = totalPairwiseScore3
        candidateNetwork.candidateEdges = candidateEdgeList

        #Using the set of edges in the candidate network, construct an Maximum Spanning Tree using pairwiseScore as the edge weight.
        maximalSpanningTree = MST.spanningTree(candidateEdgeList, lambda x:x.pairwiseScore, maximal=True)

        mstScore = 0
        for edge in maximalSpanningTree:
            mstScore += edge.pairwiseScore
        candidateNetwork.mstScore = mstScore
        candidateNetwork.mstEdges = maximalSpanningTree

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
    mstScoreY = []
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
        f.write("mstScore\t")
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
            f.write(str(row.mstScore))
            f.write("\t")
            f.write(str(row.totalEditDistance))
            f.write("\t")
            if correctCount == len(answer[0]):
                print("FILE: " + filename + " CORRECT ANSWER: ROW: " + str(rowNumber) + " SCORE: " + str(row.score) + " EDIT DIST: " + str(row.totalEditDistance) + " PAIR SCORE: " + str(row.pairwiseScore) + " MST SCORE: " + str(row.mstScore))
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
            mstScoreY.append(row.mstScore)

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

    plt.scatter(scoreX, mstScoreY)
    plt.title(filename)
    plt.xlabel("Score")
    plt.ylabel("MST Score")
    if correctRow is not None:
        plt.axvline(x=correctRow.score, c='k', linestyle='dashed')
        plt.axhline(y=correctRow.mstScore, c='k', linestyle='dashed')
        plt.axhline(y=0, c='k')
    plt.show()



    for e in correctRow.candidateEdges:
        print(str(e.A.ID) + " -> " + str(e.B.ID) + ": " + str(e.pairwiseScore))

    print("------------------------")
    for e in correctRow.mstEdges:
        print(str(e.A.ID) + " -> " + str(e.B.ID) + ": " + str(e.pairwiseScore))


