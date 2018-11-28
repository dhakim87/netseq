from collections import defaultdict
import levenshtein as Lev


"""An Amino Acid"""
class Amino:
    def __init__(self, tsvRow):
        self.tsvRow = tsvRow
        self.name = tsvRow["uc"]
        self.mass = float(tsvRow["mass"])

    def __str__(self):
        return self.name;

    def __repr__(self):
        return self.name;

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
    compositionSet = {} #Allows interning of AminoCompositions for faster lookup and equality
    
    def __init__(self, sortedAminos):
        self.sortedAminos = sortedAminos
        self.fastLookup = {} #Caches and returns results of the toCanonicalName function
        self.fastLookupAmino = {}

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
        
        self.fastLookup[mass] = s
        return s
        
    def toCanonicalAmino(self, mass):
        val = None
        if mass in self.fastLookupAmino:
            return self.fastLookupAmino[mass]

        for amino in self.sortedAminos:
            if mass > amino.mass + AminoLanguage.MASS_EPSILON:
                continue
            if mass < amino.mass - AminoLanguage.MASS_EPSILON:
                break
            if val == None:
                val = amino
            else:
                print("Warning, mass mapped to multiple amino acids, picking " + val.name)
        
        if val == None:
            raise Exception("Unknown Amino Weight: " + str(mass))
        
        self.fastLookupAmino[mass] = val
        return val

    """The canonical name of an array of masses is an array of strings that each uniquely represent amino acids within MASS_EPSILON of each mass"""
    def toCanonicalNameArr(self, massList):
        result = []
        for m in massList:
            result.append(self.toCanonicalName(m))
        return result
    
    def toCanonicalAminoArr(self, massList):
        result = []
        for m in massList:
            result.append(self.toCanonicalAmino(m))

    """The set of canonical name array spins is what you get by rotating the strings in a canonical name array.  Because we are working with cyclopeptides, all of these spins are considered identical"""
    def generateAllCanonicalNameSpins(self, canonicalNameArr):
        allSpins = []
        for i in range(len(canonicalNameArr)):
            allSpins.append(canonicalNameArr[i:] + canonicalNameArr[:i])
        return allSpins

    """To quickly identify spins of the same cyclopeptide, all spin sets can be converted to a spinInvariantCanonicalNameArr.  This is the rotation of the cyclopeptide that results in the least value string, sorted lexicographically.  For maintaining this interface, all that matters is that all rotations of a canonical name arr are assigned the same string value, but cyclopeptides that are not considered equal must generate different spin invariant canonical names"""
    def toSpinInvariant(self, canonicalNameArr):
        if len(canonicalNameArr) <= 1:
            return list(canonicalNameArr)
        
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
