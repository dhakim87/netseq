from collections import defaultdict
from globals import DEFAULT_ALPHABET
from eq import fuzzyEqual


#An immutable grab bag of Aminos.  Usable as keys in dictionaries
class Composition:
    #Factory Methods, these construct new Composition instances
    @classmethod
    def fromMassList(clazz, massList, alphabet=DEFAULT_ALPHABET):
        c = Composition()
        aminoStr = alphabet.toCanonicalAminoArr(massList)
        for amino in aminoStr:
            c.aminoCounts[amino] += 1
        c.finish()
        return c

    @classmethod
    def fromAminoArr(clazz, aminoArr):
        c = Composition()
        for amino in aminoArr:
            c.aminoCounts[amino] += 1
        c.finish()
        return c
        
    @classmethod
    def fromCountsDict(clazz, aminoCounts):
        c = Composition()
        for amino in aminoCounts:
            c.aminoCounts[amino] = aminoCounts[amino]
        c.finish()
        return c


    #Actual constructor
    def __init__(self):
        self.aminoCounts = defaultdict(int)

    #Methods
    def calcMass(self):
        return Composition.calcMassOnAminoCounts(self.aminoCounts)
    
    @staticmethod
    def calcMassOnAminoCounts(aminoCounts):
        totalMass = 0
        for amino in aminoCounts:
            totalMass += amino.mass * aminoCounts[amino]
        return totalMass

    """Returns a new composition resulting from following the given transition"""
    def followTransition(self, transition):
        newComp = Composition()
        for key in self.aminoCounts:
            newComp.aminoCounts[key] = self.aminoCounts[key]
        if transition.action == 0: #ADD
            newComp.aminoCounts[transition.amino1] += 1
        if transition.action == 1: #DEL
            if newComp.aminoCounts[transition.amino1] <= 0:
                raise Exception("Invalid Transition: Cannot delete " + str(transition.amino1))
            newComp.aminoCounts[transition.amino1] -= 1
            if newComp.aminoCounts[transition.amino1] == 0:
                del newComp.aminoCounts[transition.amino1]
        if transition.action == 2: #SWP
            if newComp.aminoCounts[transition.amino1] <= 0:
                raise Exception("Invalid Transition: Cannot swp " + str(transition.amino1))
            newComp.aminoCounts[transition.amino1] -= 1
            newComp.aminoCounts[transition.amino2] += 1
            if newComp.aminoCounts[transition.amino1] == 0:
                del newComp.aminoCounts[transition.amino1]

        newComp.finish()
        return newComp
        
    """Optimized version of followTransition, no error checking, modifies an aminoCounts dictionary as Composition instances are immutable"""
    @staticmethod
    def followTransitionInPlace(aminoCounts, transition):
        if transition.action == 0: #ADD
            aminoCounts[transition.amino1] += 1
        if transition.action == 1: #DEL
            aminoCounts[transition.amino1] -= 1
            if aminoCounts[transition.amino1] == 0:
                del aminoCounts[transition.amino1]
        if transition.action == 2: #SWP
            aminoCounts[transition.amino1] -= 1
            aminoCounts[transition.amino2] += 1
            if aminoCounts[transition.amino1] == 0:
                del aminoCounts[transition.amino1]

    """Enumerate all paths from the starting composition to compositions within tolerance of the target mass of length less than or equal to maximumLDist"""
    @staticmethod
    def generatePaths(startingComposition, targetMass, maximumLDist, massEpsilon, alphabet = DEFAULT_ALPHABET):
        startingMass = startingComposition.calcMass()
        
        outputComps = []
        outputPaths = []
        
        if fuzzyEqual(startingMass, targetMass, massEpsilon):
            outputPaths.append([])
            outputComps.append(startingComposition)

        #TODO FIXME HACK: Would this be noticeably faster as a dynamic programming implementation with some backtracing step or does it all work out basically the same?  I suppose it depends on how many duplicate ways there are to get to the same composition.  Probably doesn't matter much until we increase L dist above 2.
        def genPathsHelper(currentAminoCounts, transitionList, targetMass, remainingLDist, massEpsilon, alphabet, outputPaths, outputComps):
            if remainingLDist == 0:
                return
            validTransitions = Transition.generateValidTransitions(currentAminoCounts, alphabet)
            for transition in validTransitions:
                newComp = defaultdict(int)
                for key in currentAminoCounts:
                    newComp[key] = currentAminoCounts[key]
                Composition.followTransitionInPlace(newComp, transition)
                newTransitionList = transitionList + [transition]
                if fuzzyEqual(Composition.calcMassOnAminoCounts(newComp), targetMass, massEpsilon):
                    outputPaths.append(newTransitionList)
                    outputComps.append(Composition.fromCountsDict(newComp))
                genPathsHelper(newComp, newTransitionList, targetMass, remainingLDist - 1, massEpsilon, alphabet, outputPaths, outputComps)

        startingCounts = defaultdict(int)
        for key in startingComposition.aminoCounts:
            startingCounts[key] = startingComposition.aminoCounts[key]
        
        genPathsHelper(startingCounts, [], targetMass, maximumLDist, massEpsilon, alphabet, outputPaths, outputComps)
        
        outputComps = alphabet.removeStrDups(outputComps)
        outputPaths = startingComposition.filterPaths(outputPaths, targetMass, massEpsilon) #Annoying, but without this we can transition through amino acids that aren't in the initial composition: [SWP(Ile, Gly), SWP(Gly, Ala)] is a dumb path if Gly isn't in the initial sequence.
        return (outputPaths, outputComps)
    
    def filterPaths(self, possiblePaths, targetMass, massEpsilon):
        filtered = []
        startMass = self.calcMass()
        
        for path in possiblePaths:
            requiredCount = defaultdict(int)
            curMass = startMass
            for transition in path:
                curMass += transition.deltaMass()
                if transition.action == 0:
                    requiredCount[transition.amino1] -= 1
                if transition.action == 1:
                    requiredCount[transition.amino1] += 1
                if transition.action == 2:
                    requiredCount[transition.amino1] += 1
                    requiredCount[transition.amino2] -= 1
            
            #In order to transition through some amino acid you must actually have it in the initial sequence, otherwise its pointless.
            for key in requiredCount:
                if requiredCount[key] == 0:
                    requiredCount[key] = 1
            isValid = True
            for key in requiredCount:
                if self.aminoCounts[key] < requiredCount[key]:
                    isValid = False
            if not fuzzyEqual(curMass, targetMass, massEpsilon):
                isValid = False
            
            if isValid:
                filtered.append(path)
        return filtered
    
    """Enumerate all paths from any composition with the starting mass to compositions within tolerance of the target mass of length less than or equal to maximumLDist.  """
    @staticmethod
    def generatePossiblePaths(startingMass, targetMass, maximumLDist, massEpsilon, alphabet = DEFAULT_ALPHABET):
        
        outputPaths = []
        if fuzzyEqual(startingMass, targetMass, massEpsilon):
            outputPaths.append([])

        def genPathsHelper(validTransitions, curMass, transitionList, targetMass, remainingLDist, massEpsilon, alphabet, outputPaths):
            if remainingLDist == 0:
                return
            for transition in validTransitions:
                newMass = curMass + transition.deltaMass()
                newTransitionList = transitionList + [transition]
                if fuzzyEqual(newMass, targetMass, massEpsilon):
                    outputPaths.append(newTransitionList)
                genPathsHelper(validTransitions, newMass, newTransitionList, targetMass, remainingLDist - 1, massEpsilon, alphabet, outputPaths)

        fakeComp = Composition.fromAminoArr(alphabet.sortedAminos)
        allPossibleTransitions = Transition.generateValidTransitions(fakeComp.aminoCounts, alphabet)
        genPathsHelper(allPossibleTransitions, startingMass, [], targetMass, maximumLDist, massEpsilon, alphabet, outputPaths)
        return outputPaths

    def finish(self):
        list = []
        for key in sorted(self.aminoCounts.keys()):
            list.append((key, self.aminoCounts[key]))
        strVal = str(list)
        self.ID = strVal

    def __repr__(self):
        return self.ID
    
    def __str__(self):
        return self.ID
        
    def __hash__(self):
        return hash(self.ID)

    def __eq__(self, other):
        if not isinstance(self, Composition) or not isinstance(other, Composition):
            return False
        return self.ID == other.ID


class Transition:
    #ADD - 0
    #DEL - 1
    #SWP - 2

    def __init__(self, action, amino1, amino2 = None):
        self.action = action
        if action == 0: #ADD
            self.actionName = "ADD" #ADD(amino1)
        elif action == 1:
            self.actionName = "DEL" #DEL(amino1)
        elif action == 2:
            self.actionName = "SWP" #SWP(amino1 -> amino2)
        else:
            raise Exception("Unknown Action: " + str(action))
        self.amino1 = amino1
        self.amino2 = amino2
        
    def getInverse(self):
        if self.action == 0 or self.action == 1:
            return Transition(1-action, self.amino1, None)
        return Transition(self.action, self.amino2, self.amino1)
        
    def deltaMass(self):
        if self.action == 0:
            return self.amino1.mass
        if self.action == 1:
            return -self.amino1.mass
        if self.action == 2:
            return self.amino2.mass - self.amino1.mass

    def __str__(self):
        if self.amino2 == None:
            return self.actionName + "(" + str(self.amino1) + ")"
        return self.actionName + "(" + str(self.amino1) + ", " + str(self.amino2) + ")"
        
    def __repr__(self):
        if self.amino2 == None:
            return self.actionName + "(" + str(self.amino1) + ")"
        return self.actionName + "(" + str(self.amino1) + ", " + str(self.amino2) + ")"

    @staticmethod
    def generateValidTransitions(aminoCounts, alphabet = DEFAULT_ALPHABET):
        transitions = []
        #ADD
        for toAdd in alphabet.sortedAminos:
            transitions.append(Transition(0, toAdd))
        #DEL
        for amino in aminoCounts:
            if aminoCounts[amino] > 0:
                transitions.append(Transition(1, amino))
        #SWP
        for amino in aminoCounts:
            if aminoCounts[amino] > 0:
                for swapTo in alphabet.sortedAminos:
                    if amino == swapTo:
                        continue
                    transitions.append(Transition(2, amino, swapTo))
        return transitions



#Test
#allPaths = Composition.generatePaths(Composition.fromMassList([143.058, 165.042]), 114.0428, 2, 0.02)
#print(allPaths)
