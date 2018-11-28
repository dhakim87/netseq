from AminoLanguage import *

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
