import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

"""Given two lists of floats, pairwise add mass diffs to a dictionary"""
def computeDeltaMassDict(spectra1, spectra2, precision=None):
    deltaMassDict = defaultdict(int)
    for s1 in spectra1:
        for s2 in spectra2:
            if precision is None:
                deltaMassDict[s2-s1] += 1
            else:
                deltaMassDict[round(s2-s1,precision)] += 1
    return deltaMassDict

def computeDeltaMassDicts(spectraList1, spectraList2, precision=None):
    allDeltaMassDicts = []
    for spectra1 in spectraList1:
        for spectra2 in spectraList2:
            deltaMassDict = computeDeltaMassDict(spectra1, spectra2, precision)
            allDeltaMassDicts.append(deltaMassDict)
    return allDeltaMassDicts

"""Count by rounding - This should underestimate due to hard cutoffs at bins"""
def associateMassDeltasBin(deltaMassDict, massDeltaVal, precision):
    nodeScore = 0
    for key in deltaMassDict:
        if round(key, precision) == round(massDeltaVal, precision):
            nodeScore += deltaMassDict[key]
    return nodeScore

"""Count by differences - This should overestimate due to potential for assigning multiple separate real peaks to the single mass value."""
def associateMassDeltasEps(deltaMassDict, massDeltaVal, epsilon):
    nodeScore = 0
    for key in deltaMassDict:
        if abs(key - massDeltaVal) <= epsilon:
            nodeScore += deltaMassDict[key]
    return nodeScore

#Converts Dict<float, int> to List<float>, List<int> and calls plt.scatter
def scatterHelper(deltaMassDict, c='b'):
    x = []
    y = []
    for deltaM in deltaMassDict:
        count = deltaMassDict[deltaM]
        x.append(deltaM)
        y.append(count)

    plt.scatter(x,y,c=c)

#Plots a single deltaMassDict to a figure
def scatterSpectra(title, deltaMassDict, deltaPrecursorMass = None):
    scatterHelper(deltaMassDict)

    if deltaPrecursorMass is not None:
        plt.axvline(x=deltaPrecursorMass, label="Delta Precursor Mass", c='b')
    plt.title(title)
    plt.xlabel("Delta Mass")
    plt.ylabel("Count")
    plt.show()

#Plots multiple delta mass dicts on the same figure
def scatterSpectras(title, deltaMassDictList, deltaPrecursorMass = None, candidateDeltaMassDict = None):
    for deltaMassDict in deltaMassDictList:
        scatterHelper(deltaMassDict)
    
    if deltaPrecursorMass is not None:
        plt.axvline(x=deltaPrecursorMass, label="Delta Precursor Mass", c='b', linestyle='dashed')

    if candidateDeltaMassDict is not None:
        scatterHelper(candidateDeltaMassDict, c='r')

    plt.axvline(x=0, c='k', linestyle='dashed')
    plt.title(title)
    plt.xlabel("Delta Mass")
    plt.ylabel("Count")
    plt.show()

"""Plots a spectrum from an mgf file"""
def lineSpectrum(title, massIntensityDict):
    x = []
    y = []
    
    for m in massIntensityDict:
        x.append(m)

    x.sort()
    for m in x:
        y.append(massIntensityDict[m])
    
    plt.plot(x,y)
    plt.title(title)
    plt.xlabel("Mass")
    plt.ylabel("Intensity")
    plt.show()

"""Plots real spectra against a theoretical spectra"""
def scatterMatch(title, measuredSpectra, theoreticalSpectrum):
    for i in range(len(measuredSpectra)):
        x = measuredSpectra[i]
        y = [i]*len(measuredSpectra[i])
        plt.scatter(x,y)

    x = theoreticalSpectrum
    y = [len(measuredSpectra)] * len(theoreticalSpectrum)
    plt.scatter(x,y,c="r")


    plt.title(title)
    plt.xlabel("Mass")
    plt.ylabel("Index")
    plt.show()

def lineSpectrumMatch(title, massIntensityDict, theoreticalSpectrum):
    x = []
    y = []
    
    for m in massIntensityDict:
        x.append(m)

    x.sort()
    for m in x:
        y.append(massIntensityDict[m])
    
    plt.plot(x,y)
    plt.scatter(x,y, marker="*", c="b")

    
    x = theoreticalSpectrum
    highestIntensity = max(y)
    y = [highestIntensity] * len(theoreticalSpectrum)
    plt.scatter(x,y,c="r")
    for val in x:
        plt.axvline(x=val, c='r')

    plt.axhline(y=0, c='k', linestyle='-')


    plt.title(title)
    plt.xlabel("Mass")
    plt.ylabel("Intensity")
    plt.show()


#scatterSpectra("Hello World", {1:2,2:4,3:6,4:4,5:10})

#diffs = computeDeltaMassDict([187.0, 187.0, 195.0, 195.0], [187.0, 187.0, 195.0, 200.0])
#scatterSpectra("Hello World2", diffs)

#diffList = computeDeltaMassDicts(
#    [[187.0, 187.0, 195.0, 195.0], [187.0, 187.0, 196.0, 195.0]],
#    [[187.0, 187.0, 195.0, 200.0], [187.0, 187.0, 195.0, 214.0]],
#)
#scatterSpectras("Hello World3", diffList, 5)

#lineSpectrum("Hello", {1:2,2:4,3:6,4:4,5:10})
