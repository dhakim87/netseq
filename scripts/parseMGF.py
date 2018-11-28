def readMGF(filename):

    results = []

    headerDict = None
    massIntensityDict = None

    with open(filename, 'r') as f:
        for line in f:
            if line.strip() == "BEGIN IONS":
                headerDict = {}
                massIntensityDict = {}
                continue
            if line.strip() == "END IONS":
                results.append((headerDict, massIntensityDict))
                continue

            ss = line.split("=")
            if len(ss) == 2:
                headerDict[ss[0]] = ss[1]
                continue
            
            ss = line.split()
            if len(ss) == 2:
                massIntensityDict[float(ss[0])] = float(ss[1])
                continue

            raise Exception("OOPS" + line)

    return results

def findMasses(massIntensityDict):
    x = []
    for mass in massIntensityDict:
        x.append(mass)
    x.sort()
    return x

def genTheoreticalMasses(aminoMassList):
    x = set([])
    #Fixed point in 1/10mDa by multiplying and dividing by 10000.  Can pass in a number instead.
    for startIndex in range(len(aminoMassList)):
        sum = 10000
        for numElements in range(len(aminoMassList)):
            sum += int((aminoMassList[(startIndex + numElements)%(len(aminoMassList))]) * 10000)
            x.add(sum / 10000.0)
    x = list(x)
    x.sort()
    return x

#c785 = readMGF("../data/surugamide/spectra/cyclonovo_785.5337524.mgf")
#c799 = readMGF("../data/surugamide/spectra/cyclonovo_799.5462036.mgf")
#c856 = readMGF("../data/surugamide/spectra/cyclonovo_856.5719605.mgf")
#c870 = readMGF("../data/surugamide/spectra/cyclonovo_870.5863648.mgf")
#c884 = readMGF("../data/surugamide/spectra/cyclonovo_884.6019287.mgf")
#c898 = readMGF("../data/surugamide/spectra/cyclonovo_898.6170044.mgf")
#c912 = readMGF("../data/surugamide/spectra/cyclonovo_912.6369019.mgf")
#c914 = readMGF("../data/surugamide/spectra/cyclonovo_914.61.mgf")
#
#import spectralPlotting
#for result in [c785, c799, c856, c870, c884, c898, c912, c914]:
#    spectralPlotting.lineSpectrum(result[0][0]["PEPMASS"], result[0][1])

#cAll = readMGF("../data/surugamide/spectra/all_eightspetra.mgf")
#import spectralPlotting
#for result in cAll:
#    spectralPlotting.lineSpectrum(result[0]["PEPMASS"], result[1])
#
