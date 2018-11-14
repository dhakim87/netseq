def LevenshteinHelper(A, B, i, j, lookup):
    if (i == 0):
        return j
    if (j == 0):
        return i
    
    if (i, j) in lookup:
        return lookup[(i,j)]

    indicator = (A[i-1] == B[j-1])
    if indicator:
        indicator = 0
    else:
        indicator = 1

    result = min(
        LevenshteinHelper(A, B, i-1, j, lookup) + 1,
        LevenshteinHelper(A, B, i, j-1, lookup) + 1,
        LevenshteinHelper(A, B, i-1, j-1, lookup) + indicator)

    lookup[(i,j)] = result
    return result

def compute(A, B):
    return LevenshteinHelper(A, B, len(A), len(B), {})

def compute2(A, B, maxDist):
    w = len(A)
    h = len(B)

    #Initial row of the array
    row = []
    for x in range(w):
        row.append(x+1)
    for y in range(h):
        newRow = [0] * w
        minVal = maxDist + 1
        for x in range(w):
            indicator = (A[x] == B[y])
            if indicator:
                indicator = 0
            else:
                indicator = 1
        
            if x == 0:
                val = min((y+1) + 1, y + indicator, row[x] + 1)
            else:
                val = min(newRow[x-1] + 1, row[x-1] + indicator, row[x] + 1)
            newRow[x] = val
            minVal = min(minVal, val)

        if minVal > maxDist:
            return None
#        print row
#        raw_input("K.")
        row = newRow
    return row[-1]

def computeWithBT(A, B, maxDist):
    w = len(A)
    h = len(B)

    backTraces = []
    row = []
    backTrace = []
    for x in range(w):
        row.append(x+1)
    for y in range(h):
        newRow = [0] * w
        minVal = maxDist + 1
        for x in range(w):
            indicator = (A[x] == B[y])
            if indicator:
                indicator = 0
            else:
                indicator = 1
        
            btSet = set([])
            if x == 0:
                leftVal = (y+1) + 1
                diagVal = (y+indicator)
                downVal = row[x] + 1
                val = min(leftVal, diagVal, downVal)
                if leftVal == val:
                    btSet.add('L') #Left
                if diagVal == val:
                    if indicator == 0:
                        btSet.add('K') #Keep
                    else:
                        btSet.add('S') #Swap
                if downVal == val:
                    btSet.add('D')

            else:
                leftVal = newRow[x-1] + 1
                diagVal = row[x-1] + indicator
                downVal = row[x] + 1
                val = min(leftVal, diagVal, downVal)
                if leftVal == val:
                    btSet.add('L') #Left
                if diagVal == val:
                    if indicator == 0:
                        btSet.add('K') #Keep
                    else:
                        btSet.add('S') #Swap
                if downVal == val:
                    btSet.add('D')

            backTrace.append(btSet)
            newRow[x] = val
            minVal = min(minVal, val)

        if minVal > maxDist:
            return None
        row = newRow

        print row
        print backTrace
        backTraces.append(backTrace)
        backTrace = []
        
    def recursiveBacktrack(x, y, backTraces, scratch):
        #enumerates paths.
        if x in scratch and y in scratch[x]:
            return scratch[x][y]
        
        if x < 0 or y < 0:
            if x >= 0:
                return ["L" * (x - -1)]
            if y >= 0:
                return ["D" * (y - -1)]
            return [""]

        allPaths = []
        for backPtr in backTraces[y][x]:
            if backPtr == 'L':
                paths = recursiveBacktrack(x-1, y, backTraces, scratch)
            if backPtr == 'D':
                paths = recursiveBacktrack(x, y-1, backTraces, scratch)
            if backPtr == 'K' or backPtr == 'S':
                paths = recursiveBacktrack(x-1, y-1, backTraces, scratch)
            for p in paths:
                allPaths.append(p + backPtr)

        if x not in scratch:
            scratch[x] = {}
        scratch[x][y] = allPaths
        return allPaths
        
    followedTraces = recursiveBacktrack(x, y, backTraces, {})

    return (row[-1], followedTraces)


#print("HELLO WORLD\n")
#A = "Eval"
#B = "EvalDEvilish"
#tuple = computeWithBT(A, B, 10)
#print(A + " -- " + B + " Dist: " + str(tuple[0]))
#for trace in tuple[1]:
#    print("\t" + trace)
