from eq import fuzzyEqual
from collections import defaultdict

def score(S, T, eps):
    #Define fuzzyEqual(a,b) as abs(b-a) < eps
    #Return count of all tuples (s, t) such that:
    #s is an element of S
    #t is an element of T and
    #fuzzyEqual(s, t)
    results = []
    for s in S:
        for t in T:
            if fuzzyEqual(s, t, eps):
                results.append((s,t))
    return results

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
