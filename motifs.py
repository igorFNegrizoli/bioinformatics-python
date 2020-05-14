from oriFinding import Neighbors, HammingDistance, NumberToPattern
import random        
import math
#any(substring in string for substring in substring_list)
#It will return True if any of the substrings in substring_list is contained in string.
        
def patIsinDnaWdMism(pat, Dna, d):
    neigh=Neighbors(pat, d)
    for i in Dna:
        if not any(substring in i for substring in neigh):
            return False
    return True

def MotifEnumeration(Dna, k, d):
    allSlices=[]
    for i in Dna:                          
        for j in range(len(i)-k+1):
            if i[j:j+k] not in allSlices:
                allSlices.append(i[j:j+k])

    sliceNeighbors=[]
    for i in allSlices:
        aux = Neighbors(i, d)
        for j in aux:
            if j not in sliceNeighbors:
                sliceNeighbors.append(j)

    allSlices=[]
    for i in sliceNeighbors:
        if patIsinDnaWdMism(i, Dna, d):
            if i not in allSlices:
                allSlices.append(i)
    
    allSlices=list(dict.fromkeys(allSlices))
    return allSlices

#Receive a list os strings (genomes) and returns a dictionary with the 
#symbols as the nucleotides and the contents as frequency lists
def Count(Motifs):
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(0)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count

#Creates the profile matrix
def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = Count(Motifs)
    for symbol in "ACGT":
        for i in range(k):
            profile[symbol][i] = profile[symbol][i]/t
    return profile

#Creates a consensus string based on the frequency of nucleotides
def Consensus(Motifs):
    consensus = ""
    count = Count(Motifs)
    k = len(Motifs[0])
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus

#Froam a list of Motifs, it returns the score of te motif matrix
def Score(Motifs):
    score = 0
    consensus = Consensus(Motifs)
    lenStr = len(Motifs[0])
    lenMotifs = len(Motifs)
    for i in range(lenMotifs):
        for j in range(lenStr):
            if consensus[j] != Motifs[i][j]:
                score += 1
    return score

#Probability of a given string being the clock genome
def Pr(Text, Profile):
    prob = 1
    lenStr = len(Text)
    for i in range(lenStr):
        prob *= Profile[Text[i]][i]
    return prob

#Finds out the most probale k-mer considering the profile probability matrix
def ProfileMostProbableKmer(text, k, profile):
    bestProb = 0
    mostNotable = ""
    lenText = len(text)
    for i in range(lenText - k + 1):
        sliceStr = text[i: k+i]
        if Pr(sliceStr, profile) > bestProb:
            bestProb = Pr(sliceStr, profile)
            mostNotable = sliceStr
    if mostNotable == "": mostNotable = text[:k]
    return mostNotable

#subfunction of d(pattern, dna)
def subd(pattern, data):
    bestScore=float("inf")
    k = len(pattern)
    lenData = len(data)
    for i in range(lenData-k+1):
        scr = HammingDistance(pattern, data[i:i+k])
        if scr < bestScore:
            bestScore = scr
    return bestScore

#subfunction of medianString(dna, k)
def d(pattern, dna):
    bestScore=0
    for i in dna:
        bestScore += subd(pattern, i)
    return bestScore

#nreturns all possible kmers
def allPossKmers(k):
    allKmers=[]
    for i in range(4**k):
        allKmers.append(NumberToPattern(i, k))
    return allKmers

#finds the kmer with the better score amongst the dna
def medianString(dna, k):
    bestDist = float("inf")
    medianString = ""
    allKmers = allPossKmers(k)
    for i in allKmers:
        dist = d(i, dna)
        if dist < bestDist:
            bestDist=dist
            medianString=i
    return medianString


#Finds out whats the most probable motifs (clock genes) in the Dna patterns
def GreedyMotifSearch(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for m in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][m:m+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

#Calculates the uncertanty of a probability distribution
def Scoreentropy(p):
    score = 0
    for i in "ACGT":
        for j in range(len(p['A'])):
            if p[i][j] == 0.0:
                score+=0
            else:
                score += p[i][j]*(math.log(p[i][j],2))
    return score * -1

#Its like Counts but its pseudocounts. This is an annotation
def CountWithPseudocounts(Motifs):
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(1)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count

#Creates the profile matrix but with pseudocounts
def ProfileWithPseudocounts(Motifs):
    t = len(Motifs) + 4
    k = len(Motifs[0])
    profile = CountWithPseudocounts(Motifs)
    for symbol in "ACGT":
        for i in range(k):
            profile[symbol][i] = profile[symbol][i]/t
    return profile

#Finds out whats the most probable motifs (clock genes) in the Dna patterns
def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for m in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][m:m+k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

#Uses consensusWithPseudocounts()
def ScoreWithPseudocounts(Motifs):
    score = 0
    consensus = Consensus(Motifs)
    lenStr = len(Motifs[0])
    lenMotifs = len(Motifs)
    for i in range(lenMotifs):
        for j in range(lenStr):
            if consensus[j] != Motifs[i][j]:
                score += 1
    return score

#Create a list of 4-kmers to be motifs from a list of dna strings and a profile matrix
def Motifs(Profile, Dna):
    mostProbPats = []
    lenDna = len(Dna)
    for i in range(lenDna):
        mostProbPats.append(ProfileMostProbableKmer(Dna[i], len(Profile["A"]), Profile))
    return mostProbPats

#Select random motifs from a list of dnas. Duh
def RandomMotifs(Dna, k, t):
    randMotifs = []
    dnaLen = len(Dna[0])
    for i in range(t):
        pos = random.randint(0, dnaLen - k)
        randMotifs.append(Dna[i][pos: pos+k])
    return randMotifs

#Monte Carlo motif iteration
def RandomizedMotifSearch(Dna, k, t):
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs 

#Whole Monte carlo motif search
def RepeatedRandomizedMotifSearch(Dna, k, t, N):
    BestScore = float('inf')
    BestMotifs = []
    for i in range(N):
        Motifs = RandomizedMotifSearch(Dna, k, t)
        CurrScore = Score(Motifs)
        if CurrScore < BestScore:
            BestScore = CurrScore
            BestMotifs = Motifs
    return BestMotifs

#Noramalize probavilities in the dictionary
def Normalize(Probabilities):
    total = 0
    for i in Probabilities.keys():
        total += Probabilities[i]
    for i in Probabilities.keys():
        Probabilities[i] = Probabilities[i]/total
    return Probabilities

#Select a random key of the dictionary Probabilities
def WeightedDie(Probabilities):
    total = 0
    value = random.uniform(0, 1)
    for i in Probabilities.keys():
        total += Probabilities[i]
        if value <= total:
            return i

#Select a k-mer in Text according with the profile matrix
def ProfileGeneratedString(Text, profile, k):
    n = len(Text)
    probabilities = {}
    for i in range(0,n-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)
    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)

#N vou explicar essa merda em comentario
def GibbsSampler(Dna, k, t, N):
    Motifs = RepeatedRandomizedMotifSearch(Dna, k, t, 1500)
    bestMotifs = Motifs
    for j in range(N):
        i = random.randint(0,t-1)
        reducedMotifs = []
        for a in range(len(Motifs)):
            if a != i:
                reducedMotifs.append(Motifs[i])
        Profile = ProfileWithPseudocounts(reducedMotifs)
        Motifs[i] = ProfileGeneratedString(Dna[i], Profile, k)
        if Score(Motifs) < Score(bestMotifs):
            bestMotifs = Motifs
    return bestMotifs