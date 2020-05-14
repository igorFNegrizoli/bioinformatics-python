import math
import random


#print(*data, sep=' ')
#Printa a lista data separada por espaços

#f = open("Vibrio_cholerae.txt", "r")
#Genome = f.read()

###

#Ori finding problem

###

#Finds how many time a Pattern repeats itself in a Text
def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count 

#Returns a dictionary with all the k length substrings in text and its frequencies
def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] = PatternCount(Text, Pattern)
    return freq

#Returns a list of substrings with most occurrences in Text
def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        # add each key to words whose corresponding frequency value is equal to m
        if freq[key]==m:
            words.append(key)
    return words

#More efficient frequency counting
def ComputingFrequencies(text, k):
    frequencyArray = []
    for i in range(4**k):
        frequencyArray.append(0)
    lenText = len(text)
    for i in range(lenText-k+1):
        frequencyArray[PatternToNumber(text[i:i+k])] += 1
    return frequencyArray

#Only to use as a subfunction
def findClumps(Text, k, t):
    freq=[]
    Le = len(Text)
    for i in range(Le-k+1):
            Pattern = Text[i:i+k]
            occur = PatternCount(Text, Pattern)
            if occur >= t:
                if not Pattern in freq:
                    freq.append(Pattern)
    return freq

#Clump finding
def ClumpFinding(Genome, k, L, t):
    freq = []
    Le = len(Genome)
    for i in range(Le-L+1):
        Patterns = findClumps(Genome[i:i+L], k, t)
        for j in Patterns:
            if not j in freq:
                freq.append(j)
    return freq

#Used in the optimized clump finding problem
def PatternToNumber(Pattern):
    number=""
    for i in range(len(Pattern)):
        if Pattern[i]=="A":
            number+="0"
        elif Pattern[i]=="C":
            number+="1"
        elif Pattern[i]=="G":
            number+="2"
        else:
            number+="3"
    return int(number, 4)

#int i, b
#returns string
#used in NumberToPattern
def base_convert(i, b):
    result = ""
    while i > 0:
            result = str(i % b) + result
            i = i // b
    return result

#Used in the optimized clump finding
def NumberToPattern(Number, t):
    pat = ""
    num = str(base_convert(Number, 4))
    for i in range(len(num)):
        if num[i] == "0":
            pat += "A"
        if num[i] == "1":
            pat += "C"
        if num[i] == "2":
            pat += "G"
        if num[i] == "3":
            pat += "T"
    while len(pat)<t:
        pat = "A"+pat
    return pat

#Finds frequent words considering mutations
def FrequentWordsWithMismatches(Text, k, d):
    freq = {}
    for i in range(4**k):
        pattern = NumberToPattern(i, k)
        freq[pattern] = ApproximatePatternCount(Text, pattern, d)
    bestKmers = []
    maxOcurr = -1
    for i in freq.keys():
        if freq[i] == maxOcurr:
            bestKmers.append(i)
        elif freq[i] > maxOcurr:
            bestKmers = []
            bestKmers.append(i)
            maxOcurr = freq[i]
    return bestKmers

#Find the most frequent kmer with max d mismatches in a text, it counts its reverse complements too
def FrequentWordsWithMismatchesAndReverseComplements(Text, k, d):
    freq = {}
    for i in range(4**k):
        pattern = NumberToPattern(i, k)
        freq[pattern] = ApproximatePatternCount(Text, pattern, d)
        patRC = ReverseComplement(pattern)
        freq[pattern] += ApproximatePatternCount(Text, patRC, d)
    bestKmers = []
    maxOcurr = -1
    for i in freq.keys():
        if freq[i] == maxOcurr:
            bestKmers.append(i)
        elif freq[i] > maxOcurr:
            bestKmers = []
            bestKmers.append(i)
            maxOcurr = freq[i]
    return bestKmers

#Finds patterns with equal or less than d mismatches
def Neighbors(Pattern, d):
    Neighborhood = []
    k = len(Pattern)
    for i in range(4**k):
        pattern = NumberToPattern(i, k)
        hDist = HammingDistance(pattern, Pattern)
        if hDist <= d and pattern not in Neighborhood:
            Neighborhood.append(pattern)
    return Neighborhood

#Finds the reverse complement of Pattern
def ReverseComplement(Pattern):
    Pattern = Reverse(Pattern) # reverse all letters in a string
    Pattern = Complement(Pattern) # complement each letter in a string
    return Pattern

#Reverse given string
def Reverse(Pattern):
    slicedString=Pattern[::-1]
    return slicedString

#Transform a pattern in its complementary string of nucleotids
def Complement(Pattern):
    lista = list(Pattern)
    for i in range(len(Pattern)):
        if lista[i] == "A":
            lista[i] =  "T"
        elif lista[i] == "T":
            lista[i] =  "A"
        elif lista[i] == "G":
            lista[i] =  "C"
        else:
            lista[i] = "G"
    Pattern = ''.join(lista)
    return Pattern

#Return a list of positions in wich a given pattern starts. Funciona
def PatternMatching(Pattern, Genome):
    positions = []
    n = len(Genome)
    k = len(Pattern)
    for i in range(n-k+1):
        slice = Genome[i : i+k]
        if slice == Pattern:
            positions.append(i)
    return positions

#Count the number of a Symbol in various windows of half the lenght of the genome (considering a circular genome). This algorithm is highly inneficcient
def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(ExtendedGenome[i:i+(n//2)], symbol)
    return array

#Count the number of a Symbol in various windows of half the lenght of the genome (considering a circular genome). This algorithm is eficcient
def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]

    array[0] = PatternCount(Genome[0:n//2], symbol)

    for i in range(1, n):
        array[i] = array[i-1]
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array

#Returns an array with the differences between G and C in a Genome
def SkewArray(Genome):
    skew = []
    skew.append(0)
    for i in range(1, len(Genome)+1):
        skew.append(skew[i-1])
        if Genome[i-1]=="C":
            skew[i] -= 1
        elif Genome[i-1]=="G":
            skew[i] += 1
    return skew

#Finds the position in a genome where the skew diagram attains a minimum.
def MinimumSkew(Genome):
    positions = []
    skewArray = SkewArray(Genome)
    positions.append(0)
    for i in range(1, len(skewArray)):
        if skewArray[i] < skewArray[positions[0]]:
            positions = []
            positions.append(i)
        elif skewArray[i] == skewArray[positions[0]]:
            positions.append(i)
    return positions

#Detect how many differences (Hamming distance) exist comparing two genomes
def HammingDistance(p, q):
    hammingDistance = 0
    if len(p) != len(q):
        print ("HammingDistance(p, q) error. len(", p, ") != len(", q,")")
        return -1
    for i in range(len(p)):
        if p[i] != q[i]:
            hammingDistance += 1
    return hammingDistance

# Find all approximate occurrences of a pattern in a string. Condidering mutations
def ApproximatePatternMatching(Text, Pattern, d):
    positions = []
    lenP = len(Pattern)
    lenT = len(Text)
    for i in range(lenT - lenP + 1):
        if HammingDistance(Pattern, Text[i:lenP+i]) <= d:
            positions.append(i)
    return positions

#The name says it all, no secret here
def ApproximatePatternCount(Text, Pattern, d):
    count = len(ApproximatePatternMatching(Text, Pattern, d))
    return count
