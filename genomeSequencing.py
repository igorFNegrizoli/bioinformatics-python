from oriFinding import PatternToNumber, NumberToPattern, ReverseComplement
import itertools
#subfunction for ordered_composition
def orderedKmerInsert(list, kmer):
    numberList = []
    for i in list:
        numberList.append(PatternToNumber(i))
    kmerVal = PatternToNumber(kmer)
    for i in range(len(numberList)):
        if numberList[i]>=kmerVal:
            list.insert(i, kmer)
            return list
    list.append(kmer)
    return list

def ordered_composition(text, k):
    output = []
    lenText = len(text)
    for i in range(lenText-k+1):
        pattern = text[i:i+k]
        output = orderedKmerInsert(output, pattern)
    print("Ordered composition finished!!!")
    return output

def composition(text, k):
    output = []
    lenText = len(text)
    for i in range(lenText-k+1):
        output.append(text[i:i+k])
    return output

def pathToGenome(lista):
    genome = ""
    genome += lista[0]
    lenLista = len(lista)
    for i in range(1, lenLista):
        genome += lista[i][-1]
    return genome

def overlappers(lista, iBase):
    out = []
    lenLista = len(lista)
    for i in range(lenLista):
        if i == iBase:
            continue
        if lista[iBase][1:] == lista[i][:-1]:
            out.append(lista[i])
    return out

def overlapGraph(lista):
    out = {}
    lenLista = len(lista)
    for i in range(lenLista):
        ovlp = overlappers(lista, i)
        if len(ovlp) != 0:
            out[lista[i]] = ovlp
    return out

def listaAsString(lista):
    out = ""
    for i in lista:
        out += i
        out += ","
    #out = out[:-1]
    return out

def printOverlapGraph(lista):
    dicti = overlapGraph(lista)
    f = open("files/dataoutput.txt",'w')
    for i in dicti.keys():
        dicti[i] = list(dict.fromkeys(dicti[i]))
        txt = i
        txt += " -> "
        txt += listaAsString(dicti[i])
        txt = txt[:-1]
        txt += "\n"
        f.write(txt)
    f.close()
    print("Done!!")

def deBruijnGraph(text, k):
    graph = {}
    lenText = len(text)
    for i in range(lenText-k+1):
        if text[i:i+k-1] not in graph.keys():
            graph[text[i:i+k-1]] = []
        graph[text[i:i+k-1]].append(text[i+1:i+k])
    print("De Brujin graph is done!")
    return graph

#print graph in file "V -> V1, V2, ..., Vn"
def printGraph(dicti):
    f = open("files/dataoutput.txt",'w')
    lenKmer=0
    for x in dicti.keys():
        lenKmer = len(dicti[x][0])
        break
    maxValue = 4**lenKmer
    while(len(dicti) > 0):
        minValue = maxValue
        minKey = ""
        for i in dicti.keys():
            keyValue = PatternToNumber(i)
            if keyValue < minValue:
                minValue = keyValue
                minKey = i
        txt = minKey
        txt += " -> "
        txt += listaAsString(dicti[minKey])
        txt = txt[:-1]
        txt += "\n"
        f.write(txt)
        dicti.pop(minKey)
    f.close()
    #This chunk just for removing the last empty line ffs
    f = open("files/dataoutput.txt",'r')
    d=f.read()
    f.close()
    m=d.split("\n")
    s="\n".join(m[:-1])
    fd=open("files/dataoutput.txt","w+")
    for i in range(len(s)):
        fd.write(s[i])
    fd.close()

    print("Done printing the graph!")

def printGraphInTerminal(dicti):
    for i in dicti.keys():
        txt = i
        txt += " -> "
        txt += listaAsString(dicti[i])
        txt = txt[:-1]
        print(txt)

def deBruijnGraphWithKmers(kmers):
    graph = {}
    for i in kmers:
        if i[:-1] in graph.keys():
            graph[i[:-1]] = orderedKmerInsert(graph[i[:-1]], i[1:])
        else:
            graph[i[:-1]] = []
            graph[i[:-1]] = orderedKmerInsert(graph[i[:-1]], i[1:])
    print("De Brujin graph with kmers is done!")
    return graph
    
#does the graph has any edges?
def graphHasEdges(graph):
    for i in graph.keys():
        if graph[i] != []:
            return True
    return False

def constructGraph():
    with open("files/datainput.txt", 'r') as file:
        graph = dict((line.strip().split(' -> ') for line in file))
        for key in graph:
            graph[key] = graph[key].split(',')
    return graph

def findVertWithEdge(graph):
    for i in graph.keys():
        if graph[i] != []:
            return i

def findEulerianCycle(graph):
    curr_path = []
    circuit = []
    nextvert = ""
    while graphHasEdges(graph):
        if curr_path == []:
            vert = findVertWithEdge(graph)
        else:
            vert = nextvert
            nextvert = graph[vert][0]
            graph[vert].remove(nextvert)
            vert = nextvert
        while(graph[vert] != []):
            curr_path.append(vert)
            nextvert = graph[vert][0]
            graph[vert].remove(nextvert)
            vert = nextvert
        circuit.append(vert)
        for i in curr_path[::-1]:
            if graph[i] == []:
                curr_path.pop()
                circuit.append(i)
            else:
                nextvert = i
                break
    #f = open('files/dataoutput.txt','w')
    #f.write("->".join(circuit[::-1]))
    #f.close()
    return circuit[::-1]

def findEulerianPath(graph):
    v = ""
    w = ""
    gValues=graph.values()
    gValues = list(itertools.chain(*gValues))
    for i in graph.keys():
        if i not in gValues:
            gValues.append(i)
    gValues = list(set(gValues))
    for i in gValues:
        if i not in graph.keys():
            graph[i]=[]
        inCount = 0
        for j in graph.values():
            if i in j:
                inCount+=1
        outCount = len(graph[i])
        if outCount-inCount == 1:
            v = i
        if inCount-outCount == 1:
            w = i
        if v!="" and w!="":
            break
    if v=="" and w=="":
        print("v =",v,"w =",w)
        print("Maybe what you need is a cycle, not a path")
        return False
    curr_path = []
    circuit = []
    nextvert = ""
    while graphHasEdges(graph):
        if curr_path == []:
            vert = v
        else:
            vert = nextvert
            nextvert = graph[vert][0]
            graph[vert].remove(nextvert)
            vert = nextvert
        while(graph[vert] != []):
            if not graphHasEdges(graph):
                if vert == w:
                    graph[w].append(v)
            curr_path.append(vert)
            nextvert = graph[vert][0]
            graph[vert].remove(nextvert)
            vert = nextvert
        circuit.append(vert)
        for i in curr_path[::-1]:
            if graph[i] == []:
                curr_path.pop()
                circuit.append(i)
            else:
                nextvert = i
                break
    return(circuit[::-1])
    #f = open('files/dataoutput.txt','w')
    #f.write("->".join(circuit[::-1]))
    #f.close()

def stringReconstruction(patterns):
    graph = deBruijnGraphWithKmers(patterns)
    path = findEulerianPath(graph)
    text = pathToGenome(path)
    return text

def allKmers(k):
    output = []
    for i in range(4**k):
        output.append(NumberToPattern(i, k))
    return output

def dec_to_bin(x):
    return str(bin(x)[2:])

def allBinaryStrings(k):
    output = []
    for i in range(2**k):
        num = dec_to_bin(i)
        while(len(num)<k):
            num = "0"+num
        output.append(num)
    return output

def binUniversalCircularString(k):
    kmers = allBinaryStrings(k)
    graph = deBruijnGraphWithKmers(kmers)
    lista = findEulerianCycle(graph)
    #print(lista)
    text = lista[0]
    for i in range(1,len(lista)):
        text = text + lista[i][-1]
    text = text[:-k+1]
    return text

def reconstructionFromPairs(pairs, k, d):
    graph = {}
    for i in pairs:
        pair = i.split("|")
        node1=pair[:]
        node2=pair[:]
        node1[0]=pair[0][:-1]
        node1[1]=pair[1][:-1]
        node2[0]=pair[0][1:]
        node2[1]=pair[1][1:]
        node1="|".join(node1)
        node2="|".join(node2)
        if node1 not in graph.keys():
            graph[node1] = []
            graph[node1].append(node2)
        else:
            graph[node1].append(node2)
    path = findEulerianPath(graph)
    prefix=""
    sufix=""
    for i in range(len(path)-1):
        aux = path[i].split("|")
        prefix = prefix + aux[0][0]
        sufix = sufix + aux[1][0]
    aux = path[-1].split("|")
    prefix = prefix+aux[0]
    sufix = sufix + aux[1]
    output=""
    for i in range(d+k, len(prefix)):
        if prefix[i:] == sufix[:-i]:
            output=prefix[:i]+sufix
            return output
    print("Prefix and suffix do not overlap")

def is1in1out(graph, node):
    if node not in graph.keys():
        return False
    if len(graph[node]) != 1:
        return False
    inCount = 0
    for i in graph.keys():
        for j in graph[i]:
            if j == node:
                inCount +=1
                if inCount > 1:
                    return False
    if inCount==1:
        return True

def maximalNonBranchingPaths(graph):
    paths = []
    cyclicNodes = []
    for node in graph.keys():
        if not is1in1out(graph, node):
            if len(graph[node])>0:
                for i in graph[node]:
                    nonBranchingPath = []
                    nonBranchingPath.append(node)
                    w = i
                    nonBranchingPath.append(w)
                    while(is1in1out(graph, w)):
                        w = graph[w][0]
                        nonBranchingPath.append(w)
                    paths.append(nonBranchingPath)
        else:
            cyclicNodes.append(node)
    for i in paths:
        for j in i:
            if j in cyclicNodes:
                cyclicNodes.remove(j)
    while(cyclicNodes != []):
        v = cyclicNodes[0]
        w = graph[v][0]
        cycle = []
        cycle.append(v)
        cyclicNodes.remove(v)
        while(w!=v):
            cycle.append(w)
            cyclicNodes.remove(w)
            w = graph[w][0]
        cycle.append(v)
        paths.append(cycle)
    return paths

def generateContigs(patterns):
    graph = deBruijnGraphWithKmers(patterns)
    paths = maximalNonBranchingPaths(graph)
    contigs = []
    for i in paths:
        contig = i[0]
        for j in range(1,len(i)):
            contig += i[j][-1]
        contigs.append(contig)
    return contigs

def RNAtoProtein(rna):
    dict = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"", "UAG":"",
    "UGU":"C", "UGC":"C", "UGA":"", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"}
    
    output = ""
    lenRna = len(rna)
    for i in range(0, lenRna, 3):
        output += dict[rna[i:i+3]]
    return output

def DNAtoRNA(dna):
    rna = ""
    for i in dna:
        if i == "T":
            rna += "U"
        else:
            rna += i
    return rna

def findAminoInDNA(dna, amino):
    aminoacids = {"M" : ["ATG"],
    "I": ["ATA", "ATC", "ATT"],
    "A":["GCT", "GCA", "GCC", "GCG"],
    "S":["TCA", "TCC", "TCG", "TCT"],
    "F": ["TTC","TTT"],
    "P":["CCA", "CCC", "CCG", "CCT"],
    "C": ["TGC","TGT"],
    "K": ["AAG","AAA"],
    "H": ["CAT","CAC"],
    "D": ["GAT","GAC"],
    "V": ["GTA", "GTC", "GTG", "GTT"],
    "L": ["TTG","TTA","CTA", "CTC", "CTG", "CTT"],
    "W": ["TGG"],
    "T": ["ACA", "ACC", "ACG", "ACT"],
    "R": ["AGA", "AGG", "CGA", "CGG","CGT", "CGC"],
    "Y": ["TAT", "TAC"]}

    lenDna = len(dna)
    lenAmino = len(amino)*3
    outp = []
    for i in range(lenDna-lenAmino+1):
        slc = dna[i:i+lenAmino]
        sliceRvComp = ReverseComplement(slc)
        slc = DNAtoRNA(slc)
        sliceRvComp = DNAtoRNA(sliceRvComp)
        slc = RNAtoProtein(slc)
        sliceRvComp = RNAtoProtein(sliceRvComp)
        if slc == amino or sliceRvComp == amino:
            outp.append(dna[i:i+lenAmino])
    return outp

def subpeptides(peptide):
    l = len(peptide)
    ls = []
    looped = peptide + peptide
    for start in range(0, l):
        for length in range(1, l):
            ls.append((looped[start:start + length]))
    ls.append(peptide)
    return ls

def aminomass(amino):
    dict = {"G":57, "A":71, "S":87, "P":97, "V":99, "T":101, 
    "C":103, "I":113, "L":113, "N":114, "D":115, "K":128, "Q":128, 
    "E":129, "M":131, "H":137, "F":147, "R":156, "Y":163, "W":186}
    mass = 0
    for i in amino:
        if isinstance(i, int):
            mass += i
        else:
            mass += dict[i]
    return mass

def Cyclospectrum(peptide):
    #if 0 in peptide:
        #peptide.remove(0)
    subpeps = subpeptides(peptide)
    outp= []
    outp.append(0)
    while(subpeps != []):
        lowest = subpeps[0]
        lowestval = aminomass(lowest)
        for i in subpeps:
            if aminomass(i) < lowestval:
                lowest = i
                lowestval = aminomass(i)
        outp.append(lowestval)
        subpeps.remove(lowest)
    return outp

def linearSpectrum(peptide):
    spec = []
    spec.append(0)
    for i in range(1,len(peptide)):
        for j in range(len(peptide)-i+1):
            spec.append(aminomass(peptide[j:j+i]))
    spec.append(aminomass(peptide))
    return sorted(spec)


#Alphabet = {57: 'G', 71: 'A', 87: 'S', 97: 'P', 99: 'V', 101: 'T', 103: 'C', 113:'I/L', 114: 'N', 115: 'D', 128: 'K/Q', 129: 'E',131: 'M', 137: 'H', 147: 'F', 156: 'R', 163: 'Y', 186: 'W'}
#alphabet = {"G":57, "A":71, "S":87, "P":97, "V":99, "T":101, "C":103, "I":113, "L":113, "N":114, "D":115, "K":128, "Q":128, "E":129, "M":131, "H":137, "F":147, "R":156, "Y":163, "W":186}

def CountingMassRec(Mass, masslist):
    if Mass == 0: return 1, masslist
    if Mass < 57: return 0, masslist
    if Mass in masslist: return masslist[Mass], masslist
    n = 0
    for i in Alphabet:
        k, masslist = CountingMassRec(Mass - i, masslist)
        n += k
    masslist[Mass] = n
    return n, masslist

def CountingMass(mass):
    return CountingMassRec(mass, {})[0]

def sumList(lista):
    sum = 0
    for i in lista:
        sum+=i
    return sum


def expand(peptides):
    alphabet = [57,71,87,97,99,101,103,113,114,115,128,129,131,137,147,156,163,186]
    newpeps = []
    for i in peptides:
        for j in alphabet:
            aux = i[:]
            aux.append(j)
            newpeps.append(aux[:])
    return newpeps

def expandStr(peptides):
    alphabet = ["G","A","S","P","V","T","C","I","N","D","K","E","M","H","F","R","Y","W"]
    outp = []
    for i in peptides:
        for j in alphabet:
            outp.append(i+j)
    return outp

def isLAinLB(la, lb):
    for i in la:
        if i not in lb:
            return False
    return True


def cyclopeptideSequencing(spectrum):
    peptides = [[0]]
    outp = []
    while(peptides != []):
        peptides = expand(peptides)
        consistents = peptides[:]
        for i in peptides:
            if aminomass(i) == spectrum[-1]:
                if Cyclospectrum(i) == spectrum:
                    outp.append(i)
                consistents.remove(i)
            else:
                if not isLAinLB(linearSpectrum(i), spectrum):
                    consistents.remove(i)
        peptides = consistents[:]
    return outp
                    