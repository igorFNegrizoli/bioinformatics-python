def DPChange(money, coins):
    minNumCoins = [0]*(money+1)
    for i in range(1, money+1):
        minNumCoins[i] = max(minNumCoins) + 1
        for coin in coins[::-1]:
            if i >= coin:
                if minNumCoins[i-coin] < minNumCoins[i]:
                    minNumCoins[i] = minNumCoins[i-coin] + 1
    return minNumCoins[-1]

def manhattanTourist(n,m,matdown,matright):
    mat = []
    for i in range(n+1):
        mat.append([])
    for i in mat:
        for j in range(m+1):
            i.append(0)
    #path algorithm
    #setting the values of left and upper borders
    for i in range(1,m+1):
        mat[0][i] = mat[0][i-1] + matright[0][i-1]
    for i in range(1,n+1):
        mat[i][0] = mat[i-1][0] + matdown[i-1][0]
    #calculating the rest
    for i in range(1,n+1):
        for j in range(1,m+1):
            val1 = mat[i][j-1] + matright[i][j-1]
            val2 = mat[i-1][j] + matdown[i-1][j]
            if val1 > val2:
                mat[i][j] = val1
            else:
                mat[i][j] = val2 
    return(mat[n][m])

def manhattanExercise(n,m):
    #reads down matrix
    file = open("files/datainput.txt", "r")
    patts = file.read().strip().split('\n')
    file.close()
    matdown = []
    for i in patts:
        matdown.append(list(map(int, i.strip().split(' '))))
    #reads right matrix
    file = open("files/datainput2.txt", "r")
    patts = file.read().strip().split('\n')
    file.close()
    matright = []
    for i in patts:
        matright.append(list(map(int, i.strip().split(' '))))
    return(manhattanTourist(n,m,matdown,matright))

def longest_common_subsequence(v, w):
#Okay, I admit it, this code isnt mine, I got lazy soryyy. Copied it from this link
#https://github.com/Messi12/bioinformatics-algorithms-coursera/blob/master/Assignment_06C.py
#I was so lazy that I couldnt even download numpy (the original code uses it), so I just did those stupid loops
    # Initialize the array S and iterate through all character of v and w.
    #S = zeros((len(v)+1,len(w)+1), dtype=int)
    S = []
    for i in range(len(v)+1):
        S.append([])
    for i in S:
        for j in range(len(w)+1):
            i.append(0)
    for i in range(len(v)):
        for j in range(len(w)):
            if v[i] == w[j]:
                S[i+1][j+1] = S[i][j]+1
            else:
                S[i+1][j+1] = max(S[i+1][j],S[i][j+1])
    # Recover a maximum substring.
    longest_sseq = ''
    i,j = len(v), len(w)
    while i*j != 0:
        if S[i][j] == S[i-1][j]:
            i -= 1
        elif S[i][j] == S[i][j-1]:
            j -= 1
        else:
            longest_sseq = v[i-1] + longest_sseq
            i -= 1
            j -= 1

    return longest_sseq

def longestPathInDAG(graph, ori, dest):
    lngstDists = {}
    for i in graph.keys():
        lngstDists[i] = [0,0]
    listinha = []
    listinha.append(ori)
    while(listinha != []):
        node = listinha.pop(0)
        for i in graph[node]:
            if lngstDists[i[0]][0] < lngstDists[node][0] + i[1]:
                lngstDists[i[0]][0] = lngstDists[node][0] + i[1]
                lngstDists[i[0]][1] = node
            if i[0] not in listinha:
                if graph[i[0]] != []:
                    listinha.append(i[0])
    return lngstDists


def lngstPathDagExercise(ori, dest):
    ###in
    graph = {}
    f = open('files/datainput.txt','r')
    lines = []
    for i in f.read().strip().split('\n'):
        lines.append(i)
    for i in lines:
        aux = i.split('->')
        if int(aux[0]) not in graph.keys():
            graph[int(aux[0])] = []
        a = list(map(int, aux[1].split(':')))
        if a[0] not in graph.keys():
            graph[a[0]] = []
        graph[int(aux[0])].append(a)
    f.close()
    ###out
    lngstDists = longestPathInDAG(graph, ori, dest)
    listinha = []
    listinha.append(dest)
    while(ori not in listinha):
        listinha.append(lngstDists[listinha[-1]][1])
    print(lngstDists[dest][0])
    print("->".join(list(map(str, listinha[::-1]))))

def calc_combinations(mass, coins):
    if mass==0:
        return 1
    useful_coins = [coin for coin in coins if coin<=mass]
    if len(useful_coins)==0:
        return 0
    if len(useful_coins)==1:
        result = 1 if useful_coins[0]==mass else 0
        return result
    result = sum([calc_combinations(mass-i, coins) for i in coins])
    return result

'''
def global_alignment(v, w, scoring_matrix, sigma):

    # Initialize the matrices.
    S = [[0 for repeat_j in xrange(len(w)+1)] for repeat_i in xrange(len(v)+1)]
    backtrack = [[0 for repeat_j in xrange(len(w)+1)] for repeat_i in xrange(len(v)+1)]

    # Initialize the edges with the given penalties.
    for i in xrange(1, len(v)+1):
        S[i][0] = -i*sigma
    for j in xrange(1, len(w)+1):
        S[0][j] = -j*sigma

    # Fill in the Score and Backtrack matrices.
    for i in xrange(1, len(v)+1):
        for j in xrange(1, len(w)+1):
            scores = [S[i-1][j] - sigma, S[i][j-1] - sigma, S[i-1][j-1] + scoring_matrix[v[i-1], w[j-1]]]
            S[i][j] = max(scores)
            backtrack[i][j] = scores.index(S[i][j])

    # Quick lambda function to insert indels.
    insert_indel = lambda word, i: word[:i] + '-' + word[i:]

    # Initialize the aligned strings as the input strings.
    v_aligned, w_aligned = v, w

    # Get the position of the highest scoring cell in the matrix and the high score.
    i, j = len(v), len(w)
    max_score = str(S[i][j])

    # Backtrack to the edge of the matrix starting at the highest scoring cell.
    while i*j != 0:
        if backtrack[i][j] == 0:
            i -= 1
            w_aligned = insert_indel(w_aligned, j)
        elif backtrack[i][j] == 1:
            j -= 1
            v_aligned = insert_indel(v_aligned, i)
        else:
            i -= 1
            j -= 1

    # Prepend the necessary preceeding indels to get to (0,0).
    for repeat in xrange(i):
        w_aligned = insert_indel(w_aligned, 0)
    for repeat in xrange(j):
        v_aligned = insert_indel(v_aligned, 0)

    return max_score, v_aligned, w_aligned
'''