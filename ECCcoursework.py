# function HammingG
# input: a number r
# output: G, the generator matrix of the (2^r-1,2^r-r-1) Hamming code


import numpy as np


def hammingGeneratorMatrix(r):
    n = 2**r-1
    
    #construct permutation pi
    pi = []
    for i in range(r):
        pi.append(2**(r-i-1))
    for j in range(1,r):
        for k in range(2**j+1,2**(j+1)):
            pi.append(k)

    #construct rho = pi^(-1)
    rho = []
    for i in range(n):
        rho.append(pi.index(i+1))

    #construct H'
    H = []
    for i in range(r,n):
        H.append(decimalToVector(pi[i],r))

    #construct G'
    GG = [list(i) for i in zip(*H)]
    for i in range(n-r):
        GG.append(decimalToVector(2**(n-r-i-1),n-r))

    #apply rho to get Gtranpose
    G = []
    for i in range(n):
        G.append(GG[rho[i]])

    #transpose    
    G = [list(i) for i in zip(*G)]

    return G


def decimalToVector(i, l):
    vector = []
    while len(vector) != l: #checks if vector is required length
        vector.append(i%2) #adds 1 or 0 to vector if its odd or even
        i = i//2 #divides i by 2 and rounds down if decimal
    return vector[::-1] #flips vector to give correct binary


def repetitionEncoder(m, n):
    return m*n #multiplies vector by n


def repetitionDecoder(v):
    n0 = v.count(0) #number of 0s
    n1 = v.count(1) #number of 1s
    m = max(n0, n1) #returns the largest number of n0 and n1
    if n0 == n1:
        return [] #returns [0] if same number of 0s as 1s
    elif m == n0:
        return [0] #returns [0] if more 0s than 1s
    elif m == n1:
        return [1] #returns [1] if more 1s than 0s


def message(a):
    l = len(a)
    for i in range(2, 999): #calculates the value of some r >= 2
        if 2**i - 2*i - 1 >= l:
            r = i
            break
    k = 2**r - r - 1 #calculates the value of k
    part1 = decimalToVector(l, r) #first r terms are binary of l
    part2 = a #terms from r+1 to r+l are a
    part3 = repetitionEncoder([0], k-r-l) #terms from r+l+1 to k are [0]
    return part1 + part2 + part3


def hammingEncoder(m):
    l = len(m)
    for i in range(2, 999): #checks if m satisfies 2**r-r-1 for some rc>= 2
        temp = 2**i-i-1
        if temp == l:
            r = i
            break
        elif temp > l: #returns [] only when all feasible values of r have been checked
            return []
    G = hammingGeneratorMatrix(r) #generates hamming matrix with dimension r
    return [i%2 for i in list(np.dot(m, G))] #returns encoded message using matrix multiplication


def hammingDecoder(v): #using syndrome decoding
    l = len(v)
    for i in range(2, 999): #checks if m satisfies 2**r-1 for some r >= 2
        temp = 2**i-1
        if temp == l:
            r = i
            break
        elif temp > l: #returns [] only when all feasible values of r have been checked
            return []
    Htranspose = [decimalToVector(i, r) for i in range(1, 2**r)] #matrix where each column is each binary number
    vHtranspose = [i%2 for i in list(np.dot(v, Htranspose))] #multiplies v and Htranspose
    i = int("".join([str(i) for i in vHtranspose]), 2) #converts binary position to decimal
    if i == 0: #adds 1 to the position i-1 (since numbering starts at 0)
        return v #returns codeword
    else:
        v[i-1] += 1
        v[i-1] = v[i-1] % 2
        return v #returns codeword


def messageFromCodeword(c):
    l = len(c)
    for i in range(2, 999): #checks if m satisfies 2**r-1 for some r >= 2
        temp = 2**i-1
        if temp == l:
            r = i
            break
        elif temp > l: #returns [] only when all feasible values of r have been checked
            return []
    for i in range(r-1,-1,-1): #counts from r-1 to 0
        del c[2**i-1] #remove positions with a power of 2 (-1 since counting starts from 0)
    return c


def dataFromMessage(m):
    k = len(m)
    for i in range(2, 999): #checks if m satisfies 2**r-r-1 for some r >= 2
        temp = 2**i-i-1
        if temp == k:
            r = i
            break
        elif temp > k: #returns [] only when all feasible values of r have been checked
            return []
    lbin = [str(i) for i in m[:r]] #first r positions contain l in binary
    l = int("".join(lbin), 2) #converts l to decimal
    if k - r < l: #checks if value of l is correct (checks part 1)
        return []
    elif m[r+l:k] != [0]*(k-r-l): #checks if correct number of 0s at end (checks part 3)
        return []
    else:
        return m[r: r + l] #returns origional data (part 2)