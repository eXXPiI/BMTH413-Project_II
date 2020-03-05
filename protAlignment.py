
"""
Preamble
# Program: protAlignment.py
# Author: Jonathan Myers
# Date: Wed Jan 1 15:52:00 2020
# Purpose: Determine protein alignment for loaded data.
# Arguments: None.
# Loads: None.
# Calls: None.
# Returns: None.
"""

def dataFileOpen(stg = "protEnc.txt"):
    file = open(stg,"r")
    filedata = file.readlines()
    temp1 = [element.split(':') for element in filedata]
    family = [int(element[0]) for element in temp1]
    values = [element[1].split(",") for element in temp1]
    _ = [element.pop() for element in values]
    for seq in values:
        for i in range(len(seq)):
            seq[i] = float(seq[i])
    file.close()
    return family,values

def scoreWrite(scores,tau_possible):
    import numpy as np
    import os
    os.mkdir('EndScores')
    os.chdir('./EndScores')
    stg = "EndScores_{}.csv"
    np.savetxt("TauPossible.csv",tau_possible,fmt="%1.3f",delimiter=",")
    for i in range(len(tau_possible)):
        newstg = stg.format(i)
        np.savetxt(newstg,scores[:,:,i],fmt="%1.3f",delimiter=",")

def protDictionary(families, seqs):
    # Given: a list of sequences and corresponding family ID's
    # Returns: a dictionary with family ID as the key and a list of sequences in that family as the value
    from collections import defaultdict
    dic = defaultdict(list)
    for i in range(len(families)):
        dic[families[i]].append(seqs[i])
    return dic

def fullAlign(data):
    import numpy as np
    tau_quant = 2
    tau_possible = np.linspace(0.1,0.9,tau_quant) #loop through
    mn_scores = len(data)
    scores = np.zeros((mn_scores,mn_scores,tau_quant))
    for k in range(len(tau_possible)):
        for i in range(scores.shape[0]):
            for j in range(scores.shape[0]):
                seq1 = data[i]
                seq2 = data[j]
                bestscore = individualAlign(tau_possible[k],seq1,seq2)
                scores[i,j,k] = bestscore
    scoreWrite(scores,tau_possible)

def individualAlign(tau,seq1,seq2):
    import numpy as np
    m = len(seq1)
    n = len(seq2)
    basealign = np.zeros((m,n,2))
    for i in range(1,m):
        for j in range(1,n):
            residue1 = seq1[i-1]
            residue2 = seq2[j-1]
            propensity = np.abs(residue1-residue2)/(residue1+residue2)
            temp_pair = basealign[i-1,j-1,0] + propensity
            temp_gap1 = basealign[i,j-1,0] + tau
            temp_gap2 = basealign[i-1,j,0] + tau
            temp_set = [temp_pair,temp_gap1,temp_gap2]
            score = min(temp_set)
            backtrack = temp_set.index(min(temp_set))
            basealign[i,j,0] = score
            basealign[i,j,1] = backtrack
    return basealign[m-1,n-1,0]

def main():
    family,values = dataFileOpen("protEnc.txt")
    fullAlign(values)

#main()

# M02 End Program 