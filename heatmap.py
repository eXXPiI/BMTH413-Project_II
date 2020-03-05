"""
Preamble
# Program: heatmap.py
# Author: Jonathan Myers and Kathi Munoz-Hofmann
# Date: Wed Jan 14 2020
# Purpose: Determine heatmap for loaded data.
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
    
# def scoreDifference(stg1, stg2 = "Family"):
#     import numpy as np
#     matrix = np.loadtxt(stg1,delimiter=",",ndmin=2)
#     file = open(stg,"r")
#     seq = file.readline().split(",")
#     file.close()
#     return matrix,seq

def scoreWrite(scores):
    import numpy as np
    import os
    os.mkdir('EndScores')
    os.chdir('./EndScores')
    stg = "EndScores_{}.csv"
    newstg = stg.format(0)
    np.savetxt(newstg,scores,fmt="%1.3f",delimiter=",")

def protDictionary(families, seqs):
    # Given: a list of sequences and corresponding family ID's
    # Returns: a dictionary with family ID as the key and a list of sequences in that family as the value
    from collections import defaultdict
    dic = defaultdict(list)
    for i in range(len(families)):
        dic[families[i]].append(seqs[i])
    return dic

def newMatrix(dic):
  families = []
  values = []
  for key in dic.keys():
    families.append(key)
    for i in range(len(dic.get(key))):
      values.append(dic.get(key)[i])
  return families, values

def fullAlign(data):
    import numpy as np
    #tau_quant = 2
    #tau_possible = np.linspace(0.1,0.9,tau_quant) #loop through
    mn_scores = len(data)
    scores = np.zeros((mn_scores,mn_scores))
    #for k in range(len(tau_possible)):
    tau = 0.619
    for i in range(scores.shape[0]):
        for j in range(scores.shape[0]):
            seq1 = data[i]
            seq2 = data[j]
            bestscore = individualAlign(tau,seq1,seq2)
            scores[i,j] = bestscore
    scoreWrite(scores)

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
    dic = protDictionary(family,values)
    family,values = newMatrix(dic)
    fullAlign(values)
    

main()

# M02 End Program 