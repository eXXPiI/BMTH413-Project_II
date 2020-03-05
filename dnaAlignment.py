
"""
Preamble
# Program: dnaAlignment.py
# Author: Jonathan Myers
# Date: Wed Jan 1 15:52:00 2020
# Purpose: Determine DNA alignment for loaded data.
# Arguments: None.
# Loads: None.
# Calls: None.
# Returns: None.
"""

def dataFileOpen(stg = "DPalign1.txt"):
    file = open(stg,"r")
    filedata = file.readlines()
    data = [element.split()[0] for element in filedata]
    file.close()
    return data

def scoreFileOpen(stg = "SubstitutionScore.txt"):
    import numpy as np
    matrix = np.loadtxt(stg,delimiter=",",skiprows=1,ndmin=2)
    file = open(stg,"r")
    seq = file.readline().split(",")
    file.close()
    return matrix,seq

def scoreWrite(stg,scores):
    import numpy as np
    np.savetxt(stg,scores,fmt="%1.3f",delimiter=",")

def fullAlign(data,stg = "EndScores.csv"):
    import numpy as np
    mn_scores = len(data)
    mn_bases = len(data[0])+1
    scores = np.zeros((mn_scores,mn_scores))
    subs_matrix,subs_bases = scoreFileOpen()
    for i in range(scores.shape[0]):
        for j in range(scores.shape[0]):
            seq1 = data[i]
            seq2 = data[j]
            bestscore = individualAlign(mn_bases,subs_matrix,subs_bases,seq1,seq2)
            scores[i,j] = bestscore
    scoreWrite(stg,scores)

def individualAlign(mn_bases,subs_matrix,subs_bases,seq1,seq2):
    import numpy as np
    basealign = np.zeros((mn_bases,mn_bases,2))
    lastvalue = len(subs_bases)-1
    for k in range(1,mn_bases):
        basealign[0,k] = k
        basealign[k,0] = k
    for i in range(1,mn_bases):
        for j in range(1,mn_bases):
            base1 = seq1[i-1].capitalize()
            base2 = seq2[j-1].capitalize()
            index1 = subs_bases.index(base1)
            index2 = subs_bases.index(base2)
            temp_pair = basealign[i-1,j-1,0] + subs_matrix[index1,index2]
            temp_gap1 = basealign[i,j-1,0] + subs_matrix[lastvalue,index2]
            temp_gap2 = basealign[i-1,j,0] + subs_matrix[index1,lastvalue]
            temp_set = [temp_pair,temp_gap1,temp_gap2]
            score = min(temp_set)
            backtrack = temp_set.index(min(temp_set))
            basealign[i,j,0] = score
            basealign[i,j,1] = backtrack
    return basealign[mn_bases-1,mn_bases-1,0]

def main():
    data = dataFileOpen("DPalign1.txt")
    fullAlign(data,"EndScoresTrue.csv")

main()

# M02 End Program 