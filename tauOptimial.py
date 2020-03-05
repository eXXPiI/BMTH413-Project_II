
"""
Preamble
# Program: tauOptimial.py
# Author: Jonathan Myers
# Date: Mon Jan 13 19:34:45 2020
# Purpose: None.
# Arguments: None.
# Loads: None.
# Calls: None.
# Returns: None.
"""

def tauOptimize(stg1 = "./EndScores_200_Options"):
    import numpy as np
    import os
    stg2 = "EndScores_{}.csv"
    families = np.loadtxt("Family.csv",delimiter=",")
    uniqueFamilies,uniqueCount = np.unique(families,return_counts=True)
    os.chdir(stg1)
    num = len(os.listdir())-1
    tauPossible = np.loadtxt("TauPossible.csv",delimiter=",")
    tracking = np.zeros((num))
    for i in range(num):
        stg3 = stg2.format(i)
        matrix = np.loadtxt(stg3,delimiter=",",ndmin=2)
        for j in range(matrix.shape[0]):
            family = families[j]
            uniqueIndex = np.where(uniqueFamilies == family)[0][0]
            familyCount = uniqueCount[uniqueIndex]
            row = matrix[j,:]
            minimumIndices = row.argsort()[1:familyCount]
            for k in minimumIndices:
                if family == families[k]:
                    tracking[i] += 1
    maximalSuccessIndices = np.where(tracking == np.amax(tracking))[0]
    bestTau = list(tauPossible[maximalSuccessIndices])
    return bestTau

def main():
    import numpy as np
    bestTau = tauOptimize("./EndScores_200_Options")
    np.savetxt("BestTau.csv",bestTau,fmt="%1.3f",delimiter=",")

main()

# M02 End Program