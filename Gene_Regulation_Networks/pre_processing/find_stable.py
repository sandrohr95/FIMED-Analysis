# -*- coding: utf-8 -*-
"""
@author: Sandro Hurtado
"""

import numpy as np
import math as mt
#values is a numpy array (genes as rows, columns as patients)
#n is the number of the most stable genes you want to get

def findStable(values, n):
    totalGenes=values.shape[0]
    
    V=np.zeros((totalGenes,totalGenes))
    M=np.zeros(totalGenes)

    for i in range(0, totalGenes):
        for j in range(0, totalGenes):
            if i!=j:
                val=[]
                for colVect in values.T:
                    val.append(mt.log(colVect[i]/colVect[j], 2))
                V[i,j]=np.std(val)
        M[i]=np.mean(V[i,:])
    
    return M.argsort()[:n]

def findUnstable(values, n):
    totalGenes=values.shape[0]
    
    V=np.zeros((totalGenes,totalGenes))
    M=np.zeros(totalGenes)
    
#    print(values.T)

    for i in range(0, totalGenes):
        for j in range(0, totalGenes):
            if i!=j:
                val=[]
                for colVect in values.T:
                    val.append(mt.log(colVect[i]/colVect[j], 2)) #We look how variate de value of this gene
                V[i,j]=np.std(val)
        M[i]=np.mean(V[i,:])
    
    return (-M).argsort()[:n]
