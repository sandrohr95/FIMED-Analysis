import numpy as np
import math as mt
from scipy import stats as sst

def process(geneLabels, classLabels, sampleIds, values, flags, hkLabel='Housekeeping', endoLabel='Endogenous', negLabel='Negative', eps=.15, bc=1):
    factors=[]
    endoLists=[]
    hkeepers=[]
    
    #background correction (thresholding)
    if bc==1:
        for vector in values:
            vector[vector < eps] = eps
            negMean=sst.gmean(vector[np.where(classLabels==negLabel)])
            vector[vector < negMean] = negMean

    #normalization (houskeepers only)
    for vector in values:
        mask = np.where(classLabels==hkLabel)
        hkeeping=vector[mask]
        hkLabels=geneLabels[mask]
        mask = np.where(classLabels==endoLabel)
        endo=vector[mask]
        endoLabels=geneLabels[mask]
        endoLists.append(endo)
        hkeepers.append(hkeeping)

    A=np.zeros((len(hkLabels),len(hkLabels)))
    M=np.zeros(len(hkLabels))

    # for each gene a series of pairwise ratios
    for i in range(0, len(hkLabels)):
        for j in range(0, len(hkLabels)):
            if i!=j:
                val=[]
                for hkVect in hkeepers:
                    val.append(mt.log(hkVect[i]/hkVect[j], 2))
                A[i,j]=np.std(val)
        M[i]=np.mean(A[i,:])
    
    oldFact=[]

    for n in range(3,10):
        v=100
        factors=[]
        indexes=M.argsort()[:n]
        for hkVect in hkeepers:
            reduced=hkVect[indexes]
            gm=sst.gmean(reduced) 
            factors.append(gm)
        av=np.mean(factors)
        factors=av/factors
        if len(oldFact)>0:
            v=np.std(np.log2(oldFact/factors))
        if v<eps:
            #print(n)
            #print(factors)
            #print(hkLabels[indexes])
            break
        oldFact=factors

    results = []
    results.append(endoLabels)
    for i in range(0,len(factors)):
        f=factors[i]
        endo=endoLists[i]
        results.append(endo*f)
    results.append(flags)
    results.append(sampleIds)
    #print (flags)
    return results


