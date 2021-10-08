import numpy as np
import sklearn.metrics as skm
from scipy.stats import pearsonr

def parse(input, posLabel='Positive', negLabel='Negative'):
    results=[]
    dt=[]
    flags=[]

    read=0
    for line in input.splitlines():
        if line.strip().startswith('</Code_Summary>'):
            read=0;
        if read==1 and not line.strip().startswith('CodeClass'):
            className,gene,tax,value=line.strip().split(",")
            dt.append(np.array([gene, className, np.float32(value)]))
            #print(line) 
        if line.strip().startswith('<Code_Summary>'):
            read=1
        if line.strip().startswith('FovCount,'):
            fovCount=line.strip().split(",")[1]
        if line.strip().startswith('FovCounted,'):
            fovCounted=line.strip().split(",")[1]
        if line.strip().startswith('BindingDensity,'):
            bd=float(line.strip().split(",")[1])

    fov=float(fovCounted)/float(fovCount)

    dt = np.array(dt)
   
    flags=[]
    # Quality control flags
    #FOV
    if fov<0.75:
        flags.append('FOV')
    #binding density
    if bd<0.1 or bd>1.8:
        flags.append('BindingDensity_SPRINT')
    if bd<0.1 or bd>2.25:
        flags.append('BindingDensity_MAX/FLEX')
    #linearity of positive controls
    positives=dt[np.any(dt==posLabel, axis=1)]
    positives=positives[~np.any(positives=='POS_F(0.125)', axis=1)]
    known=np.log2(np.array([item[item.find("(")+1:item.find(")")] for item in positives[:,0]]).astype(np.float32))
    measured = np.log2(positives[:,2].astype(np.float32))
    r2=pearsonr(known, measured)[0]
    r2=r2*r2
    if r2<0.95:
        flags.append('PC_Linearity')
    # limit of detection
    posE=dt[np.any(dt=='POS_E(0.5)', axis=1)]
    negatives=dt[np.any(dt==negLabel, axis=1)]
    negatives=negatives[:,2].astype(np.float32)
    mean=np.mean(negatives)
    #print(np.mean(negatives)+(2*np.std(negatives)))
    if posE[0,2].astype(np.float32) < (np.mean(negatives)+(2*np.std(negatives))):
        flags.append('DetectionLimit')
#[positives.find("(")+1:positives.find(")")]
    results.append(dt[:,0]) # at 0 gene name
    results.append(dt[:,1]) # at 1 class name
    results.append(dt[:,2].astype(np.float32)) # at 2 raw value
    results.append(flags) # at 3 flags if any
    return results

        
