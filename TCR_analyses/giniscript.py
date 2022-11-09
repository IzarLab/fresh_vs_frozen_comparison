import numpy as np

### title: Assign TCR clonality information to T cells to rds files for a particular dataset
### author: Yiping Wang date: 11/08/2022

#code from https://github.com/oliviaguest/gini/blob/master/gini.py
def gini(array):
    """Calculate the Gini coefficient of a numpy array."""
    # based on bottom eq:
    # http://www.statsdirect.com/help/generatedimages/equations/equation154.svg
    # from:
    # http://www.statsdirect.com/help/default.htm#nonparametric_methods/gini.htm
    # All values are treated equally, arrays must be 1d:
    array = array.flatten()
    if np.amin(array) < 0:
        # Values cannot be negative:
        array -= np.amin(array)
    # Values cannot be 0:
    array += 0.0000001
    # Values must be sorted:
    array = np.sort(array)
    # Index per array element:
    index = np.arange(1,array.shape[0]+1)
    # Number of array elements:
    n = array.shape[0]
    # Gini coefficient:
    return ((np.sum((2 * index - n  - 1) * array)) / (n * np.sum(array)))

if __name__=="__main__":
    inFI = open("/mnt/vdb/home/ubuntu2/rlelist_noncum.txt")
    header = inFI.readline()
    occurrences = []
    for line in inFI:
        words = line.strip().split('\t')
        occurrences.append(int(words[1]))

    gini_simul=[]
    for i in range(len(occurrences)):
        gini_simul+=[float(i) for j in range(occurrences[i])]
    gini_simul=np.asarray(gini_simul)
    gini_val = gini(gini_simul)
    outFI = open('/mnt/vdb/home/ubuntu2/gini_val.txt','w')
    outFI.write('gini\n')
    outFI.write(str(gini_val)+'\n')
    outFI.close()
