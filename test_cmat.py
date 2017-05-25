import numpy
import sys,os

#end of contact matrix construction
cmat_ones = np.nonzero(cmat)
T_cmat_ones = np.transpose(cmat_ones)
#print(cmat_ones)
cluster = []
print(cmat_ones[0])
for idx,i in enumerate(T_cmat_ones):
    in_cluster = []
    in_cluster.append(i[0])
    next_indices = np.where(T_cmat_ones[0] == i[1])
    
