import numpy as np
import sys,os
import matplotlib.pyplot as plt
#sample contact matrix:
#place points down
coord = []
box_size = 10.0
part = 10
cut = 1
for i in range(part):
    coord.append([np.random.rand()*box_size,np.random.rand()*box_size])
coord = np.asarray(coord)
ax = plt.scatter(coord[:,0], coord[:,1])
plt.show()
#make contact matrix:
cmat = np.zeros(part*part).reshape(part,part)
for idxi, i in enumerate(coord):
    for idxj, j in enumerate(coord):
        if idxj > idxi:
            #calculate distance
            dist = i - j
            #apply pbc
            pbc_dist = [adj - float(int(adj/(0.5*box_size)))*box_size for adj in dist]
            if np.linalg.norm(pbc_dist) < cut:
                cmat[idxi][idxj] = 1
                cmat[idxj][idxi] = 1
print(cmat)
cmat_ones = np.nonzero(cmat)
T_cmat_ones = np.transpose(cmat_ones)
print(cmat_ones)
cluster = []
print(cmat_ones[0])
for idx,i in enumerate(T_cmat_ones):
    in_cluster = []
    in_cluster.append(i[0])
    next_indices = np.where(T_cmat_ones[0] == i[1])
    
