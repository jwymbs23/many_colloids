import sys, glob, os, re
import numpy as np
import matplotlib.pyplot as plt
#import seaborn as sns
from sklearn.cluster import DBSCAN
import pprint
flist = []
#for file in glob.glob("traj_rho_0.30_k_50.00_eps_0.00.lammpstrj"):#*.lammpstrj"):
#    flist.append(file)
flist = ["traj_n20_r0.20_k20.00_c5.00_200.lammpstrj"]
#"20_frames.lammpstrj"
#file = "traj_rho_0.30_k_50.00_eps_2.00.lammpstrj"
number_flag = 0
number_of_atoms = 0
data_flag = 0
box_flag = 0
atom_count = 0
data_list = []
box_size = 0
rcut_big2 = 5*5
rcut_small2 = 1.1*1.1
#print(flist)
#sns.set_context('poster')
#sns.set_color_codes()
#plot_kwds = {'alpha' : 0.25, 's' : 80, 'linewidths':0}
#def plot_cluster(data, algorithm, args, kwds):
    #labels = algorithm(*args, **kwds).fit_predict(data)
    #palette = sns.color_palette('deep', np.unique(labels).max() + 1) #dunno what this does
    #colors = [palette[x] if x >= 0 else (0.0,0.0,0.0) for x in labels]
    #plt.scatter(data.T[0], data.T[1], c = colors, **plot_kwds)
    #frame = plt.gca()
    #plt.title('Clusters found by {}'.format(str(algorithm.__name__)), fontsize=24)

def get_box_neighbors(b):
    neighbors=np.asarray([0,0,0,0,0,0,0,0,0])
    pm = [-1,0,1]
    row_number = np.floor(b/(float(bin_num)))
    column_number = b-bin_num*row_number
    count = 0
    for i in pm:
        for j in pm:
            neighbors[count] = ((column_number+j)%bin_num) + bin_num*((row_number+i)%bin_num)
            count += 1
    return neighbors
#    print(b, neighbors, neighbors_test)
#    print(np.floor(b/float(bin_num)+0.0000001), row_number, column_number)


#test_box_neighbors:
#bin_num = 4
#for i in range(bin_num*bin_num):
#    get_box_neighbors(i)
#sys.exit()

#read in init info
for fname in flist:
    print(fname)
    number_flag = 0
    number_of_atoms = 0
    data_flag = 0
    box_flag = 0
    atom_count = 0
    data_list = []
    box_size = 0
    #read_in_file
    with open(fname) as f:
        for line in f:
            no_return = line.strip('\n')
            line_list = no_return.split()
            if number_flag == 1:
                number_of_atoms = int(line_list[0])
                number_flag = 0
            if "NUMBER" in line_list:
                number_flag = 1
            if box_flag == 1:
                box_size = float(line_list[1]) - float(line_list[0])
                box_flag = 0
            if "BOX" in line_list:
                box_flag = 1
            if box_size != 0:
                break
    print(box_size,number_of_atoms)
    bin_num = 15
    bin_size = box_size/float(bin_num)
    print(bin_size)
    frame_count = 0
    big_particle = [ [] for i in range(bin_num*bin_num)]
    #big_particle[0].append([1,2,3])
    #print(big_particle)
    #print(big_particle[0])
    #sys.exit()
    small_particle = [ [] for i in range(bin_num*bin_num)]
    #np.array([]).reshape(0,4)
    small_partition = []#np.array([])
    big_partition = []#np.array([])
    small_ind = 0
    big_ind = 0
    with open(fname) as f:
        for line in f:
            #no_return = line.strip('\n')
            line_list = line.split()
            if data_flag == 1:
                atom_count = atom_count + 1
                if int(line_list[1]) == 5:
                    #print(line_list)
                    xy_big = np.asarray([float(i) for i in line_list[2:-1] ])
                    #print(big_particle)
                    xy_big += -box_size*(xy_big >= box_size) + box_size * (xy_big <= 0.)
                    pbin = int(xy_big[0]/bin_size) + bin_num*int(xy_big[1]/bin_size)
                    big_particle[pbin].append([ big_ind,xy_big[0],xy_big[1] ])
                    big_ind += 1
                    #big_partition.append(int(xy_big[0]/bin_size) + bin_num*int(xy_big[1]/bin_size))
                if int(line_list[1]) == 3:
                    xy_small = np.asarray([float(i) for i in line_list[2:-1] ])
                    xy_small += - box_size*(xy_small >= box_size) + box_size * (xy_small <= 0.)
                    #get next line 
                    line2 = next(f)
                    atom_count += 1
                    line_list2 = line2.split()
                    xy_small_dir = np.asarray([float(i) for i in line_list2[2:-1] ])
                    xy_small_dir += -box_size*(xy_small_dir >= box_size) + box_size * (xy_small_dir <= 0.)
                    #get orientation
                    diff = xy_small_dir - xy_small
                    pbc_diff = [i - float(int(i/(0.5*box_size)))*box_size for i in diff]
                    
                    pbc_diff = pbc_diff/(np.hypot(pbc_diff[0],pbc_diff[1]))#np.sqrt(np.sum(diff*diff)))
                    pbin = int(xy_small[0]/(bin_size)) + bin_num*int(xy_small[1]/(bin_size))
                    if(pbin > bin_num*bin_num):
                        print(pbin, xy_small, 'huh')
                    small_particle[pbin].append([ small_ind, xy_small[0], xy_small[1], pbc_diff[0], pbc_diff[1] ])
                    small_ind += 1
            if atom_count == number_of_atoms:
                data_flag = 0
                atom_count = 0
#                plt.scatter(small_particle.T[0],small_particle.T[1],s=10.,c=small_partition, cmap = 'flag')#, cmap=matplotlib.colors.ListedColormap(colors)
#                plt.scatter(big_particle.T[0],big_particle.T[1],s=100,c=big_partition, cmap = 'flag')
#                if frame_count%10 == 0:
#                    plt.show()
                #sys.exit()
                #
                # 
                #analyze data here:
                #distance matrix:
                cmat = np.zeros((small_ind + big_ind)*(small_ind + big_ind)).reshape(small_ind + big_ind, small_ind + big_ind)
                for box in range(bin_num*bin_num):
                    #print(big_particle[box])
                    neighbors = get_box_neighbors(box)
                    for is1, spart in enumerate(small_particle[box]):
                        for nbox in neighbors:
                            #loop over big particles
                            spartarray = np.asarray(spart[1:3])
                            for bpart in big_particle[nbox]:
                                #test distance cutoff
                                dist = np.asarray(bpart[1:]) - spartarray
                                pbc_dist = [i - float(int(i/(0.5*box_size)))*box_size for i in dist]
                                if pbc_dist[0]*pbc_dist[0] + pbc_dist[1]*pbc_dist[1] < rcut_big2:
                                    #test direction condition
                                    angle = np.arccos(np.dot(pbc_dist,spart[3:])/(np.linalg.norm(pbc_dist)))
                                    #print('big', angle)
                                    if angle < np.pi*0.5:
                                        cmat[bpart[0]][big_ind + spart[0]] = 1
                                        cmat[big_ind + spart[0]][bpart[0]] = 1
                                        
                            #loop over other small particles
                            for is2, spart2 in enumerate(small_particle[nbox]):
                                if(is2 > is1):
                                    dist = np.asarray(spart2[1:3]) - spartarray
                                    pbc_dist = [i - float(int(i/(0.5*box_size)))*box_size for i in dist]
                                    if pbc_dist[0]*pbc_dist[0] + pbc_dist[1]*pbc_dist[1] < rcut_small2:
                                        #test direction condition
                                        angle = np.arccos(np.clip(np.dot(spart2[3:],spart[3:]),-1,1))
                                        #print('small', angle)
                                        if angle < np.pi*0.5:
                                            cmat[big_ind + spart2[0]][big_ind + spart[0]] = 1
                                            cmat[big_ind + spart[0]][big_ind + spart2[0]] = 1
                #end of contact matrix construction
                cmat_ones = np.nonzero(cmat)
                T_cmat_ones = np.transpose(cmat_ones)
                #print(cmat_ones)
                cluster = []
                print(cmat_ones[0])
                for idx,i in enumerate(T_cmat_ones):
                    in_cluster = []
                    in_cluster.append(i[0])
                    next_indices = np.where(cmat_ones[0] == i[1])
                    
                    print(idx, i, in_cluster, next_indices)
                    #walk_through_cluster(i)
                    #cluster.append(walk_through_cluster)
                #def walk_through_cluster(i):
                    #in_cluster.append(i[1])
                    #walk_through_cluster(i[1])
                #for i in cmat_ones[0]:
                #    if 
                #DBSCAN(min_samples = 1).fit_predict(cmat)
                                    #print(bpart[0], spart[0], angle*180/3.14159265)
                    #for ind_s,small in enumerate(small_particle):
                        
                #partition box into smaller boxes, and assign box number to all big and small particles
                #for each big particle, 
                if frame_count == 100:
                    #plt.scatter(small_particle.T[0], small_particle.T[1])
                    #plt.show()
                    plot_cluster(small_particle, cluster.DBSCAN, (), {'eps':5})
                    plt.show()
                    sys.exit()
                #
                #
                #
                #print(frame_count, small_particle,len(big_particle))
                big_particle = [[] for i in range(bin_num*bin_num)]
                small_particle = [[] for i in range(bin_num*bin_num)]
                small_ind = 0
                big_ind = 0
                print(frame_count)
            #include this to skip first 9 non-data lines in each frame
            if "id" in line_list:
                data_flag = 1
                atom_count = 0
                frame_count = frame_count + 1
                #print(number_of_atoms)
                #        print(data_list)
        sys.exit()
    dmat_t4 = np.delete(np.delete(dmat_t4.astype(float)/float(frame_count),-1,0),-1,1)
    dmat_t2 = np.delete(np.delete(dmat_t2.astype(float)/float(frame_count),-1,0),-1,1)
    print(dmat_t2, dmat_t4)
    global_max = max(np.amax(dmat_t2),np.amax(dmat_t4))
    global_min = min(np.amax(dmat_t2),np.amax(dmat_t4))
    fig = plt.figure()
    ax = fig.add_subplot(1,2,1)
    ax.set_aspect('equal')
    im = ax.imshow(dmat_t2, vmin = 0, vmax = global_max)
    #plt.colorbar()
    ax = fig.add_subplot(1,2,2)
    ax.set_aspect('equal')
    im = ax.imshow(dmat_t4,vmin = 0, vmax = global_max)
    #plt.colorbar()
    #ax = fig.add_subplot(1,3,3)
    #ax.set_aspect('equal')
    #im = ax.imshow(dmat_t2 + dmat_t4,vmin = 0, vmax = global_max)
    #plt.colorbar()
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(im, cax=cbar_ax, cmap = 'Greys')
    #plt.show()
    dens_name = fname.replace("traj","dens")
    plt.savefig(dens_name.replace(".lammpstrj",".png"))
