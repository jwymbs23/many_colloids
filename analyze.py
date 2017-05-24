import sys, glob, os, re
import numpy as np
import matplotlib.pyplot as plt
#import seaborn as sns
import sklearn.cluster as cluster

flist = []
#for file in glob.glob("traj_rho_0.30_k_50.00_eps_0.00.lammpstrj"):#*.lammpstrj"):
#    flist.append(file)
flist = ["20_frames.lammpstrj"]
#file = "traj_rho_0.30_k_50.00_eps_2.00.lammpstrj"
number_flag = 0
number_of_atoms = 0
data_flag = 0
box_flag = 0
atom_count = 0
data_list = []
box_size = 0
rcut_big2 = 5*5
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
    big_particle = []#np.array([]).reshape(0,2)
    small_particle = []#np.array([]).reshape(0,4)
    small_partition = []#np.array([])
    big_partition = []#np.array([])
    with open(fname) as f:
        for line in f:
            #no_return = line.strip('\n')
            line_list = line.split()
            if data_flag == 1:
                atom_count = atom_count + 1
                if int(line_list[1]) == 5:
                    #print(line_list)
                    xy_big = [float(i) for i in line_list[2:-1] ]
                    #print(big_particle)
                    big_particle.append([xy_big[0],xy_big[1]])# = np.r_[big_particle, [[xy_big[0], xy_big[1] ]] ]
                    big_partition.append(int(xy_big[0]/bin_size) + bin_num*int(xy_big[1]/bin_size))# = np.r_[big_partition, [int(xy_big[0]/bin_size) + bin_num*int(xy_big[1]/bin_size)] ]
                if int(line_list[1]) == 3:
                    xy_small = [float(i) for i in line_list[2:-1] ]
                    #get next line 
                    line2 = next(f)
                    atom_count += 1
                    line_list2 = line2.split()
                    xy_small_dir = [float(i) for i in line_list2[2:-1] ]
                    small_particle.append([xy_small[0], xy_small[1], xy_small_dir[0], xy_small_dir[1]])# = np.r_[ small_particle,[[xy_small[0], xy_small[1], diff[0], diff[1] ]] ]
                    small_partition.append(int(xy_small[0]/bin_size) + bin_num*int(xy_small[1]/bin_size))
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
                bound_small = []
                for ind_b,big in enumerate(big_particle):
                    box = big_partition[ind_b]
                    neighbors = get_box_neighbors(box)
                    for ind_s,small in enumerate(small_particle):
                        if small_partition[ind_s] in neighbors:
                            #test distance cutoff
                            dist = np.asarray(big) - np.asarray(small[:2])
                            pbc_dist = [i - float(int(i/(0.5*box_size)))*box_size for i in dist]
                            if pbc_dist[0]*pbc_dist[0] + pbc_dist[1]*pbc_dist[1] < rcut_big2:
                                #test direction condition
                                angle = np.math.atan2(np.linalg.det([dist,small[2:]]),np.dot(dist,small[2:]))#arccos(np.clip(np.dot(dist,small[2:])/np.norm(dist)/np.norm(small[2:]),-1,1))
                                if angle < np.pi*0.5:
                                    bound_small.append([small,ind_s])
                        
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
                big_particle = []#np.array([]).reshape(0,2)
                small_particle = []#np.array([]).reshape(0,4)
                small_partition = []#np.array([])
                big_partition = []#np.array([])
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
