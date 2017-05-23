import sys, glob, os, re
import numpy as np
import matplotlib.pyplot as plt

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
#print(flist)

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
    frame_count = 0
    with open(fname) as f:
        for line in f:
            frame_count = frame_count + 1
            no_return = line.strip('\n')
            line_list = no_return.split()
            if data_flag == 1 and atom_count < number_of_atoms:
                if int(line_list[1]) == 5:
                    print(line_list)
                if int(line_list[1]) == 3:
                    print(line_list)
                    line = next(f)
                    no_return = line.strip('\n')
                    line_list = no_return.split()
                    print(line_list)
                    #sys.exit()
                atom_count = atom_count + 1
                if atom_count == 220:
                    sys.exit()
            #include this to skip first 9 non-data lines in each frame
            if "id" in line_list:
                data_flag = 1
                atom_count = 0
                #print(number_of_atoms)
                #        print(data_list)
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
