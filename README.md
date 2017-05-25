# many_colloids
analyze.py reads in a trajectory frame by frame
1) read in particle coordinates, and calculate particles' orientations frame by frame
2) for each frame, make a contact matrix depending on particles' relative distance and orientation
3) determine particle clusters from the contact matrix
4) explore how cluster size distribution changes
