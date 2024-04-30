import numpy as np 

def getRibosomeSites(cytoplasm, N_edges,riboNum=500):
    
    ribosome_centers = np.full((N_edges[0], N_edges[1], N_edges[2]), False)
    ribosomes = np.full((N_edges[0], N_edges[1], N_edges[2]), False)

    # riboNum = 500
    
    cyto_coords = np.argwhere(cytoplasm==True)

    ribo_centers = []

    for i in range(riboNum):
        rand_Idx = np.random.randint(len(cyto_coords), size=1)[0]
        ribo_centers.append(cyto_coords[rand_Idx])
        cyto_coords = np.delete(cyto_coords, rand_Idx, 0)

    for center_point in ribo_centers:
        
        x_int = center_point[0]
        y_int = center_point[1]
        z_int = center_point[2]
        
        ribosomes[x_int,y_int,z_int] = True

        ribosomes[x_int+1,y_int,z_int] = True
        ribosomes[x_int-1,y_int,z_int] = True
        ribosomes[x_int,y_int+1,z_int] = True
        ribosomes[x_int,y_int-1,z_int] = True
        ribosomes[x_int,y_int,z_int+1] = True
        ribosomes[x_int,y_int,z_int-1] = True
    
    return ribosomes


def getDNAsites(DNAfile, N_edges,lattice_spacing,N_2_x,N_2_y,N_2_z):
    
#     DNAfile needs to be bin
    
    DNAsites = np.full((N_edges[0], N_edges[1], N_edges[2]), False)
        
    with open(DNAfile,'rb') as f:

        DNAbin = np.fromfile(f,dtype=np.float64,count=-1)

    DNAcoords = DNAbin.reshape((3,DNAbin.shape[0]//3),order='F').T

    # print(DNAcoords.shape)
    pos = []
    for DNAparticle in DNAcoords:

        x = DNAparticle[0]
        y = DNAparticle[1]
        z = DNAparticle[2]

        x_lattice = (x*1e-9)//(10*lattice_spacing)+N_2_x
        y_lattice = (y*1e-9)//(10*lattice_spacing)+N_2_y
        z_lattice = (z*1e-9)//(10*lattice_spacing)+N_2_z

        DNAsites[int(x_lattice),int(y_lattice),int(z_lattice)] = True
        pos.append([int(x_lattice),int(y_lattice),int(z_lattice)])
    return DNAsites, pos


