from  pythtb import *
import numpy as np # numerics for matrices

# This is an extension to PythTB, which allows to define TB-models
# using the Slater-Koster parameters.
# See also J. C. Slater and G. F. Koster, 
# Simplified LCAO method for the periodic potential problem, Phys. Rev. 94, 1498 (1954).
# PythTB is availabe at http://www.physics.rutgers.edu/pythtb/

def init_model_SK(dim_k,dim_r,lat,atoms):
    r"""
    This function initializes the TB model with orbitals assigned to
atoms. The Slater-Koster parametrization is used to calculate the 
hopping amplitudes.

    :param dim_k: Dimensionality of reciprocal space, i.e., specifies how
      many directions are considered to be periodic.

    :param dim_r: Dimensionality of real space, i.e., specifies how many
      real space lattice vectors there are and how many coordinates are
      needed to specify the orbital coordinates.
    .. note:: Parameter *dim_r* can be larger than *dim_k*! For example,
      a polymer is a three-dimensional molecule (one needs three
      coordinates to specify orbital positions), but it is periodic
      along only one direction. For a polymer, therefore, we should
      have *dim_k* equal to 1 and *dim_r* equal to 3. See similar example
      here: :ref:`trestle-example`.

    :param lat: Array containing lattice vectors in Cartesian
      coordinates (in arbitrary units). In example the below, the first
      lattice vector has coordinates [1.0,0.5] while the second
      one has coordinates [0.0,2.0].  By default, lattice vectors
      are an identity matrix.

    :param atoms: Array containing reduced coordinates of all atoms in a unit cell, the types of orbitals localized at the atom position, and their onsite energies. The orbitals are: 0=s, 1=p-x, 2=p-y, 3=p-z, 4=d-xy, 5=d-yz, 6=d-zx, 7=d-x^2-y^2, 8=d-3r^2-z^2.

    Example usage::

       # Creates sp3 model of rocksalt SnTe
       tb = init_model_SK(3, 3,
                lat=[[0.5,0.5,0.0],[0.0,0.5,0.5],[0.5,0.0,0.5]],
                atoms = [[[0.,0.,0.],[0,1,2,3],[-6.578,1.659,1.659,1.659]],[[0.5,-0.5,0.5],[0,1,2,3],[-12.067,-0.167,-0.167,-0.167]]])

    """

    # number of atoms
    n_at = len(atoms)
    print('Number of atoms: ',n_at)
    print('Atom')
    # lists of orbitals and energies
    orb=[]
    eons=[]
    for i in range(n_at):
        print(i)
        print('coordinates:', atoms[i][0])
        # number of orbitals at each atom
        n_orb = len(atoms[i][1])
        for j in range(n_orb):
            print('orbital:',atoms[i][1][j],'energy:',atoms[i][2][j])
            orb.append(atoms[i][0])
            eons.append(atoms[i][2][j])

    # create TB model
    model = tb_model(dim_k,dim_r,lat,orb,nspin=2)

    # set onsite energies
    model.set_onsite(eons)

    return model

def set_hop_SK(model,atoms,n_at_i,n_at_j,lv,sss=0.,sps=0.,pss=0.,pps=0.,ppp=0.,sds=0.,pds=0.,pdp=0.,dss=0.,dps=0.,dpp=0.,dds=0.,ddp=0.,ddd=0.):
    # set the hoppings according to the Slater-Koster parametrization
    # model...TB model to add hoppings to
    # atoms......the list of atoms and their Orbitals
    # n_at_i.....number of atom_i in the list
    # n_at_j.....number of atom_j in the list
    # lv.........lattice vector connecting the atoms
    # parameters:
    #   sss......ss-sigma
    #   sps......sp-sigma
    #   pss......ps-sigma
    #   pps......pp-sigma
    #   ppp......pp-pi
    #   sds......sd-sigma
    #   pds......pd-sigma
    #   pdp......pd-pi
    #   dss......ds-sigma
    #   dps......dp-sigma
    #   dpp......dp-pi
    #   dds......dd-sigma
    #   ddp......dd-pi
    #   ddd......dd-delta

    # dimension
    dim = model._dim_r
    # lattice vectors
    lat=model.get_lat()

    # atoms
    at_i = atoms[n_at_i]
    at_j = atoms[n_at_j]

    # number of orbitals stored below atom_i
    noba_i = 0
    for i in range(n_at_i):
        n_orb = len(atoms[n_at_i][1])
        noba_i += n_orb
    # number of orbitals stored below atom_j
    noba_j = 0
    for j in range(n_at_j):
        n_orb = len(atoms[n_at_j][1])
        noba_j += n_orb

    # distance in r.c.
    dist_tmp = np.array(lv)+np.array(at_j[0])-np.array(at_i[0])
    # distance in cartesian coordinates
    dist = np.zeros(dim)
    for i in range(dim):
        dist += dist_tmp[i]*lat[i]
    # direction cosines
    dc = np.zeros(3)
    for i in range(dim):
        dc[i] = dist[i]/np.sqrt(np.dot(dist,dist))
    print()
    print('Generating the hoppings from the SK parameters')
    print('Direction cosines:',dc)

    # all connections between the orbitals
    n_orb_i = len(at_i[1])
    n_orb_j = len(at_j[1])
    for i in range(n_orb_i):
        for j in range(n_orb_j):
            orb_i = at_i[1][i]
            orb_j = at_j[1][j]

            if orb_i==0:
                if orb_j==0:
                    # s-s
                    amp = sss
                elif orb_j==1:
                    # s-x
                    amp = dc[0]*sps
                elif orb_j==2:
                    # s-y
                    amp = dc[1]*sps
                elif orb_j==3:
                    # s-z
                    amp = dc[2]*sps
                elif orb_j==4:
                    # s-xy
                    amp = np.sqrt(3.0) * dc[0] * dc[1] * sds
                elif orb_j==5:
                    # s-yz                    
                    amp = np.sqrt(3.0) * dc[1] * dc[2] * sds
                elif orb_j==6: 
                    # s-zx                
                    amp = np.sqrt(3.0) * dc[2] * dc[0] * sds
                elif orb_j==7: 
                    # s-x2y2
                    amp = 0.5 * np.sqrt(3.0) * (dc[0]**2 - dc[1]**2) * sds
                elif orb_j==8:
                    # s-z2r2
                    amp = (dc[2]**2 - 0.5 * (dc[0]**2 + dc[1]**2)) * sds
            elif orb_i==1:
                if orb_j==0:
                    # x-s
                    amp = -dc[0]*pss
                elif orb_j==1:
                    # x-x
                    amp = (dc[0]**2)*pps+(1.-dc[0]**2)*ppp
                elif orb_j==2:
                    # x-y
                    amp = dc[0]*dc[1]*(pps-ppp)
                elif orb_j==3:
                    # x-z
                    amp = dc[0]*dc[2]*(pps-ppp)
                elif orb_j==4: 
                    # x-xy
                    amp = np.sqrt(3.0) * dc[0]**2 * dc[1] * pds + dc[1] * (1.0 - 2.0 * dc[0]**2)* pdp
                elif orb_j==5: 
                    # x_yz      
                    amp = np.sqrt(3.0) * dc[0] * dc[1] * dc[2] * pds - 2.0 * dc[0] * dc[1] * dc[2] * pdp
                elif orb_j==6: 
                    # x_zx      
                    amp = np.sqrt(3.0) * dc[0]**2 * dc[2] * pds + dc[2] * (1.0 - 2.0 * dc[0]**2)* pdp
                elif orb_j==7: 
                    # x_x2y2    
                    amp = 0.5 * np.sqrt(3.0) * dc[0] * (dc[0]**2 - dc[1]**2) * pds + dc[0] * (1.0 - dc[0]**2 + dc[1]**2)* pdp
                elif orb_j==8: 
                    # x_z2r2    
                    amp = dc[0] * (dc[2]**2 - 0.5 * (dc[0]**2 + dc[1]**2)) * pds - np.sqrt(3.0) * dc[0] * dc[2]**2 * pdp                    
            elif orb_i==2:
                if orb_j==0:
                    # y-s
                    amp = -dc[1]*pss
                elif orb_j==1:
                    # y-x
                    amp = dc[1]*dc[0]*(pps-ppp)
                elif orb_j==2:
                    # y-y
                    amp = (dc[1]**2)*pps+(1.-dc[1]**2)*ppp
                elif orb_j==3:
                    # y-z
                    amp = dc[1]*dc[2]*(pps-ppp)
                elif orb_j==4: 
                    # y_xy      
                    amp = np.sqrt(3.0) * dc[1]**2 * dc[0] * pds + dc[0] * (1.0 - 2.0 * dc[1]**2)* pdp
                elif orb_j==5: 
                    # y_yz      
                    amp = np.sqrt(3.0) * dc[1]**2 * dc[2] * pds + dc[2] * (1.0 - 2.0 * dc[1]**2)* pdp
                elif orb_j==6: 
                    # y_zx      
                    amp = np.sqrt(3.0) * dc[1] * dc[2] * dc[0] * pds - 2.0 * dc[1] * dc[2] * dc[0] * pdp
                elif orb_j==7: 
                    # y_x2y2    
                    amp = 0.5 * np.sqrt(3.0) * dc[1] * (dc[0]**2 - dc[1]**2) * pds - dc[1] * (1.0 + dc[0]**2 - dc[1]**2)* pdp
                elif orb_j==8: 
                    # y_z2r2    
                    amp = dc[1] * (dc[2]**2 - 0.5 * (dc[0]**2 + dc[1]**2)) * pds - np.sqrt(3.0) * dc[1] * dc[2]**2 * pdp
            elif orb_i==3:
                if orb_j==0:
                    # z-s
                    amp = -dc[2]*pss
                elif orb_j==1:
                    # z-x
                    amp = dc[2]*dc[0]*(pps-ppp)
                elif orb_j==2:
                    # z-y
                    amp = dc[2]*dc[1]*(pps-ppp)
                elif orb_j==3:
                    # z-z
                    amp = (dc[2]**2)*pps+(1.-dc[2]**2)*ppp
                elif orb_j==4: 
                    # z_xy      
                    amp = np.sqrt(3.0) * dc[2] * dc[0] * dc[1] * pds - 2.0 * dc[2] * dc[0] * dc[1] * pdp
                elif orb_j==5: 
                    # z_yz      
                    amp = np.sqrt(3.0) * dc[2]**2 * dc[1] * pds + dc[1] * (1.0 - 2.0 * dc[2]**2)* pdp
                elif orb_j==6: 
                    # z_zx      
                    amp = np.sqrt(3.0) * dc[2]**2 * dc[0] * pds + dc[0] * (1.0 - 2.0 * dc[2]**2)* pdp
                elif orb_j==7: 
                    # z_x2y2    
                    amp = 0.5 * np.sqrt(3.0) * dc[2] * (dc[0]**2 - dc[1]**2) * pds - dc[2] * (dc[0]**2 - dc[1]**2)* pdp
                elif orb_j==8: 
                    # z_z2r2    
                    amp = dc[2] * (dc[2]**2 - 0.5 * (dc[0]**2 + dc[1]**2)) * pds + np.sqrt(3.0) * dc[2] * (dc[0]**2 + dc[1]**2) * pdp
            elif orb_i==4:
                if orb_j==0: 
                    # xy_s      
                    amp = np.sqrt(3.0) * (-dc[0]) * (-dc[1]) * dss
                elif orb_j==1: 
                    # xy_x      
                    amp = np.sqrt(3.0) * dc[0]**2 * (-dc[1]) * dps + (-dc[1]) * (1.0 - 2.0 * dc[0]**2) * dpp
                elif orb_j==2: 
                    # xy_y      
                    amp = np.sqrt(3.0) * dc[1]**2 * (-dc[0]) * dps + (-dc[0]) * (1.0 - 2.0 * dc[1]**2) * dpp
                elif orb_j==3: 
                    # xy_z      
                    amp = np.sqrt(3.0) * (-dc[2]) * (-dc[0]) * (-dc[1]) * dps - 2.0 * (-dc[2]) * (-dc[0]) * (-dc[1]) * dpp
                elif orb_j==4: 
                    # xy_xy     
                    amp = 3.0 * dc[0]**2 * dc[1]**2 * dds + (dc[0]**2 + dc[1]**2 - 4.0 * dc[0]**2 * dc[1]**2) * ddp + (dc[2]**2  + dc[0]**2 * dc[1]**2) * ddd
                elif orb_j==5: 
                    # xy_yz     
                    amp = 3.0 * dc[0] * dc[1]**2 * dc[2] * dds + dc[0] * dc[2] * (1.0 - 4.0 * dc[1]**2) * ddp + dc[0] * dc[2] * (dc[1]**2 - 1.0) * ddd
                elif orb_j==6: 
                    # xy_zx     
                    amp = 3.0 * dc[0]**2 * dc[1] * dc[2] * dds + dc[1] * dc[2] * (1.0 - 4.0 * dc[0]**2) * ddp + dc[1] * dc[2] * (dc[0]**2 - 1.0) * ddd
                elif orb_j==7: 
                    # xy_x2y2   
                    amp = 1.50 * dc[0] * dc[1] * (dc[0]**2 - dc[1]**2) * dds + 2.0 * dc[0] * dc[1] * (dc[1]**2 - dc[0]**2) * ddp + 0.5 * dc[0] * dc[1] * (dc[0]**2 - dc[1]**2) * ddd
                elif orb_j==8: 
                    # xy_z2r2   
                    amp = np.sqrt(3.0)* dc[0] * dc[1] * (dc[2]**2 - 0.5 * (dc[0]**2 + dc[1]**2)) * dds - 2.0 * np.sqrt(3.0) * dc[0] * dc[1] * dc[2]**2 * ddp + 0.5 * np.sqrt(3.0) * dc[0] * dc[1] * (1.0 + dc[2]**2) * ddd
            elif orb_i==5:
                if orb_j==0: 
                    # yz_s      
                    amp = np.sqrt(3.0) * (-dc[1]) * (-dc[2]) * dss
                elif orb_j==1: 
                    # yz_x      
                    amp = np.sqrt(3.0) * (-dc[0]) * (-dc[1]) * (-dc[2]) * dps - 2.0 * (-dc[0]) * (-dc[1]) * (-dc[2]) * dpp
                elif orb_j==2: 
                    # yz_y      
                    amp = np.sqrt(3.0) * dc[1]**2 * (-dc[2]) * dps + (-dc[2]) * (1.0 - 2.0 * dc[1]**2) * dpp
                elif orb_j==3: 
                    # yz_z      
                    amp = np.sqrt(3.0) * dc[2]**2 * (-dc[1]) * dps + (-dc[1]) * (1.0 - 2.0 * dc[2]**2) * dpp
                elif orb_j==4: 
                    # yz_xy     
                    amp = 3.0 * (-dc[0]) * dc[1]**2 * (-dc[2]) * dds + (-dc[0]) * (-dc[2]) * (1.0 - 4.0 * dc[1]**2) * ddp + (-dc[0]) * (-dc[2]) * (dc[1]**2 - 1.0) * ddd
                elif orb_j==5: 
                    # yz_yz     
                    amp = 3.0 * dc[1]**2 * dc[2]**2 * dds + (dc[1]**2 + dc[2]**2 - 4.0 * dc[1]**2 * dc[2]**2) * ddp + (dc[0]**2  + dc[1]**2 * dc[2]**2) * ddd
                elif orb_j==6: 
                    # yz_zx     
                    amp = 3.0 * dc[1] * dc[2]**2 * dc[0] * dds + dc[1] * dc[0] * (1.0 - 4.0 * dc[2]**2) * ddp + dc[1] * dc[0] * (dc[2]**2 - 1.0) * ddd
                elif orb_j==7: 
                    # yz_x2y2   
                    amp = 1.50 * dc[1] * dc[2] * (dc[0]**2 - dc[1]**2) * dds - dc[1] * dc[2] * (1.0 + 2.0 * (dc[0]**2 - dc[1]**2)) * ddp + dc[1] * dc[2] * (1.0 + 0.5 * (dc[0]**2 - dc[1]**2)) * ddd
                elif orb_j==8: 
                    # yz_z2r2   
                    amp = np.sqrt(3.0)* dc[1] * dc[2] * (dc[2]**2 - 0.5 * (dc[0]**2 + dc[1]**2)) * dds + np.sqrt(3.0) * dc[1] * dc[2] * (dc[0]**2 + dc[1]**2 - dc[2]**2) * ddp - 0.5 * np.sqrt(3.0) * dc[1] * dc[2] * (dc[0]**2 + dc[1]**2) * ddd
            elif orb_i==6:
                if orb_j==0: 
                    # zx_s      
                    amp = np.sqrt(3.0) * (-dc[2]) * (-dc[0]) * dss
                elif orb_j==1: 
                    # zx_x      
                    amp = np.sqrt(3.0) * dc[0]**2 * (-dc[2]) * dps + (-dc[2]) * (1.0 - 2.0 * dc[0]**2)* dpp
                elif orb_j==2: 
                    # zx_y      
                    amp = np.sqrt(3.0) * (-dc[1]) * (-dc[2]) * (-dc[0]) * dps - 2.0 * (-dc[1]) * (-dc[2]) * (-dc[0]) * dpp
                elif orb_j==3: 
                    # zx_z      
                    amp = np.sqrt(3.0) * dc[2]**2 * (-dc[0]) * dps + (-dc[0]) * (1.0 - 2.0 * dc[2]**2)* dpp
                elif orb_j==4: 
                    # zx_xy     
                    amp = 3.0 * dc[0]**2 * (-dc[1]) * (-dc[2]) * dds + (-dc[1]) * (-dc[2]) * (1.0 - 4.0 * dc[0]**2) * ddp + (-dc[1]) * (-dc[2]) * (dc[0]**2 - 1.0) * ddd
                elif orb_j==5: 
                    # zx_yz     
                    amp = 3.0 * (-dc[1]) * dc[2]**2 * (-dc[0]) * dds + (-dc[1]) * (-dc[0]) * (1.0 - 4.0 * dc[2]**2) * ddp + (-dc[1]) * (-dc[0]) * (dc[2]**2 - 1.0) * ddd
                elif orb_j==6: 
                    # zx_zx     
                    amp = 3.0 * dc[2]**2 * dc[0]**2 * dds + (dc[2]**2 + dc[0]**2 - 4.0 * dc[2]**2 * dc[0]**2) * ddp + (dc[1]**2  + dc[2]**2 * dc[0]**2) * ddd
                elif orb_j==7: 
                    # zx_x2y2   
                    amp = 1.50 * dc[2] * dc[0] * (dc[0]**2 - dc[1]**2) * dds + dc[2] * dc[0] * (1.0 - 2.0 * (dc[0]**2 - dc[1]**2)) * ddp - dc[2] * dc[0] * (1.0 - 0.5 * (dc[0]**2 - dc[1]**2)) * ddd
                elif orb_j==8: 
                    # zx_z2r2   
                    amp = np.sqrt(3.0)* dc[0] * dc[2] * (dc[2]**2 - 0.5 * (dc[0]**2 + dc[1]**2)) * dds + np.sqrt(3.0) * dc[0] * dc[2] * (dc[0]**2 + dc[1]**2 - dc[2]**2) * ddp - 0.5 * np.sqrt(3.0) * dc[0] * dc[2] * (dc[0]**2 + dc[1]**2) * ddd
            elif orb_i==7:
                if orb_j==0: 
                    # x2y2_s    
                    amp = 0.5 * np.sqrt(3.0) * (dc[0]**2 - dc[1]**2) * dss
                elif orb_j==1: 
                    # x2y2_x    
                    amp = 0.5 * np.sqrt(3.0) * (-dc[0]) * (dc[0]**2 - dc[1]**2) * dps + (-dc[0]) * (1.0 - dc[0]**2 + dc[1]**2) * dpp
                elif orb_j==2: 
                    # x2y2_y    
                    amp = 0.5 * np.sqrt(3.0) * (-dc[1]) * (dc[0]**2 - dc[1]**2) * dps - (-dc[1]) * (1.0 + dc[0]**2 - dc[1]**2) * dpp
                elif orb_j==3: 
                    # x2y2_z    
                    amp = 0.5 * np.sqrt(3.0) * (-dc[2]) * (dc[0]**2 - dc[1]**2) * dps - (-dc[2]) * (dc[0]**2 - dc[1]**2) * dpp
                elif orb_j==4: 
                    # x2y2_xy   
                    amp = 1.50 * (-dc[0]) * (-dc[1]) * (dc[0]**2 - dc[1]**2) * dds + 2.0 * (-dc[0]) * (-dc[1]) * (dc[1]**2 - dc[0]**2) * ddp + 0.5 * (-dc[0]) * (-dc[1]) * (dc[0]**2 - dc[1]**2) * ddd
                elif orb_j==5: 
                    # x2y2_yz   
                    amp = 1.50 * (-dc[1]) * (-dc[2]) * (dc[0]**2 - dc[1]**2) * dds - (-dc[1]) * (-dc[2]) * (1.0 + 2.0 * (dc[0]**2 - dc[1]**2)) * ddp + (-dc[1]) * (-dc[2]) * (1.0 + 0.5 * (dc[0]**2 - dc[1]**2)) * ddd
                elif orb_j==6: 
                    # x2y2_zx   
                    amp = 1.50 * (-dc[2]) * (-dc[0]) * (dc[0]**2 - dc[1]**2) * dds + (-dc[2]) * (-dc[0]) * (1.0 - 2.0 * (dc[0]**2 - dc[1]**2)) * ddp - (-dc[2]) * (-dc[0]) * (1.0 - 0.5 * (dc[0]**2 - dc[1]**2)) * ddd
                elif orb_j==7: 
                    # x2y2_x2y2 
                    amp = 0.750 * (dc[0]**2 - dc[1]**2) * (dc[0]**2 - dc[1]**2) * dds + (dc[0]**2 + dc[1]**2 - (dc[0]**2 - dc[1]**2) * (dc[0]**2 - dc[1]**2)) * ddp + (dc[2]**2 + 0.250 * (dc[0]**2 - dc[1]**2) * (dc[0]**2 - dc[1]**2)) * ddd
                elif orb_j==8: 
                    # x2y2_z2r2 
                    amp = 0.5 * np.sqrt(3.0) * (dc[0]**2 - dc[1]**2) * (dc[2]**2 - 0.5 * (dc[0]**2 + dc[1]**2)) * dds + np.sqrt(3.0) * dc[2]**2 * (dc[1]**2 - dc[0]**2) * ddp + 0.250 * np.sqrt(3.0) * (1.0 + dc[2]**2) * (dc[0]**2 - dc[1]**2) * ddd
            elif orb_i==8:
                if orb_j==0: 
                    # z2r2_s    
                    amp = (dc[2]**2 - 0.5 * (dc[0]**2 + dc[1]**2)) * dss
                elif orb_j==1: 
                    # z2r2_x    
                    amp = (-dc[0]) * (dc[2]**2 - 0.5 * (dc[0]**2 + dc[1]**2)) * dps - np.sqrt(3.0) * (-dc[0]) * dc[2]**2 * dpp
                elif orb_j==2: 
                    # z2r2_y    
                    amp = (-dc[1]) * (dc[2]**2 - 0.5 * (dc[0]**2 + dc[1]**2)) * dps - np.sqrt(3.0) * (-dc[1]) * dc[2]**2 * dpp
                elif orb_j==3: 
                    # z2r2_z    
                    amp = (-dc[2]) * (dc[2]**2 - 0.5 * (dc[0]**2 + dc[1]**2)) * dps + np.sqrt(3.0) * (-dc[2]) * (dc[0]**2 + dc[1]**2) * dpp
                elif orb_j==4: 
                    # z2r2_xy   
                    amp = np.sqrt(3.0)* (-dc[0]) * (-dc[1]) * (dc[2]**2 - 0.5 * (dc[0]**2 + dc[1]**2)) * dds - 2.0 * np.sqrt(3.0) * (-dc[0]) * (-dc[1]) * dc[2]**2 * ddp + 0.5 * np.sqrt(3.0) * (-dc[0]) * (-dc[1]) * (1.0 + dc[2]**2) * ddd
                elif orb_j==5: 
                    # z2r2_yz   
                    amp = np.sqrt(3.0)* (-dc[1]) * (-dc[2]) * (dc[2]**2 - 0.5 * (dc[0]**2 + dc[1]**2)) * dds + np.sqrt(3.0) * (-dc[1]) * (-dc[2]) * (dc[0]**2 + dc[1]**2 - dc[2]**2) * ddp - 0.5 * np.sqrt(3.0) * (-dc[1]) * (-dc[2]) * (dc[0]**2 + dc[1]**2) * ddd
                elif orb_j==6: 
                    # z2r2_zx   
                    amp = np.sqrt(3.0)* (-dc[0]) * (-dc[2]) * (dc[2]**2 - 0.5 * (dc[0]**2 + dc[1]**2)) * dds + np.sqrt(3.0) * (-dc[0]) * (-dc[2]) * (dc[0]**2 + dc[1]**2 - dc[2]**2) * ddp - 0.5 * np.sqrt(3.0) * (-dc[0]) * (-dc[2]) * (dc[0]**2 + dc[1]**2) * ddd
                elif orb_j==7: 
                    # z2r2_x2y2 
                    amp = 0.5 * np.sqrt(3.0) * (dc[0]**2 - dc[1]**2) * (dc[2]**2 - 0.5 * (dc[0]**2 + dc[1]**2)) * dds + np.sqrt(3.0) * dc[2]**2 * (dc[1]**2 - dc[0]**2) * ddp + 0.250 * np.sqrt(3.0) * (1.0 + dc[2]**2) * (dc[0]**2 - dc[1]**2) * ddd
                elif orb_j==8: 
                    # z2r2_z2r2 
                    amp = (dc[2]**2 - 0.5 * (dc[0]**2 + dc[1]**2)) * (dc[2]**2 - 0.5 * (dc[0]**2 + dc[1]**2)) * dds + 3.0 * dc[2]**2 * (dc[0]**2 + dc[1]**2) * ddp + 0.750 * (dc[0]**2 + dc[1]**2) * (dc[0]**2 + dc[1]**2) * ddd
                    
            # set the hopping if not onsite
            abs_d = np.sqrt(np.dot(dist,dist))
            if abs_d>0.0000001:
                model.set_hop(amp, noba_i+i, noba_j+j, lv)
            elif orb_i != orb_j:
                model.set_hop(amp, noba_i+i, noba_j+j, lv)
    
    # return the TB model with the hoppings added
    return model

def set_SOC_onsite_p(model,atoms,soc):
    # set the onsite SOC of p-orbitals for all atoms
    # works only with spin
    # my_model...TB model to add hoppings to
    # atoms......the list of atoms and their Orbitals
    # soc........the list of onsite SOC values
    
    # number of atoms
    n_at = len(atoms)
    
    # dimension
    dim = model._dim_r
    lv = [0]*dim

    # loop over all atoms
    for at in range(n_at):
        # number of orbitals stored below atom at
        noba = 0
        for i in range(at):
            n_orb = len(atoms[i][1])
            noba += n_orb
        
        # all combinations of orbitals
        n_orb = len(atoms[at][1])
        for i in range(n_orb):
            orb_1=atoms[at][1][i]
            for j in range(i+1,n_orb):
                orb_2=atoms[at][1][j]
                
                amp = [0.,0.,0.,0.]
                if orb_1 ==1:
                    if orb_2 ==2:
                        amp = [0,0,0,-1j*soc[at]]
                    elif orb_2 ==3:
                        amp = [0,0,1j*soc[at],0]
                elif orb_1 ==2:
                    if orb_2 ==1:
                        amp = [0,0,0,1j*soc[at]]
                    elif orb_2 ==3:
                        amp = [0,-1j*soc[at],0,0]
                elif orb_1 ==3:
                    if orb_2 ==1:
                        amp = [0,0,-1j*soc[at],0.]
                    elif orb_2 ==2:
                        amp = [0,1j*soc[at],0,0]

                model.set_hop(amp, noba+i, noba+j, lv)
    
    # return the TB model with SOC added
    return model

def set_SOC_onsite_d(model,atoms,soc):
    # set the onsite SOC of d-orbitals for all atoms
    # works only with spin
    # my_model...TB model to add hoppings to
    # atoms......the list of atoms and their orbitals
    # soc........the list of onsite SOC values
    
    # number of atoms
    n_at = len(atoms)
    
    # dimension
    dim = model._dim_r
    lv = [0]*dim

    # loop over all atoms
    for at in range(n_at):
        # number of orbitals stored below atom at
        noba = 0
        for i in range(at):
            n_orb = len(atoms[i][1])
            noba += n_orb
        
        # all combinations of orbitals
        n_orb = len(atoms[at][1])
        for i in range(n_orb):
            orb_1=atoms[at][1][i]
            for j in range(i+1,n_orb):
                orb_2=atoms[at][1][j]
                
                amp = [0.,0.,0.,0.]
                if orb_1 ==4:
                    if orb_2 ==5:
                        amp = [0,0,1j*soc[at],0]
                    elif orb_2 ==6:
                        amp = [0,-1j*soc[at],0,0]
                    elif orb_2 ==7:
                        amp = [0,0,0,2.*1j*soc[at]]
                elif orb_1 ==5:
                    if orb_2 ==4:
                        amp = [0,0,-1j*soc[at],0]
                    elif orb_2 ==6:
                        amp = [0,0,0,1j*soc[at]]
                    elif orb_2 ==7:
                        amp = [0,-1j*soc[at],0,0]
                    elif orb_2 ==8:
                        amp = [0,-1j*np.sqrt(3.0)*soc[at],0,0]
                elif orb_1 ==6:
                    if orb_2 ==4:
                        amp = [0,1j*soc[at],0,0]
                    elif orb_2 ==5:
                        amp = [0,0,0,-1j*soc[at]]
                    if orb_2 ==7:
                        amp = [0,0,-1j*soc[at],0]
                    elif orb_2 ==8:
                        amp = [0,0,1j*np.sqrt(3.0)*soc[at],0]
                elif orb_1 ==7:
                    if orb_2 ==4:
                        amp = [0,0,0,-2.0*1j*soc[at]]
                    elif orb_2 ==5:
                        amp = [0,1j*soc[at],0,0]
                    if orb_2 ==6:
                        amp = [0,0,1j*soc[at],0]
                elif orb_1 ==8:
                    if orb_2 ==5:
                        amp = [0,1j*np.sqrt(3.0)*soc[at],0,0]
                    elif orb_2 ==6:
                        amp = [0,0,-1j*np.sqrt(3.0)*soc[at],0]

                model.set_hop(amp, noba+i, noba+j, lv, mode="add")
    
    # return the TB model with SOC added
    return model
