# SKext2PythTB
Extension to PythTB enabling to define TB models with Slater-Koster parameters

PythTB is a package for electronic structure calculations in the tight-binfing approximation.
See http://www.physics.rutgers.edu/pythtb

In the Slater-Koster approximation, the hopping parameters of the tight-binding approximation
are expressed in terms of orbital overlaps (e.g., (sp-sigma), or (pp-pi)) and the direction 
cosine relating the respective positions of the orbitals and their orientation in the Cartesian
coordinates.
See J. C. Slater and G. F. Koster, Simplified LCAO method for the periodic potential problem, 
Phys. Rev. 94, 1498 (1954).

## sk_extension.py
Includes the routines necessary to define the Hamilton matrix with the Slater-Koster parameters.

### init_model_SK(dim_k,dim_r,lat,atoms)
Initializes and returns the TB model defined with pythtb.tb_model
**dim_k**: Dimensionality of reciprocal space, same as in pythtb_tb_model.
**dim_r**: Dimensionality of real space, same as in pythtb_tb_model.
**lat**: Array containing lattice vectors in Cartesian coordinates, same as in pythtb_tb_model
**atoms**: List containing reduced coordinates of all atoms in the unit cell, the types of orbitals localized at the atom position, and their onsite energies. The orbitals are: 0=s, 1=p-x, 2=p-y, 3=p-z, 4=d-xy, 5=d-yz, 6=d-zx, 7=d-x^2-y^2, 8=d-3r^2-z^2.

### set_hop_SK(model,atoms,n_at_i,n_at_j,lv,sss=0.,sps=0.,pss=0.,pps=0.,ppp=0.,sds=0.,pds=0.,pdp=0.,dss=0.,dps=0.,dpp=0.,dds=0.,ddp=0.,ddd=0.)
Set the hoppings according to the Slater-Koster parametrization.
**model**: TB model to add hoppings to
**atoms**: the list of atoms and their Orbitals
**n_at_i**: number of atom_i in the list
**n_at_j**: number of atom_j in the list
**lv**: lattice vector connecting the atoms
**SK parameters**:
  sss......ss-sigma
  sps......sp-sigma
  pss......ps-sigma
  pps......pp-sigma
  ppp......pp-pi
  sds......sd-sigma
  pds......pd-sigma
  pdp......pd-pi
  dss......ds-sigma
  dps......dp-sigma
  dpp......dp-pi
  dds......dd-sigma
  ddp......dd-pi
  ddd......dd-delta

### set_SOC_onsite_p(model,atoms,soc)
Set the onsite SOC of p-orbitals for all atoms.
Works only with spin.
**model**: TB model to add hoppings to
**atoms**: the list of atoms and their orbitals
**soc**: the list of onsite SOC values

### set_SOC_onsite_d(model,atoms,soc)
Set the onsite SOC of d-orbitals for all atoms.
Works only with spin.
**model**: TB model to add hoppings to
**atoms**: the list of atoms and their orbitals
**soc**: the list of onsite SOC values

## Examples

Includes four examples for the usage of SKext2PythTB. PythTB must be installed or pythtb.py must be present in the directory together with sk_extension.py.

### 1d chain

test_1d_chain.py

1d chain, single *s* orbital in the unit cell. 

On-site energy $E_s=0.0$, nearest-neighbor hopping $(ss\sigma)=-1.0$.

Band structure shown in 1d_bands.pdf.

### FCC Cu

test_Cu.py

Cu in the FCC structure. 

$sp^3d^5$ parametrization.

Nearest and next-nearest neighbors for the hoppings.

Band structure shown in Cu_bands.pdf.

### HgTe

test_HgTe.py

HgTe in the zincblende structure.

$sp^3d^5$ parametrization.

Nearest and next-nearest neighbors for the hoppings.

Parameters from T. Rauch *et al.*, PRL **114**, 236805 (2015).

Band structure shown in HgTe_bands.pdf.

### SnTe

test_SnTe.py

SnTe in the rocksalt structure. Topological crystalline insulator.

$sp^3$ parametrization.

Nearest neighbors for the hoppings.

Parameters from C. Lent *et al.*, Superlattices and Microstructures **2**, 491 (1986).

Band structure shown in SnTe_bands.pdf.

### 

### 