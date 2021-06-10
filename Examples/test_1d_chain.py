from sk_extension import *
import matplotlib.pyplot as plt

# initialize the 1d chain
    
# set up the TB model using the Slater-Koster parametrization
lat=[[1.0]]
atoms = [[[0.],[0],[0.0]]]

my_model = init_model_SK(1, 1, lat, atoms)

# hoppings
# sss - (ss-sigma)
# sps - (sp-sigma)
# pdp - (pd-pi)
# ddd - (dd-delta)
# etc.

sss = -1.0     

# set hoppings
# 1st neighbours
my_model = set_hop_SK(my_model,atoms, 0, 0, [ 1], sss)

my_model.display()

# Calculate the band structure
# generate list of k-points following a segmented path in the BZ
# list of nodes (high-symmetry points) that will be connected
path=[[-0.5],[0.0],[0.5]]
# labels of the nodes
label=(r'$-\pi/a$',r'$0$',r'$\pi/a$')
# total number of interpolated k-points along the path
nk=100

# call function k_path to construct the actual path
(k_vec,k_dist,k_node)=my_model.k_path(path,nk)

print('---------------------------------------')
print('starting calculation')
print('---------------------------------------')
print('Calculating bands...')

# obtain eigenvalues to be plotted
evals=my_model.solve_all(k_vec)

# figure for bandstructure

fig, ax = plt.subplots()
# specify horizontal axis details
# set range of horizontal axis
ax.set_xlim(k_node[0],k_node[-1])
# put tickmarks and labels at node positions
ax.set_xticks(k_node)
ax.set_xticklabels(label)
# add vertical lines at node positions
for n in range(len(k_node)):
  ax.axvline(x=k_node[n],linewidth=0.5, color='k')
# put title
ax.set_title("1d chain band structure")
ax.set_xlabel("Path in k-space")
ax.set_ylabel("Band energy")

ax.set_ylim(-3,3)
# plot bands
for j in range(len(evals)):
    ax.plot(k_dist,evals[j],color='k')

# make an PDF figure of a plot
fig.tight_layout()
fig.savefig("1d_bands.pdf")

print('Done.\n')
