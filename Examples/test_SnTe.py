from sk_extension import *
import matplotlib.pyplot as plt

# initialize the SnTe sp3 model from Lent et al
    
# set up the TB model using the Slater-Koster parametrization
lat=[[0.5,0.5,0.0],[0.0,0.5,0.5],[0.5,0.0,0.5]]
atoms = [[[0.,0.,0.],[0,1,2,3],[-6.578,1.659,1.659,1.659]],[[0.5,-0.5,0.5],[0,1,2,3],[-12.067,-0.167,-0.167,-0.167]]]

my_model = init_model_SK(3, 3, lat, atoms)

# hoppings
# sss - (ss-sigma)
# sps - (sp-sigma)
# pdp - (pd-pi)
# ddd - (dd-delta)
# etc.

sss = -0.510
sps = 0.949
pss = -0.198
pps = 2.218
ppp = -0.446

# set hoppings
# 1st neighbours
# x
my_model = set_hop_SK(my_model,atoms, 0, 1, [ 0, 0, 0], sss,sps,pss,pps,ppp)
my_model = set_hop_SK(my_model,atoms, 0, 1, [ -1, 1, -1], sss,sps,pss,pps,ppp)
# y
my_model = set_hop_SK(my_model,atoms, 0, 1, [ -1, 0, 0], sss,sps,pss,pps,ppp)
my_model = set_hop_SK(my_model,atoms, 0, 1, [ 0, 1, -1], sss,sps,pss,pps,ppp)
# z
my_model = set_hop_SK(my_model,atoms, 0, 1, [ -1, 1, 0], sss,sps,pss,pps,ppp)
my_model = set_hop_SK(my_model,atoms, 0, 1, [ 0, 0, -1], sss,sps,pss,pps,ppp)

# add SOC
soc = [0.296,0.282]
my_model = set_SOC_onsite_p(my_model,atoms,soc)

my_model.display()

# Calculate the band structure
# generate list of k-points following a segmented path in the BZ
# list of nodes (high-symmetry points) that will be connected
path=[[0.,0.,0.],[.5,.5,0.],[3./4.,.5,0.25],[.5,.5,.5],[0.,0.,0.],[3./4.,3./8.,3./8.]]
# labels of the nodes
label=(r'$\Gamma $', r'$X$', r'$W$',r'$L$', r'$\Gamma $',r'$K$')
# total number of interpolated k-points along the path
nk=500

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
ax.set_title("SnTe band structure")
ax.set_xlabel("Path in k-space")
ax.set_ylabel("Band energy")

ax.set_ylim(-2,2)
# plot bands
for j in range(12):
    ax.plot(k_dist,evals[j],color='k')

# make an PDF figure of a plot
fig.tight_layout()
fig.savefig("SnTe_bands.pdf")

print('Done.\n')
