from sk_extension import *
import matplotlib.pyplot as plt

# initialize the Cu sp3d5 model fitted by Tomas Rauch
    
# set up the TB model using the Slater-Koster parametrization
lat=[[0.5,0.5,0.0],[0.0,0.5,0.5],[0.5,0.0,0.5]]
atoms = [[[0.,0.,0.],[0,1,2,3,4,5,6,7,8],[3.06,10.66,10.66,10.66,-2.648,-2.648,-2.648,-2.648,-2.648]]]

my_model = init_model_SK(3, 3, lat, atoms)

# hoppings
# sss - (ss-sigma)
# sps - (sp-sigma)
# pdp - (pd-pi)
# ddd - (dd-delta)
# etc.

sss1 = -1.0228800000000000     
sps1 =  1.5743100000000001     
pss1 =  1.5743100000000001     
pps1 =  2.6760999999999999     
ppp1 =  0.26395000000000002     
sds1 = -0.42272888180000001     
pds1 = -0.44749124309999999     
pdp1 =  0.23850779840000000     
dss1 = -0.42272888180000001     
dps1 = -0.44749124309999999     
dpp1 =  0.23850779840000000     
dds1 = -0.34912208259999999     
ddp1 =  0.24490247420000000     
ddd1 = -5.55112275000000030E-002

sss2 = -1.25172375999999998E-002
sps2 =  0.16612551160000000     
pss2 =  0.16612551160000000     
pps2 =  0.73321079629999997     
ppp2 =  0.11510416290000000     
sds2 = -0.11592000000000000     
pds2 = -7.29299999999999948E-002
pdp2 =  4.36700000000000005E-002
dss2 = -0.11592000000000000     
dps2 = -7.29299999999999948E-002
dpp2 =  4.36700000000000005E-002
dds2 = -6.13599999999999979E-002
ddp2 =  3.27897201999999971E-002
ddd2 = -3.94565099999999967E-003

# set hoppings
# 1st neighbours
my_model = set_hop_SK(my_model,atoms, 0, 0, [ 1, 0, 0], sss1,sps1,pss1,pps1,ppp1,sds=sds1,pds=pds1,pdp=pdp1,dss=dss1,dps=dps1,dpp=dpp1,dds=dds1,ddp=ddp1,ddd=ddd1)
my_model = set_hop_SK(my_model,atoms, 0, 0, [ 0, 1, 0], sss1,sps1,pss1,pps1,ppp1,sds=sds1,pds=pds1,pdp=pdp1,dss=dss1,dps=dps1,dpp=dpp1,dds=dds1,ddp=ddp1,ddd=ddd1)
my_model = set_hop_SK(my_model,atoms, 0, 0, [ 0, 0, 1], sss1,sps1,pss1,pps1,ppp1,sds=sds1,pds=pds1,pdp=pdp1,dss=dss1,dps=dps1,dpp=dpp1,dds=dds1,ddp=ddp1,ddd=ddd1)
my_model = set_hop_SK(my_model,atoms, 0, 0, [ 1,-1, 0], sss1,sps1,pss1,pps1,ppp1,sds=sds1,pds=pds1,pdp=pdp1,dss=dss1,dps=dps1,dpp=dpp1,dds=dds1,ddp=ddp1,ddd=ddd1)
my_model = set_hop_SK(my_model,atoms, 0, 0, [ 0, 1,-1], sss1,sps1,pss1,pps1,ppp1,sds=sds1,pds=pds1,pdp=pdp1,dss=dss1,dps=dps1,dpp=dpp1,dds=dds1,ddp=ddp1,ddd=ddd1)
my_model = set_hop_SK(my_model,atoms, 0, 0, [-1, 0, 1], sss1,sps1,pss1,pps1,ppp1,sds=sds1,pds=pds1,pdp=pdp1,dss=dss1,dps=dps1,dpp=dpp1,dds=dds1,ddp=ddp1,ddd=ddd1)

# 2nd neighbours
my_model = set_hop_SK(my_model,atoms, 0, 0, [ 1,-1, 1], sss2,sps2,pss2,pps2,ppp2,sds=sds2,pds=pds2,pdp=pdp2,dss=dss2,dps=dps2,dpp=dpp2,dds=dds2,ddp=ddp2,ddd=ddd2)
my_model = set_hop_SK(my_model,atoms, 0, 0, [-1, 1, 1], sss2,sps2,pss2,pps2,ppp2,sds=sds2,pds=pds2,pdp=pdp2,dss=dss2,dps=dps2,dpp=dpp2,dds=dds2,ddp=ddp2,ddd=ddd2)
my_model = set_hop_SK(my_model,atoms, 0, 0, [ 1, 1,-1], sss2,sps2,pss2,pps2,ppp2,sds=sds2,pds=pds2,pdp=pdp2,dss=dss2,dps=dps2,dpp=dpp2,dds=dds2,ddp=ddp2,ddd=ddd2)


# add SOC
soc_d = [0.0579]
my_model = set_SOC_onsite_d(my_model,atoms,soc_d)

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
ax.set_title("Cu band structure")
ax.set_xlabel("Path in k-space")
ax.set_ylabel("Band energy")

ax.set_ylim(-10,8)
# plot bands
for j in range(18):
    ax.plot(k_dist,evals[j],color='k')

# make an PDF figure of a plot
fig.tight_layout()
fig.savefig("Cu_bands.pdf")

print('Done.\n')
