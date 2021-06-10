from sk_extension import *
import matplotlib.pyplot as plt

# initialize the HgTe sp3d5 model fitted by Tomas Rauch
    
# set up the TB model using the Slater-Koster parametrization
lat=[[0.5,0.5,0.0],[0.0,0.5,0.5],[0.5,0.0,0.5]]
atoms = [[[0.,0.,0.],[0,1,2,3,4,5,6,7,8],[0.1441,4.7989,4.7989,4.7989,-6.2819,-6.2819,-6.2819,-6.2819,-6.2819]],[[0.25,0.25,0.25],[0,1,2,3,4,5,6,7,8],[-10.1191,-1.2874,-1.2874,-1.2874,9.7816,9.7816,9.7816,9.7816,9.7816]]]

my_model = init_model_SK(3, 3, lat, atoms)

# hoppings
# sss - (ss-sigma)
# sps - (sp-sigma)
# pdp - (pd-pi)
# ddd - (dd-delta)
# etc.

# 1st neighbours Hg->Te
sss1 = 0.652875363768169
sps1 =  1.58154985033517
pss1 = 0.138917060362976
pps1 = 1.5302
ppp1 =-0.28968
sds1 =  1.71851381114842     
pds1 =-8.874398770073622E-003
pdp1 =  1.01816156748847     
dss1 =-0.678772036896918     
dps1 = 0.864895327149847     
dpp1 = 1.383208024574967E-002
dds1 =-0.903921527244890     
ddp1 =-0.985108201035476     
ddd1 =-5.190033671149452E-003

# 2nd neighbours Hg->Hg
sss2_Hg_Hg =-9.006034776345326E-002
sps2_Hg_Hg =-0.100012318930304    
pss2_Hg_Hg =-0.100012318930304    
pps2_Hg_Hg = 0.434689044462726    
ppp2_Hg_Hg =-0.279487808587277    
sds2_Hg_Hg =-0.111597972915523     
pds2_Hg_Hg =-0.404462077133189     
pdp2_Hg_Hg =-0.156317681440761     
dss2_Hg_Hg =-0.111597972915523     
dps2_Hg_Hg =-0.404462077133189     
dpp2_Hg_Hg =-0.156317681440761     
dds2_Hg_Hg =-9.360294973299493E-002
ddp2_Hg_Hg = 5.000878217200200E-002
ddd2_Hg_Hg = 3.114273894912767E-002

# 2nd neighbours Te->Te
sss2_Te_Te =-3.927452645188628E-002
sps2_Te_Te = 0.380944817503874     
pss2_Te_Te = 0.380944817503874     
pps2_Te_Te = 0.279890091293738     
ppp2_Te_Te = 1.045782905484532E-002
sds2_Te_Te =-0.605144766828820     
pds2_Te_Te = 0.703366721880547     
pdp2_Te_Te = 9.667287823598020E-002
dss2_Te_Te =-0.605144766828820     
dps2_Te_Te = 0.703366721880547     
dpp2_Te_Te = 9.667287823598020E-002
dds2_Te_Te = 0.658708072029976     
ddp2_Te_Te =-0.971938310011725     
ddd2_Te_Te = 0.114760749393371     

# set hoppings
# 1st neighbours
my_model = set_hop_SK(my_model,atoms, 0, 1, [ 0, 0, 0], sss1,sps1,pss1,pps1,ppp1,sds=sds1,pds=pds1,pdp=pdp1,dss=dss1,dps=dps1,dpp=dpp1,dds=dds1,ddp=ddp1,ddd=ddd1)
my_model = set_hop_SK(my_model,atoms, 0, 1, [ 0, 0,-1], sss1,sps1,pss1,pps1,ppp1,sds=sds1,pds=pds1,pdp=pdp1,dss=dss1,dps=dps1,dpp=dpp1,dds=dds1,ddp=ddp1,ddd=ddd1)
my_model = set_hop_SK(my_model,atoms, 0, 1, [-1, 0, 0], sss1,sps1,pss1,pps1,ppp1,sds=sds1,pds=pds1,pdp=pdp1,dss=dss1,dps=dps1,dpp=dpp1,dds=dds1,ddp=ddp1,ddd=ddd1)
my_model = set_hop_SK(my_model,atoms, 0, 1, [ 0,-1, 0], sss1,sps1,pss1,pps1,ppp1,sds=sds1,pds=pds1,pdp=pdp1,dss=dss1,dps=dps1,dpp=dpp1,dds=dds1,ddp=ddp1,ddd=ddd1)

# 2nd neighbours Hg->Hg
my_model = set_hop_SK(my_model,atoms, 0, 0, [ 0, 1, 0], sss2_Hg_Hg,sps2_Hg_Hg,pss2_Hg_Hg,pps2_Hg_Hg,ppp2_Hg_Hg,sds=sds2_Hg_Hg,pds=pds2_Hg_Hg,pdp=pdp2_Hg_Hg,dss=dss2_Hg_Hg,dps=dps2_Hg_Hg,dpp=dpp2_Hg_Hg,dds=dds2_Hg_Hg,ddp=ddp2_Hg_Hg,ddd=ddd2_Hg_Hg)
my_model = set_hop_SK(my_model,atoms, 0, 0, [ 0, 0, 1], sss2_Hg_Hg,sps2_Hg_Hg,pss2_Hg_Hg,pps2_Hg_Hg,ppp2_Hg_Hg,sds=sds2_Hg_Hg,pds=pds2_Hg_Hg,pdp=pdp2_Hg_Hg,dss=dss2_Hg_Hg,dps=dps2_Hg_Hg,dpp=dpp2_Hg_Hg,dds=dds2_Hg_Hg,ddp=ddp2_Hg_Hg,ddd=ddd2_Hg_Hg)
my_model = set_hop_SK(my_model,atoms, 0, 0, [ 1, 0, 0], sss2_Hg_Hg,sps2_Hg_Hg,pss2_Hg_Hg,pps2_Hg_Hg,ppp2_Hg_Hg,sds=sds2_Hg_Hg,pds=pds2_Hg_Hg,pdp=pdp2_Hg_Hg,dss=dss2_Hg_Hg,dps=dps2_Hg_Hg,dpp=dpp2_Hg_Hg,dds=dds2_Hg_Hg,ddp=ddp2_Hg_Hg,ddd=ddd2_Hg_Hg)
my_model = set_hop_SK(my_model,atoms, 0, 0, [ 1, 0,-1], sss2_Hg_Hg,sps2_Hg_Hg,pss2_Hg_Hg,pps2_Hg_Hg,ppp2_Hg_Hg,sds=sds2_Hg_Hg,pds=pds2_Hg_Hg,pdp=pdp2_Hg_Hg,dss=dss2_Hg_Hg,dps=dps2_Hg_Hg,dpp=dpp2_Hg_Hg,dds=dds2_Hg_Hg,ddp=ddp2_Hg_Hg,ddd=ddd2_Hg_Hg)
my_model = set_hop_SK(my_model,atoms, 0, 0, [-1, 1, 0], sss2_Hg_Hg,sps2_Hg_Hg,pss2_Hg_Hg,pps2_Hg_Hg,ppp2_Hg_Hg,sds=sds2_Hg_Hg,pds=pds2_Hg_Hg,pdp=pdp2_Hg_Hg,dss=dss2_Hg_Hg,dps=dps2_Hg_Hg,dpp=dpp2_Hg_Hg,dds=dds2_Hg_Hg,ddp=ddp2_Hg_Hg,ddd=ddd2_Hg_Hg)
my_model = set_hop_SK(my_model,atoms, 0, 0, [ 0, 1,-1], sss2_Hg_Hg,sps2_Hg_Hg,pss2_Hg_Hg,pps2_Hg_Hg,ppp2_Hg_Hg,sds=sds2_Hg_Hg,pds=pds2_Hg_Hg,pdp=pdp2_Hg_Hg,dss=dss2_Hg_Hg,dps=dps2_Hg_Hg,dpp=dpp2_Hg_Hg,dds=dds2_Hg_Hg,ddp=ddp2_Hg_Hg,ddd=ddd2_Hg_Hg)


# 2nd neighbours Te->Te
my_model = set_hop_SK(my_model,atoms, 1, 1, [ 0, 1, 0], sss2_Te_Te,sps2_Te_Te,pss2_Te_Te,pps2_Te_Te,ppp2_Te_Te,sds=sds2_Te_Te,pds=pds2_Te_Te,pdp=pdp2_Te_Te,dss=dss2_Te_Te,dps=dps2_Te_Te,dpp=dpp2_Te_Te,dds=dds2_Te_Te,ddp=ddp2_Te_Te,ddd=ddd2_Te_Te)
my_model = set_hop_SK(my_model,atoms, 1, 1, [ 0, 0, 1], sss2_Te_Te,sps2_Te_Te,pss2_Te_Te,pps2_Te_Te,ppp2_Te_Te,sds=sds2_Te_Te,pds=pds2_Te_Te,pdp=pdp2_Te_Te,dss=dss2_Te_Te,dps=dps2_Te_Te,dpp=dpp2_Te_Te,dds=dds2_Te_Te,ddp=ddp2_Te_Te,ddd=ddd2_Te_Te)
my_model = set_hop_SK(my_model,atoms, 1, 1, [ 1, 0, 0], sss2_Te_Te,sps2_Te_Te,pss2_Te_Te,pps2_Te_Te,ppp2_Te_Te,sds=sds2_Te_Te,pds=pds2_Te_Te,pdp=pdp2_Te_Te,dss=dss2_Te_Te,dps=dps2_Te_Te,dpp=dpp2_Te_Te,dds=dds2_Te_Te,ddp=ddp2_Te_Te,ddd=ddd2_Te_Te)
my_model = set_hop_SK(my_model,atoms, 1, 1, [ 1, 0,-1], sss2_Te_Te,sps2_Te_Te,pss2_Te_Te,pps2_Te_Te,ppp2_Te_Te,sds=sds2_Te_Te,pds=pds2_Te_Te,pdp=pdp2_Te_Te,dss=dss2_Te_Te,dps=dps2_Te_Te,dpp=dpp2_Te_Te,dds=dds2_Te_Te,ddp=ddp2_Te_Te,ddd=ddd2_Te_Te)
my_model = set_hop_SK(my_model,atoms, 1, 1, [-1, 1, 0], sss2_Te_Te,sps2_Te_Te,pss2_Te_Te,pps2_Te_Te,ppp2_Te_Te,sds=sds2_Te_Te,pds=pds2_Te_Te,pdp=pdp2_Te_Te,dss=dss2_Te_Te,dps=dps2_Te_Te,dpp=dpp2_Te_Te,dds=dds2_Te_Te,ddp=ddp2_Te_Te,ddd=ddd2_Te_Te)
my_model = set_hop_SK(my_model,atoms, 1, 1, [ 0, 1,-1], sss2_Te_Te,sps2_Te_Te,pss2_Te_Te,pps2_Te_Te,ppp2_Te_Te,sds=sds2_Te_Te,pds=pds2_Te_Te,pdp=pdp2_Te_Te,dss=dss2_Te_Te,dps=dps2_Te_Te,dpp=dpp2_Te_Te,dds=dds2_Te_Te,ddp=ddp2_Te_Te,ddd=ddd2_Te_Te)

# add SOC
soc_p = [0.52, 0.34]
soc_d = [0.31, 0.22]
my_model = set_SOC_onsite_p(my_model,atoms,soc_p)
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
ax.set_title("HgTe band structure")
ax.set_xlabel("Path in k-space")
ax.set_ylabel("Band energy")

ax.set_ylim(-2,2)
# plot bands
for j in range(int(len(evals[:,0]))):
    ax.plot(k_dist,evals[j],color='k')

# make an PDF figure of a plot
fig.tight_layout()
fig.savefig("HgTe_bands.pdf")

print('Done.\n')
