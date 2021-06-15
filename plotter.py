# Import useful packages
import uproot
import pandas
import awkward as ak
import numpy as np
import matplotlib as mpl
#matplotlib.use('Agg')  # need for not displaying plots when running batch jobs
import matplotlib.pyplot as plt


# ---------- DATA FRAMES ----------- #

# Open root file and its trees
file = uproot.open('/data/jchishol/mntuple_ttbar_6_parton_ejets.root')
tree_truth = file['parton'].arrays()
tree_reco = file['reco'].arrays()

# Create pandas dataframe
df = ak.to_pandas(tree_reco['klfitter_bestPerm_topHad_pt'])          # Add reco pt to df
df.rename(columns={'values':'reco pt'},inplace=True)                 # Rename first column appropriately
df['truth pt'] = ak.to_pandas(tree_truth['MC_thad_afterFSR_pt'])     # Add truth pt to df
df['truth pt'] = df['truth pt']/1000				     # Convert from MeV to GeV for truth pt
df['reco match'] = ak.to_pandas(tree_reco['isMatched'])		     # Add reco isMatched
df['truth match'] = ak.to_pandas(tree_truth['isMatched'])            # Add truth isMatched (might not really need this one)
df['likelihood'] = ak.to_pandas(tree_reco['klfitter_logLikelihood'])     # Add logLikelihood of klfitter

# Cut rows that are not matched between truth and reco
indexNames = df[df['reco match'] == 0 ].index
df.drop(indexNames,inplace=True)

# Calculate the resolution
df['resolution'] = df['truth pt']/df['reco pt']-1

# Cut rows where likelihood <= some number (currently saved as a copy and not the original df)
cutA = df[df['likelihood']>-52]
cutB = df[df['likelihood']>-50]
cutC = df[df['likelihood']>-48]

# ------------- PLOTTING -------------- #

# Choose ticks/binning (may be better choices, this was just in the paper)
ticks = [0,50,100,160,225,300,360,475,1000]
ticks_labels = ['0','50','100','160','225','300','360','475','1000']
n = len(ticks)
ran = ticks[::n-1]

# Create 2D array of truth vs reco pt (which can be plotted also)
#fig,ax = plt.subplots()
#ax.set_xticks(ticks)
#ax.set_yticks(ticks)
H, xedges, yedges = np.histogram2d(df['reco pt'],df['truth pt'],bins=ticks,range=[ran,ran])
#plt.pcolormesh(xedges,yedges,cmap=plt.cm.Blues)
#plt.xlabel('Reco p$_T^{t,had}$ [GeV]')
#plt.ylabel('Truth p$_T^{t,had}$ [GeV]')
#plt.colorbar()

# Normalize the rows out of 100 and round to integers
norm = (np.rint((H.T/np.sum(H,axis=1))*100)).T.astype(int)

# Plot truth vs reco pt with normalized rows
plt.figure('Normalized 2D Plot')
masked_norm = np.ma.masked_where(norm==0,norm)  # Needed to make the zero bins whitee
plt.imshow(masked_norm,extent=[0,8,0,8],cmap=plt.cm.Blues,origin='lower')
plt.xticks(np.arange(n),ticks_labels,fontsize=9,rotation=-25)
plt.yticks(np.arange(n),ticks_labels,fontsize=9)
plt.xlabel('Reco p$_T^{t,had}$ [GeV]',fontsize=12)
plt.ylabel('Truth p$_T^{t,had}$ [GeV]',fontsize=12)
plt.clim(0,100)
plt.colorbar()

# Label the content of each bin
for i in range (8):
       for j in range(8):
		if masked_norm.T[i,j] != 0:   # Don't label empty bins
                	plt.text(i+0.5,j+0.5,masked_norm.T[i,j],color="k",fontsize=6,weight="bold",ha="center",va="center")
plt.savefig('plots/Normalized_2D_thad_pt')


# Create plot of resolution with and without likelihood cut
plt.figure('Resolution')
plt.hist(df['resolution'],bins=100,range=(-1,1),histtype='step')
plt.hist(cutA['resolution'],bins=100,range=(-1,1),histtype='step')
plt.hist(cutB['resolution'],bins=100,range=(-1,1),histtype='step')
plt.hist(cutC['resolution'],bins=100,range=(-1,1),histtype='step')
plt.xlabel('p$_T^{t,had}$ Resolution')
plt.ylabel('Counts')
plt.legend(['No Cuts','logLikelihood > -52','logLikelihood > -50','logLikelihood > -48'])
plt.savefig('plots/Resolution_thad_pt')

# Show plots
plt.show()




print("done :)")
