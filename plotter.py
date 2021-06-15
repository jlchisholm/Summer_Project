# Import useful packages
import uproot
import pandas
import awkward as ak
import matplotlib
#matplotlib.use('Agg')  # need for not displaying plots when running batch jobs
import matplotlib.pyplot as plt

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

# Cut rows where likelihood <= -50 (currently saved as a copy and not the original df)
cut = df[df['likelihood']>-50]

# Printing for testing
# print df
# print cut

# Create 2D plot of truth vs reco pt
plt.figure('2D')
plt.hist2d(df['reco pt'],df['truth pt'],bins=10,cmap=plt.cm.Blues)
plt.xlabel('Reco p$_T^{t,had}$ [GeV]')
plt.ylabel('Truth p$_T^{t,had}$ [GeV]')
plt.colorbar()

# Create plot of resolution with and without likelihood cut
plt.figure('Resolution')
plt.hist(df['resolution'],bins=100,range=(-1,1),histtype='step')
plt.hist(cut['resolution'],bins=100,range=(-1,1),histtype='step')
plt.xlabel('p$_T^{t,had}$ Resolution')
plt.ylabel('Counts')
plt.legend(['No Cuts','logLikelihood > -50'])
plt.show()



# SCRAP WORK

#plt.subplot(1,3,1)
#plt.hist(branches1['MC_thad_afterFSR_pt'],bins=10)
#plt.subplot(1,3,2)
#plt.hist(branches2['klfitter_bestPerm_topHad_pt'],bins=10)
#plt.subplot(1,3,3)
#plt.hist([branches1['MC_thad_afterFSR_pt'],branches2['klfitter_bestPerm_topHad_pt']],bins=10)

#plt.figure()
#(n, bins, patches) = plt.hist(branches1['MC_thad_afterFSR_pt'],bins=10)
#print n

#plt.figure()
#plt.hist2d(x,y,bins=8,range=[[0,400],[0,400]],cmap=plt.cm.Blues)

#ticks = [0,50,100,160,225,300,360,475,1000]
#(fig, ax) = plt.subplots()
#ax.set_aspect("equal")
#(hist, xbins, ybins, im) = plt.hist2d(x,y,bins=8,cmap=plt.cm.Blues) # returns 2D array of values (I think?), x bin edges, y bin edges, and image
#print hist
#plt.xticks(x,ticks)
#plt.yticks(y,ticks)
#plt.xlabel('klfitter_bestPerm_topHad_pt')
#plt.ylabel('PseudoTop_Reco_top_had_pt')
#plt.colorbar()
#for i in range (len(xbins)-1):
#	for j in range(len(ybins)-1):
		#ax.text(xbins[i]+25,ybins[j]+25,round(hist.T[i,j]/60000,2),color="k",ha="center",va="center") # x position, y position, string of text, etc
#		ax.text(xbins[i]+25,ybins[j]+25,i,color="k",ha="center",va="center")
#plt.show()
#plt.savefig('/home/jchishol/Summer_Project/test_plot.png')


print("done :)")
