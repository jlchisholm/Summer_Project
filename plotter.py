# Import useful packages
import uproot
import pandas
import awkward as ak
import numpy as np
import matplotlib
#matplotlib.use('Agg')  # need for not displaying plots when running batch jobs
import matplotlib.pyplot as plt

# List some useful strings
var = ['pt','eta','y','phi','m','E','pout']
labels = ['$p_T^{t,had}$','$\eta^{t,had}$','$y^{t,had}$','$\phi^{t,had}$','$m^{t,had}$','$E^{t,had}$','$p_{out}^{t,had}$']
units = [' [GeV]','','','',' [GeV]',' [GeV]',' [GeV]'] 
truth = 'MC_thad_afterFSR_'
reco = 'klfitter_bestPerm_topHad_'


# ---------- IMPORTING DATA ---------- #

# Open first root file and its trees
file0 = uproot.open('/data/jchishol/mntuple_ttbar_0_parton_ejets.root')
tree_truth0 = file0['parton'].arrays()
tree_reco0 = file0['reco'].arrays()

# Create pandas dataframe
df = ak.to_pandas(tree_reco0['isMatched'])                              # Add reco isMatched to df
df.rename(columns={'values':'reco isMatched'},inplace=True)             # Rename first column appropriately
#df['truth isMatched'] = ak.to_pandas(tree_truth['isMatched'])          # Add truth isMatched (might not really need this one)
df['likelihood'] = ak.to_pandas(tree_reco0['klfitter_logLikelihood'])   # Add logLikelihood of klfitter

# Add variables to be plotted
for i in range(len(var)):
	df['reco '+var[i]] = ak.to_pandas(tree_reco0[reco+var[i]])
	df['truth '+var[i]] = ak.to_pandas(tree_truth0[truth+var[i]])

# Close root file
file0.close()

# Make a function that will append data to the frame
# Input: Data frame and file number (in string form)
def Append_Data(df,n):
	
	# Open root file and its trees
	filen = uproot.open('/data/jchishol/mntuple_ttbar_'+n+'_parton_ejets.root')
	tree_truth_n = filen['parton'].arrays()
	tree_reco_n = filen['reco'].arrays()

	# Create pandas data frame
	df_addon = ak.to_pandas(tree_reco_n['isMatched'])
	df_addon.rename(columns={'values':'reco isMatched'},inplace=True)
	df_addon['likelihood'] = ak.to_pandas(tree_reco_n['klfitter_logLikelihood'])

	# Add variables to be plotted
	for i in range(len(var)):
		df_addon['reco '+var[i]] = ak.to_pandas(tree_reco_n[reco+var[i]])
		df_addon['truth '+var[i]] = ak.to_pandas(tree_truth_n[truth+var[i]])

	# Append data to the main data frame
	df = df.append(df_addon,ignore_index=True)

	# Close root file
	filen.close()

	return df

# Append data from each data file to the main data frame
for n in range(1,8):
	df = Append_Data(df,str(n))

# Cut rows that are not matched between truth and reco
indexNames = df[df['reco isMatched'] == 0 ].index
df.drop(indexNames,inplace=True)

# Convert some truth values from MeV to GeV to match reco units
df['truth pt'] = df['truth pt']/1000
df['truth m'] = df['truth m']/1000
df['truth E'] = df['truth E']/1000
df['truth pout'] = df['truth pout']/1000

df['truth eta'] =df['truth eta']/1000000000    # No idea what's going on here

# Calculate the resolution
for i in range(len(var)):
	df['resolution '+var[i]] = df['truth '+var[i]]/df['reco '+var[i]]-1

# Cut rows where likelihood <= some number (currently saved as a copy and not the original df)
cutA = df[df['likelihood']>-52]
cutB = df[df['likelihood']>-50]
cutC = df[df['likelihood']>-48]


# ------------- 2D PLOTS ------------- #

# Make a function that will make the 2D plots
# Input: index of variable to be plotted and axes ticks
def plot_2D(i,ticks):

	# Define useful constants
	ticks_labels = map(str,ticks)
	n = len(ticks)
	ran = ticks[::n-1]

	# Create 2D array of truth vs reco variable (which can be plotted also)
	#fig,ax = plt.subplots()
	#ax.set_xticks(ticks)
	#ax.set_yticks(ticks)
	H, xedges, yedges = np.histogram2d(cutB['reco '+var[i]],cutB['truth '+var[i]],bins=ticks,range=[ran,ran])
	#plt.pcolormesh(xedges,yedges,cmap=plt.cm.Blues)
	#plt.xlabel('Reco '+labels[i]+units[i])
	#plt.ylabel('Truth '+labels[i]+units[i])
	#plt.colorbar()

	# Normalize the rows out of 100 and round to integers
	norm = (np.rint((H.T/np.sum(H,axis=1))*100)).T.astype(int)

	# Plot truth vs reco pt with normalized rows
	plt.figure(var[i]+' Normalized 2D Plot')
	masked_norm = np.ma.masked_where(norm==0,norm)  # Needed to make the zero bins whitee
	plt.imshow(masked_norm,extent=[0,n-1,0,n-1],cmap=plt.cm.Blues,origin='lower')
	plt.xticks(np.arange(n),ticks_labels,fontsize=9,rotation=-25)
	plt.yticks(np.arange(n),ticks_labels,fontsize=9)
	plt.xlabel('Reco '+labels[i]+units[i])
	plt.ylabel('Truth '+labels[i]+units[i])
	plt.clim(0,100)
	plt.colorbar()

	# Label the content of each bin
	for j in range (n-1):
	       for k in range(n-1):
			if masked_norm.T[j,k] != 0:   # Don't label empty bins
	                	plt.text(j+0.5,k+0.5,masked_norm.T[j,k],color="k",fontsize=6,weight="bold",ha="center",va="center")
	plt.savefig('plots/Normalized_2D_thad_'+var[i],bbox_inches='tight')

# Plot pt
#ticks = [0,50,100,160,225,300,360,475,1000]   # Choose ticks/binning (may be better choices, this was just in the paper)
#ticks_labels = ['0','50','100','160','225','300','360','475','1000']
ticks_pt = range(0,450,50)
plot_2D(0,ticks_pt)

# Plot eta
ticks_eta = np.arange(-2,2.5,0.5)
plot_2D(1,ticks_eta)

# Plot y
ticks_y = np.arange(-2,2.5,0.5) 
plot_2D(2,ticks_y)

# Plot phi
ticks_phi = np.arange(-3,3.5,0.5)
plot_2D(3,ticks_phi)

# Plot m
ticks_m = range(80,260,20)
plot_2D(4,ticks_m)

# Plot E
ticks_E = range(100,1000,100)
plot_2D(5,ticks_E)

# Plot pout
ticks_pout = range(-200,250,50)
plot_2D(6,ticks_pout)


# ---------- RESOLUTION PLOTS ---------- #

# Create plots of resolution with and without likelihood cut
for i in range(len(var)):
	plt.figure(var[i]+' Resolution')
	plt.hist(df['resolution '+var[i]],bins=100,range=(-1,1),histtype='step')
	plt.hist(cutA['resolution '+var[i]],bins=100,range=(-1,1),histtype='step')
	plt.hist(cutB['resolution '+var[i]],bins=100,range=(-1,1),histtype='step')
	plt.hist(cutC['resolution '+var[i]],bins=100,range=(-1,1),histtype='step')
	plt.xlabel(labels[i]+' Resolution')
	plt.ylabel('Counts')
	plt.legend(['No Cuts','logLikelihood > -52','logLikelihood > -50','logLikelihood > -48'])
	plt.savefig('plots/Resolution_thad_'+var[i],bbox_inches='tight')

# Show plots
#plt.show()


print("done :)")
