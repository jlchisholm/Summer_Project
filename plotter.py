# Import useful packages
import uproot
import pandas
import awkward as ak
import numpy as np
import matplotlib
#matplotlib.use('Agg')  # need for not displaying plots when running batch jobs
from matplotlib import pyplot as plt
from matplotlib import colors

# List some useful strings
var = ['pt','eta','y','phi','m','E','pout']
labels = ['$p_T^{t,had}$','$\eta^{t,had}$','$y^{t,had}$','$\phi^{t,had}$','$m^{t,had}$','$E^{t,had}$','$p_{out}^{t,had}$']
units = [' [GeV]','','','',' [GeV]',' [GeV]',' [GeV]'] 
truth = 'MC_thad_afterFSR_'
reco = 'klfitter_bestPerm_topHad_'


# ---------- FUNCTION DEFINITIONS ---------- #

# Imports data from 0th file and creates a data frame
# Input: type of data (e.g. parton_ejets)
# Output: data frame
def Create_DF(name): 
	
	# Open first root file and its trees
	file0 = uproot.open('/data/jchishol/mntuple_ttbar_0_'+name+'.root')
	tree_truth0 = file0['parton'].arrays()
	tree_reco0 = file0['reco'].arrays()

	# Create pandas dataframe
	df = ak.to_pandas(tree_reco0['isMatched'])                              # Add reco isMatched to df
	df.rename(columns={'values':'reco isMatched'},inplace=True)             # Rename first column appropriately
	#df['truth isMatched'] = ak.to_pandas(tree_truth0['isMatched'])          # Add truth isMatched (might not really need this one)
	df['likelihood'] = ak.to_pandas(tree_reco0['klfitter_logLikelihood'])   # Add logLikelihood of klfitter

	# Add variables to be plotted
	for i in range(len(var)):
        	df['reco '+var[i]] = ak.to_pandas(tree_reco0[reco+var[i]])
        	df['truth '+var[i]] = ak.to_pandas(tree_truth0[truth+var[i]])

	# Close root file
	file0.close()
	
	# Return main data frame
	return df


# Appends data from the rest of the files to the existing data frame
# Input: data frame, file number (in string form), and type of data
# Output data frame
def Append_Data(df,n,name):

        # Open root file and its trees
        filen = uproot.open('/data/jchishol/mntuple_ttbar_'+n+'_'+name+'.root')
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

	print 'Appended file '+name+'_'+n

        return df


# Drops unmatched data, converts units, and calculates resolution
# Input: data frame
# Output: data frame
def Edit_Data(df):

	# Cut rows that are not matched between truth and reco
	indexNames = df[df['reco isMatched'] == 0 ].index
	df.drop(indexNames,inplace=True)

	# Convert some truth values from MeV to GeV to match reco units
	df['truth pt'] = df['truth pt']/1000
	df['truth m'] = df['truth m']/1000
	df['truth E'] = df['truth E']/1000
	df['truth pout'] = df['truth pout']/1000

	# Calculate the resolution
	#for i in range(len(var)):
        #	df['resolution '+var[i]] = df['truth '+var[i]]/df['reco '+var[i]]-1
	df['resolution pout'] = (df['truth pout']/df['reco pout'])-1

	return df


# Creates and saves the 2D histogram normalized by rows
# Input: cutB of df, variable index, chosen bins/ticks, and data type
# Output: saves histogram as png
def Plot_2D(cutB,i,ticks,ticks_labels,name):

	# Define useful constants
        n = len(ticks)
        ran = ticks[::n-1]

        # Create 2D array of truth vs reco variable (which can be plotted also)
        #plt.figure('2D')
	#fig,ax = plt.subplots()
        #ax.set_xticks(ticks)
        #ax.set_yticks(ticks)
        H, xedges, yedges, im = plt.hist2d(np.clip(cutB['reco '+var[i]],ticks[0],ticks[-1]),np.clip(cutB['truth '+var[i]],ticks[0],ticks[-1]),bins=ticks,range=[ran,ran])
        #plt.pcolormesh(xedges,yedges,cmap=plt.cm.Blues)
        #plt.xlabel('Reco '+labels[i]+units[i])
        #plt.ylabel('Truth '+labels[i]+units[i])

	# Normalize the rows out of 100 and round to integers
	norm = (np.rint((H/np.sum(H.T,axis=1))*100)).T.astype(int)

        # Plot truth vs reco pt with normalized rows
        plt.figure(name+' '+var[i]+' Normalized 2D Plot')
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
        
	# Save the figure as a png
	plt.savefig('plots/Normalized_2D_'+name+'_thad_'+var[i],bbox_inches='tight')
	print 'Saved Figure: Normalized_2D_'+name+'_thad_'+var[i]
	

# Defines bins/ticks and implements Plot_2D for each variable
# Input: cutB of df, and data type
def Create_2D_Plots(cutB, name):
	
	# pt
	#ticks = [0,50,100,160,225,300,360,475,1000]   # Choose ticks/binning (may be better choices, this was just in the paper)
	#ticks_labels = ['0','50','100','160','225','300','360','475','1000']
	ticks_pt = range(0,500,50)
	ticks_labels_pt = map(str,ticks_pt)
        ticks_labels_pt[-1] = r'$\infty$'
	Plot_2D(cutB,0,ticks_pt,ticks_labels_pt,name)

	# eta
	ticks_eta = np.arange(-6,6.5,1)
	ticks_labels_eta = map(str,ticks_eta)
	ticks_labels_eta[0] = '-6.28'
	ticks_labels_eta[-1] = '6.28'
	Plot_2D(cutB,1,ticks_eta,ticks_labels_eta,name)
	
	# y
	ticks_y = np.arange(-2.5,3,0.5)
	ticks_labels_y = map(str,ticks_y)
	ticks_labels_y[0] = '-'+r'$\infty$'
	ticks_labels_y[-1] = r'$\infty$'
	Plot_2D(cutB,2,ticks_y,ticks_labels_y,name)

	# phi
	ticks_phi = np.arange(-3,3.5,0.5)
	ticks_labels_phi = map(str,ticks_phi)
	ticks_labels_phi[0] = '-3.14'
	ticks_labels_phi[-1] = '3.14'
	Plot_2D(cutB,3,ticks_phi,ticks_labels_phi,name)

	# m
	ticks_m = range(100,260,20)
	ticks_labels_m = map(str,ticks_m)
	ticks_labels_m[0] = '-'+r'$\infty$'
	ticks_labels_m[-1] = r'$\infty$'
	Plot_2D(cutB,4,ticks_m,ticks_labels_m,name)

	# E
	ticks_E = range(100,1000,100)
	ticks_labels_E = map(str,ticks_E)
	ticks_labels_E[-1] = r'$\infty$'
	Plot_2D(cutB,5,ticks_E,ticks_labels_E,name)

	# pout
	ticks_pout = range(-250,300,50)
	ticks_labels_pout = map(str,ticks_pout)
	ticks_labels_pout[0] = '-'+r'$\infty$'
	ticks_labels_pout[-1] = r'$\infty$'
	Plot_2D(cutB,6,ticks_pout,ticks_labels_pout,name)


# Creates and saves the resolution plots for each variable
# Input: data frame, cutA of df, cutB, cutC, and data type
def Plot_Res(df,cutA,cutB,cutC,name):
	
	# Create each resolution plot
	for i in range(len(var)):
        	plt.figure(name+' '+var[i]+' Resolution')
        	plt.hist(df['resolution '+var[i]],bins=100,range=(-1,1),histtype='step')
        	plt.hist(cutA['resolution '+var[i]],bins=100,range=(-1,1),histtype='step')
        	plt.hist(cutB['resolution '+var[i]],bins=100,range=(-1,1),histtype='step')
        	plt.hist(cutC['resolution '+var[i]],bins=100,range=(-1,1),histtype='step')
        	plt.xlabel(labels[i]+' Resolution')
        	plt.ylabel('Counts')
        	plt.legend(['No Cuts','logLikelihood > -52','logLikelihood > -50','logLikelihood > -48'])
        	plt.savefig('plots/Resolution_'+name+'_thad_'+var[i],bbox_inches='tight')
		
		print 'Saved Figure: Resolution_'+name+'_thad_'+var[i]


# ------------- PARTON EJETS ------------- #

# Create the main data frame
df_ejets = Create_DF('parton_ejets')

# Append all data to the main frame
for n in range(1,8):
        df_ejets = Append_Data(df_ejets,str(n),'parton_ejets')

# Edit(?) Data (remove unmatched, convert units, calculate resolution, etc.)
df_ejets = Edit_Data(df_ejets)

# Cut rows where likelihood <= some number (currently saved as a copy and not the original df)
cutA_ejets = df_ejets[df_ejets['likelihood']>-52]
cutB_ejets = df_ejets[df_ejets['likelihood']>-50]
cutC_ejets = df_ejets[df_ejets['likelihood']>-48]

# Create 2D plots
Create_2D_Plots(cutB_ejets,'parton_ejets')

# Create Resolution Plots
Plot_Res(df_ejets,cutA_ejets,cutB_ejets,cutC_ejets,'parton_ejets')



# ------------- PARTON MJETS ------------- #

# Create the main data frame
df_mjets = Create_DF('parton_mjets')

# Append all data to the main frame
for n in range(1,8):
        df_mjets = Append_Data(df_mjets,str(n),'parton_mjets')

# Edit(?) Data (remove unmatched, convert units, calculate resolution, etc.)
df_mjets = Edit_Data(df_mjets)

# Cut rows where likelihood <= some number (currently saved as a copy and not the original df)
cutA_mjets = df_mjets[df_mjets['likelihood']>-52]
cutB_mjets = df_mjets[df_mjets['likelihood']>-50]
cutC_mjets = df_mjets[df_mjets['likelihood']>-48]

# Create 2D plots
Create_2D_Plots(cutB_mjets,'parton_mjets')

# Create Resolution Plots
Plot_Res(df_mjets,cutA_mjets,cutB_mjets,cutC_mjets,'parton_mjets')






print("done :)")
