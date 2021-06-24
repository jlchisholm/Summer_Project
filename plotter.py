# Import useful packages
import uproot
import pandas
import awkward as ak
import numpy as np
import matplotlib
#matplotlib.use('Agg')  # need for not displaying plots when running batch jobs
from matplotlib import pyplot as plt
from matplotlib import colors

# List some useful strings, organized by particle
par = {'thad' : ['klfitter_bestPerm_topHad_','MC_thad_afterFSR_'],
       'tlep' : ['klfitter_bestPerm_topLep_','MC_tlep_afterFSR_'],
       'ttbar': ['klfitter_bestPerm_ttbar_','MC_ttbar_afterFSR_']}
var = {'thad' : ['pt','eta','y','phi','m','E','pout'], 
       'tlep': ['pt','eta','y','phi','m','E','pout'], 
       'ttbar' : ['pt','eta','y','phi','m','E','dphi','Ht','yboost']}
labels = {'thad' : ['$p_T^{t,had}$','$\eta^{t,had}$','$y^{t,had}$','$\phi^{t,had}$','$m^{t,had}$','$E^{t,had}$','$p_{out}^{t,had}$'], 
	  'tlep' : ['$p_T^{t,lep}$','$\eta^{t,lep}$','$y^{t,lep}$','$\phi^{t,lep}$','$m^{t,lep}$','$E^{t,lep}$','$p_{out}^{t,lep}$'],
	  'ttbar' : ['$p_T^{t\overline{t}}$','$\eta^{t\overline{t}}$','$y^{t\overline{t}}$','$\phi^{t\overline{t}}$','$m^{t\overline{t}}$','$E^{t\overline{t}}$','$\Delta\phi^{t\overline{t}}$','$H_T^{t\overline{t}}$','$y_{boost}^{t\overline{t}}$']}
units = {'thad' : [' [GeV]','','','',' [GeV]',' [GeV]',' [GeV]'], 
	 'tlep' : [' [GeV]','','','',' [GeV]',' [GeV]',' [GeV]'], 
	 'ttbar': [' [GeV]','','','',' [GeV]',' [GeV]','',' [GeV]','']} 


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
	for p in par:                     # For each particle
		for v in var[p]:          # For each variable of that particle
			df['reco '+p+' '+v] = ak.to_pandas(tree_reco0[par[p][0]+v])      # Append the reco values
        		df['truth '+p+' '+v] = ak.to_pandas(tree_truth0[par[p][1]+v])    # Append the true values
		
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
        for p in par:                     # For each particle
                for v in var[p]:          # For each variable of that particle
                        df_addon['reco '+p+' '+v] = ak.to_pandas(tree_reco_n[par[p][0]+v])      # Append the reco values
                        df_addon['truth '+p+' '+v] = ak.to_pandas(tree_truth_n[par[p][1]+v])    # Append the true values
        
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

	for p in par:                 # For each particle
		for v in var[p]:      # For each variable of that particle
        		
			# Convert from MeV to GeV for some truth values
			if v in ['pt','m','E','pout','Ht']:             
                                df['truth '+p+' '+v] = df['truth '+p+' '+v]/1000

			# Calculate the resolutions
			df['resolution '+p+' '+v] = df['truth '+p+' '+v]/df['reco '+p+' '+v]-1

	return df


# Creates and saves the 2D histogram normalized by rows
# Input: cutB of df, variable index, chosen bins/ticks, tick labels, data type, and particle
# Output: saves histogram as png
def Plot_2D(cutB,i,ticks,ticks_labels,name,particle):

	# Define useful constants
        n = len(ticks)
        ran = ticks[::n-1]

        # Create 2D array of truth vs reco variable (which can be plotted also)
        #plt.figure('2D')
	#fig,ax = plt.subplots()
        #ax.set_xticks(ticks)
        #ax.set_yticks(ticks)
        H, xedges, yedges, im = plt.hist2d(np.clip(cutB['reco '+particle+' '+var[particle][i]],ticks[0],ticks[-1]),np.clip(cutB['truth '+particle+' '+var[particle][i]],ticks[0],ticks[-1]),bins=ticks,range=[ran,ran])
        #plt.pcolormesh(xedges,yedges,cmap=plt.cm.Blues)
        #plt.xlabel('Reco '+labels[i]+units[i])
        #plt.ylabel('Truth '+labels[i]+units[i])

	# Normalize the rows out of 100 and round to integers
	norm = (np.rint((H/np.sum(H.T,axis=1))*100)).T.astype(int)

        # Plot truth vs reco pt with normalized rows
        plt.figure(name+' '+particle+' '+var[particle][i]+' Normalized 2D Plot')
        masked_norm = np.ma.masked_where(norm==0,norm)  # Needed to make the zero bins white
        plt.imshow(masked_norm,extent=[0,n-1,0,n-1],cmap=plt.cm.Blues,origin='lower')
        plt.xticks(np.arange(n),ticks_labels,fontsize=9,rotation=-25)
        plt.yticks(np.arange(n),ticks_labels,fontsize=9)
        plt.xlabel('Reco '+labels[particle][i]+units[particle][i])
        plt.ylabel('Truth '+labels[particle][i]+units[particle][i])
        plt.clim(0,100)
        plt.colorbar()

        # Label the content of each bin
        for j in range (n-1):
               for k in range(n-1):
                        if masked_norm.T[j,k] != 0:   # Don't label empty bins
                                plt.text(j+0.5,k+0.5,masked_norm.T[j,k],color="k",fontsize=6,weight="bold",ha="center",va="center")
        
	# Save the figure as a png
	plt.savefig('plots/'+name+'/Normalized_2D_'+name+'_'+particle+'_'+var[particle][i],bbox_inches='tight')
	print 'Saved Figure: Normalized_2D_'+name+'_'+particle+'_'+var[particle][i]
	
	plt.close()	


# Defines bins/ticks and implements Plot_2D for each variable
# Input: cutB of df, data type, and particle
def Create_2D_Plots(cutB, name, particle):
	
	# eta
	ticks_eta = np.arange(-6,7.5,1.5)
	ticks_labels_eta = map(str,ticks_eta)
	ticks_labels_eta[0] = '-6.28'
	ticks_labels_eta[-1] = '6.28'
	Plot_2D(cutB,1,ticks_eta,ticks_labels_eta,name,particle)
	
	# y
	ticks_y = np.arange(-2.5,3,0.5)
	ticks_labels_y = map(str,ticks_y)
	ticks_labels_y[0] = '-'+r'$\infty$'
	ticks_labels_y[-1] = r'$\infty$'
	Plot_2D(cutB,2,ticks_y,ticks_labels_y,name,particle)

	# phi
	ticks_phi = np.arange(-3,3.75,0.75)
	ticks_labels_phi = map(str,ticks_phi)
	ticks_labels_phi[0] = '-3.14'
	ticks_labels_phi[-1] = '3.14'
	Plot_2D(cutB,3,ticks_phi,ticks_labels_phi,name,particle)

	if particle == 'thad' or particle == 'tlep':

		# pt
        	#ticks = [0,50,100,160,225,300,360,475,1000]   # Choose ticks/binning (may be better choices, this was just in the paper)
        	#ticks_labels = ['0','50','100','160','225','300','360','475','1000']
        	ticks_pt = range(0,500,50)
        	ticks_labels_pt = map(str,ticks_pt)
        	ticks_labels_pt[-1] = r'$\infty$'
        	Plot_2D(cutB,0,ticks_pt,ticks_labels_pt,name,particle)
	
		# m
        	ticks_m = range(100,260,20)
        	ticks_labels_m = map(str,ticks_m)
        	ticks_labels_m[0] = '-'+r'$\infty$'
        	ticks_labels_m[-1] = r'$\infty$'
        	Plot_2D(cutB,4,ticks_m,ticks_labels_m,name,particle)
	
		# E
        	ticks_E = range(100,1000,100)
        	ticks_labels_E = map(str,ticks_E)
        	ticks_labels_E[-1] = r'$\infty$'
        	Plot_2D(cutB,5,ticks_E,ticks_labels_E,name,particle)

		# pout
		ticks_pout = range(-250,300,50)
		ticks_labels_pout = map(str,ticks_pout)
		ticks_labels_pout[0] = '-'+r'$\infty$'
		ticks_labels_pout[-1] = r'$\infty$'
		Plot_2D(cutB,6,ticks_pout,ticks_labels_pout,name,particle)

	elif particle == 'ttbar':
		
		# pt
        	ticks_pt = range(0,500,50)
        	ticks_labels_pt = map(str,ticks_pt)
        	ticks_labels_pt[-1] = r'$\infty$'
        	Plot_2D(cutB,0,ticks_pt,ticks_labels_pt,name,particle)

		# m
        	ticks_m = range(200,1100,100)
        	ticks_labels_m = map(str,ticks_m)
        	ticks_labels_m[-1] = r'$\infty$'
        	Plot_2D(cutB,4,ticks_m,ticks_labels_m,name,particle)
		
		# E
        	ticks_E = range(200,2200,200)
        	ticks_labels_E = map(str,ticks_E)
        	ticks_labels_E[-1] = r'$\infty$'
        	Plot_2D(cutB,5,ticks_E,ticks_labels_E,name,particle)

		# dphi
		ticks_dphi = np.arange(0,3.5,0.5)
		ticks_labels_dphi = map(str,ticks_dphi)
		ticks_labels_dphi[-1] = '3.14'
		Plot_2D(cutB,6,ticks_dphi,ticks_labels_dphi,name,particle)		

		# Ht
		ticks_Ht = range(0,1100,100)
		ticks_labels_Ht = map(str,ticks_Ht)
		ticks_labels_Ht[-1] = r'$\infty$'
		Plot_2D(cutB,7,ticks_Ht,ticks_labels_Ht,name,particle)

		# yboost
		ticks_yboost = np.arange(-3,3.5,0.75)
		ticks_labels_yboost = map(str,ticks_yboost)
		ticks_labels_yboost[0] = '-3.14'
		ticks_labels_yboost[-1] = '3.14'
		Plot_2D(cutB,8,ticks_yboost,ticks_labels_yboost,name,particle)


# Creates and saves the resolution plots for each variable
# Input: data frame, cutA of df, cutB, cutC, data type, and particle
def Plot_Res(df,cutA,cutB,cutC,name,particle):
	
	# Create each resolution plot
	for i in range(len(var[particle])):
        	plt.figure(name+' '+particle+' '+var[particle][i]+' Resolution')
        	plt.hist(df['resolution '+particle+' '+var[particle][i]],bins=100,range=(-1,1),histtype='step')
        	plt.hist(cutA['resolution '+particle+' '+var[particle][i]],bins=100,range=(-1,1),histtype='step')
        	plt.hist(cutB['resolution '+particle+' '+var[particle][i]],bins=100,range=(-1,1),histtype='step')
        	plt.hist(cutC['resolution '+particle+' '+var[particle][i]],bins=100,range=(-1,1),histtype='step')
        	plt.xlabel(labels[particle][i]+' Resolution')
        	plt.ylabel('Counts')
        	plt.legend(['No Cuts','logLikelihood > -52','logLikelihood > -50','logLikelihood > -48'])
        	plt.savefig('plots/'+name+'/Resolution_'+name+'_'+particle+'_'+var[particle][i],bbox_inches='tight')
		print 'Saved Figure: Resolution_'+name+'_'+particle+'_'+var[particle][i]
		
		plt.close()


# Creates the data frame and plots for a given set of data files
# Input: data type
def Run_Plotting(name):

	# Create the main data frame
	df = Create_DF(name)

	# Append all data to the main frame
	for n in range(1,8):
        	df = Append_Data(df,str(n),name)

	# Edit data (remove unmatched, convert units, calculate resolution, etc.)
	df = Edit_Data(df)

	# Cut rows where likelihood <= some number (currently saved as a copy and not the original df)
	cutA = df[df['likelihood']>-52]
	cutB = df[df['likelihood']>-50]
	cutC = df[df['likelihood']>-48]

	# Create 2D plots
	for p in par:
		Create_2D_Plots(cutB,name,p)

	# Create Resolution Plots
	for p in par:
		Plot_Res(df,cutA,cutB,cutC,name,p)


# ------------- MAIN CODE ------------- #

Run_Plotting('parton_ejets')
#Run_Plotting('parton_mjets')

print("done :)")
