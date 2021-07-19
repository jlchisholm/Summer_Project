# Import useful packages
import uproot
import pandas as pd
import awkward as ak
import numpy as np
import matplotlib
matplotlib.use('Agg')  # need for not displaying plots when running batch jobs
from matplotlib import pyplot as plt
from matplotlib import colors
from scipy.optimize import curve_fit


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

# Defines bins/ticks for each variable

# eta
ticks_eta = np.arange(-6,7.5,1.5)
ticks_labels_eta = map(str,ticks_eta)
ticks_labels_eta[0] = '-'+r'$\infty$'
ticks_labels_eta[-1] = r'$\infty$'

# y
ticks_y = np.arange(-2.5,3,0.5)
ticks_labels_y = map(str,ticks_y)
ticks_labels_y[0] = '-'+r'$\infty$'
ticks_labels_y[-1] = r'$\infty$'

# phi
ticks_phi = np.arange(-3,3.75,0.75)
ticks_labels_phi = map(str,ticks_phi)
ticks_labels_phi[0] = '-3.14'
ticks_labels_phi[-1] = '3.14'

# thad or tlep:

# pt
ticks_pt_t = range(0,500,50)
ticks_labels_pt_t = map(str,ticks_pt_t)
ticks_labels_pt_t[-1] = r'$\infty$'

# m
ticks_m_t = range(100,260,20)
ticks_labels_m_t = map(str,ticks_m_t)
ticks_labels_m_t[0] = '-'+r'$\infty$'
ticks_labels_m_t[-1] = r'$\infty$'

# E
ticks_E_t = range(100,1000,100)
ticks_labels_E_t = map(str,ticks_E_t)
ticks_labels_E_t[-1] = r'$\infty$'

# pout
ticks_pout_t = range(-250,300,50)
ticks_labels_pout_t = map(str,ticks_pout_t)
ticks_labels_pout_t[0] = '-'+r'$\infty$'
ticks_labels_pout_t[-1] = r'$\infty$'

# ttbar:

# pt
ticks_pt_ttbar = range(0,500,50)
ticks_labels_pt_ttbar = map(str,ticks_pt_ttbar)
ticks_labels_pt_ttbar[-1] = r'$\infty$'

# m
ticks_m_ttbar = range(200,1100,100)
ticks_labels_m_ttbar = map(str,ticks_m_ttbar)
ticks_labels_m_ttbar[-1] = r'$\infty$'

# E
ticks_E_ttbar = range(200,2200,200)
ticks_labels_E_ttbar = map(str,ticks_E_ttbar)
ticks_labels_E_ttbar[-1] = r'$\infty$'

# dphi
ticks_dphi_ttbar = np.arange(0,4,0.5)
ticks_labels_dphi_ttbar = map(str,ticks_dphi_ttbar)
ticks_labels_dphi_ttbar[-1] = r'$\infty$'

# Ht
ticks_Ht_ttbar = range(0,1100,100)
ticks_labels_Ht_ttbar = map(str,ticks_Ht_ttbar)
ticks_labels_Ht_ttbar[-1] = r'$\infty$'

# yboost
ticks_yboost_ttbar = np.arange(-3,3.5,0.75)
ticks_labels_yboost_ttbar = map(str,ticks_yboost_ttbar)
ticks_labels_yboost_ttbar[0] = '-'+r'$\infty$'
ticks_labels_yboost_ttbar[-1] = r'$\infty$'


ticks = {'thad' : [ticks_pt_t, ticks_eta, ticks_y, ticks_phi, ticks_m_t, ticks_E_t, ticks_pout_t],
         'tlep' : [ticks_pt_t, ticks_eta, ticks_y, ticks_phi, ticks_m_t, ticks_E_t, ticks_pout_t],
         'ttbar' : [ticks_pt_ttbar, ticks_eta, ticks_y, ticks_phi, ticks_m_ttbar, ticks_E_ttbar, ticks_dphi_ttbar, ticks_Ht_ttbar, ticks_yboost_ttbar]}
ticks_labels = {'thad' : [ticks_labels_pt_t, ticks_labels_eta, ticks_labels_y, ticks_labels_phi, ticks_labels_m_t, ticks_labels_E_t, ticks_labels_pout_t],
                'tlep' : [ticks_labels_pt_t, ticks_labels_eta, ticks_labels_y, ticks_labels_phi, ticks_labels_m_t, ticks_labels_E_t, ticks_labels_pout_t],
                'ttbar' : [ticks_labels_pt_ttbar, ticks_labels_eta, ticks_labels_y, ticks_labels_phi, ticks_labels_m_ttbar, ticks_labels_E_ttbar, ticks_labels_dphi_ttbar, ticks_labels_Ht_ttbar, ticks_labels_yboost_ttbar]}



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

	# Need to drop rows with inf or NaN values (a fix for now)
	df.replace([np.inf, -np.inf], np.nan, inplace=True)
	df.dropna(subset=['truth thad y'],inplace=True)

	for p in par:                 # For each particle
		for v in var[p]:      # For each variable of that particle
        		
			# Convert from MeV to GeV for some truth values
			if v in ['pt','m','E','pout','Ht']:             
                                df['truth '+p+' '+v] = df['truth '+p+' '+v]/1000

			# Calculate the resolutions
			if v in ['eta','y','phi']:    # Don't want to divide by 0, so we calculate residual resolution here
				df['resolution '+p+' '+v] = df['reco '+p+' '+v]-df['truth '+p+' '+v]
			else:
				df['resolution '+p+' '+v] = (df['reco '+p+' '+v]-df['truth '+p+' '+v])/df['truth '+p+' '+v]

	return df


# Creates and saves the 2D histogram normalized by rows
# Input: cutB of df, variable index, chosen bins/ticks, tick labels, data type, and particle
# Output: saves histogram as png
def Plot_2D(cutB,i,name,particle):

	# Define useful constants
        tk = ticks[particle][i]
	tkls = ticks_labels[particle][i]
	n = len(tk)
        ran = tk[::n-1]

        # Create 2D array of truth vs reco variable (which can be plotted also)
        H, xedges, yedges, im = plt.hist2d(np.clip(cutB['reco '+particle+' '+var[particle][i]],tk[0],tk[-1]),np.clip(cutB['truth '+particle+' '+var[particle][i]],tk[0],tk[-1]),bins=tk,range=[ran,ran])

	# Normalize the rows out of 100 and round to integers
	norm = (np.rint((H/np.sum(H.T,axis=1))*100)).T.astype(int)

        # Plot truth vs reco pt with normalized rows
        plt.figure(name+' '+particle+' '+var[particle][i]+' Normalized 2D Plot')
        masked_norm = np.ma.masked_where(norm==0,norm)  # Needed to make the zero bins white
        plt.imshow(masked_norm,extent=[0,n-1,0,n-1],cmap=plt.cm.Blues,origin='lower')
        plt.xticks(np.arange(n),tkls,fontsize=9,rotation=-25)
        plt.yticks(np.arange(n),tkls,fontsize=9)
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
	plt.savefig('plots/'+name+'/'+particle+'/Normalized_2D_'+name+'_'+particle+'_'+var[particle][i],bbox_inches='tight')
	print 'Saved Figure: Normalized_2D_'+name+'_'+particle+'_'+var[particle][i]
	
	plt.close()	


# Creates and saves the resolution plots for each variable
# Input: data frame, cutA of df, cutB, cutC, data type, and particle
def Plot_Res(df,cutA,cutB,cutC,name,particle):
	
	# Create each resolution plot
	for i in range(len(var[particle])):
        	
		# Define the absolute range for the resolution plot
		v = var[particle][i]
		ran = 5 if v=='pout' else 1		

		# Create 1D resolution plot
		plt.figure(name+' '+particle+' '+var[particle][i]+' Resolution')
        	plt.hist(df['resolution '+particle+' '+var[particle][i]],bins=100,range=(-ran,ran),histtype='step')
        	plt.hist(cutA['resolution '+particle+' '+var[particle][i]],bins=100,range=(-ran,ran),histtype='step')
        	plt.hist(cutB['resolution '+particle+' '+var[particle][i]],bins=100,range=(-ran,ran),histtype='step')
        	plt.hist(cutC['resolution '+particle+' '+var[particle][i]],bins=100,range=(-ran,ran),histtype='step')
        	plt.xlabel(labels[particle][i]+' Resolution')
        	plt.ylabel('Counts')
        	plt.legend(['No Cuts','logLikelihood > -52','logLikelihood > -50','logLikelihood > -48'])
        	
		# Save and close figure
		plt.savefig('plots/'+name+'/'+particle+'/Resolution_'+name+'_'+particle+'_'+var[particle][i],bbox_inches='tight')
		print 'Saved Figure: Resolution_'+name+'_'+particle+'_'+var[particle][i]
		plt.close()


# Creates and saves the resolution vs variable plots
# Input: df, cutA, cutB, cutC, data type, particle, resolution (y) variable, (x) variable, index of y variable, index of x variable, and [first bin edge, last bin edge, width of bins]
def Plot_Res_vs_Var(df,cutA,cutB,cutC,name,particle,y_var,x_var,iy,ix,bins):

	data = []   # Array to hold var vs fwhm values
	for bottom_edge in np.arange(bins[0],bins[1],bins[2]):	

		# Set some helpful variables
		top_edge = bottom_edge+bins[2]
		middle = bottom_edge+(bins[2]/2)

		# Look at resolution at a particular value of var
		cut_temp = df[df['truth '+particle+' '+x_var]>=bottom_edge]      # Should I fold in edges of first and last?
		cut_temp = cut_temp[cut_temp['truth '+particle+' '+x_var]<top_edge]
		cut_tempA = cutA[cutA['truth '+particle+' '+x_var]>=bottom_edge]      # Should I fold in edges of first and last?
                cut_tempA = cut_tempA[cut_tempA['truth '+particle+' '+x_var]<top_edge]
		cut_tempB = cutB[cutB['truth '+particle+' '+x_var]>=bottom_edge]      # Should I fold in edges of first and last?
                cut_tempB = cut_tempB[cut_tempB['truth '+particle+' '+x_var]<top_edge]
		cut_tempC = cutC[cutC['truth '+particle+' '+x_var]>=bottom_edge]      # Should I fold in edges of first and last?
                cut_tempC = cut_tempC[cut_tempC['truth '+particle+' '+x_var]<top_edge]

		# Get standard deviations with dataframe method
                sigma = cut_temp['resolution '+particle+' '+y_var].std()
		sigmaA = cut_tempA['resolution '+particle+' '+y_var].std()
		sigmaB = cut_tempB['resolution '+particle+' '+y_var].std()
		sigmaC = cut_tempC['resolution '+particle+' '+y_var].std()

		# Append the data from this particular var
		data.append([middle,bins[2]/2,sigma,sigmaA,sigmaB,sigmaC])		

	# Convert the array to a pandas data frame for easy reading	
	df_res = pd.DataFrame(data,columns=(x_var,'p/m','sigma','sigmaA','sigmaB','sigmaC'))

	# Create a scatter plot
	plt.figure(y_var+' Resolution vs '+x_var)
	plt.scatter(df_res[x_var], df_res['sigma'])
	plt.scatter(df_res[x_var], df_res['sigmaA'])
	plt.scatter(df_res[x_var], df_res['sigmaB'])
	plt.scatter(df_res[x_var], df_res['sigmaC'])
	plt.legend(['No Cuts','logLikelihood > -52','logLikelihood > -50','logLikelihood > -48'])
	plt.xlabel('Truth '+labels[particle][ix]+units[particle][ix])
	plt.ylabel(labels[particle][iy]+' Resolution Width')

	# Save and close figure
        plt.savefig('plots/'+name+'/'+particle+'/'+y_var+'_Resolution_vs_'+x_var+'_'+name+'_'+particle,bbox_inches='tight')
        print 'Saved Figure: '+y_var+'_Resolution_vs_'+x_var+'_'+name+'_'+particle
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
		for i in range(len(var[p])):
			Plot_2D(cutB,i,name,p)

	# Create resolution plots
	for p in par:
		Plot_Res(df,cutA,cutB,cutC,name,p)
		Plot_Res_vs_Var(df,cutA,cutB,cutC,name,p,'eta','eta',1,1,[-2.4,2.4,0.1])
		Plot_Res_vs_Var(df,cutA,cutB,cutC,name,p,'pt','pt',0,0,[90,510,20])	
		Plot_Res_vs_Var(df,cutA,cutB,cutC,name,p,'pt','eta',0,1,[-2.4,2.4,0.1])
		Plot_Res_vs_Var(df,cutA,cutB,cutC,name,p,'eta','pt',1,0,[90,510,20])



# ------------- MAIN CODE ------------- #

Run_Plotting('parton_ejets')
Run_Plotting('parton_mjets')

print("done :)")
