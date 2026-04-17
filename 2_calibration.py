#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 25 15:45:24 2025

@author: mrembauv
"""
#%% ===========================================================================
# Import libraries
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt
plt.style.use('style.mplstyle') # Comment if the Arial font is not installed
plt.close('all')
import pandas as pd
inch = 2.54 # inch to centimeter

# Numerical
import numpy as np
import pandas as pd
import netCDF4 as nc

# Statistics
import statsmodels.formula.api as smf
from sklearn.cross_decomposition import PLSRegression
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import root_mean_squared_error
from sklearn.metrics import r2_score,mean_squared_error
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram,set_link_color_palette

#%% ===========================================================================
# Define functions
# =============================================================================

#Define function to plot data according to the Southern Ocean zone
def zone_scatter(ax,x,y,zone,zone_list,zone_col,sta,add_labels):
	for i in range(len(zone_list)):
		ax.scatter(x[zone==zone_list[i]],y[zone==zone_list[i]],color=zone_col[i],s=7,label=zone_list[i],edgecolor='k',lw=0.5)
		if add_labels == True :
			for j in range(len(sta)):
				ax.text(x[j],y[j],sta[j])
				

def clean_diatom_species(data) :
	# merge Chaetoceros spores and vegetative stages
	data['Chaetoceros_sub._Hyalochaete'] = data['Chaetoceros_sub._Hyalochaete'] + data['Chaetoceros_resting_spore']
	data=data.drop(['Chaetoceros_resting_spore'],axis=1)
	
	# merge Eucampia spores and vegetative stages
	data['Eucampia_antarctica'] = data['Eucampia_antarctica'] + data['Eucampia_resting_spore']
	data=data.drop(['Eucampia_resting_spore'],axis=1)
	
	# Merge Thalassiosira resting spore with Thalassiosira_antarctica.
	data['Thalassiosira_antarctica'] = data['Thalassiosira_antarctica'] + data['Thalassiosira_resting_spore']
	data=data.drop(['Thalassiosira_resting_spore'],axis=1)
	
	# Merge the unresolved F. separanda/rhombica with Fragilariopsis spp.
	data['Fragilariopsis_spp.'] = data['Fragilariopsis_spp.'] + data['Fragilariopsis_separanda_rhombica']
	data=data.drop(['Fragilariopsis_separanda_rhombica'],axis=1)
	
	# Delete unidentified species (not informative)
	data=data.drop(['Centric_<_20_µm','Centric_>_20_µm','Pennate_unid.','Thalassiosira_spp.'],axis=1)
	
	# Return the cleaned data
	return(data)


#%% ===========================================================================
# Open the sediment core and sediment trap datasets
# =============================================================================

# Open the sediment core dataset
data_core = pd.read_excel('data/sediment_cores/sed_core_data.ods', engine='odf', sheet_name='data')
# Apply the species cleaning function on the sediment core data
data_core = clean_diatom_species(data_core)
diat_sel_core  = np.array(np.sum(data_core.iloc[:,6:],axis=0) != 0)
print('Number of diatom species/taxa groups in the sediment core data :')
print(np.sum(data_core.iloc[:,6:].sum(axis=0)!=0))

# Open sediment trap dataset
data_trap = pd.read_excel('data/sediment_traps/sed_trap_data.ods', engine='odf', sheet_name='cup_scale')
print('Number of diatom species/taxa groups in the sediment trap data :')
print(np.sum(data_trap.iloc[:,15:].sum(axis=0)!=0))

# Apply the species cleaning function on the sediment trap data
data_trap = clean_diatom_species(data_trap)

# Remove unreliable data (unresolved diatom taxonomy or too shallow sediment trap deployments)
data_trap = data_trap[(data_trap['reference']!='Salter_2012')] # Unresolved diatom taxonomy
data_trap = data_trap[(data_trap['reference']!='Blain_2021_2022')] # A3 neritic station (<300 m)
data_trap = data_trap[(data_trap['reference']!='Rembauville_2015')]# A3 neritic station (<300 m)
data_trap = data_trap[(data_trap['reference']!='Rembauville_2017')] # KERFIX neritic station (<300 m)

# Remove samples where POC or PIC or BSI fluxes = 0
data_trap=data_trap[(data_trap['poc_µmol_m2_d'] > 0) & (data_trap['pic_µmol_m2_d'] > 0) & (data_trap['bsi_µmol_m2_d'] > 0)]

# Calculate BSi:POC and PIC:POC ratios
data_trap.insert(11,'pic_poc',data_trap['pic_µmol_m2_d']/data_trap['poc_µmol_m2_d'])
data_trap.insert(12,'bsi_poc',data_trap['bsi_µmol_m2_d']/data_trap['poc_µmol_m2_d'])

# Average data by season
data_season = pd.DataFrame(columns = data_trap.columns[np.hstack([np.arange(7),np.arange(11,data_trap.shape[1])])])

sta_list = np.unique(data_trap['station'])

for i in np.arange(len(sta_list)):
	# select each station
	data_sub = data_trap[data_trap['station']==sta_list[i]]
	# identify the month corresponding to the midpoint of the sampling period
	midpoint = data_sub['start_date'] + (data_sub['stop_date']-data_sub['start_date'])/2
	month = np.array([i.month for i in midpoint])
	#identify the season
	spring = (month==9) | (month==10) | (month==11) # SON
	summer = (month==12) | (month==1) | (month==2) # DJF
	fall = (month==3) | (month==4) | (month==5) # MAM
	winter = (month==6) | (month==7) | (month==8)# JJA

	# add new line with seasonnaly averaged data to data frame
	new_line = pd.DataFrame(pd.concat([data_sub.iloc[0,:7],data_sub.iloc[spring,11:].mean(axis=0)])).T
	data_season = pd.concat([data_season,new_line],axis=0,ignore_index=True)
	new_line = pd.DataFrame(pd.concat([data_sub.iloc[0,:7],data_sub.iloc[summer,11:].mean(axis=0)])).T
	data_season = pd.concat([data_season,new_line],axis=0,ignore_index=True)
	new_line = pd.DataFrame(pd.concat([data_sub.iloc[0,:7],data_sub.iloc[fall,11:].mean(axis=0)])).T
	data_season = pd.concat([data_season,new_line],axis=0,ignore_index=True)
	new_line = pd.DataFrame(pd.concat([data_sub.iloc[0,:7],data_sub.iloc[winter,11:].mean(axis=0)])).T
	data_season = pd.concat([data_season,new_line],axis=0,ignore_index=True)

# Remove empty lines (no data available for a given season)
line_sel = np.sum(data_season.iloc[:,7:],axis=1) != 0
data_season=data_season[line_sel]

# Save the seasonnaly-averaged sediment trap data
data_season.to_csv('data/output/transfer_functions/data_trap_season.csv', sep=';')

#%% ===========================================================================
# Extract clean vectors at seasonal scale from the entire dataset
# =============================================================================

# Work on seasonnal data : clean diatom species and select diatom data
data_trap = data_season
diat_sel_trap = np.array(np.sum(data_trap.iloc[:,13:],axis=0) != 0)

# Total selector (merging diatom species present in both trap & core data)
diat_sel_total = (diat_sel_core) & (diat_sel_trap)
diat_trap = data_trap.iloc[:,13:].to_numpy().astype('float64')[:,diat_sel_total]
diat_core = data_core.iloc[:,6:].to_numpy().astype('float64')[:,diat_sel_total]

# Extract data_trap into numpy vectors for staistics
lat = data_trap['lat_degN'].to_numpy().astype('float64')
lon = data_trap['lon_degE'].to_numpy().astype('float64')
depth = data_trap['depth_m'].to_numpy().astype('float64')
poc = data_trap['poc_µmol_m2_d'].to_numpy().astype('float64')/1000 # Convert µmol to mmol/m2/d
pic_poc = data_trap['pic_poc'].to_numpy().astype('float64')
ref = data_trap['reference'].to_numpy()
sta = data_trap['station'].to_numpy()
zone = data_trap['zone'].to_numpy()

zone_col = ['red','lime','deepskyblue','white']
zone_list = ['SAZ','PFZ','POOZ','SIZ']

# Define chemical matrix
chem = np.vstack([poc,pic_poc]).T

# Clean species names (remove '_')
species = np.array(list(data_core.iloc[:,6:]))[diat_sel_total]
species_clean = []
for i in species:
	spe = i.replace('_',' ')
	species_clean.append(spe)
species_clean=np.array(species_clean)

# Abbreviate genus name
species_clean_abbr = np.empty(len(species_clean),dtype='<U42')
for i in range(len(species_clean_abbr)):
	words = species_clean[i].split(' ')
	if words[1] != 'spp.':
		cap = words[0][0]+'.'
		words[0] = cap
		species_clean_abbr[i] = ' '.join(words)
	else :
		species_clean_abbr[i] = species_clean[i]
	if species_clean[i] =='Thalassiothrix antarctica':
		species_clean_abbr[i] = 'Tx. antarctica'
	
	
#%% ===========================================================================
# Plot the relative abundance and seclect the most abundant species (Fig. 2)
# =============================================================================

criterion = 0.1 # 0.1 % after Esper and Gersonde (2014)

# Sort sediment core species
diat_core_mean = np.nanmean(diat_core,axis=0)
diat_core_max = np.nanmax(diat_core,axis=0)
diat_core_min = np.nanmin(diat_core,axis=0)
species_core_sorted  = [x for _, x in sorted(zip(diat_core_mean, species_clean_abbr))]
diat_core_max_sorted = [x for _, x in sorted(zip(diat_core_mean, diat_core_max))]
diat_core_min_sorted = [x for _, x in sorted(zip(diat_core_mean, diat_core_min))]
diat_core_mean_sorted = np.sort(diat_core_mean)

# Sort sediment trap species
diat_trap_mean = np.nanmean(diat_trap,axis=0)
diat_trap_max = np.nanmax(diat_trap,axis=0)
diat_trap_min = np.nanmin(diat_trap,axis=0)
species_trap_sorted  = [x for _, x in sorted(zip(diat_trap_mean, species_clean_abbr))]
diat_trap_max_sorted = [x for _, x in sorted(zip(diat_trap_mean, diat_trap_max))]
diat_trap_min_sorted = [x for _, x in sorted(zip(diat_trap_mean, diat_trap_min))]
diat_trap_mean_sorted = np.sort(diat_trap_mean)

# Identify indicative species
diat_sel_percent = (diat_core_mean>criterion) & (diat_trap_mean>criterion)
species_selected = species_clean_abbr[diat_sel_percent]

# Plot sediment core and sediment trap relative abundances
fig, ax = plt.subplots(1,2,figsize=(18/inch,18/inch))

# Plot sediment core species relative abundace
y = np.arange(len(species_core_sorted))
ax[0].barh(y,diat_core_mean_sorted,facecolor='w',ec='k',lw=0.5)
for i in range(len(species_core_sorted)):
	ax[0].plot([diat_core_mean_sorted[i],diat_core_max_sorted[i]],[y[i],y[i]],'-k',lw=0.5)
	if species_core_sorted[i] in np.array(species_selected):
		ax[0].barh(y[i],diat_core_mean_sorted[i],facecolor='lightgrey',ec='k',lw=0.5)
ax[0].set_yticklabels(species_core_sorted,fontstyle='italic',fontsize=6)

# Plot sediment trap species relative abundace
ax[1].barh(y,diat_trap_mean_sorted,facecolor='w',ec='k',lw=0.5)
for i in range(len(species_trap_sorted)):
	ax[1].plot([diat_trap_mean_sorted[i],diat_trap_max_sorted[i]],[y[i],y[i]],'-k',lw=0.5)
	if species_trap_sorted[i] in np.array(species_selected):
		ax[1].barh(y[i],diat_trap_mean_sorted[i],facecolor='lightgrey',ec='k',lw=0.5)
ax[1].set_yticklabels(species_trap_sorted,fontstyle='italic',fontsize=6)

# Cosmetics
lab = ['a','b']
for i in [0,1]:
	ax[i].set_xlim([1e-4,1e2])
	ax[i].set_xscale('log')
	ax[i].set_yticks(y)
	ax[i].set_ylim([-0.75,len(species_clean)-0.25])
	ax[i].set_xlabel('Mean relative abundance (%)')
	ax[i].grid(axis='x',lw=0.5)
	ax[i].set_axisbelow(True)
	ax[i].annotate(lab[i],(-0.7,1.02),xycoords='axes fraction',fontweight='bold')

ax[0].set_title('Sediment cores',fontsize=8,fontweight='bold')
ax[1].set_title('Sediment traps',fontsize=8,fontweight='bold')

plt.subplots_adjust(left=0.2,right=0.97,wspace=1,top=0.95,bottom=0.08)

plt.savefig('fig/fig_03_diatom_species.png',dpi=300)
plt.savefig('fig/fig_03_diatom_species.svg')
plt.savefig('fig/fig_03_diatom_species.pdf')


#%% Select species > criterion % relative abundace
species_clean = species_clean[diat_sel_percent]
species_clean_abbr = species_clean_abbr[diat_sel_percent]
diat_core_full = diat_core
diat_trap_full = diat_trap
diat_core = diat_core[:,diat_sel_percent]
diat_trap = diat_trap[:,diat_sel_percent]

# Print abundaces ranges of informative species
min_core = str(min(np.sum(diat_core,axis=1)))
max_core = str(max(np.sum(diat_core,axis=1)))
min_trap = str(min(np.sum(diat_trap,axis=1)))
max_trap = str(max(np.sum(diat_trap,axis=1)))

print('The sum of informative species account for',min_core,'to',max_core,'% of the sediment core data')
print('The sum of informative species account for',min_trap,'to',max_trap,'% of the sediment trap data')

# normalise indicative species matrices to 100%
diat_core = (diat_core.T/np.sum(diat_core,axis=1)*100).T
diat_trap = (diat_trap.T/np.sum(diat_trap,axis=1)*100).T

#%% ===========================================================================
#  Plot sediment trap fluxes and diatom community composition (Fig. 3)
# =============================================================================

# Sort stations according to hydrological zones
sta_list = ['47S','MS2','54S','P3','P2','61S','MS3','MS4','WSC','MS5','PZB-1']
sta_lab = np.array(['47S', '', '', '', 'MS2', '', '', '', '54S', '', '', '','P3', '', '', '',
      'P2', '', '', '','61S', '', '', '', 'MS3', '', '', '','MS4', '', '', '','WSC', '', '', '',
       'MS5', '', '','','PZB-1', '', '', ''])

fig,ax = plt.subplots(1,3,figsize=(18/inch,18/inch))

# Plot the diatom species and the POC and PIC:POC ratio
# Define key species for the sed trap plot
key_species = np.hstack([diat_trap_full[:,species=='Azpeitia_tabularis'],
diat_trap_full[:,species=='Nitzschia_bicapitata']+diat_trap_full[:,species=='Nitzschia_sicula'],
diat_trap_full[:,species=='Pseudo-nitzschia_heimii']+diat_trap_full[:,species=='Pseudo-nitzschia_lineola']+diat_trap_full[:,species=='Pseudo-nitzschia_turgiduloides'],
diat_trap_full[:,species=='Thalassionema_nitzschioides']+diat_trap_full[:,species=='Thalassionema_nitzschioides_var_capitulata']+diat_trap_full[:,species=='Thalassionema_nitzschioides_var_lanceolata']+diat_trap_full[:,species=='Thalassionema_nitzschioides_var_parva'],
diat_trap_full[:,species=='Fragilariopsis_kerguelensis'],
						diat_trap_full[:,species=='Fragilariopsis_curta'],
						diat_trap_full[:,species=='Fragilariopsis_cylindrus'],
						diat_trap_full[:,species=='Fragilariopsis_spp.']+diat_trap_full[:,species=='Fragilariopsis_doliolus']+diat_trap_full[:,species=='Fragilariopsis_obliquecostata']+diat_trap_full[:,species=='Fragilariopsis_pseudonana']+diat_trap_full[:,species=='Fragilariopsis_rhombica']+diat_trap_full[:,species=='Fragilariopsis_ritscheri']+diat_trap_full[:,species=='Fragilariopsis_separanda']+diat_trap_full[:,species=='Fragilariopsis_vanheurckii'],
						diat_trap_full[:,species=='Thalassiosira_gracilis']+diat_trap_full[:,species=='Thalassiosira_gracilis_var_expecta'],
						diat_trap_full[:,species=='Thalassiosira_lentiginosa']+diat_trap_full[:,species=='Thalassiosira_antarctica']+diat_trap_full[:,species=='Thalassiosira_eccentrica']+diat_trap_full[:,species=='Thalassiosira_gravida']+diat_trap_full[:,species=='Thalassiosira_lineata']+diat_trap_full[:,species=='Thalassiosira_maculata']+diat_trap_full[:,species=='Thalassiosira_oestrupii']+diat_trap_full[:,species=='Thalassiosira_oliverana']+diat_trap_full[:,species=='Thalassiosira_tumida']+diat_trap_full[:,species=='Thalassiosira_spp.'],
						diat_trap_full[:,species=='Corethron_pennatum'],
						diat_trap_full[:,species=='Chaetoceros_sub._Hyalochaete']])

colors = ['gold','orange','chocolate','saddlebrown','blue','deepskyblue','lightskyblue','lightblue','lightgreen','yellowgreen','green','darkgreen']

lab_legend = ['Azpeitia tabularis', 'Nitzschia spp.','Pseudo-nitzschia spp.', 'Thalassionema nitzschioides','Fragilariopsis kerguelensis','Fragilariopsis curta', 'Fragilariopsis cylindrus','Fragilariopsis spp.','Thalassiosira gracilis','Thalassiosira spp.','Corethron pennatum','Chaetoceros sub. Hyalochaete']

# Stacked relative abundance plot
start = 0
for i in sta_list :
	# select data
	sel = (sta==i)
	if i =='MS5': # There is no winter data for MS5
		y = start + np.arange(sum(sel)+1)
		poc_var = np.hstack([poc[sel],np.nan])
		pic_poc_var = np.hstack([pic_poc[sel],np.nan])
	else :
		y = start + np.arange(sum(sel))
		poc_var = poc[sel]
		pic_poc_var=pic_poc[sel]
	
	# Plot diatom relative abundance
	data = key_species[sel,:]
	for j in range(data.shape[0]): # for each seasonal average
		left = 0
		for k in range(data.shape[1]): # for each species
			ax[0].barh(y[j],data[j,k],left = left,color=colors[k],height=0.75,label=lab_legend[k])
			left += data[j,k]
		ax[0].barh(y[j],100,left = left,color='lightgrey',height=0.75)
			
	# Plot POC fluc and PIC:POC ratio
	ax[1].plot(poc_var,y,'-ok',markersize=2,lw=0.5)
	ax[2].plot(pic_poc_var,y,'-ok',markersize=2,lw=0.5)
	start = max(y)+1

# Cosmetics
ticks = [-0.5,3.5,11.5,27.5,43.5]
for i in [0,1,2] :
	ax[i].set_ylim([max(ticks),min(ticks)])
	for j in range(len(ticks)):
		ax[i].axhline(ticks[j],color='k',linewidth=0.5)
	if i>0 :
		ax[i].set_yticks([])

ax[0].set_yticks(range(len(sta_lab)))
ax[0].set_yticklabels(sta_lab)
ax[0].set_xlim([0,100])
ax[0].set_xlabel('Relative abundance (%)')

ax[1].set_xlim([0,2.5])
ax[1].set_xticks([0,.5,1,1.5,2,2.5])
ax[1].set_xlabel('POC flux\n$(mmol\ m^{-2}\ d^{-1})$')

ax[2].set_xlim([0,2])
ax[2].set_xlabel('PIC:POC\n (mol:mol)')

plt.subplots_adjust(left=0.08,right=0.92,top=0.85,bottom=0.1)

# Legend with diatom species names
handles, labels = ax[0].get_legend_handles_labels()
handles = handles[:len(lab_legend)]
legend = fig.legend(handles = handles,fontsize=8,ncol=3,bbox_to_anchor=(0.91, 0.985),frameon=False)
for text in legend.get_texts():
    text.set_fontstyle("italic")

# Add SO zones
axtwin = ax[2].twinx()
axtwin.set_ylim(ax[2].get_ylim())
axtwin.set_yticks(ticks)
axtwin.tick_params('y',length=5)
axtwin.set_yticklabels([])
y_pos = [1.75,7.75,20,35]
for i in range(len(zone_list)):
	axtwin.text(2.1,y_pos[i],zone_list[i],va='center')
	
# Add labels to pannels
lab = ['a','b','c']
for i in [0,1,2]:
	ax[i].annotate(lab[i],(0,1.01),xycoords='axes fraction',fontweight='bold')

'''
plt.savefig('fig/fig_02_biogeochem_fluxes.png',dpi=300)
plt.savefig('fig/fig_02_biogeochem_fluxes.pdf')
'''


#%% ===========================================================================
# Factor analysis with varimax rotation and multiple linear regresion (MLR)
# =============================================================================

from factor_analyzer import FactorAnalyzer

# Perform factor analysis (keeping 4 factor according to the scree plot)
X = np.log(diat_core+1)
X = StandardScaler().fit(X).transform(X)
fa = FactorAnalyzer(n_factors=4,rotation='varimax').fit(X)
X = np.log(diat_trap+1)
X = StandardScaler().fit(X).transform(X)
X_transformed = fa.transform(X)

# Scree plot for factor selection (Supplementary figure 1)
loadings = fa.loadings_
eig = fa.get_eigenvalues()
fig, ax = plt.subplots(1,1,figsize=(9/inch,7.5/inch))
ax.bar(np.arange(len(eig[1][:10]))+1,eig[1][:10],fc='lightgrey',ec='k',width=0.7,lw=0.5)
ax.axhline(1,color='k',linestyle='--',lw=0.5,zorder=-10)
ax.set_ylim([0,6])
ax.set_xlim([0.4,10.6])
ax.set_xticks(np.arange(10)+1)
ax.set_ylabel('Eigenvalue')
ax.set_xlabel('Factor')
plt.tight_layout()
plt.savefig('fig/fig_S1_scree_plot.png', dpi=300)
plt.savefig('fig/fig_S1_scree_plot.pdf')

# Save loadings for Table 3
FA = FactorAnalyzer(n_factors=4,rotation='varimax').fit(X)
data = np.vstack([species_clean,np.round(FA.loadings_.T,2)]).T
np.savetxt('data/output/transfer_functions/FA_loadings.txt', data, delimiter=";", fmt="%s") 

# Select final MLR models explicitely (Table 4)
group_names = ['g'+str(i+1) for i in range(X_transformed.shape[1])]
df = pd.DataFrame(data=np.hstack([chem,X_transformed]),columns=np.hstack([['poc','pic_poc'],group_names]))
chem_pred_mlr = np.empty(chem.shape)
# MLR for POC
mlr1 = smf.ols(formula='poc ~ g1+g2+g3+g4',data=df).fit()
print(mlr1.summary())
chem_pred_mlr[:,0] = mlr1.predict()
# MLR for PIC
mlr2 = smf.ols(formula='pic_poc ~ g1+g2+g3+g4',data=df).fit()
print(mlr2.summary())
chem_pred_mlr[:,1] = mlr2.predict()



#%% ===========================================================================
# Partial least square regression (PLSR)
# =============================================================================

plsr = PLSRegression(n_components=2) # keep n components accoring to the broken stick model (see below)
plsr.fit(X, chem)
chem_pred_plsr = plsr.predict(X)

#Check PLS prediction capability as function of number of components
comp = np.arange(1,X.shape[1])
r2_pls = []
for i in comp :
	plsr_test = PLSRegression(n_components=i)
	plsr_test.fit(X, chem)
	r2_pls.append(plsr_test.score(X, chem))
r2_pls = np.array(r2_pls)
r2_pls = r2_pls/max(r2_pls)*100
var_comp = np.hstack([r2_pls[0],np.diff(r2_pls)])

# Broken stick model (Jackson 1993)
broken_stick = np.flip(np.cumsum(1/np.arange(X.shape[1]-1,0,-1)))
broken_stick = broken_stick/np.sum(broken_stick)*100

# Compare variance with broken stick model to select the number of components (Fig. S2)
fig,ax = plt.subplots(figsize=(9/inch,7.5/inch))
# Plot PLSR variance
bar = ax.bar(comp,np.sort(var_comp)[::-1],fc='lightgrey',ec='k',label='PLSR',width=0.7,lw=0.5)
brok, = ax.plot(comp,broken_stick,'-ok',lw=1,ms=3,label='Broken stick')

# Cosmetics
ax.set_xlim([0.4,10.6])
ax.set_xticks(np.arange(10)+1)
ax.set_ylim([0,60])
ax.set_ylabel('Variance (%)')
ax.set_xlabel('Component')
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1],labels[::-1],frameon=False)
plt.tight_layout()

# Save PLSR components figure
plt.savefig('fig/fig_S2_PLSR_components.png',dpi=300)
plt.savefig('fig/fig_S2_PLSR_components.pdf')

# PLSR beta coefficients figure
plsr.fit(X, StandardScaler().fit(chem).transform(chem))
beta = plsr.coef_
# save PLSR beta coefficients
data_export = np.vstack([species_clean_abbr,np.round(beta,3)]).T.astype(str)
np.savetxt('data/output/transfer_functions/PLSR_beta.txt',data_export,delimiter='\t',
		   fmt='%s\t%s\t%s')
fig,ax = plt.subplots(1,2,figsize=(18/inch,10/inch))
plt.subplots_adjust(left=0.25,right=0.97,top=0.9,bottom=0.15,wspace=1.2)
species_poc  = [x for _, x in sorted(zip(beta[0,:], species_clean_abbr))]
beta_poc_sorted = np.sort(beta[0,:])
species_pic_poc  = [x for _, x in sorted(zip(beta[1,:], species_clean_abbr))]
beta_pic_poc_sorted = np.sort(beta[1,:])
var = [beta_poc_sorted,beta_pic_poc_sorted]
var_lab =['POC','PIC:POC']
species_lab = [species_poc,species_pic_poc]
lab = ['a','b']

x = np.arange(len(species_clean_abbr))
for i in [0,1]:
	ax[i].barh(x,var[i],fc='lightgrey',ec='k',lw=0.5,height=0.75)
	ax[i].axvline(0,color='k',lw=0.5)
	ax[i].set_yticks(x)
	ax[i].set_yticklabels(species_lab[i],fontstyle='italic')
	ax[i].set_xlim([-.2,.2])
	ax[i].set_ylim([-0.5,21.5])
	ax[i].set_xlabel('$\\beta$ coefficient')
	ax[i].set_title(var_lab[i],fontsize=8,fontweight='bold')
	ax[i].annotate(lab[i],(-0.9,0.97),xycoords='axes fraction',fontweight='bold',fontsize='8')
	ax[i].grid(axis='x')
	ax[i].set_axisbelow(True)

plt.savefig('fig/fig_S4_beta_PLSR.png',dpi=300)
plt.savefig('fig/fig_S4_beta_PLSR.pdf')


#%% ===========================================================================
#  Gradient boosting regression (GBR)
# =============================================================================

chem_pred_gbr = np.empty(chem.shape)

# See below for the selection tree depth (2) and number of estimation (30)
gbr1 = GradientBoostingRegressor(max_depth=2, n_estimators=30)
gbr2 = GradientBoostingRegressor(max_depth=2, n_estimators=30)

# GBR for POC
gbr1.fit(X, chem[:,0])
chem_pred_gbr[:,0] = gbr1.predict(X)
gbr1_features = gbr1.feature_importances_*100
# GBR for PIC
gbr2.fit(X, chem[:,1])
chem_pred_gbr[:,1] = gbr2.predict(X)
gbr2_features = gbr2.feature_importances_*100
# Save features importance
data = np.vstack([species_clean,np.round(np.vstack([gbr1_features,gbr2_features]),2)]).T
np.savetxt('data/output/transfer_functions/GBR_features_importance.txt', data, delimiter=";", fmt="%s") 

# Sort species names according to their contribution to the GBR
species_poc  = [x for _, x in sorted(zip(gbr1_features, species_clean_abbr))]
gbr1_features_sorted = np.sort(gbr1_features)
species_pic_poc  = [x for _, x in sorted(zip(gbr2_features, species_clean_abbr))]
gbr2_features_sorted = np.sort(gbr2_features)

# Plot feature importances for each GBR model
fig,ax = plt.subplots(1,2,figsize=(18/inch,10/inch))
plt.subplots_adjust(left=0.25,right=0.97,top=0.9,bottom=0.15,wspace=1.2)
x = np.arange(X.shape[1])+1
ax[0].barh(x,gbr1_features_sorted,fc='lightgrey',ec='k',lw=0.5,height=0.75)
ax[1].barh(x,gbr2_features_sorted,fc='lightgrey',ec='k',lw=0.5,height=0.75)
for i in [0,1] :
	ax[i].set_yticks(x)
	ax[i].set_xscale('log')
	ax[i].set_xlim([1e-2,1e2])
	ax[i].set_ylim([0.5,22.5])
	ax[i].grid(axis='x',lw=0.5,zorder=-10)
	ax[i].set_axisbelow(True)
	ax[i].set_xlabel('Species importance (%)')
ax[0].set_yticklabels(species_poc,fontstyle='italic',fontsize=8)
ax[1].set_yticklabels(species_pic_poc,fontstyle='italic',fontsize=8)
ax[0].set_title('POC',fontweight='bold',fontsize=8)
ax[1].set_title('PIC:POC',fontweight='bold',fontsize=8)
lab = ['a','b']
for i in [0,1]:
	ax[i].annotate(lab[i],(-0.9,0.97),xycoords='axes fraction',fontweight='bold',fontsize='8')

# Save figure
plt.savefig('fig/fig_S5_GBR_features_importance.png',dpi=300)
plt.savefig('fig/fig_S5_GBR_features_importance.pdf')


#%% ===========================================================================
# Transfer function calibration plot (Fig. 4)
# =============================================================================

# Create figure
fig, ax = plt.subplots(3,2,figsize = (9/inch,12/inch))
x_labels = ['Observed $(mmol\ m^{-2}\ d^{-1})$',
		  'Observed (mol:mol)']
cols = ['POC flux','PIC:POC ratio']
rows = ['MLR','PLSR','GBR']
lab = ['a','b','c','d','e','f']
c=0

# Y-axis ranges
poc_min = -0.2
poc_max = 2.2
pic_poc_min = -0.2
pic_poc_max = 2.2

for i in [0,1]:
	# Plot MLR results
	zone_scatter(ax[0,i],chem[:,i],chem_pred_mlr[:,i],zone,zone_list,zone_col,sta,add_labels=False)
	r2 = np.round(r2_score(chem[:,i],chem_pred_mlr[:,i]),2)
	ax[0,i].annotate('$R^{2}=$'+str(r2),(0.5,0.84),xycoords='axes fraction',ha='center')
	# Plot PLSR results
	zone_scatter(ax[1,i],chem[:,i],chem_pred_plsr[:,i],zone,zone_list,zone_col,sta,add_labels=False)
	r2 = np.round(r2_score(chem[:,i],chem_pred_plsr[:,i]),2)
	ax[1,i].annotate('$R^{2}=$'+str(r2),(0.5,0.84),xycoords='axes fraction',ha='center')
	# PLot GBR results
	zone_scatter(ax[2,i],chem[:,i],chem_pred_gbr[:,i],zone,zone_list,zone_col,sta,add_labels=False)
	r2 = np.round(r2_score(chem[:,i],chem_pred_gbr[:,i]),2)
	ax[2,i].annotate('$R^{2}=$'+str(r2),(0.5,0.84),xycoords='axes fraction',ha='center')
		
	# Axes labels
	ax[2,i].set_xlabel(x_labels[i])
	ax[0,i].set_title(cols[i],fontsize=8,fontweight='bold')

	# Add 1:1 relationship
	for j in [0,1,2]:
		# Add 1:1 relationship
		ax[j,i].set_xlim([-0.2,2.2])
		ax[j,i].set_ylim([-0.2,2.2])
		ax[j,i].plot(ax[j,i].get_xlim(),ax[j,i].get_xlim(),'-k',lw=0.5,zorder=-10)

# Add labels
for i in [0,1,2]:
	ax[i,0].annotate(rows[i],(-0.42,0.5),xycoords='axes fraction',rotation=90,va='center',fontweight='bold')
	for j in [0,1]:
		ax[i,j].annotate(lab[c],(0.07,0.85),xycoords='axes fraction',fontweight='bold')
		c+= 1
		
ax[1,0].set_ylabel('Estimated')

plt.subplots_adjust(left=0.17,right=0.95,wspace=0.5,hspace=0.4,top=0.93)

# Save figure
plt.savefig('fig/fig_04_proxy_calibration.png',dpi=300)
plt.savefig('fig/fig_04_proxy_calibration.pdf')


#%% =============================================================================
#  Bootstrap to choose GBR tree depth and number of estimators 
#	-- WARNING : takes ~20 min with standard CPU --
# =============================================================================
'''
estim = 50 # Start with 50 estimation as a basis
train_fraction = 2/3

def gbr_iterator(X,var,estim,depth,train_fraction):
	# Apply GBR
	reg = GradientBoostingRegressor(n_estimators=estim,max_depth=depth)
	X_train, X_test, y_train, y_test = train_test_split(X, var, train_size=train_fraction)
	y_pred = reg.fit(X_train, y_train).predict(X_test)
	# Calculate RMSE
	return root_mean_squared_error(y_test, y_pred)

depth = [1,2,3,4,5]
n = 10000 # Number of bootstrapping iterations
rmse_poc_depth = np.empty((n,len(depth)))
rmse_pic_poc_depth = np.empty((n,len(depth)))

for i in range(n):
	#Print iteration
	if i % 100 == 0:
		print(i)
	for j in range(len(depth)):
		rmse_poc_depth[i,j] = gbr_iterator(X,chem[:,0],estim,depth[j],train_fraction)
		rmse_pic_poc_depth[i,j] = gbr_iterator(X,chem[:,1],estim,depth[j],train_fraction)

# Save RMSEs as output		
np.savetxt('data/output/transfer_functions/rmse_poc_depth.txt', rmse_poc_depth, delimiter='\t',fmt='%1.4e')
np.savetxt('data/output/transfer_functions/rmse_pic_poc_depth.txt', rmse_pic_poc_depth, delimiter='\t',fmt='%1.4e')

# Bootstrap to choose the number of GBR estimators

estim = 50
train_fraction = 2/3

def gbr_iterator(X,var,depth,estim,train_fraction,random_state):
	# Apply GBR
	reg = GradientBoostingRegressor(n_estimators=estim,max_depth=depth)
	X_train, X_test, y_train, y_test = train_test_split(X, var, train_size=train_fraction, random_state=random_state)
	reg.fit(X_train, y_train)
	# Calculate RMSE
	rmse = np.empty(estim)
	for k, y_pred in enumerate(reg.staged_predict(X_test)):
	    rmse[k] = mean_squared_error(y_test, y_pred)
	# Return the RMSE
	return rmse
	
n = 10000 # Number of bootstrapping iterations
rmse_poc_estim = np.empty((n,estim))
rmse_pic_poc_estim = np.empty((n,estim))

for i in range(n):
	#Print iteration
	if i % 100 == 0:
		print(i)
	rmse_poc_estim[i,:] = gbr_iterator(X,chem[:,0],2,estim,train_fraction,i)
	rmse_pic_poc_estim[i,:] = gbr_iterator(X,chem[:,1],2,estim,train_fraction,i)

# Save RMSEs as output
#np.savetxt('data/output/transfer_functions/rmse_poc_estim.txt', rmse_poc_estim, delimiter='\t',fmt='%1.4e')
#np.savetxt('data/output/transfer_functions/rmse_pic_poc_estim.txt', rmse_pic_poc_estim, delimiter='\t',fmt='%1.4e')
'''

# Plot RMSE for different tree depth and estimators numbers

# Open RMSEs saved as output
rmse_poc_depth = np.loadtxt('data/output/transfer_functions/rmse_poc_depth.txt',delimiter='\t')
rmse_pic_poc_depth = np.loadtxt('data/output/transfer_functions/rmse_pic_poc_depth.txt',delimiter='\t')
rmse_poc_estim = np.loadtxt('data/output/transfer_functions/rmse_poc_estim.txt',delimiter='\t')
rmse_pic_poc_estim = np.loadtxt('data/output/transfer_functions/rmse_pic_poc_estim.txt',delimiter='\t')

def plot_rmse_depth(ax,var,style,lw):
	mean_rmse = np.nanmean(var,axis=0)
	mean_rmse_scaled = (mean_rmse-np.mean(mean_rmse))/np.std(mean_rmse)
	x = np.arange(len(mean_rmse))+1
	ax.plot(x,mean_rmse_scaled,style,markersize=3,lw=lw)

def plot_rmse_estim(ax,var,style,handle):
	mean_rmse = np.nanmean(var,axis=0)
	mean_rmse_scaled = (mean_rmse-np.mean(mean_rmse))/np.std(mean_rmse)
	x = np.arange(len(mean_rmse))+1
	ax.plot(x,mean_rmse_scaled,style,lw=1,label=handle)

fig,ax = plt.subplots(1,2,figsize=(12/inch,5.5/inch))
plot_rmse_depth(ax[0],rmse_poc_depth,'-ok',1)
plot_rmse_depth(ax[0],rmse_pic_poc_depth,'--ok',1)
poc = plot_rmse_estim(ax[1],rmse_poc_estim,'-k','POC flux')
pic = plot_rmse_estim(ax[1],rmse_pic_poc_estim,'--k','PIC:POC')
ax[1].legend(frameon=False,handlelength=1.8)

# Cosmetics
ax[0].set_xlim([0.5,5.5])
ax[0].set_xticks([1,2,3,4,5])
ax[0].set_yticks([-2,-1,0,1,2])

ax[1].set_xlim([0,50])
ax[1].set_xticks([0,10,20,30,40,50])
ax[1].set_yticks([-1,0,1,2,3,4,5])

ax[0].set_xlabel("Tree depth")
ax[0].set_ylabel("Standardized RMSE")
ax[1].set_xlabel("Estimators")
ax[1].set_ylabel("Standardized RMSE")
plt.tight_layout(w_pad=5)

lab = ['a','b']
for i in [0,1] :
	ax[i].annotate(lab[i],(0.05,0.9),xycoords='axes fraction',fontweight='bold')
# Save figure
plt.savefig('fig/fig_S3_gbr_depth_estimators.png',dpi=300)
plt.savefig('fig/fig_S3_gbr_depth_estimators.pdf')



#%% ===========================================================================
#  Evaluation of Transfer function precision using bootstraping (Figure 5)
#	-- WARNING : takes ~10 min with standard CPU --
# =============================================================================
'''
range_poc = np.max(chem[:,0])-np.min(chem[:,0])
range_pic_poc = np.max(chem[:,1])-np.min(chem[:,1])
range_chem = [range_poc,range_pic_poc]

n = 10000
train_fraction = 2/3

# MLR evaluation
print('MLR evaluation (please wait) \n')
X = np.log(diat_trap+1)
X = StandardScaler().fit(X).transform(X)
X = fa.transform(X)
rmse_mlr = np.empty((n,2))
for i in range(len(rmse_mlr)):
	#Print iteration
	if i % 100 == 0:
		print(i)
	# Split dataset into train/test subsets
	X_train, X_test, Y_train, Y_test = train_test_split(X, chem, train_size=train_fraction)
	df_train = pd.DataFrame(data=np.hstack([Y_train,X_train]),columns=['poc','pic_poc']+group_names)
	df_test = pd.DataFrame(X_test,columns=group_names)
	Y_pred = np.empty(Y_test.shape)
	# Predict with MLR
	mlr1 = smf.ols(formula='poc ~ g1 + g2 + g3 + g4',data=df_train).fit()
	Y_pred[:,0] = mlr1.predict(df_test)
	mlr2 = smf.ols(formula='pic_poc ~ g1 + g2 + g3+ g4',data=df_train).fit()
	Y_pred[:,1] = mlr2.predict(df_test)
	# Store RMSE for each variable
	for j in [0,1]:
		RMSE = root_mean_squared_error(Y_test[:,j], Y_pred[:,j])
		rmse_mlr[i,j] = RMSE/range_chem[j]*100


# PLSR evaluation
print('PLSR evaluation (please wait)\n')
X = np.log(diat_trap+1)
X = StandardScaler().fit(X).transform(X)
rmse_plsr = np.empty((n,2))
for i in range(len(rmse_plsr)):
	# Split dataset into train/test subsets
	X_train, X_test, Y_train, Y_test = train_test_split(X, chem, train_size=train_fraction)
	# Predict with PLSR
	plsr.fit(X_train,Y_train)
	Y_pred = plsr.predict(X_test)
	# Store RMSE for each variable
	for j in [0,1]:
		RMSE = root_mean_squared_error(Y_test[:,j], Y_pred[:,j])
		rmse_plsr[i,j] = RMSE/range_chem[j]*100

print('GBR evaluation (please wait) \n')
# GBR evaluation
rmse_gbr = np.empty((n,2))
for i in range(len(rmse_gbr)):
	#Print iteration
	if i % 100 == 0:
		print(i)
	# Split dataset into train/test subsets
	X_train, X_test, Y_train, Y_test = train_test_split(X, chem, train_size=train_fraction)
	Y_pred = np.empty(Y_test.shape)
	# Predict with GBR
	gbr1.fit(X_train, Y_train[:,0])
	Y_pred[:,0] = gbr1.predict(X_test)
	gbr2.fit(X_train, Y_train[:,1])
	Y_pred[:,1] = gbr2.predict(X_test)
	# Store RMSE for each variable (POC PIC:POC)
	for j in [0,1]:
		RMSE = root_mean_squared_error(Y_test[:,j], Y_pred[:,j])
		rmse_gbr[i,j] = RMSE/range_chem[j]*100


# np.savetxt('data/output/transfer_functions/rmse_mlr.txt', rmse_mlr, delimiter='\t',fmt='%1.4e')
# np.savetxt('data/output/transfer_functions/rmse_plsr.txt', rmse_plsr, delimiter='\t',fmt='%1.4e')
# np.savetxt('data/output/transfer_functions/rmse_gbr.txt', rmse_gbr, delimiter='\t',fmt='%1.4e')
'''

# Plot RMSE distribution for each transfer function and variable

# Load rmse from the bootstrapping  saved as output
rmse_mlr = np.loadtxt('data/output/transfer_functions/rmse_mlr.txt',delimiter='\t')
rmse_plsr = np.loadtxt('data/output/transfer_functions/rmse_plsr.txt',delimiter='\t')
rmse_gbr = np.loadtxt('data/output/transfer_functions/rmse_gbr.txt',delimiter='\t')

var = [rmse_mlr,rmse_plsr,rmse_gbr]

bins=range(0,40)
lab = ['a','b','c','d','e','f']
cols = ['POC flux','PIC:POC ratio']
rows = ['MLR','PLSR','GBR']

fig, ax = plt.subplots(3,2,figsize = (9/inch,12/inch))
for i in [0,1,2]:
	for j in [0,1]:
		data = var[i][:,j]
		mean = str(np.round(np.nanmean(data),1))
		std = str(np.round(np.nanstd(data),1))
		weights = np.zeros_like(data)+1 / data.size
		ax[i,j].hist(data,weights =weights,bins=bins,edgecolor='k',facecolor='lightgrey',lw=0.5)
		ax[i,j].axvline(np.nanmean(data),color='r',lw=0.5)
		ax[i,j].axvline(np.nanmean(data)+np.nanstd(data),color='r',lw=0.5,linestyle='--')
		ax[i,j].axvline(np.nanmean(data)-np.nanstd(data),color='r',lw=0.5,linestyle='--')
		value = mean+' %\n'+'± '+std
		ax[i,j].annotate(value,(0.65,0.72),xycoords='axes fraction',color='r')

c=0
for i in [0,1,2]:
	ax[i,0].annotate(rows[i],(-0.65,0.5),xycoords='axes fraction',rotation=90,va='center',fontweight='bold')
	for j in [0,1]:
		ax[i,j].set_xlim([0,40])
		ax[i,j].set_ylim([0,0.2])
		ax[i,j].set_xticks([0,10,20,30,40])
		ax[0,j].set_title(cols[j],fontsize=8,fontweight='bold')
		ax[2,j].set_xlabel('RMSE (%)')
		ax[i,j].annotate(lab[c],(0.07,0.85),xycoords='axes fraction',fontweight='bold')
		c+=1
	
ax[2,1].set_xlabel('RMSE (%)')
ax[1,0].set_ylabel('Relative frequency')

plt.subplots_adjust(left=0.22,right=0.95,wspace=0.5,hspace=0.4,top=0.93)
# Save figure
plt.savefig('fig/fig_05_precision.png', dpi=300)
plt.savefig('fig/fig_05_precision.pdf')

