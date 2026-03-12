#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 18 11:18:35 2025

@author: mrembauv
"""

#%% ===========================================================================
#  Execute the transfer function calibration
# =============================================================================
exec(open("./2_calibration.py").read())


# =============================================================================
#  Import libraries
# =============================================================================
# Plotting libraries
import matplotlib.pyplot as plt
plt.style.use('style.mplstyle') # Comment if the Arial font is not installed
plt.close('all')
inch = 2.54 # inch to centimeter
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

# Data analysis and statistics
import numpy as np
import pandas as pd
import netCDF4 as nc
import statsmodels.formula.api as smf
from sklearn.metrics import r2_score,mean_squared_error
from scipy.ndimage import uniform_filter1d as mv_avg

# mapping lirbaries
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.path as mpath

# =============================================================================
# Define useful functions
# =============================================================================
# Plot age model for each individual core
def plot_age_model(d_cal,age_cal,d,t,core_name):
	fig_age, ax_age = plt.subplots(1,1)
	ax_age.plot(d_cal,age_cal,'or')
	ax_age.plot(d,t,'-k')
	ax_age.set_xlabel('Depth (m)')
	ax_age.set_ylabel('Age (ka BP)')
	ax_age.set_title('Core ' + core_name)
	
# Figure for each core
window = 3 # window for moving average
def plot_core(core_name,ax_global,start,stop,window,lab_left,lab_right,lab):
	global mlr_handle,plsr_handle,gbr_handle
	data = df[df['Core']==core_name]
	data = data.sort_values(by=['Depth_m'])
	
	# Select diatom species with important mean abundace
	diat_core = data.iloc[:,6:].to_numpy().astype('float64')
	diat_core = diat_core[:,diat_sel_total]
	diat_core = diat_core[:,diat_sel_percent]
	
	# Normalise the selected diatom data to 100 %
	diat_core = (diat_core.T/np.sum(diat_core,axis=1)*100).T

	# Apply the MLR model to diatom groups
	mlr_pred = np.empty([diat_core.shape[0],2])
	X = np.log(diat_core+1)
	X = StandardScaler().fit(X).transform(X)
	X = fa.transform(X)
	df_mlr = pd.DataFrame(X,columns=group_names)
	mlr_pred[:,0] = mlr1.predict(df_mlr)
	mlr_pred[:,1]= mlr2.predict(df_mlr)

	# Apply the PLS, GBR to the full species list
	X = np.log(diat_core+1)
	X = StandardScaler().fit(X).transform(X)
	plsr_pred = plsr.predict(X)
	gbr_pred = np.empty([diat_core.shape[0],2])
	gbr_pred[:,0] = gbr1.predict(X)
	gbr_pred[:,1] = gbr2.predict(X)
	
	
	# Convert predicted POC flux to annual flux (molC/m2/y)
	mlr_pred[:,0] = mlr_pred[:,0]*365/1000
	plsr_pred[:,0] = plsr_pred[:,0]*365/1000
	gbr_pred[:,0] = gbr_pred[:,0]*365/1000
	
	# Create mutliproxy dataset to be saved
	multiproxy = np.vstack([d,t,mlr_pred.T,plsr_pred.T,gbr_pred.T]).T
	fname = core_name.replace('/','_')+'.txt'
	header = 'depth_[cm]\ttime_[kyBP]\tMLR_POC_flux_[mol/m2/y]\tMLR_PIC:POC_[mol:mol]\tPLSR_POC_flux_[mol/m2/y]\tPLSR_PIC:POC_[mol:mol]\tGBR_POC_flux_[mol/m2/y]\tGBR_PIC:POC_[mol:mol]'
	np.savetxt('data/output/reconstructions/'+fname,multiproxy,header=header,delimiter='\t',fmt='%1.3f')
	
	# axes properties
	ax_global.set_xlabel('Age (ka)')
	ax_global.set_xlim([start,stop])
	ax_global.patch.set_alpha(0)
	ax_global.set_yticks([])
	ax_global.set_zorder(-10)
	
	pos = ax_global.get_position()
	ax_height = pos.height*0.6
	ax1 = fig.add_axes([pos.x0,pos.y1-ax_height,pos.width,ax_height])
	ax2 = fig.add_axes([pos.x0,pos.y0,pos.width,ax_height])
	for side in ['top','right','bottom']:
		ax1.spines[side].set_visible(False)
	for side in ['top','left','bottom']:
		ax2.spines[side].set_visible(False)
	ax2.yaxis.tick_right()
	ax2.yaxis.set_label_position("right")
	
	axes = [ax1,ax2]
	for i in [0,1]:
	
		# Plot MLR result
		mlr_handle, = axes[i].plot(t,mlr_pred[:,i],'-b',lw=0.5,alpha=0.5,label='MLR')
		# PLot PLSR results
		plsr_handle, = axes[i].plot(t,plsr_pred[:,i],'-g',lw=0.5,alpha=0.5,label='PLSR')
		# Plot GBR results
		gbr_handle, = axes[i].plot(t,gbr_pred[:,i],'-r',lw=0.5,alpha=0.5,label='GBR')
		
		# Plot multi-proxy average
		multi_proxy = np.nanmean(np.vstack([mlr_pred[:,i],plsr_pred[:,i],gbr_pred[:,i]]),axis=0)
		axes[i].plot(t,mv_avg(multi_proxy,size=window),'-k',lw=0.5)
		
		# Tme range and remove x ticks
		axes[i].set_xlim(start,stop)
		axes[i].set_xticks([])
		axes[i].patch.set_alpha(0)
		
	# Y ranges
	ax1.set_ylim([-0.2,0.6])
	ax1.set_yticks([0,0.2,0.4,0.6])
	ax2.set_ylim([-0,1.6])
	
	#X minor locator
	ax_global.xaxis.set_minor_locator(AutoMinorLocator())
	
	# ACR and MIS
	ax_global.axvspan(13,14.7, color='lightblue') # ACR
	ax_global.axvspan(14,29, color='lightgrey') # MIS2
	ax_global.axvspan(57,71, color='lightgrey') # MIS4
	ax_global.axvspan(130,190, color='lightgrey') # MIS6
	
	# Axes labels
	if lab_left == True :
		ax1.set_ylabel('POC flux\n$(mol\ m^{-2}\ y^{-1})$')
	if lab_right == True :
		ax2.set_ylabel('PIC:POC\n(mol:mol)')
	
	# Core name and panel label
	ax_global.annotate(core_name,(0.5,0.93),xycoords='axes fraction',fontsize=8,horizontalalignment='center',fontweight='bold')
	ax_global.annotate(lab,(0.03,0.93),xycoords='axes fraction',fontsize=8,fontweight='bold')
	
	return mlr_pred,plsr_pred,gbr_pred


# Define box plot function
def box_plot(ax,data,pos,fc,lw):
	boxprops = dict(linewidth=lw)
	whiskerprops = dict(linewidth=lw)
	capprops = dict(linewidth=lw)
	medianprops = dict(linewidth=lw)
	meanprops = dict(linewidth=lw)
	bplot = ax.boxplot(data, positions=pos,widths=0.15,whis=1,showfliers=False,patch_artist=True,
					showmeans=True,meanline=True,
					boxprops=boxprops,whiskerprops=whiskerprops,capprops=capprops,medianprops=medianprops,meanprops=meanprops)
	for patch in bplot['boxes'] :
	    patch.set_facecolor(fc)
	for patch in bplot['medians']:
		patch.set_color('k')
	for patch in bplot['means']:
		patch.set_color('k')
	return bplot


# Paired U-test (MIS 2 vs. Holo)
def stat_test(data,ax,pos,offset):
	from scipy.stats import mannwhitneyu
	y0 = np.mean(data[0])+(np.quantile(data[0],0.75)-np.quantile(data[0],0.25))
	y1 = np.mean(data[1])+(np.quantile(data[1],0.75)-np.quantile(data[1],0.25))
	y = max([y0,y1])+offset
	p = mannwhitneyu(data[0],data[1], alternative='two-sided').pvalue
	ax.plot([pos-0.1,pos+0.1],[y,y],'-k',lw=0.5)
	if p<0.01:	
		ax.text(pos,y,'**',ha='center')
	elif (p>=0.01) & (p<0.05) :
		ax.text(pos,y,'*',ha='center')
	else:
		ax.text(pos,y+0.02,'ns',ha='center')

# Open reconstructions and extract data by sector and period
def extract_data(file,sector) :
	global atl_poc_holo,atl_poc_mis,atl_pic_poc_holo,atl_pic_poc_mis,ind_poc_holo,ind_poc_mis,ind_pic_poc_holo,ind_pic_poc_mis,pac_poc_holo,pac_poc_mis,pac_pic_poc_holo,pac_pic_poc_mis
	# Define periods
	holo_start = 11.7
	mis_start = 29
	mis_stop = 14
	
	# Extract data for each period
	proxy = np.loadtxt(file,delimiter='\t',skiprows=1)
	t = proxy[:,1]
	poc_mean = np.nanmean(proxy[:,[2,4,6]],axis=1)
	pic_poc_mean = np.nanmean(proxy[:,[3,5,7]],axis=1)
	holo = (t<holo_start)
	mis = (t<mis_start) & (t>mis_stop)
	
	# Store data in the appropriate array
	if sector=='atl' :
		atl_poc_holo = np.hstack([atl_poc_holo,poc_mean[holo]])
		atl_pic_poc_holo = np.hstack([atl_pic_poc_holo,pic_poc_mean[holo]])
		atl_poc_mis = np.hstack([atl_poc_mis,poc_mean[mis]])
		atl_pic_poc_mis = np.hstack([atl_pic_poc_mis,pic_poc_mean[mis]])
	elif sector=='ind' :
		ind_poc_holo = np.hstack([ind_poc_holo,poc_mean[holo]])
		ind_pic_poc_holo = np.hstack([ind_pic_poc_holo,pic_poc_mean[holo]])
		ind_poc_mis = np.hstack([ind_poc_mis,poc_mean[mis]])
		ind_pic_poc_mis = np.hstack([ind_pic_poc_mis,pic_poc_mean[mis]])
	elif sector=='pac' :
		pac_poc_holo = np.hstack([pac_poc_holo,poc_mean[holo]])
		pac_pic_poc_holo = np.hstack([pac_pic_poc_holo,pic_poc_mean[holo]])
		pac_poc_mis = np.hstack([pac_poc_mis,poc_mean[mis]])
		pac_pic_poc_mis = np.hstack([pac_pic_poc_mis,pic_poc_mean[mis]])

#%% =============================================================================
#  Open sediment dataset
# =============================================================================

# Open sediment core dataset
df = pd.read_excel('data/sediment_cores/sed_core_data.ods', engine='odf', sheet_name='data')
df = clean_diatom_species(df)
df = df.drop(0) # error from core PS58_270-5 : first line is an outlier
lat_core = df.groupby('Core')['Lat_degN'].mean()
lon_core = df.groupby('Core')['Lon_degE'].mean()

# Open Fronts derived from altimetry (Park 2019)
data_fronts = nc.Dataset('data/mapping/park_fronts_2019.nc')
lon_SAF,lat_SAF = data_fronts['LonSAF'][:].filled(np.nan),data_fronts['LatSAF'][:].filled(np.nan)
lon_PF,lat_PF = data_fronts['LonPF'][:],data_fronts['LatPF'][:]
# Clean outliers to avoid projection mistake
lon_SAF[lon_SAF == max(lon_SAF)]=np.nan
lon_PF[lon_PF == max(lon_PF)]=np.nan

# Open Sea ice edge location
lon_ice,lat_ice = np.loadtxt('data/mapping/sea_ice_extent_september_1980_2010.txt',skiprows=1,unpack=True)
lon_ice[lon_ice == max(lon_ice)]=np.nan
	
#%% ===========================================================================
# Apply transfer functions and plot reconstructed fluxes (Figure 6)
# =============================================================================

# Create figure	
plt.close('all')
fig = plt.figure(figsize=(20/inch,25/inch))
gs = GridSpec(4, 3, figure=fig, left=0.1,bottom=0.045, right=0.92,top=0.98,wspace=0.35,hspace=0.25)

# Map in the center
data_proj = ccrs.PlateCarree()
ax_map = fig.add_subplot(gs[1:3,1:2],projection=ccrs.SouthPolarStereo())
ax_map.set_extent([-180, 180, -90, -40],data_proj)

# Add coastline and land
land = cfeature.NaturalEarthFeature(name='coastline', category='physical',
    scale='50m',facecolor='lightgrey',edgecolor='k',linewidth=0.5)
ax_map.add_feature(land)

# Add seasonal ice zone
sie, = ax_map.plot(lon_ice,lat_ice,'-b',lw=0.5,transform=data_proj,label='SIE')
pf, = ax_map.plot(lon_PF,lat_PF,'-g',transform=data_proj,lw=0.5,label='PF')
saf, = ax_map.plot(lon_SAF,lat_SAF,'-r',transform=data_proj,lw=0.5,label='SAF')

# Add sediment cores location
cores = ax_map.scatter(lon_core,lat_core,marker='s',transform=data_proj,color='k',linewidth=0.5,s=10,fc='yellow',zorder=15)

# Set circular axis map boundary
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)
ax_map.set_boundary(circle, transform=ax_map.transAxes)

# Add oceans delimitations
ax_map.plot([-70.8,-70.6],[-80,-52.6],'-k',transform=data_proj,lw=0.5,zorder=-15)
ax_map.plot([20,20.5],[-80,-32],'-k',transform=data_proj,lw=0.5,zorder=-15)
ax_map.plot([148,148],[-80,-42.2],'-k',transform=data_proj,lw=0.5,zorder=-15)
	

# =============================================================================
# Atlantic sector
# =============================================================================

# Core PS97/72-1 (Vorrath et al., 2023) ==============================
core_name = 'PS97/72-1'
data = df[df['Core']==core_name]
data = data.sort_values(by=['Depth_m'])
d = data['Depth_m']
age_file = 'data/sediment_cores/vorrath/PS97_72-1_age.tab'
d_cal,age_cal = np.loadtxt(age_file,skiprows=25,usecols=[0,9],unpack=True)
d_cal = np.hstack([0,d_cal,12])
age_cal = np.hstack([0,age_cal,16.5])
t = np.interp(d,d_cal,age_cal)

#plot_age_model(d_cal,age_cal,d,t,core_name)

ax = fig.add_subplot(gs[0,0])
mlr_pred,plsr_pred,gbr_pred = plot_core(core_name,ax,0,15,window,lab_left=True, lab_right=False,lab='a')
ax.legend(handles = [mlr_handle,plsr_handle,gbr_handle],
		  frameon=False,loc='center left',bbox_to_anchor=[0.5,0.45],
		  fontsize=8,labelspacing=0.2)


# Core PS1786-1 (Esper and Gersonde, 2008) ==============================
core_name = 'PS1786-1'
# Age model from Jacot Des Combes et al., 2008
data = df[df['Core']==core_name]
data = data.sort_values(by=['Depth_m'])
d = data['Depth_m']
d_cal = np.array([0.0275,0.7325,1.23,1.825,2.35])
age_cal = np.array([1.895,5.581,9.534,18.327,27.5])
d_cal = np.hstack([d_cal])
age_cal = np.hstack([age_cal])
t = np.interp(d,d_cal,age_cal)

#plot_age_model(d_cal,age_cal,d,t,core_name)

ax = fig.add_subplot(gs[0,1])
mlr_pred,plsr_pred,gbr_pred = plot_core(core_name,ax,0,28,window,lab_left=False, lab_right=False,lab='b')


# Core PS1768-8 (Zielinski et al., 1998) ==============================
core_name = 'PS1768-8'
data = df[df['Core']==core_name]
data = data.sort_values(by=['Depth_m'])
d = data['Depth_m']
age_file = 'data/sediment_cores/zielinski/'+core_name+'_age_model.tab'
d_cal,age_cal = np.loadtxt(age_file,skiprows=13,usecols=[0,1],unpack=True)
d_cal = np.hstack([0,d_cal])
age_cal = np.hstack([0,age_cal])
t = np.interp(d,d_cal,age_cal)

#plot_age_model(d_cal,age_cal,d,t,core_name)

ax = fig.add_subplot(gs[0,2])
mlr_pred,plsr_pred,gbr_pred = plot_core(core_name,ax,0,150,window,lab_left=False, lab_right=True,lab='c')


# =============================================================================
# Indian sector
# =============================================================================

# Core KH-10-7 COR1GC (Orme et al. 2020) =================================
core_name = 'COR1GC'
data = df[df['Core']==core_name]
data = data.sort_values(by=['Depth_m'])
t = data['Age_kyBP'] # Age already in database for this core
d = data['Depth_m']

ax = fig.add_subplot(gs[1,2])
mlr_pred,plsr_pred,gbr_pred = plot_core(core_name,ax,0,15,window,lab_left=False, lab_right=True,lab='d')

# Core PS2606-6 (Jacot des Combes et al., 2008) ===========================
core_name = 'PS2606-6'
data = df[df['Core']==core_name]
data = data.sort_values(by=['Depth_m'])
d = data['Depth_m']
# Extended Age model from Civel-Mazens et al. 2024
age_file = 'data/sediment_cores/esper/'+core_name+'_SSTd.tab'
d_cal,age_cal = np.loadtxt(age_file,skiprows=18,usecols=[0,1],unpack=True)
t = np.interp(d,d_cal,age_cal)
#plot_age_model(d_cal,age_cal,d,t,core_name)

ax = fig.add_subplot(gs[2,2])
mlr_pred,plsr_pred,gbr_pred = plot_core(core_name,ax,0,50,window,lab_left=False, lab_right=True,lab='e')



# =============================================================================
# Pacific sector
# =============================================================================

# Core PS75/072-4 (Benz et al., 2016) ============================
core_name = 'PS75/072-4'
data = df[df['Core']==core_name]
data = data.sort_values(by=['Depth_m'])
#data = data[data['Depth_m']<1.7] # Keep only the well dated part of the core
d = data['Depth_m']
age_file = 'data/sediment_cores/benz/PS75_072-4_age_model.tab'
d_cal,age_cal = np.loadtxt(age_file,skiprows=15,usecols=[0,1],unpack=True)

d_cal = np.hstack([0,d_cal])
age_cal = np.hstack([0,age_cal])
t = np.interp(d,d_cal,age_cal)

#plot_age_model(d_cal,age_cal,d,t,core_name)

ax = fig.add_subplot(gs[3,2])
mlr_pred,plsr_pred,gbr_pred = plot_core(core_name,ax,0,150,window,lab_left=False, lab_right=True,lab='f')


# Core PS58/270-5 (Benz et al., 2005) ==============================
core_name = 'PS58/270-5'
data = df[df['Core']==core_name]
data = data.sort_values(by=['Depth_m'])
d = data['Depth_m']
age_file = 'data/sediment_cores/benz/PS58_270-5_age_model.tab'
d_cal,age_cal = np.loadtxt(age_file,skiprows=15,usecols=[0,1],unpack=True)
d_cal = np.hstack([d_cal])
age_cal = np.hstack([age_cal])
t = np.interp(d,d_cal,age_cal)

#plot_age_model(d_cal,age_cal,d,t,core_name)

ax = fig.add_subplot(gs[3,1])
mlr_pred,plsr_pred,gbr_pred = plot_core(core_name,ax,0,150,window,lab_left=False, lab_right=False,lab='g')


# Core PS58/271-1 (Benz et al., 2005) =============================
core_name = 'PS58/271-1'
data = df[df['Core']==core_name]
data = data.sort_values(by=['Depth_m'])
d = data['Depth_m']
age_file = 'data/sediment_cores/benz/PS58_271-1_age_model.tab'
d_cal,age_cal = np.loadtxt(age_file,skiprows=15,usecols=[0,1],unpack=True)
d_cal = np.hstack([d_cal,25])
age_cal = np.hstack([age_cal,150])
t = np.interp(d,d_cal,age_cal)

#plot_age_model(d_cal,age_cal,d,t,core_name)

ax = fig.add_subplot(gs[3,0])
mlr_pred,plsr_pred,gbr_pred = plot_core(core_name,ax,0,150,window,lab_left=True, lab_right=False,lab='h')


# Core PS75/054-1 (Benz et al., 2016) =============================
core_name = 'PS75/054-1'
data = df[df['Core']==core_name]
data = data.sort_values(by=['Depth_m'])
d = data['Depth_m']
age_file = 'data/sediment_cores/benz/PS75_054-1_age_model.tab'
d_cal,age_cal = np.loadtxt(age_file,skiprows=15,usecols=[0,1],unpack=True)
d_cal = d_cal[[0,1,2,3,7,8,9,10,11,12,13]] # Remove age data comming from core PS75/056-1 (other core)
age_cal = age_cal[[0,1,2,3,7,8,9,10,11,12,13]] # Remove age data comming from core PS75/056-1 (other core)
t = np.interp(d,d_cal,age_cal)

#plot_age_model(d_cal,age_cal,d,t,core_name)

ax = fig.add_subplot(gs[2,0])
mlr_pred,plsr_pred,gbr_pred = plot_core(core_name,ax,6,40,window,lab_left=True, lab_right=False,lab='i')


# Core PS58/274-1 (Benz et al., 2016) ===========================
core_name = 'PS58/274-1'
data = df[df['Core']==core_name]
data = data.sort_values(by=['Depth_m'])
d = data['Depth_m']
age_file = 'data/sediment_cores/benz/PS58_274-1_age_model.tab'
d_cal,age_cal = np.loadtxt(age_file,skiprows=15,usecols=[0,1],unpack=True)
d_cal = np.hstack([d_cal,25])
age_cal = np.hstack([age_cal,100])
t = np.interp(d,d_cal,age_cal)

#plot_age_model(d_cal,age_cal,d,t,core_name)

ax = fig.add_subplot(gs[1,0])
mlr_pred,plsr_pred,gbr_pred = plot_core(core_name,ax,0,100,window,lab_left=True, lab_right=False,lab='j')


# Save figure
plt.savefig('fig/fig_06_cores.png',dpi=300)
plt.savefig('fig/fig_06_cores.svg')
plt.savefig('fig/fig_06_cores.pdf')



#%% ===========================================================================
#  Boxplots comparing biological pump across basins during MIS 2
# =============================================================================

# Define empty lists to store extracted data

atl_poc_holo = np.array([])
atl_pic_poc_holo =  np.array([])
atl_poc_mis = np.array([])
atl_pic_poc_mis =  np.array([])

ind_poc_holo = np.array([])
ind_pic_poc_holo = np.array([])
ind_poc_mis = np.array([])
ind_pic_poc_mis = np.array([])

pac_poc_holo = np.array([])
pac_pic_poc_holo = np.array([])
pac_poc_mis = np.array([])
pac_pic_poc_mis = np.array([])


# Atlantic sector
extract_data('data/output/reconstructions/PS1786-1.txt',sector='atl')
extract_data('data/output/reconstructions/PS1768-8.txt',sector='atl')

# Indian sector
extract_data('data/output/reconstructions/COR1GC.txt',sector='ind')
extract_data('data/output/reconstructions/PS2606-6.txt',sector='ind')

# Pacfic sector
extract_data('data/output/reconstructions/PS75_072-4.txt',sector='pac')
extract_data('data/output/reconstructions/PS58_270-5.txt',sector='pac')
extract_data('data/output/reconstructions/PS58_271-1.txt',sector='pac')
extract_data('data/output/reconstructions/PS75_054-1.txt',sector='pac')
extract_data('data/output/reconstructions/PS58_274-1.txt',sector='pac')


# Modern annual flux data : average trap data by station to get annual flux data
glob_poc_mod = np.array([])
glob_pic_poc_mod = np.array([])
stations = np.unique(sta)
for i in range(len(stations)) :
	sel = (sta==stations[i])
	glob_poc_mod = np.append(glob_poc_mod, np.nanmean(chem[sel,0]))
	glob_pic_poc_mod = np.append(glob_pic_poc_mod, np.nanmean(chem[sel,1]))
glob_poc_mod = glob_poc_mod*365/1000

# Merge all basins for global SO
glob_poc_mis = np.hstack([atl_poc_mis,ind_poc_mis,pac_poc_mis])
glob_pic_poc_mis = np.hstack([atl_pic_poc_mis,ind_pic_poc_mis,pac_pic_poc_mis])
glob_poc_holo = np.hstack([atl_poc_holo,ind_poc_holo,pac_poc_holo])
glob_pic_poc_holo = np.hstack([atl_pic_poc_holo,ind_pic_poc_holo,pac_pic_poc_holo])


# Create figure
plt.close('all')
fig, ax = plt.subplots(2,1,figsize=(9/inch,12/inch))

# Add modern POC values from sed trap
mean = np.mean(glob_poc_mod)
std = np.std(glob_poc_mod)
ax[0].fill_between([0,5],mean-std,mean+std,color='lightgrey',zorder=-20)
ax[0].axhline(mean,lw=1,color='k',linestyle='-',zorder=-10)

# Add modern PIC:POC values from sed trap
mean = np.mean(glob_pic_poc_mod)
std = np.std(glob_pic_poc_mod)
ax[1].fill_between([0,5],mean-std,mean+std,color='lightgrey',zorder=-20)
ax[1].axhline(mean,lw=1,color='k',linestyle='-',zorder=-10)

# Boxplots
pos1 = np.array([0.85,1.85,2.85])
pos2 = np.array([1.15,2.15,3.15])
bp_holo = box_plot(ax[0],[atl_poc_holo,ind_poc_holo,pac_poc_holo],pos=pos1,fc='salmon',lw=0.5)
bp_mis = box_plot(ax[0],[atl_poc_mis,ind_poc_mis,pac_poc_mis],pos=pos2,fc='lightblue',lw=0.5)
box_plot(ax[1],[atl_pic_poc_holo,ind_pic_poc_holo,pac_pic_poc_holo],pos=pos1,fc='salmon',lw=0.5)
box_plot(ax[1],[atl_pic_poc_mis,ind_pic_poc_mis,pac_pic_poc_mis],pos=pos2,fc='lightblue',lw=0.5)
# GLobal SO in bold
box_plot(ax[0],glob_poc_holo,pos=np.array([3.85]),fc='salmon',lw=1)
box_plot(ax[0],glob_poc_mis,pos=np.array([4.15]),fc='lightblue',lw=1)
box_plot(ax[1],glob_pic_poc_holo,pos=np.array([3.85]),fc='salmon',lw=1)
box_plot(ax[1],glob_pic_poc_mis,pos=np.array([4.15]),fc='lightblue',lw=1)

# Add U-test significance levles
offset=0.1
stat_test([atl_poc_holo,atl_poc_mis],ax[0],1,offset-0.07)
stat_test([ind_poc_holo,ind_poc_mis],ax[0],2,offset-0.1)
stat_test([pac_poc_holo,pac_poc_mis],ax[0],3,offset-0.03)
stat_test([glob_poc_holo,glob_poc_mis],ax[0],4,offset-0.025)
offset=0.1
stat_test([atl_pic_poc_holo,atl_pic_poc_mis],ax[1],1,offset+0.05)
stat_test([ind_pic_poc_holo,ind_pic_poc_mis],ax[1],2,offset+0.12)
stat_test([pac_pic_poc_holo,pac_pic_poc_mis],ax[1],3,offset+0.05)
stat_test([glob_pic_poc_holo,glob_pic_poc_mis],ax[1],4,offset+0.1)
		
# Add legend
ax[0].legend([bp_holo["boxes"][0], bp_mis["boxes"][0]], ['Holocene', 'MIS 2'],
			 loc='upper center',bbox_to_anchor=(0.5, 1.15),frameon=False,ncols=2)

# Set axis labels and limits
for i in [0,1]:
	ax[i].set_xlim([0.5,4.5])
	ax[i].set_xticks([1,2,3,4])
	ax[i].set_xticklabels([])
ax[0].set_ylim([-0.1,0.5])
ax[0].set_yticks([0,.1,.2,.3,.4,.5])
ax[1].set_ylim([0,1.4])

# Addd samples number
pos1 = np.array([0.85,1.85,2.85,3.85])
pos2 = np.array([1.15,2.15,3.15,4.15])
obs = [str(len(atl_poc_holo)),
	   str(len(atl_poc_mis)),
	   str(len(ind_poc_holo)),
	   str(len(ind_poc_mis)),
	   str(len(pac_poc_holo)),
	   str(len(pac_poc_mis)),
	   str(len(glob_poc_holo)),
	   str(len(glob_poc_mis))]
for i in np.arange(len(obs)):
	ax[0].text(np.sort(np.hstack([pos1,pos2]))[i],-0.09,obs[i],ha='center')
	ax[1].text(np.sort(np.hstack([pos1,pos2]))[i],0.02,obs[i],ha='center')
ax[0].set_ylabel('POC flux $(mol\ m^{-2}\ y^{-1})$')
ax[1].set_ylabel('PIC:POC (mol:mol)')
ax[1].set_xticklabels(['Atlantic','Indian','Pacific','Global SO'])
ax[1].get_xticklabels()[-1].set_fontweight('bold')

ax[0].annotate('a',(0.025,0.93),xycoords='axes fraction',fontsize=8,fontweight='bold')
ax[1].annotate('b',(0.025,0.93),xycoords='axes fraction',fontsize=8,fontweight='bold')

plt.tight_layout()
# Save figure 
plt.savefig('fig/fig_07_MIS2_H_comparison.png',dpi=300)
plt.savefig('fig/fig_07_MIS2_H_comparison.pdf')

# =============================================================================
# Print MIS 2 vs. Holocene changes
# =============================================================================
print('Modern sediment trap')
print('Mean POC flux modern :',str(np.mean(glob_poc_mod)))
print('Mean PIC:POC modern :',str(np.mean(glob_pic_poc_mod)))
print('\n')

print('Global SO')
print('Mean POC flux MIS 2 :',str(np.mean(glob_poc_mis)))
print('Mean POC flux holo :',str(np.mean(glob_poc_holo)))
print('Mean PIC:POC MIS 2 :',str(np.mean(glob_pic_poc_mis)))
print('Mean PIC:POC holo :',str(np.mean(glob_pic_poc_holo)))
print('\n')

print('Atlantic')
print('Mean POC flux MIS 2 :',str(np.mean(atl_poc_mis)))
print('Mean POC flux holo :',str(np.mean(atl_poc_holo)))
print('Mean PIC:POC MIS 2 :',str(np.mean(atl_pic_poc_mis)))
print('Mean PIC:POC holo :',str(np.mean(atl_pic_poc_holo)))
print('\n')

print('Indian')
print('Mean POC flux MIS 2 :',str(np.mean(ind_poc_mis)))
print('Mean POC flux holo :',str(np.mean(ind_poc_holo)))
print('Mean PIC:POC MIS2 :',str(np.mean(ind_pic_poc_mis)))
print('Mean PIC:POC holo :',str(np.mean(ind_pic_poc_holo)))
print('\n')

print('Pacific')
print('Mean POC flux MIS 2 :',str(np.mean(pac_poc_mis)))
print('Mean POC flux holo :',str(np.mean(pac_poc_holo)))
print('Mean PIC:POC MIS 2 :',str(np.mean(pac_pic_poc_mis)))
print('Mean PIC:POC holo :',str(np.mean(pac_pic_poc_holo)))
print('\n')


print('POC flux decrease from MIS 2 to holo :')
print((np.mean(glob_poc_mis)-np.mean(glob_poc_holo))/np.mean(glob_poc_mis))
print('POC flux decrease from MIS 2 to modern :')
print((np.mean(glob_poc_mis)-np.mean(glob_poc_mod))/np.mean(glob_poc_mis))
print('PIC:POC ratio increase from MIS 2 to holo :')
print((np.mean(glob_pic_poc_holo)-np.mean(glob_pic_poc_mis))/np.mean(glob_pic_poc_mis))
print('PIC:POC ratio increase from MIS 2 to modern :')
print((np.mean(glob_pic_poc_mod)-np.mean(glob_pic_poc_mis))/np.mean(glob_pic_poc_mis))
