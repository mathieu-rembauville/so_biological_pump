#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 23 11:18:12 2025

@author: mrembauv
"""
# Plotting libraries
import matplotlib.pyplot as plt
plt.style.use('style.mplstyle') # Comment if the Arial font is not installed
plt.close('all')
inch = 2.54 # inch to centimeter

# Data analysis and statistics
import numpy as np
import pandas as pd
import netCDF4 as nc

# mapping lirbaries
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.path as mpath


#%% ===========================================================================
# Open datasets
# =============================================================================

# Open core locations
data_core = pd.read_excel('data/sediment_cores/sed_core_data.ods', engine='odf', sheet_name='data')
lat_core = data_core.groupby('Core')['Lat_degN'].mean()
lon_core = data_core.groupby('Core')['Lon_degE'].mean()
name_core = np.unique(data_core['Core'])

# Open trap locations
data_trap = pd.read_excel('data/sediment_traps/sed_trap_data.ods', engine='odf', sheet_name='cup_scale')
# Remove unreliable data (diatom taxonomy unresolved or too sahllow sed trap deployments)
data_trap = data_trap[(data_trap['reference']!='Salter_2012')]
data_trap = data_trap[(data_trap['reference']!='Blain_2021_2022')]
data_trap = data_trap[(data_trap['reference']!='Rembauville_2015')]
data_trap = data_trap[(data_trap['reference']!='Rembauville_2017')]

# Extract data
lat_trap = data_trap.groupby('station')['lat_degN'].mean()
lon_trap = data_trap.groupby('station')['lon_degE'].mean()
name_trap = np.unique(data_trap['station'])
zone = data_trap.groupby('station')['zone'].unique()

# Declare colors
zone_col = ['red','lime','deepskyblue','white']
zone_list = ['SAZ','PFZ','POOZ','SIZ']

#%% ===========================================================================
# Map sediment trap and sediment core locations
# =============================================================================

# Create figure
plt.close('all')
fig = plt.figure(figsize=[12/inch,10/inch])
data_proj = ccrs.PlateCarree()
ax = plt.subplot(1, 1, 1, projection=ccrs.SouthPolarStereo())
ax.set_extent([-180, 180, -90, -40],data_proj)

# Add ETOPO bathymetry in the background
img = plt.imread('data/mapping/etopo.jpg')
img_extent = (-180, 180, -90, 90)
ax.imshow(img, origin='upper', extent=img_extent, transform=data_proj,alpha=0.75)

# Add coastline and land
land = cfeature.NaturalEarthFeature(name='coastline', category='physical',
    scale='50m',facecolor='lightgrey',edgecolor='k',linewidth=0.3)
ax.add_feature(land)

# Open fronts derived from satellite altimetry (Park 2019)
data_fronts = nc.Dataset('data/mapping/park_fronts_2019.nc')
lon_SAF,lat_SAF = data_fronts['LonSAF'][:].filled(np.nan),data_fronts['LatSAF'][:].filled(np.nan)
lon_PF,lat_PF = data_fronts['LonPF'][:],data_fronts['LatPF'][:]
# Clean data to avoid projection mistake
lon_SAF[lon_SAF == max(lon_SAF)]=np.nan
lon_PF[lon_PF == max(lon_PF)]=np.nan
saf, = ax.plot(lon_SAF,lat_SAF,'-r',transform=data_proj,lw=1,label='SAF')
pf, = ax.plot(lon_PF,lat_PF,'-g',transform=data_proj,lw=1,label='PF')

# Sea ice edge
lon_ice,lat_ice = np.loadtxt('data/mapping/sea_ice_extent_september_1980_2010.txt',skiprows=1,unpack=True)
lon_ice[lon_ice == max(lon_ice)]=np.nan
sie, = ax.plot(lon_ice,lat_ice,'b',lw=1,transform=data_proj,label='SIE')

fs = 8
# Add sediment trap location
traps = []
for i in range(len(zone_list)):
	h = ax.scatter(lon_trap[zone==zone_list[i]],lat_trap[zone==zone_list[i]],transform=data_proj,	
			facecolor=zone_col[i],edgecolors='k',s=15,lw=0.5,label=zone_list[i],zorder=10)
	traps.append(h)

# Add sediment cores location
cores = ax.scatter(lon_core,lat_core,marker='s',transform=data_proj,lw=0.5,color='k',s=15,fc='yellow',label='Core name',zorder=20)

# Set circular axis boundary
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)
ax.set_boundary(circle, transform=ax.transAxes)

# add legend
leg1 = plt.legend(handles=[saf,pf,sie],frameon=False,bbox_to_anchor=(-0.06,0.6),alignment='left',fontsize=8,
				  title='Fronts',title_fontproperties={'weight':'bold'})
leg2 = plt.legend(handles=traps,frameon=False,bbox_to_anchor=(0.02,0.38),fontsize=8,
				  title='Sediment traps',alignment='left',title_fontproperties={'weight':'bold'})
leg3 = plt.legend(handles=[cores],frameon=False,bbox_to_anchor=(0.035,0.13),fontsize=8,
				  title='Sediment cores',alignment='left',title_fontproperties={'weight':'bold'})
plt.gca().add_artist(leg1)
plt.gca().add_artist(leg2)

plt.subplots_adjust(bottom=0.05,top=0.95,right=1.05)

# Save figure
plt.savefig('fig/fig_01_map_raw.png',dpi=300)
plt.savefig('fig/fig_01_map_raw.pdf',dpi=300)
plt.savefig('fig/fig_01_map_raw.svg',dpi=300)


