#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 16:12:29 2025

@author: mrembauv
"""

# =============================================================================
#  Import Libraries
# =============================================================================
import matplotlib.pyplot as plt
plt.style.use('style.mplstyle') # Comment if the Arial font is not installed
plt.close('all')
inch = 2.54 # inch to centimeter
from matplotlib.ticker import (AutoMinorLocator,FixedLocator)

# Data analysis and statistics
import numpy as np

# =============================================================================
# Define sub_axis function (stacked paleo-plots)
# =============================================================================
def sub_axis(ax,loc):
	if loc=='left':
		for side in ['top','right','bottom']:
			ax.spines[side].set_visible(False)
	elif loc=='right':
		for side in ['top','left','bottom']:
			ax.spines[side].set_visible(False)
		ax.yaxis.tick_right()
		ax.yaxis.set_label_position("right")

# Define colors
lightblue = 'cornflowerblue'
lightred = 'lightcoral'
ice = 'w'
alpha = 0.5


#%% ===========================================================================
# Core 1768-8 (Atlantic) - comparison with %TOC and % CaCO3
# =============================================================================

# Declare figure and age range
plt.close('all')
start = 0
stop = 150
fig,ax_global = plt.subplots(1,1,figsize=(9/inch,15/inch))
ax_global.set_xlabel('Age (ka)')
ax_global.set_xlim([start,stop])
ax_global.patch.set_alpha(0)
ax_global.set_yticks([])
ax_global.set_zorder(-10)
plt.subplots_adjust(left=0.16,right=0.82,bottom=0.07,top=0.95)

# Add axes
pos = ax_global.get_position()
ax_height = pos.height*0.3
ax1 = fig.add_axes([pos.x0,pos.y1-ax_height,pos.width,ax_height]) # pCO2
sub_axis(ax1,'left')
ax2 = fig.add_axes([pos.x0,pos.y1-ax_height*0.8,pos.width,ax_height*0.8]) # dust
sub_axis(ax2,'right')
ax3 = fig.add_axes([pos.x0,pos.y0+1.4*ax_height,pos.width,ax_height*1.1]) # POC export flux proxy
sub_axis(ax3,'left')
ax4 = fig.add_axes([pos.x0,pos.y0+1.4*ax_height,pos.width,ax_height*1.1]) # TOC sediment core
sub_axis(ax4,'right')
ax5 = fig.add_axes([pos.x0,pos.y0+0.5*ax_height,pos.width,ax_height*1.1]) # PIC:POC export proxy
sub_axis(ax5,'left')
ax6 = fig.add_axes([pos.x0,pos.y0+0.5*ax_height,pos.width,ax_height]) # PIC sediment core
sub_axis(ax6,'right')
ax7 = fig.add_axes([pos.x0,pos.y0,pos.width,ax_height*0.5]) # IRD
sub_axis(ax7,'left')

axes = [ax1,ax2,ax3,ax4,ax5,ax6,ax7]
for i in range(len(axes)):
	# Tme range and remove x ticks
	axes[i].set_xlim(start,stop)
	axes[i].set_xticks([])
	axes[i].patch.set_alpha(0)
	
# Add MIS and ACR
ax_global.xaxis.set_minor_locator(AutoMinorLocator())
ax_global.axvspan(13,14.7, color='lightblue') # ACR
ax_global.axvspan(14,29, color='lightgrey') # MIS2
ax_global.axvspan(57,71, color='lightgrey') # MIS4
ax_global.axvspan(130,190, color='lightgrey') # MIS6

# Plot pCO2
t,pco2,sigma = np.loadtxt('data/sediment_cores/bereiter/bereiter_co2.csv',delimiter=';',skiprows=15,unpack=True)
ax1.plot(t/1000,pco2,'-k',lw=0.5,zorder=10)
ax1.set_ylabel('\t $CO_2$ (ppmv)')
ax1.set_ylim([150,300])
ax1.set_yticks([180,220,260,300])
ax1.yaxis.set_minor_locator(FixedLocator(np.arange(180,300,20)))

# Plot dust
d,t,dust = np.loadtxt('data/sediment_cores/lambert/EDC_DustFlux_1cm.tab',delimiter='\t',skiprows=18,unpack=True)
ax2.fill_between(t[::10],dust[::10],y2=0,color='orange',zorder=5,lw=0)
ax2.plot(t,dust,c='darkgoldenrod',lw=0.1,zorder=5)
ax2.axhline(0,color='k',lw=0.5,zorder=10)
ax2.axvline(start,color='k',lw=0.5,zorder=10)
ax2.axvline(stop,color='k',lw=0.5,zorder=10)
ax2.set_ylabel('\t$\leftarrow$ EDC dust flux\n\t\t(mg $m^{-2}$ $y^{-1}$)',color='darkgoldenrod')
ax2.set_ylim([45,0])
ax2.set_yticks([0,10,20,30])
ax2.yaxis.set_minor_locator(FixedLocator(np.arange(0,30,5)))

# Plot reconstructed POC flux and PIC:POC export ratio
proxy = np.loadtxt('data/output/reconstructions/PS1768-8.txt',delimiter='\t',skiprows=1)
t = proxy[:,1]
poc_mean = np.nanmean(proxy[:,[2,4,6]],axis=1)
poc_min = np.nanmin(proxy[:,[2,4,6]],axis=1)
poc_max = np.nanmax(proxy[:,[2,4,6]],axis=1)
pic_poc_mean = np.nanmean(proxy[:,[3,5,7]],axis=1)
pic_poc_min = np.nanmin(proxy[:,[3,5,7]],axis=1)
pic_poc_max = np.nanmax(proxy[:,[3,5,7]],axis=1)

ax3.fill_between(t,poc_min,poc_max,color=lightblue,alpha=alpha,edgecolor='none')
ax3.plot(t,poc_mean,'-b',lw=0.5)
ax3.set_ylim([-0.15,0.48])
ax3.set_yticks([0,.1,.2,.3,.4])
ax3.set_ylabel('\tPOC flux (mol $m^{-2}$ $y^{-1}$)',color='b')
ax3.yaxis.set_minor_locator(FixedLocator(np.arange(0,0.4,0.05)))

ax5.fill_between(t,pic_poc_min,pic_poc_max,color=lightred,alpha=alpha,edgecolor='none')
ax5.plot(t,pic_poc_mean,'-r',lw=0.5)
ax5.set_ylim([0.1,1.3])
ax5.set_yticks([.2,.4,.6,.8,1])
ax5.set_ylabel('PIC:POC (mol:mol)        ',color='r')
ax5.yaxis.set_minor_locator(FixedLocator(np.arange(0.2,1,0.1)))

# Plot TOC and CaCO3 content
bulk = np.loadtxt('data/sediment_cores/bulk/PS1768-8_bulk_data.tab',delimiter='\t',skiprows=32)
d = bulk[:,0]
PIC = bulk[:,6]
TOC = bulk[:,7]
age_file = 'data/sediment_cores/zielinski/PS1768-8_age_model.tab'
d_cal,age_cal = np.loadtxt(age_file,skiprows=13,usecols=[0,1],unpack=True)
t = np.interp(d,d_cal,age_cal)+1

ax4.plot(t,TOC,'-k',lw=0.5)
ax4.set_ylabel('      TOC (%)')	
ax4.set_ylim([0,0.5])
ax4.set_yticks([.1,.2,.3,.4,.5])
ax4.yaxis.set_minor_locator(FixedLocator(np.arange(0.1,0.5,0.05)))

ax6.plot(t,PIC,'-k',lw=0.5)
ax6.set_ylabel('\t     CaCO$_3$ (%)')
ax6.set_ylim([-7,18])
ax6.set_yticks([0,5,10,15])
ax6.yaxis.set_minor_locator(FixedLocator(np.arange(0,15,2.5)))

# Plot IRD
d,ird = np.loadtxt('data/sediment_cores/bulk/PS1768-8_IRD.tab',delimiter='\t',skiprows=32,unpack=True)
t = np.interp(d,d_cal,age_cal)
ax7.fill_between(t,y1=0,y2=ird,fc=ice)
ax7.plot(t,ird,'-k',lw=0.2)
ax7.set_ylim([0,16])
ax7.set_yticks([0,5,10])
ax7.yaxis.set_minor_locator(FixedLocator(np.arange(0,10,2.5)))
ax7.set_ylabel('IRD (#/10 $cm^{3}$)')
ax7.axhline(0,color='k',lw=0.5,zorder=10)

# Pannels labels
ax_global.annotate('a',(0.03,0.97),xycoords='axes fraction',fontsize=8,fontweight='bold')
ax_global.annotate('b',(0.03,0.7),xycoords='axes fraction',fontsize=8,fontweight='bold')
ax_global.annotate('c',(0.03,0.41),xycoords='axes fraction',fontsize=8,fontweight='bold')
ax_global.annotate('d',(0.03,0.13),xycoords='axes fraction',fontsize=8,fontweight='bold')
ax_global.set_title('Core PS1768-8',fontweight='bold',fontsize=8)

# Save figure
plt.savefig('fig/fig_08_comparison_PS1768-8.pdf')
plt.savefig('fig/fig_08_comparison_PS1768-8.png',dpi=300)

#%% ===========================================================================
# Core 1786-1 (Atlantic) - comparison with 15N, %TOC and % CaCO3
# =============================================================================

# Declare figure and age range
plt.close('all')
start = 0
stop = 30
fig,ax_global = plt.subplots(1,1,figsize=(9/inch,17/inch))
ax_global.set_xlabel('Age (ka)')
ax_global.set_xlim([start,stop])
ax_global.patch.set_alpha(0)
ax_global.set_yticks([])
plt.subplots_adjust(left=0.16,right=0.82,bottom=0.07,top=0.95)

# Add axes
pos = ax_global.get_position()
ax_height = pos.height*0.3
ax1 = fig.add_axes([pos.x0,pos.y1-ax_height*0.9,pos.width,ax_height*0.9]) # pCO2
sub_axis(ax1,'left')
ax2 = fig.add_axes([pos.x0,pos.y1-ax_height*0.8,pos.width,ax_height*0.8]) # dust
sub_axis(ax2,'right')
ax3 = fig.add_axes([pos.x0,pos.y0+1.67*ax_height,pos.width,ax_height*0.9]) # POC export flux proxy
sub_axis(ax3,'left')
ax4 = fig.add_axes([pos.x0,pos.y0+1.67*ax_height,pos.width,ax_height*0.9]) # d15N sediment core
sub_axis(ax4,'right')
ax5 = fig.add_axes([pos.x0,pos.y0+0.35*ax_height,pos.width,ax_height*0.9]) # PIC:POC export proxy
sub_axis(ax5,'left')
ax6 = fig.add_axes([pos.x0,pos.y0+0.5*ax_height,pos.width,ax_height*0.5]) # PIC sediment core
sub_axis(ax6,'right')
ax7 = fig.add_axes([pos.x0,pos.y0,pos.width,ax_height*0.5]) # IRD
sub_axis(ax7,'left')
axtoc = fig.add_axes([pos.x0,pos.y0+1.15*ax_height,pos.width,ax_height*0.5]) # TOC sediment core
sub_axis(axtoc,'right')

axes = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,axtoc]
for i in range(len(axes)):
	# Tme range and remove x ticks
	axes[i].set_xlim(start,stop)
	axes[i].set_xticks([])
	axes[i].patch.set_alpha(0)
	
# Add MIS and ACR
ax_global.xaxis.set_minor_locator(AutoMinorLocator())
ax_global.axvspan(13,14.7, color='lightblue',zorder=-40) # ACR
ax_global.axvspan(14,29, color='lightgrey',zorder=-40) # MIS2
ax_global.axvspan(57,71, color='lightgrey',zorder=-40) # MIS4
ax_global.axvspan(130,190, color='lightgrey',zorder=-40) # MIS6

# Plot pCO2
t,pco2,sigma = np.loadtxt('data/sediment_cores/bereiter/bereiter_co2.csv',delimiter=';',skiprows=15,unpack=True)
ax1.plot(t/1000,pco2,'-k',lw=0.5,zorder=10)
ax1.set_ylabel('\t $CO_2$ (ppmv)')
ax1.set_ylim([150,300])
ax1.set_yticks([180,220,260,300])
ax1.yaxis.set_minor_locator(FixedLocator(np.arange(180,300,20)))

# Plot dust
d,t,dust = np.loadtxt('data/sediment_cores/lambert/EDC_DustFlux_1cm.tab',delimiter='\t',skiprows=18,unpack=True)
ax2.fill_between(t[::10],dust[::10],y2=0,color='orange',zorder=5,lw=0)
ax2.plot(t,dust,c='darkgoldenrod',lw=0.1,zorder=5)
ax2.plot(0,color='k',lw=0.5,zorder=10)
ax2.axvline(start,color='k',lw=0.5,zorder=10)
ax2.axvline(stop,color='k',lw=0.5,zorder=10)
ax2.set_ylabel('\t$\leftarrow$ EDC dust flux\n\t\t(mg $m^{-2}$ $y^{-1}$)',color='darkgoldenrod')
ax2.set_ylim([45,0])
ax2.set_yticks([0,10,20,30])
ax2.yaxis.set_minor_locator(FixedLocator(np.arange(0,30,5)))

# Plot reconstructed POC flux and PIC:POC export ratio
proxy = np.loadtxt('data/output/reconstructions/PS1786-1.txt',delimiter='\t',skiprows=1)
t = proxy[:,1]
poc_mean = np.nanmean(proxy[:,[2,4,6]],axis=1)
poc_min = np.nanmin(proxy[:,[2,4,6]],axis=1)
poc_max = np.nanmax(proxy[:,[2,4,6]],axis=1)
pic_poc_mean = np.nanmean(proxy[:,[3,5,7]],axis=1)
pic_poc_min = np.nanmin(proxy[:,[3,5,7]],axis=1)
pic_poc_max = np.nanmax(proxy[:,[3,5,7]],axis=1)

ax3.fill_between(t,poc_min,poc_max,color=lightblue,alpha=alpha,edgecolor='none')
ax3.plot(t,poc_mean,'-b',lw=0.5)
ax3.set_ylim([0,0.5])
ax3.set_yticks([0,.1,.2,.3,.4])
ax3.set_ylabel('POC flux (mol $m^{-2}$ $y^{-1}$)',color='b')
ax3.yaxis.set_minor_locator(FixedLocator(np.arange(0,0.4,0.05)))

ax5.fill_between(t,pic_poc_min,pic_poc_max,color=lightred,alpha=alpha,edgecolor='none')
ax5.plot(t,pic_poc_mean,'-r',lw=0.5)
ax5.set_ylim([0.1,1.4])
ax5.set_yticks([.4,.6,.8,1,1.2])
ax5.set_ylabel('   PIC:POC (mol:mol)',color='r')
ax5.yaxis.set_minor_locator(FixedLocator(np.arange(0.4,1.2,0.1)))

# Plot d15N from Jacot des Combes et al., 2008
d,t,d13,d15,cn = np.loadtxt('data/sediment_cores/bulk/PS1786-1_d13C_d15N_opal.tab',delimiter='\t',skiprows=17,unpack=True)
ax4.plot(t,d15,'-k',lw=0.5)
ax4.set_ylim([1,6])
ax4.set_yticks([2,3,4,5])
ax4.yaxis.set_minor_locator(FixedLocator(np.arange(2,5,0.5)))
ax4.set_ylabel('${\delta}^{15}$N diatom-bound\n(‰ vs. air) ')	

# Plot CaCO3 content
d_cal = np.array([0.0275,0.7325,1.23,1.825,2.35])
age_cal = np.array([1.895,5.581,9.534,18.327,27.5])
d,w,po,dg,de, caco3,toc,ts,ird= np.loadtxt('data/sediment_cores/bulk/PS1786-1_sedimentology.tab',delimiter='\t',skiprows=20,unpack=True)
t = np.interp(d,d_cal,age_cal)
ax6.plot(t,caco3,'-k',lw=0.5)
ax6.set_ylim([0,0.8])
ax6.set_yticks([0,.2,.4,.6,.8])
ax6.yaxis.set_minor_locator(FixedLocator(np.arange(0,.8,0.1)))
ax6.set_ylabel('CaCO$_3$ (%)')

# Plot TOC content
axtoc.plot(t[1:],toc[1:],'-k',lw=0.5)
axtoc.set_ylim([0.1,0.5])
axtoc.set_yticks([0.1,0.2,0.3,0.4,0.5])
axtoc.yaxis.set_minor_locator(FixedLocator(np.arange(0.1,0.5,0.05)))
axtoc.set_ylabel('TOC (%)')

# Plot IRD
ax7.fill_between(t,y1=0,y2=ird,fc=ice)
ax7.plot(t,ird,'-k',lw=0.5)
ax7.set_ylim([0,5])
ax7.set_yticks([0,2,4])
ax7.yaxis.set_minor_locator(FixedLocator(np.arange(0,4,1)))
ax7.set_ylabel('IRD (#/10 $cm^{3}$)     ')
ax7.axhline(0,color='k',lw=0.5,zorder=10)

# Pannels labels
ax_global.annotate('a',(0.03,0.975),xycoords='axes fraction',fontsize=8,fontweight='bold')
ax_global.annotate('b',(0.03,0.69),xycoords='axes fraction',fontsize=8,fontweight='bold')
ax_global.annotate('c',(0.03,0.47),xycoords='axes fraction',fontsize=8,fontweight='bold')
ax_global.annotate('d',(0.03,0.31),xycoords='axes fraction',fontsize=8,fontweight='bold')
ax_global.annotate('e',(0.03,0.1),xycoords='axes fraction',fontsize=8,fontweight='bold')
ax_global.set_title('Core PS1786-1',fontweight='bold',fontsize=8)

# Save figure
plt.savefig('fig/fig_09_comparison_PS1786-1.pdf')
plt.savefig('fig/fig_09_comparison_PS1786-1.png',dpi=300)

#%% ===========================================================================
# Core 2606-6 (Indian) - comparison with d15N and % CaCO3
# =============================================================================

# Declare figure and age range
plt.close('all')
start = 0
stop = 50
fig,ax_global = plt.subplots(1,1,figsize=(9/inch,15/inch))
ax_global.set_xlabel('Age (ka)')
ax_global.set_xlim([start,stop])
ax_global.patch.set_alpha(0)
ax_global.set_yticks([])
ax_global.set_zorder(-10)
plt.subplots_adjust(left=0.16,right=0.82,bottom=0.07,top=0.95)
	
# Add axes
pos = ax_global.get_position()
ax_height = pos.height*0.3
ax1 = fig.add_axes([pos.x0,pos.y1-ax_height,pos.width,ax_height]) # pCO2
sub_axis(ax1,'left')
ax2 = fig.add_axes([pos.x0,pos.y1-ax_height*0.8,pos.width,ax_height*0.8]) # dust
sub_axis(ax2,'right')
ax3 = fig.add_axes([pos.x0,pos.y0+1.4*ax_height,pos.width,ax_height*1.1]) # POC export flux proxy
sub_axis(ax3,'left')
ax4 = fig.add_axes([pos.x0,pos.y0+1.6*ax_height,pos.width,ax_height*1.1]) # diatombound d15N
sub_axis(ax4,'right')
ax5 = fig.add_axes([pos.x0,pos.y0+0.5*ax_height,pos.width,ax_height]) # PIC:POC export proxy
sub_axis(ax5,'left')
ax6 = fig.add_axes([pos.x0,pos.y0+0.5*ax_height,pos.width,ax_height*0.9]) # PIC sediment core
sub_axis(ax6,'right')
ax7 = fig.add_axes([pos.x0,pos.y0,pos.width,ax_height*0.5]) # IRD
sub_axis(ax7,'left')

axes = [ax1,ax2,ax3,ax4,ax5,ax6,ax7]
for i in range(len(axes)):
	# Tme range and remove x ticks
	axes[i].set_xlim(start,stop)
	axes[i].set_xticks([])
	axes[i].patch.set_alpha(0)

# Add MIS and ACR
ax_global.xaxis.set_minor_locator(AutoMinorLocator())
ax_global.axvspan(13,14.7, color='lightblue') # ACR
ax_global.axvspan(14,29, color='lightgrey') # MIS2
ax_global.axvspan(57,71, color='lightgrey') # MIS4
ax_global.axvspan(130,190, color='lightgrey') # MIS6

# Plot pCO2
t,pco2,sigma = np.loadtxt('data/sediment_cores/bereiter/bereiter_co2.csv',delimiter=';',skiprows=15,unpack=True)
ax1.plot(t/1000,pco2,'-k',lw=0.5,zorder=10)
ax1.set_ylabel('\t $CO_2$ (ppmv)')
ax1.set_ylim([150,300])
ax1.set_yticks([180,220,260,300])
ax1.yaxis.set_minor_locator(FixedLocator(np.arange(180,300,20)))

# Plot dust
d,t,dust = np.loadtxt('data/sediment_cores/lambert/EDC_DustFlux_1cm.tab',delimiter='\t',skiprows=18,unpack=True)
ax2.fill_between(t[::10],dust[::10],y2=0,color='orange',zorder=5,lw=0)
ax2.plot(t,dust,c='darkgoldenrod',lw=0.1,zorder=5)
ax2.axhline(0,color='k',lw=0.5,zorder=10)
ax2.axvline(start,color='k',lw=0.5,zorder=10)
ax2.axvline(stop,color='k',lw=0.5,zorder=10)
ax2.set_ylabel('\t$\leftarrow$ EDC dust flux\n\t\t(mg $m^{-2}$ $y^{-1}$)',color='darkgoldenrod')
ax2.set_ylim([45,0])
ax2.set_yticks([0,10,20,30])
ax2.yaxis.set_minor_locator(FixedLocator(np.arange(0,30,5)))

# Plot reconstructed POC flux and PIC:POC export ratio
proxy = np.loadtxt('data/output/reconstructions/PS2606-6.txt',delimiter='\t',skiprows=1)
t = proxy[:,1]
poc_mean = np.nanmean(proxy[:,[2,4,6]],axis=1)
poc_min = np.nanmin(proxy[:,[2,4,6]],axis=1)
poc_max = np.nanmax(proxy[:,[2,4,6]],axis=1)
pic_poc_mean = np.nanmean(proxy[:,[3,5,7]],axis=1)
pic_poc_min = np.nanmin(proxy[:,[3,5,7]],axis=1)
pic_poc_max = np.nanmax(proxy[:,[3,5,7]],axis=1)

ax3.fill_between(t,poc_min,poc_max,color=lightblue,alpha=alpha,edgecolor='none')
ax3.plot(t,poc_mean,'-b',lw=0.5)
ax3.set_ylim([-0.15,0.5])
ax3.set_yticks([0,.1,.2,.3,.4,])
ax3.set_ylabel('\tPOC flux (mol $m^{-2}$ $y^{-1}$)',color='b')
ax3.yaxis.set_minor_locator(FixedLocator(np.arange(0,0.4,0.05)))

ax5.fill_between(t,pic_poc_min,pic_poc_max,color=lightred,alpha=alpha,edgecolor='none')
ax5.plot(t,pic_poc_mean,'-r',lw=0.5)
ax5.set_ylim([0.1,1.4])
ax5.set_yticks([.2,.4,.6,.8,1,1.2])
ax5.set_ylabel('PIC:POC (mol:mol)',color='r')
ax5.yaxis.set_minor_locator(FixedLocator(np.arange(0.2,1.2,0.1)))

# Plot diatom-bound 15N (Jacot Des Combres, 2008)
# Extended Age model from Civel-Mazens et al. 2024
age_file = 'data/sediment_cores/esper/PS2606-6_SSTd.tab'
d_cal,age_cal = np.loadtxt(age_file,skiprows=18,usecols=[0,1],unpack=True)
d,t,d13,d15,c_n = np.loadtxt('data/sediment_cores/bulk/PS2606-6_d13C_d15N_opal.tab',delimiter='\t',skiprows=17,unpack=True)
t = np.interp(d,d_cal,age_cal)
ax4.plot(t,d15,'-k',lw=0.5)
ax4.set_ylim([0,6])
ax4.set_yticks([0,1,2,3,4])
ax4.yaxis.set_minor_locator(FixedLocator(np.arange(0,4,0.5)))
ax4.set_ylabel('${\delta}^{15}$N diatom-bound\t\t\t\n(‰ vs. air)              ')

# Plot CaCO3 from Kuhn, Gerhard (2011): Sedimentology on various cores from the South Atlantic
data  = np.loadtxt('data/sediment_cores/bulk/PS2606-6_sedimentology.tab',delimiter='\t',skiprows=21)
d = data[:,0]
bsi = data[:,2]
qz = data[:,3]
caco3 = data[:,5]
t = np.interp(d,d_cal,age_cal)
sel = (t<np.max(age_cal))
t = np.interp(d,d_cal,age_cal)

ax6.plot(t[sel],caco3[sel],'-k',lw=0.5)
ax6.set_ylim([-2,13])
ax6.set_yticks([0,4,8,12])
ax6.set_ylabel('CaCO$_3$ (%)')
ax6.yaxis.set_minor_locator(FixedLocator(np.arange(0,12,2)))

# Plot IRD
d,ird = np.loadtxt('data/sediment_cores/bulk/PS2606-6_IRD.tab',delimiter='\t',skiprows=32,unpack=True)
t = np.interp(d,d_cal,age_cal)
ax7.fill_between(t,y1=0,y2=ird,fc=ice)
ax7.plot(t,ird,'-k',lw=0.2)
ax7.set_ylim([0,6])	
ax7.set_yticks([0,2,4])
ax7.yaxis.set_minor_locator(FixedLocator(np.arange(0,4,1)))
ax7.set_ylabel('IRD (#/10 $cm^{3}$)        ')
ax7.axhline(0,color='k',lw=0.5,zorder=10)

# Pannels labels
ax_global.annotate('a',(0.03,0.97),xycoords='axes fraction',fontsize=8,fontweight='bold')
ax_global.annotate('b',(0.03,0.68),xycoords='axes fraction',fontsize=8,fontweight='bold')
ax_global.annotate('c',(0.03,0.4),xycoords='axes fraction',fontsize=8,fontweight='bold')
ax_global.annotate('d',(0.03,0.08),xycoords='axes fraction',fontsize=8,fontweight='bold')
ax_global.set_title('Core PS2606-6',fontweight='bold',fontsize=8)

# Save figure
plt.savefig('fig/fig_10_comparison_PS2606-6.pdf')
plt.savefig('fig/fig_10_comparison_PS2606-6.png',dpi=300)

#%% ===========================================================================
# Core PS75/072-4 (Pacific) - comparison with d15N and Ba/Fe
# =============================================================================

# Declare figure and age range
plt.close('all')
start = 0
stop = 200
fig,ax_global = plt.subplots(1,1,figsize=(15/inch,15/inch))
ax_global.set_xlabel('Age (ka)')
ax_global.set_xlim([start,stop])
ax_global.patch.set_alpha(0)
ax_global.set_yticks([])
ax_global.set_zorder(-10)
plt.subplots_adjust(left=0.1,right=0.88,bottom=0.07,top=0.95)

# Add axes
pos = ax_global.get_position()
ax_height = pos.height*0.3
ax1 = fig.add_axes([pos.x0,pos.y1-ax_height,pos.width,ax_height]) # pCO2
sub_axis(ax1,'left')
ax2 = fig.add_axes([pos.x0,pos.y1-ax_height*0.8,pos.width,ax_height*0.8]) # dust
sub_axis(ax2,'right')
ax3 = fig.add_axes([pos.x0,pos.y0+1.4*ax_height,pos.width,ax_height*1.1]) # POC export flux proxy
sub_axis(ax3,'left')
ax4 = fig.add_axes([pos.x0,pos.y0+1.4*ax_height,pos.width,ax_height*1.1]) # TOC sediment core
sub_axis(ax4,'right')
ax5 = fig.add_axes([pos.x0,pos.y0+0.5*ax_height,pos.width,ax_height]) # PIC:POC export proxy
sub_axis(ax5,'left')
ax6 = fig.add_axes([pos.x0,pos.y0+0.4*ax_height,pos.width,ax_height*1.1]) # Ba/Fe sediment core
sub_axis(ax6,'right')
ax7 = fig.add_axes([pos.x0,pos.y0,pos.width,ax_height*0.5]) # IRD
sub_axis(ax7,'left')

axes = [ax1,ax2,ax3,ax4,ax5,ax6,ax7]
for i in range(len(axes)):
	# Tme range and remove x ticks
	axes[i].set_xlim(start,stop)
	axes[i].set_xticks([])
	axes[i].patch.set_alpha(0)
	
# Add MIS and ACR
ax_global.xaxis.set_minor_locator(AutoMinorLocator())
ax_global.axvspan(13,14.7, color='lightblue') # ACR
ax_global.axvspan(14,29, color='lightgrey') # MIS2
ax_global.axvspan(57,71, color='lightgrey') # MIS4
ax_global.axvspan(130,190, color='lightgrey') # MIS6

# Plot pCO2
t,pco2,sigma = np.loadtxt('data/sediment_cores/bereiter/bereiter_co2.csv',delimiter=';',skiprows=15,unpack=True)
ax1.plot(t/1000,pco2,'-k',lw=0.5,zorder=10)
ax1.set_ylabel('\t $CO_2$ (ppmv)')
ax1.set_ylim([150,300])
ax1.set_yticks([180,220,260,300])
ax1.yaxis.set_minor_locator(FixedLocator(np.arange(180,300,20)))

# Plot dust
d,t,dust = np.loadtxt('data/sediment_cores/lambert/EDC_DustFlux_1cm.tab',delimiter='\t',skiprows=18,unpack=True)
ax2.fill_between(t[::10],dust[::10],y2=0,color='orange',zorder=5,lw=0)
ax2.plot(t,dust,c='darkgoldenrod',lw=0.1,zorder=5)
ax2.axhline(0,color='k',lw=0.5,zorder=10)
ax2.axvline(start,color='k',lw=0.5,zorder=10)
ax2.axvline(stop,color='k',lw=0.5,zorder=10)
ax2.set_ylabel('\t$\leftarrow$ EDC dust flux\n\t\t(mg $m^{-2}$ $y^{-1}$)',color='darkgoldenrod')
ax2.set_ylim([45,0])
ax2.set_yticks([0,10,20,30])
ax2.yaxis.set_minor_locator(FixedLocator(np.arange(0,30,5)))

# Plot reconstructed POC flux and PIC:POC export ratio
proxy = np.loadtxt('data/output/reconstructions/PS75_072-4.txt',delimiter='\t',skiprows=1)
t = proxy[:,1]
poc_mean = np.nanmean(proxy[:,[2,4,6]],axis=1)
poc_min = np.nanmin(proxy[:,[2,4,6]],axis=1)
poc_max = np.nanmax(proxy[:,[2,4,6]],axis=1)
pic_poc_mean = np.nanmean(proxy[:,[3,5,7]],axis=1)
pic_poc_min = np.nanmin(proxy[:,[3,5,7]],axis=1)
pic_poc_max = np.nanmax(proxy[:,[3,5,7]],axis=1)

ax3.fill_between(t,poc_min,poc_max,color=lightblue,alpha=alpha,edgecolor='none')
ax3.plot(t,poc_mean,'-b',lw=0.5)
ax3.set_ylim([-0.15,0.65])
ax3.set_yticks([0,.1,.2,.3,.4,.5])
ax3.set_ylabel('POC flux (mol $m^{-2}$ $y^{-1}$)',color='b')
ax3.yaxis.set_minor_locator(FixedLocator(np.arange(0,0.5,0.05)))

ax5.fill_between(t,pic_poc_min,pic_poc_max,color=lightred,alpha=alpha,edgecolor='none')
ax5.plot(t,pic_poc_mean,'-r',lw=0.5)
ax5.set_ylim([0.1,1.4])
ax5.set_yticks([.2,.4,.6,.8,1,1.2])
ax5.set_ylabel('PIC:POC (mol:mol)',color='r')
ax5.yaxis.set_minor_locator(FixedLocator(np.arange(0.2,1.2,0.1)))

# Plot diatom-bound 15N (Studer et al., 2015)
d,t,d15,d15p,d15c,d15a = np.loadtxt('data/sediment_cores/bulk/PS75_072-4_d15N.tab',delimiter='\t',skiprows=18,unpack=True)

ax4.plot(t,d15,'-k',mfc=lightblue,ms=3,lw=1,label='Total')
ax4.plot(t,d15c,'--k',ms=3,lw=0.5,label='Centric')
ax4.plot(t,d15p,'-k',ms=3,lw=0.5,label='Pennate')
ax4.set_ylim([1,8])
ax4.set_yticks([2,3,4,5,6,7])
ax4.yaxis.set_minor_locator(FixedLocator(np.arange(2,7,0.5)))
ax4.legend(frameon=False,loc='lower right',handletextpad=0.5,handlelength=1.2,fontsize=6,
		   bbox_to_anchor=[0.96, 0.1])
ax4.set_ylabel('${\delta}^{15}$N diatom-bound\n (‰ vs. air)')

# Plot Ba/Fe export proxy
d,t,bafe = np.loadtxt('data/sediment_cores/bulk/PS75_072-4_Ba-Fe.tab',delimiter='\t',skiprows=15,unpack=True)

ax6.plot(t,bafe,'-k',ms=3,lw=0.5)
ax6.set_ylim([-0.1,0.8])
ax6.set_yticks([0,.2,.4,.6])
ax6.yaxis.set_minor_locator(FixedLocator(np.arange(0,0.6,0.1)))
ax6.set_ylabel('Ba/Fe')

# Plot IRD from two other cores in the Pacific sector
# IRD in core RS15-LC48 (Somov sea) from McKay et al. (2022)
d,ird = np.loadtxt('data/sediment_cores/mckay/RS15-LC48_IRD_mckay2022.csv',delimiter='\t',skiprows=1,unpack=True)
ird = ird*10/16 # convert #/16cm 3 to #/10cm3
t=d*100/0.8+8 #Age model for core RS15-LC48 from Ohneiser et al., 2019
ax7.fill_between(t,y1=0,y2=ird,fc=ice,lw=0.2)
ax7.plot(t,ird,'--k',lw=0.5)

# IRD in core V1439 (Campbell plateau) from Carter et al. (2002)
t,ird = np.loadtxt('data/sediment_cores/carter/V1439_IRD_carter2002.csv',delimiter='\t',skiprows=1,unpack=True)
ax7.fill_between(t,y1=0,y2=ird,fc=ice)
ax7.plot(t,ird,'-k',lw=0.5)

ax7.set_ylim([0,12])	
ax7.set_yticks([0,5,10])
ax7.yaxis.set_minor_locator(FixedLocator(np.arange(0,10,2.5)))
ax7.set_ylabel('IRD (#/10 $cm^{3}$)  ')
ax7.axhline(0,color='k',lw=0.5,zorder=10)
ax7.axvline(200,color='k',lw=0.5,zorder=10)

ax_global.annotate('RS15-LC48\nSomov sea',(0.76,0.09),xycoords='axes fraction',fontsize=6)
ax_global.annotate('V1439\nCampbell plateau',(0.16,0.09),xycoords='axes fraction',fontsize=6)	

# Pannels labels
ax_global.annotate('a',(0.01,0.97),xycoords='axes fraction',fontsize=8,fontweight='bold')
ax_global.annotate('b',(0.01,0.67),xycoords='axes fraction',fontsize=8,fontweight='bold')
ax_global.annotate('c',(0.01,0.39),xycoords='axes fraction',fontsize=8,fontweight='bold')
ax_global.annotate('d',(0.01,0.11),xycoords='axes fraction',fontsize=8,fontweight='bold')

ax_global.set_title('Core PS75/072-4',fontweight='bold',fontsize=8)

# Save figure
plt.savefig('fig/fig_11_comparison_PS75_072-4.pdf')
plt.savefig('fig/fig_11_comparison_PS75_072-4.png',dpi=300)

