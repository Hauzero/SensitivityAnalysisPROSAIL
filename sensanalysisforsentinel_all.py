# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 18:46:26 2017

@author: hauserlt
"""

####################### LOADING LIBRARIES ###################
#############################################################

import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
#import prosail
import random
random.seed()
import pyprosail

################### IMPORT MODES AND SPECTRAL BANDS SENTINEL ########################
#####################################################################################

traitrange = pd.read_csv('/media/leon/FREECOM HDD/W-Documents/Sensitivity Analysis/Prosailranges.csv')
spectralsensitivityfile = '/media/leon/FREECOM HDD/Data/Scaling-Borneo/Forwardmodelling/S2-SRF_COPE-GSEG-EOPG-TN-15-0007_3.0.csv'
s2sens = pd.read_csv(spectralsensitivityfile)
traitslist = ['N','Cab','Car','Cbp','Cw','Cm','Psoil','LAI','hot','solar_zenith','solar_azimuth','view_zenith','view_azimuth','ALA']
straitslist = ['N', 'LAI', 'ALA','Cab','Car','Cm','Cw','Cbp','Psoil','hot']             
bandslist = ['B2','B3','B4', 'B5', 'B6', 'B7', 'B8A', 'B11', 'B12']
s2min = np.array( [ 439, 537, 646, 694, 731, 768, 848, 1539, 2072])
s2max = np.array( [ 535, 582, 685, 714, 749, 796, 881, 1681, 2312] )

#pyprosail.run(N, chloro, caroten, brown, EWT, LMA, psoil, LAI, hot_spot, solar_zenith, solar_azimuth, view_zenith, view_azimuth, LIDF)
Nth               = traitrange[traitrange['Symbol'] == 'N']
LAI               = traitrange[traitrange['Symbol'] == 'LAI']
ALA               = traitrange[traitrange['Symbol'] == 'ALA']
CAB               = traitrange[traitrange['Symbol'] == 'Cab']
CAR               = traitrange[traitrange['Symbol'] == 'Car']
CDM               = traitrange[traitrange['Symbol'] == 'Cm']
CW                = traitrange[traitrange['Symbol'] == 'Cw']
CBP               = traitrange[traitrange['Symbol'] == 'Cbp']
Bsoil             = traitrange[traitrange['Symbol'] == 'Bsoil']
Psoil             = traitrange[traitrange['Symbol'] == 'Psoil']
HOT             = traitrange[traitrange['Symbol'] == 'hot']
solar_zenith    = 26.945400199999998
solar_azimuth   = 50.350898700000002
view_zenith     = 3.9065937999999996
view_azimuth    = 302.40621949999996

x = float(Nth['Fix']), float(CAB['Fix']), float(CAR['Fix']), float(CBP['Fix']), float(CW['Fix']), float(CDM['Fix']), float(Psoil['Fix']), float(LAI['Fix']), float(HOT['Fix']), solar_zenith, solar_azimuth, view_zenith, view_azimuth, float(ALA['Fix'])
x = list(x)
################### RUN MODAL CURVATURE ########################
################################################################

N, chloro, caroten, brown, EWT, LMA, psoil, lai, hot_spot, solar_zenith, solar_azimuth, view_zenith, view_azimuth, LIDF = x
pro_out = pyprosail.run(N, chloro, caroten, brown, EWT, LMA, psoil, lai, hot_spot, solar_zenith, solar_azimuth, view_zenith, view_azimuth, LIDF)
s2out    = pd.DataFrame(columns={'B2','B3','B4', 'B5', 'B6', 'B7', 'B8A', 'B11', 'B12'})
s2outt   = np.array([[490,560,665,705,740,783,865,1610,2190],np.zeros(9)])
for ib, band in enumerate(bandslist):
    weigthedvalue = 0
    bandname = 'S2A_SR_AV_' + band
    sumvalue        = np.sum(s2sens[bandname])
    for ix in np.where(s2sens[bandname] != 0)[0]:
        iwl = np.where((pro_out[:,0]*1000) == int(s2sens['SR_WL'].iloc[ix]))[0][0]
        weigthedvalue   = weigthedvalue + (float(s2sens[bandname].iloc[ix]) * pro_out[iwl,1])
    s2out['band'] = (weigthedvalue/sumvalue)     
    s2outt[1,ib] =  (weigthedvalue/sumvalue)   
plt.plot(s2outt[0], s2outt[1])
avgline = s2outt

################### RUN MIN&MAX CURVES ########################
###################################################################
minlines = np.empty((len(traitslist),9))
maxlines = np.empty((len(traitslist),9))

for iv, variable in enumerate(traitslist):
        if variable in straitslist:           
            x = float(Nth['Fix']), float(CAB['Fix']), float(CAR['Fix']), float(CBP['Fix']), float(CW['Fix']), float(CDM['Fix']), float(Psoil['Fix']), float(LAI['Fix']), float(HOT['Fix']), solar_zenith, solar_azimuth, view_zenith, view_azimuth, float(ALA['Fix'])
            x = list(x)
            minvar = float(traitrange[traitrange['Symbol'] == variable]['Minimum'])
            maxvar = float(traitrange[traitrange['Symbol'] == variable]['Maximum'])
            x[iv]  = minvar
            
            s2out    = pd.DataFrame(columns={'B2','B3','B4', 'B5', 'B6', 'B7', 'B8A', 'B11', 'B12'})
            N, chloro, caroten, brown, EWT, LMA, psoil, lai, hot_spot, solar_zenith, solar_azimuth, view_zenith, view_azimuth, LIDF = x
            pro_out = pyprosail.run(N, chloro, caroten, brown, EWT, LMA, psoil, lai, hot_spot, solar_zenith, solar_azimuth, view_zenith, view_azimuth, LIDF)
            s2outt   = np.array([[490,560,665,705,740,783,865,1610,2190],np.zeros(9)])
            for ib, band in enumerate(bandslist):
                weigthedvalue = 0
                bandname = 'S2A_SR_AV_' + band
                sumvalue        = np.sum(s2sens[bandname])
                for ix in np.where(s2sens[bandname] != 0)[0]:
                    iwl = np.where((pro_out[:,0]*1000) == int(s2sens['SR_WL'].iloc[ix]))[0][0]
                    weigthedvalue   = weigthedvalue + (float(s2sens[bandname].iloc[ix]) * pro_out[iwl,1])
                s2out['band'] = (weigthedvalue/sumvalue)     
                s2outt[1,ib] =  (weigthedvalue/sumvalue)   
            plt.plot(s2outt[0], s2outt[1])
            minlines[iv] = s2outt[1]
            
            x[iv]  = maxvar
            s2out    = pd.DataFrame(columns={'B2','B3','B4', 'B5', 'B6', 'B7', 'B8A', 'B11', 'B12'})
            N, chloro, caroten, brown, EWT, LMA, psoil, lai, hot_spot, solar_zenith, solar_azimuth, view_zenith, view_azimuth, LIDF = x
            pro_out = pyprosail.run(N, chloro, caroten, brown, EWT, LMA, psoil, lai, hot_spot, solar_zenith, solar_azimuth, view_zenith, view_azimuth, LIDF)
            s2outt   = np.array([[490,560,665,705,740,783,865,1610,2190],np.zeros(9)])
            for ib, band in enumerate(bandslist):
                weigthedvalue = 0
                bandname = 'S2A_SR_AV_' + band
                sumvalue        = np.sum(s2sens[bandname])
                for ix in np.where(s2sens[bandname] != 0)[0]:
                    iwl = np.where((pro_out[:,0]*1000) == int(s2sens['SR_WL'].iloc[ix]))[0][0]
                    weigthedvalue   = weigthedvalue + (float(s2sens[bandname].iloc[ix]) * pro_out[iwl,1])
                s2out['band'] = (weigthedvalue/sumvalue)     
                s2outt[1,ib] =  (weigthedvalue/sumvalue)   
            plt.plot(s2outt[0], s2outt[1])
            print str(variable)
            maxlines[iv] = s2outt[1]          

################## SET & RUN ITERATIONS ########################
#################################################################

iterat      = 2000
total = pd.DataFrame(columns=('N', 'LAI', 'ALA','Cab','Car','Cm','Cw','Cbp','Psoil','hot', 'B2','B3','B4', 'B5', 'B6', 'B7', 'B8A', 'B11', 'B12'))
totalnp = np.empty((iterat,19))
totalnp[:] = np.nan

row = 0
while len(total) < iterat:
        N           = np.random.uniform(Nth['Minimum'],Nth['Maximum'])
        chloro      = np.random.uniform(CAB['Minimum'],CAB['Maximum'])
        EWT         = np.random.uniform(CW['Minimum'],CW['Maximum'])
        LMA         = np.random.uniform(CDM['Minimum'],CDM['Maximum'])
        brown       = np.random.uniform(CBP['Minimum'],CBP['Maximum'])
        caroten     = np.random.uniform(CAR['Minimum'],CAR['Maximum'])
        psoil       = np.random.uniform(Psoil['Minimum'],Psoil['Maximum'])
        lai         = np.random.uniform(LAI['Minimum'],LAI['Maximum'])
        hot_spot    = np.random.uniform(HOT['Minimum'],HOT['Maximum'])
        LIDF        = np.random.uniform(ALA['Minimum'],ALA['Maximum'])

        pro_out = pyprosail.run(N, chloro, caroten, brown, EWT, LMA, psoil, lai, hot_spot, solar_zenith, solar_azimuth, view_zenith, view_azimuth, LIDF)
        s2out    = pd.DataFrame(columns={'B2','B3','B4', 'B5', 'B6', 'B7', 'B8A', 'B11', 'B12'})
        s2outt   = np.array([[520,560,665,705,740,783,865,1610,2190],np.zeros(9)])
        for ib, band in enumerate(bandslist):
            weigthedvalue = 0
            bandname = 'S2A_SR_AV_' + band
            sumvalue        = np.sum(s2sens[bandname])
            for ix in np.where(s2sens[bandname] != 0)[0]:
                iwl = np.where((pro_out[:,0]*1000) == int(s2sens['SR_WL'].iloc[ix]))[0][0]
                weigthedvalue   = weigthedvalue + (float(s2sens[bandname].iloc[ix]) * pro_out[iwl,1])
            s2out['band'] = (weigthedvalue/sumvalue)     
            s2outt[1,ib] =  (weigthedvalue/sumvalue)
        plt.plot(s2outt[0], s2outt[1])
        spectra = pd.DataFrame([s2outt[1]], columns=('B2','B3','B4', 'B5', 'B6', 'B7', 'B8A', 'B11', 'B12'))
        traits  = pd.DataFrame([np.array([N, chloro, caroten, brown, EWT, LMA, psoil, lai, hot_spot, LIDF])], columns=('N', 'Cab','Car','Cbp','Cw','Cm','Psoil','LAI', 'hot', 'ALA'))
        total = total.append(pd.concat((traits,spectra),axis=1), ignore_index=True)
        totalnp[row,10:19] = s2outt[1]
        totalnp[row,0:10] = N, chloro, caroten, brown, EWT, LMA, psoil, lai, hot_spot, LIDF
        row = row +1

pvalue = 0.01
xiv = 0
fig, ax = plt.subplots(5, 2, figsize=(20,10), sharey='col', sharex='all')
fig2, ax2 = plt.subplots(5, 2, figsize=(20,10), sharey='col', sharex='all')

for iv, variable in enumerate(traitslist):
        response = np.empty((1, 2*len(bandslist)))
        if variable in straitslist:  
            if xiv > 4:
                riv = 1
            else:
                riv = 0
            if xiv > 4:
                civ = (xiv-6)
            else:
                civ = xiv
            xiv = xiv + 1
            bothlines = np.vstack((minlines[iv], maxlines[iv]))
            bothlines = np.sort(bothlines, axis=0)
            for ib, band in enumerate(bandslist):      
                response[0,ib], response[0,ib+len(bandslist)] = ((scipy.stats.pearsonr(total[variable], total[band])[0]**2)**0.5), scipy.stats.pearsonr(total[variable], total[band])[1]
            ax[civ,riv].plot(s2outt[0], response[0,0:9], marker='D', markersize=5, linestyle='--')
            ax[civ,riv].title.set_text(str(variable))
            ax[civ,riv].set_ylim(0,0.8)
            ax[civ,riv].set_ylabel('Pearson`s R' )

            for ib, band in enumerate(bandslist):
                if response[0,(9+ib)] < pvalue:
                    ax[civ,riv].axvspan(s2min[ib], s2max[ib], color='k',alpha=0.15)
            #        ax[0,1].fill_betweenx(1, s2min[ib], s2max[ib])    
                    
            #ax[0,0].plot(s2outt[0], laiminline[1], marker='D', markersize=5, linestyle='--')
            #ax[0,0].plot(s2outt[0], laimaxline[1], marker='D', markersize=5, linestyle='--')
            ax2[civ,riv].fill_between(s2outt[0], bothlines[0], bothlines[1], alpha=0.2, interpolate=True)
            ax2[civ,riv].plot(s2outt[0], avgline[1], marker='D', markersize=5, linestyle='-')
            ax2[civ,riv].title.set_text(str(variable))
            ax2[civ,riv].set_ylim(0,0.4)
            ax2[civ,riv].set_ylabel('Reflectance')