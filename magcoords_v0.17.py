#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 14:08:42 2019

@author: smrak
"""
import numpy as np
from cartomap import geogmap as gm
from datetime import datetime
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

import apexpy as ap

latlim = [-0,60]
lonlim= [-140,0]
date = datetime(2017, 8, 21, 6)

fig = gm.plotCartoMap(projection='plate', title='Geomagnetic coordinates: MLAT/MLT',
                      latlim=latlim, lonlim=lonlim, 
                      parallels = [0,10,20, 40, 60, 80, 90],
                      meridians = [-220, -180, -160,-140,-120,-100, -80,-60, -40, 0],
                      grid_linewidth=1,
                      figure=True,
                      states=False)
A = ap.Apex(date=date)

#glon = np.arange(lonlim[0]-40, lonlim[1] + 40.1, 1)
#glat = np.arange(latlim[0], latlim[1] + 0.1, 1)

#longrid, latgrid = np.meshgrid(glon, glat)

mlat_levels = np.arange(-90, 90.1, 10)
#mlat_levels = np.array([40,50,60,70])
# mlon
#mlat, mlon = A.convert(latgrid, longrid, 'geo', 'apex')
#mlon_levels = np.arange(-180,180,20)
# mlt
#mlat, mlon = A.convert(latgrid, longrid, 'geo', 'mlt', datetime=date)
mlon_levels = np.arange(0,24.2,2)

#ay = plt.contour(glon,glat, mlat, levels = mlat_levels, colors='red', transform=ccrs.PlateCarree())
#ax = plt.contour(glon,glat, mlon, levels = mlon_levels, colors='blue', linestyles ='solid',  transform=ccrs.PlateCarree())
#ax.clabel(inline=True, fmt = '%d', fontsize=12, colors='blue')
#ay.clabel(inline=True, fmt = '%d', fontsize=12, colors='red')

# MLATS
mlat_range = np.arange(mlat_levels[0], mlat_levels[-1]+0.1, 0.1)
mlon_range = np.arange(mlon_levels[0], 24.3, 0.1)
for mlon in mlon_levels:
    MLON = mlon * np.ones(mlat_range.size)
    y, x = A.convert(mlat_range,MLON, 'mlt', 'geo', datetime=date)
    if int(mlon) == 0:# or int(mlon) == 2:
        continue
    inmap = np.logical_and(x >= lonlim[0], x <= lonlim[1])
    if np.sum(inmap) > 10:
        plt.plot(np.unwrap(x,180), np.unwrap(y,90), 'b', lw=2, transform=ccrs.PlateCarree())
        ix = abs(y-np.mean(latlim)).argmin()
        mx = x[ix]-4
        my = np.mean(latlim)
        if np.logical_and(mx >= lonlim[0], mx <= lonlim[1]) and int(mlon) is not 0:
            plt.text(mx, my, str(int(mlon)), color='k', 
                     fontsize=14, backgroundcolor='white',transform=ccrs.PlateCarree())
for mlat in mlat_levels:
    MLAT = mlat * np.ones(mlon_range.size)
    gy,gx = A.convert(MLAT, mlon_range, 'mlt', 'geo', datetime=date)
    inmap = np.logical_and(gy >= latlim[0], gy <= latlim[1])
    if np.sum(inmap) > 10:
        plt.plot(np.unwrap(gx, 180), np.unwrap(gy, 90), 'b', transform=ccrs.PlateCarree())
        ix = abs(gx-np.mean(lonlim)).argmin()
        mx = np.mean(lonlim)
        my = gy[ix]-0.5
    if np.logical_and(mx >= lonlim[0], mx <= lonlim[1]) and \
    np.logical_and(my >= latlim[0], my <= latlim[1]):
        ix = abs(gx-np.mean(lonlim)).argmin()
        plt.text(mx, my, str(int(mlat)), color='k', 
             fontsize=14, backgroundcolor='white',transform=ccrs.PlateCarree())
        
#Functional
fig = gm.plotCartoMap(projection='plate', 
                      title='Geomagnetic coordinates: MLAT/MLT',
                      latlim=latlim, lonlim=lonlim, 
                      date=date,
                      #parallels = [0,10,20, 40, 60, 80, 90],
                      #meridians = [-220, -180, -160,-140,-120,-100, -80,-60, -40, 0],
                      grid_linewidth = 1,
                      figure = True,
                      states = False,
                      geomag = True,
                      gmagtype = 'apex',
                      mlon_cs = 'mlt',
                      mlon_levels = mlon_levels,
                      mlat_levels = mlat_levels,
                      mlon_colors='k',
                      mlat_colors='k',
                      mlat_labels=False)