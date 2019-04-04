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

latlim = [30,85]
lonlim= [-160,-60]
date = datetime(2017, 8, 21, 0)

fig = gm.plotCartoMap(projection='stereo', title='Geomagnetic coordinates: MLAT/MLT',
                      latlim=latlim, lonlim=lonlim, 
                      parallels = [20, 40, 60, 80, 90],
                      meridians = [-220, -180, -160,-140,-120,-100, -80,-60, -40, 0],
                      grid_linewidth=1,
                      figure=True,
                      states=False)
A = ap.Apex(date=date)

glon = np.arange(lonlim[0]-40, lonlim[1] + 40.1, 0.1)
glat = np.arange(latlim[0], latlim[1] + 0.1, 0.11)

longrid, latgrid = np.meshgrid(glon, glat)

mlat_levels = np.arange(-90, 90.1, 5)
# mlon
mlat, mlon = A.convert(latgrid, longrid, 'geo', 'apex')
mlon_levels = np.arange(-180,180,20)
# mlt
mlat, mlon = A.convert(latgrid, longrid, 'geo', 'mlt', datetime=date)
mlon_levels = np.arange(0,24.1,1)
ay = plt.contour(glon,glat, mlat, levels = mlat_levels, colors='red', transform=ccrs.PlateCarree())
ax = plt.contour(glon,glat, mlon, levels = mlon_levels, colors='blue', linestyles ='solid',  transform=ccrs.PlateCarree())
ax.clabel(inline=True, fmt = '%d', fontsize=12, colors='blue')
ay.clabel(inline=True, fmt = '%d', fontsize=12, colors='red')