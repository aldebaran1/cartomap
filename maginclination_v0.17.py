#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 19:06:18 2019

@author: smrak
"""
from datetime import datetime
import numpy as np
from cartomap import geogmap as gm
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

import igrf12

latlim = [30,85]
lonlim= [-160,-60]
date = datetime(2017, 8, 21, 11)

fig = gm.plotCartoMap(projection='stereo', title='Magnetic declination angle. Ground (solid), 150 km (dashed)',
                      latlim=latlim, lonlim=lonlim, 
                      parallels = [20, 40, 60, 80, 90],
                      meridians = [-220, -180, -160,-140,-120,-100, -80,-60, -40, 0],
                      grid_linewidth=1,
                      figure=True,
                      states=False)

glon = np.arange(0, 360, 1)
glat = np.arange(-90, 90.1, 1)


longrid, latgrid = np.meshgrid(glon, glat)
mag0 = igrf12.gridigrf12(t=date, glat=latgrid, glon=longrid, alt_km=0.0)
mag150 = igrf12.gridigrf12(t=date, glat=latgrid, glon=longrid, alt_km=150.0)

# Declination
#ai = plt.contour(longrid, latgrid, mag0.decl.values, levels=np.arange(-30,30.1,1), 
#                 cmap='bwr', transform=ccrs.PlateCarree())
#ai = plt.contour(longrid, latgrid, mag150.decl.values, levels=np.arange(-30,30.1,1), 
#                 cmap='bwr', transform=ccrs.PlateCarree(), linestyles='dashed')
#plt.contour(longrid, latgrid, mag150.decl.values, levels=np.arange(0,1,1), 
#                 colors='green', transform=ccrs.PlateCarree(), linestyles='solid')
#ai.clabel(inline=True, fmt = '%d', fontsize=12, colors='b')
#Inclination

ai = plt.contour(longrid, latgrid, mag0.incl.values, levels=np.arange(-90,90.1,1), 
                 colors='b', transform=ccrs.PlateCarree())
ai = plt.contour(longrid, latgrid, mag150.incl.values, levels=np.arange(-90,90.1,1), 
                 colors='b', transform=ccrs.PlateCarree(), linestyles='dashed')
ai.clabel(inline=True, fmt = '%d', fontsize=12, colors='b')