#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 14:08:42 2019

@author: smrak
"""

import datetime
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature.nightshade import Nightshade


fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal())
ax.coastlines(color='black', resolution='110m')
# ns
date = datetime.datetime(2017, 8, 21, 11)
ax.add_feature(Nightshade(date, alpha=0.2))


ax.set_title('Night time shading for {}'.format(date))
#ax.stock_img()

gl = ax.gridlines(crs=ccrs.PlateCarree(), color='black', draw_labels=True,
                          linestyle='--', linewidth=1)
gl.xlabels_top = False
gl.ylabels_left = False


plt.show()
