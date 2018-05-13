#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 10 10:32:04 2018

@author: smrak
"""

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
#import matplotlib.colors as colors
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

def plotCartoMap(latlim=[0,75],lonlim=[-40,40],parallels=[],meridians=[],
                 figsize=(12,8),projection='stereo',title='',resolution='110m',
                 states=True,grid_linewidth=0.5,grid_color='black',terrain=False,
                 grid_linestyle='--', background_color=None,border_color='k'):

    STATES = cfeature.NaturalEarthFeature(
            category='cultural',
            name='admin_1_states_provinces_lines',
            scale='50m',
            facecolor='none')
    if figsize is None:
        figsize = (12,8)
    plt.figure(figsize=figsize)
    if projection == 'stereo':
        ax = plt.axes(projection=ccrs.Stereographic(central_longitude=(sum(lonlim)/2)))
    if projection == 'merc':
        ax = plt.axes(projection=ccrs.Mercator())
    if projection == 'plate':
        ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=(sum(lonlim)/2)))
    if background_color is not None:
        ax.background_patch.set_facecolor(background_color)
    ax.set_title(title)
    ax.coastlines(color=border_color,resolution=resolution) # 110m, 50m or 10m
    if states:
        ax.add_feature(STATES, edgecolor=border_color)
    ax.add_feature(cfeature.BORDERS,edgecolor=border_color)
    if terrain:
        ax.stock_img()
    ax.set_extent([lonlim[0], lonlim[1], latlim[0], latlim[1]])
    
    if projection != 'merc':
        gl = ax.gridlines(crs=ccrs.PlateCarree(),color=grid_color,
                          linestyle=grid_linestyle,linewidth=grid_linewidth)
    else:
        gl = ax.gridlines(crs=ccrs.PlateCarree(),color=grid_color,draw_labels=True,
                          linestyle=grid_linestyle,linewidth=grid_linewidth)
        gl.xlabels_top=False
        gl.ylabels_right=False
    gl.xlocator = mticker.FixedLocator(meridians)
    gl.ylocator = mticker.FixedLocator(parallels)
    gl.xlabels_top = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    
    return ax