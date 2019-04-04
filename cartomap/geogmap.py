#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 10 10:32:04 2018

@author: smrak
"""
import numpy as np
from shapely import geometry as sgeom
from copy import copy
from datetime import datetime
import apexpy as ap
import igrf12

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.feature.nightshade import Nightshade
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

def find_side(ls, side):
    """
    Given a shapely LineString which is assumed to be rectangular, return the
    line corresponding to a given side of the rectangle.

    """
    minx, miny, maxx, maxy = ls.bounds
    points = {'left': [(minx, miny), (minx, maxy)],
              'right': [(maxx, miny), (maxx, maxy)],
              'bottom': [(minx, miny), (maxx, miny)],
              'top': [(minx, maxy), (maxx, maxy)], }
    return sgeom.LineString(points[side])


def lambert_xticks(ax, ticks):
    """Draw ticks on the bottom x-axis of a Lambert Conformal projection."""
    def te(xy): return xy[0]

    def lc(t, n, b): return np.vstack((np.zeros(n) + t, np.linspace(b[2], b[3], n))).T
    xticks, xticklabels = _lambert_ticks(ax, ticks, 'bottom', lc, te)
    ax.xaxis.tick_bottom()
    ax.set_xticks(xticks)
    ax.set_xticklabels([ax.xaxis.get_major_formatter()(xtick) for xtick in xticklabels])


def lambert_yticks(ax, ticks):
    """Draw ricks on the left y-axis of a Lamber Conformal projection."""
    def te(xy): return xy[1]

    def lc(t, n, b): return np.vstack((np.linspace(b[0], b[1], n), np.zeros(n) + t)).T
    yticks, yticklabels = _lambert_ticks(ax, ticks, 'left', lc, te)
    ax.yaxis.tick_left()
    ax.set_yticks(yticks)
    ax.set_yticklabels([ax.yaxis.get_major_formatter()(ytick) for ytick in yticklabels])


def _lambert_ticks(ax, ticks, tick_location, line_constructor, tick_extractor):
    """Get the tick locations and labels for an axis of a Lambert Conformal projection."""
    outline_patch = sgeom.LineString(ax.outline_patch.get_path().vertices.tolist())
    axis = find_side(outline_patch, tick_location)
    n_steps = 30
    extent = ax.get_extent(ccrs.PlateCarree())
    _ticks = []
    for t in ticks:
        xy = line_constructor(t, n_steps, extent)
        proj_xyz = ax.projection.transform_points(ccrs.PlateCarree(), xy[:, 0], xy[:, 1])
        xyt = proj_xyz[..., :2]
        ls = sgeom.LineString(xyt.tolist())
        locs = axis.intersection(ls)
        if not locs:
            tick = [None]
        else:
            tick = tick_extractor(locs.xy)
        _ticks.append(tick[0])
    # Remove ticks that aren't visible:
    ticklabels = copy(ticks)
    while True:
        try:
            index = _ticks.index(None)
        except ValueError:
            break
        _ticks.pop(index)
        ticklabels.pop(index)
    return _ticks, ticklabels


def plotCartoMap(latlim=[0, 75], lonlim=[-40, 40], parallels=[], meridians=[],
                 figsize=(12, 8), projection='stereo', title='', resolution='110m',
                 states=True, grid_linewidth=0.5, grid_color='black', terrain=False,
                 grid_linestyle='--', background_color=None, border_color='k',
                 figure=False, nightshade=False, ns_dt=None, ns_alpha=0.1,
                 geomag=False, gmagtype='apex', date=None, 
                 mlat_levels=None, mlon_levels=None, alt_km=0.0,
                 mlon_colors='blue', mlat_colors='red', mgrid_width=1,
                 mgrid_labels=True, mgrid_fontsize=12, mlon_cs='mlon',
                 incl_levels=None, decl_levels=None, igrf_param='incl'):

    STATES = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none')

    if figsize is None:
        fig = plt.figure()
    else:
        fig = plt.figure(figsize=figsize)

#    fig.gca(projection=ccrs.PlateCarree())
    if projection == 'stereo':
        ax = plt.axes(projection=ccrs.Stereographic(central_longitude=(sum(lonlim)/2)))
    elif projection == 'merc':
        ax = plt.axes(projection=ccrs.Mercator())
    elif projection == 'plate':
        ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=(sum(lonlim)/2)))
    elif projection == 'lambert':
        ax = plt.axes(projection=ccrs.LambertConformal(central_longitude=(sum(lonlim)/2),
                                                       central_latitude=(sum(latlim)/2)))
    if background_color is not None:
        ax.background_patch.set_facecolor(background_color)
    ax.set_title(title)
    ax.coastlines(color=border_color, resolution=resolution)  # 110m, 50m or 10m
    if states:
        ax.add_feature(STATES, edgecolor=border_color)
    ax.add_feature(cfeature.BORDERS, edgecolor=border_color)
    if terrain:
        ax.stock_img()
    if nightshade:
        assert ns_dt is not None
        assert ns_alpha is not None
        ax.add_feature(Nightshade(ns_dt, ns_alpha))
    

    if projection == 'merc' or projection == 'plate':
        if isinstance(meridians, np.ndarray):
            meridians = list(meridians)
        if isinstance(parallels, np.ndarray):
            parallels = list(parallels)
        
        if len(meridians) > 0 or len(parallels) > 0:
            gl = ax.gridlines(crs=ccrs.PlateCarree(), color=grid_color, draw_labels=False,
                              linestyle=grid_linestyle, linewidth=grid_linewidth)
        if len(meridians) > 0:
            gl.xlocator = mticker.FixedLocator(meridians)
            gl.xlabels_bottom = True
        else:
            gl.xlines = False
        if len(parallels) > 0:
            gl.ylocator = mticker.FixedLocator(parallels)
            ax.yaxis.set_major_formatter(LONGITUDE_FORMATTER)
            gl.ylabels_left = True
        else:
            gl.ylines = False
#        
    else:
        gl = ax.gridlines(crs=ccrs.PlateCarree(), color=grid_color,
                          linestyle=grid_linestyle, linewidth=grid_linewidth)
        ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
        ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)
        if isinstance(meridians, np.ndarray):
            meridians = list(meridians)
        if isinstance(parallels, np.ndarray):
            parallels = list(parallels)
        gl.xlocator = mticker.FixedLocator(meridians)
        gl.ylocator = mticker.FixedLocator(parallels)
        lambert_xticks(ax, meridians)
        lambert_yticks(ax, parallels)
    # Geomagnetic coordinates @ Apex
    if geomag:
        if date is None:
            date = datetime(2017, 12, 31, 0, 0, 0)
        assert isinstance(date, datetime)
        glon = np.arange(lonlim[0]-40, lonlim[1] + 40.1, 0.5)
        glat = np.arange(latlim[0], latlim[1] + 0.1, 0.1)
        longrid, latgrid = np.meshgrid(glon, glat)
            
        if gmagtype == 'apex':
            A = ap.Apex(date = date)
            
            if mlon_cs == 'mlt':
                mlat, mlon = A.convert(latgrid, longrid, 'geo', 'mlt', datetime=date)
                if mlon_levels is None:
                    mlon_levels = np.arange(0,24.1,1)
            else:
                mlat, mlon = A.convert(latgrid, longrid, 'geo', 'apex')
                if mlon_levels is None:
                    mlon_levels = np.arange(-180,180,20)
            if mlat_levels is None:
                mlat_levels = np.arange(-90,90.1,5)
            
            axy = plt.contour(glon,glat, mlat, levels = mlat_levels, colors = mlat_colors, 
                             linewidths=mgrid_width, linestyles ='solid', 
                             transform=ccrs.PlateCarree())
            axx = plt.contour(glon,glat, mlon, levels = mlon_levels, colors = mlon_colors, 
                             linewidths=mgrid_width, linestyles ='solid', 
                             transform=ccrs.PlateCarree())
            axx.clabel(inline=True, fmt = '%d', fontsize = mgrid_fontsize, colors = mlon_colors)
            axy.clabel(inline=True, fmt = '%d', fontsize = mgrid_fontsize, colors = mlat_colors)
        elif gmagtype == 'igrf':
            mag = igrf12.gridigrf12(t=date, glat=latgrid, glon=longrid, alt_km=alt_km)
            if incl_levels is None:
                incl_levels = np.arange(-90, 90.1, 2)
            if decl_levels is None:
                decl_levels = np.arange(-30, 30.1, 2)
            assert igrf_param in ('incl', 'decl')
            if igrf_param == 'incl':
                z = mag.incl.values
                lvl = incl_levels
            else:
                z = mag.decl.values
                lvl = decl_levels
            ai = plt.contour(longrid, latgrid, z, levels=lvl, 
                             colors='b', transform=ccrs.PlateCarree())
            ai.clabel(inline=True, fmt = '%d', fontsize=12, colors='b')
        else:
            print ('Wrong input argument')
    # Set Extent
    ax.set_extent([lonlim[0], lonlim[1], latlim[0], latlim[1]])
    
    return fig, ax
