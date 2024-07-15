#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 10 10:32:04 2018

@author: smrak
"""
import numpy as np
from datetime import datetime
try:
    import apexpy as ap
except: 
    pass
try:
    import igrf
except:
    pass
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.feature.nightshade import Nightshade
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
try:
    from sunrise import terminator as ter
except:
    pass

projection_dict = {'stereo': ccrs.Stereographic(), 
                  'merc': ccrs.Mercator(),
                  'plate': ccrs.PlateCarree(),
                  'lambert': ccrs.LambertConformal(),
                  'mollweide': ccrs.Mollweide(),
                  'north': ccrs.NorthPolarStereo(),
                  'south': ccrs.NorthPolarStereo(),
                  'ortographic': ccrs.Orthographic(),
                  'robinson': ccrs.Robinson(),
                  'miller': ccrs.Miller(),
                  'geostationary': ccrs.Geostationary(),
                  }

def carto_quiver(ax, x, y, vx, vy, width=1, color='b', scale=1, figure=False):
    vxx = np.divide(vx, np.cos(np.radians(y)))
    vyy = vy
    magn = np.sqrt(vx**2 + vy**2)
    magn_src_crs = np.sqrt(vxx**2 + vyy**2)
    
    Q = ax.quiver(x,y, vxx*magn/magn_src_crs, vyy*magn/magn_src_crs, 
              width=width, scale=scale, color=color, 
              transform=ccrs.PlateCarree())
    if figure:
        return Q, ax
    else:
        return ax

def plotCartoMap(latlim=[0, 75], lonlim=[-40, 40], parallels=None, meridians=None,
                 pole_center_lon=0,figsize=(12, 8), terrain=False, ax=False,
                 projection='stereo', title='', resolution='110m', lon0=None,lat0=None,
                 states=True, grid_linewidth=0.5, grid_color='black', 
                 grid_linestyle='--', background_color=None, border_color='k',
                 figure=False, nightshade=False, ns_alpha=0.1,
                 apex=False, date=None, 
                 mlat_levels=None, mlon_levels=None, alt_km=0.0,
                 mlon_colors='blue', mlat_colors='red', mgrid_width=1,
                 mgrid_labels=True, mgrid_fontsize=12, mlon_cs='mlon',
                 incl_levels=None, decl_levels=None, igrf_param='incl',
                 mlon_labels=True, mlat_labels=True, mgrid_style='--',
                 label_colors='k', apex_alt=0,
                 decl_colors='k', incl_colors='k',
                 terminator=False, terminator_altkm=350,
                 polarization_terminator=False, pt_hemisphere='north',
                 ter_color='red', ter_style='--', ter_width=2,
                 midnight=False, midnight_colors='m', midnight_width=2,
                 midnight_style='--'):
    if lonlim is None or lonlim == []:
        lonlim = [-180, 180]
    STATES = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none')
    
    
    if not ax:
        if figsize is None:
            fig = plt.figure()
        else:
            fig = plt.figure(figsize=figsize)
        
        if projection == 'stereo':
            ax = plt.axes(projection=ccrs.Stereographic(central_longitude=(sum(lonlim)/2)))
        elif projection == 'merc':
            ax = plt.axes(projection=ccrs.Mercator())
        elif projection == 'plate':
            ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=(sum(lonlim)/2)))
        elif projection == 'lambert':
            ax = plt.axes(projection=ccrs.LambertConformal(central_longitude=(sum(lonlim)/2),
                                                           central_latitude=(sum(latlim)/2)))
        elif projection == 'mollweide':
            if lon0 is not None:
                ax = plt.axes(projection=ccrs.Mollweide(central_longitude=lon0))
            else:
                ax = plt.axes(projection=ccrs.Mollweide(central_longitude=(sum(lonlim)/2)))
        elif projection == 'robinson':
            if lon0 is not None:
                ax = plt.axes(projection=ccrs.Robinson(central_longitude=lon0))
            else:
                ax = plt.axes(projection=ccrs.Robinson(central_longitude=(sum(lonlim)/2)))
        elif projection == 'north':
            if lon0 is None:
                lon0 = 0
            ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=lon0))
        elif projection == 'south':
            if lon0 is None:
                lon0 = 0
            ax = plt.axes(projection=ccrs.SouthPolarStereo(central_longitude=lon0))
        elif projection == 'ortographic':
            if lon0 is None:
                lon0 = 0
            if lat0 is None:
                lat0 = 0
            ax = plt.axes(projection=ccrs.Orthographic(central_longitude=lon0, central_latitude=lat0))
        else:
            print ("Projection is invalid. Please enter the right one. \n \
                   'stereo', 'merc', 'plate', 'lambret', mollweide', 'north, 'south', 'ortographic'")
            return 0
    # if background_color is not None:
        ax.patch.set_facecolor(background_color)
        # ax.background_img.set_facecolor(background_color)
    ax.set_title(title)
    ax.coastlines(color=border_color, resolution=resolution)  # 110m, 50m or 10m
    if states:
        ax.add_feature(STATES, edgecolor=border_color)
    ax.add_feature(cfeature.BORDERS, edgecolor=border_color)
    if terrain:
        ax.stock_img()
    if nightshade:
        assert date is not None
        assert ns_alpha is not None
        ax.add_feature(Nightshade(date, ns_alpha))
    # Draw Parralels
    if projection == 'merc' or projection == 'plate':
        if isinstance(meridians, np.ndarray):
            meridians = list(meridians)
        if isinstance(parallels, np.ndarray):
            parallels = list(parallels)
        
        gl = ax.gridlines(crs=ccrs.PlateCarree(), color=grid_color, draw_labels=False,
                          linestyle=grid_linestyle, linewidth=grid_linewidth)
        if meridians is None:
            gl.xlines = False
        else:
            if len(meridians) > 0:
                gl.xlocator = mticker.FixedLocator(meridians)
                gl.xlabels_bottom = True
            else:
                gl.ylines = False
                
        if parallels is None:
            gl.ylines = False
        else:
            if len(parallels) > 0:
                gl.ylocator = mticker.FixedLocator(parallels)
                gl.ylabels_left = True
            else:
                gl.ylines = False
    else:
        gl = ax.gridlines(crs=ccrs.PlateCarree(), color=grid_color, draw_labels=False,
                          linestyle=grid_linestyle, linewidth=grid_linewidth)
        if meridians is None:
            gl.xlines = False
        else:
            gl.xlocator = mticker.FixedLocator(meridians)
        if parallels is not None: 
            if isinstance(parallels, np.ndarray):
                parallels = list(parallels)
            gl.ylocator = mticker.FixedLocator(parallels)
        else:
            gl.ylines = False
        
    # Geomagnetic coordinates @ Apex
    if apex:
        if date is None:
            date = datetime(2017, 12, 31, 0, 0, 0)
        assert isinstance(date, datetime)
        
        A = ap.Apex(date = date)
        # Define levels and ranges for conversion
        if mlon_cs == 'mlt':
            if mlon_levels is None:
                mlon_levels = np.array([])
                mlon_range = np.arange(0, 24.01, 0.01)
            elif isinstance(mlon_levels, bool):
                if mlon_levels == False:
                    mlon_levels = np.array([])
                    mlon_range = np.arange(0, 24.01, 0.01)
            else:
                mlon_range = np.arange(mlon_levels[0], mlon_levels[-1]+0.1, 0.01)
        else:
            if mlon_levels is None:
                mlon_levels = np.array([])
                mlon_range = np.arange(-180,180,0.5)
            elif isinstance(mlon_levels, bool):
                if mlon_levels == False:
                    mlon_levels = np.array([])
                    mlon_range = np.arange(-180, 181, 0.1)
            else:
                mlon_range = np.arange(mlon_levels[0], mlon_levels[0]+362, 0.1)
        if mlat_levels is None:
            mlat_levels = np.arange(-90, 90.1, 1)
        mlat_range = np.arange(mlat_levels[0], mlat_levels[-1]+0.1, 0.1)
        
        # Do meridans
        for mlon in mlon_levels:
            MLON = mlon * np.ones(mlat_range.size)
            if mlon_cs == 'mlt':
                y, x = A.convert(mlat_range, MLON, 'mlt', 'geo', datetime=date, height=apex_alt)
            else:
                y, x  = A.convert(mlat_range, MLON, 'apex', 'geo', height=apex_alt)
            mlat_mask_extent = (y >=  latlim[0] + 0.2) & (y <=  latlim[1] - 0.2)
            # Plot meridian
            inmap = np.logical_and(x >= lonlim[0], x <= lonlim[1])
            if np.sum(inmap) > 10:
                ax.plot(np.unwrap(x[mlat_mask_extent], 180), 
                        np.unwrap(y[mlat_mask_extent], 90), c=mlon_colors, 
                         lw=mgrid_width, linestyle=mgrid_style, zorder=90,
                         transform=ccrs.PlateCarree())
                
            # Labels
            if mlon_labels:
                ix = abs(y-np.mean(latlim)).argmin()
                mx = x[ix] - 1 if mlon >=10 else x[ix] - 0.5
                my = np.mean(latlim)
                if np.logical_and(mx >= lonlim[0], mx <= lonlim[1]):
                    if mlon_cs == 'mlt' and mlon != 0:
                        ax.text(mx, my, str(int(mlon)), color=label_colors, 
                                 fontsize=14, backgroundcolor='white',
                                 transform=ccrs.PlateCarree())
                    elif mlon_cs != 'mlt' and mlon != 360:
                        ax.text(mx, my, str(int(mlon)), color=label_colors, 
                                 fontsize=14, backgroundcolor='white',
                                 transform=ccrs.PlateCarree())
        # Do parallels
        for mlat in mlat_levels:
            MLAT = mlat * np.ones(mlon_range.size)
            if mlon_cs == 'mlt':
                gy, gx = A.convert(MLAT, mlon_range, 'mlt', 'geo', datetime=date, height=apex_alt)
            else:
                
                gy, gx = A.convert(MLAT, mlon_range, 'apex', 'geo', datetime=date, height=apex_alt)
            inmap = np.logical_and(gy >= latlim[0], gy <= latlim[1])
            if np.sum(inmap) > 20:
                ax.plot(np.unwrap(gx, 180), np.unwrap(gy, 90), c=mlat_colors,
                         lw=mgrid_width, linestyle=mgrid_style, zorder=90,
                         transform=ccrs.PlateCarree())
                
            # Labels
            if mlat_labels:
                ix = abs(gx-np.mean(lonlim)).argmin()
                mx = np.mean(lonlim)
                my = gy[ix] - 0.5
                if np.logical_and(mx >= lonlim[0], mx <= lonlim[1]) and \
                np.logical_and(my >= latlim[0], my <= latlim[1]):
                    ax.text(mx, my, str(int(mlat)), color=label_colors, 
                             fontsize=14, backgroundcolor='white',
                             transform=ccrs.PlateCarree())
    
    if mlat_levels is not None:
        glon = np.arange(-180, 180.1, 4)
        glat = np.arange(-90, 90.1, 2)
        longrid, latgrid = np.meshgrid(glon, glat)
        
        incl = igrf.grid(date, glat=latgrid, glon=longrid, alt_km=alt_km).incl.values
        ax.contour(longrid, latgrid, incl, levels=mlat_levels,  zorder=90,
                   colors=mlat_colors, transform=ccrs.PlateCarree())
        # if incl_levels is not None:
        #     z = mag.incl.values
        #     for inclination in incl_levels:
        #         ax.contour(longrid, latgrid, z, levels=declination,  zorder=90,
        #                          colors=incl_colors, transform=ccrs.PlateCarree())
    # Terminators
    if terminator:
        assert date is not None
        if not isinstance(terminator_altkm, list):
            terminator_altkm = [terminator_altkm]
        for takm in terminator_altkm:
            try:
                glon_ter, glat_ter = ter.get_terminator(date, alt_km = takm)
                if glon_ter is not None and glat_ter is not None:
                    if isinstance(glon_ter, list):
                        for i in range(len(glon_ter)):
                            ax.plot(np.unwrap(glon_ter[i], 180), np.unwrap(glat_ter[i], 90),
                                c=ter_color, lw=ter_width, ls=ter_style, zorder=90,
                                transform=ccrs.PlateCarree())
                    else:
                        ax.plot(np.unwrap(glon_ter, 180), np.unwrap(glat_ter, 90),
                                c=ter_color, lw=ter_width, ls=ter_style, zorder=90,
                                transform=ccrs.PlateCarree())
            except:
                pass
    if midnight:
        mlat_range = np.arange(-89.9, 90.1, 0.1)
        MLON = 0 * np.ones(mlat_range.size)
        if mlon_cs == 'mlt':
            y, x = A.convert(mlat_range, MLON, 'mlt', 'geo', datetime=date, height=apex_alt)
        else:
            y, x  = A.convert(mlat_range, MLON, 'apex', 'geo', height=apex_alt)
        mlat_mask_extent = (y >=  latlim[0] + 0.2) & (y <=  latlim[1] - 0.2)
        # Plot meridian
        inmap = np.logical_and(x >= lonlim[0], x <= lonlim[1])
        if np.sum(inmap) > 10:
            ax.plot(np.unwrap(x[mlat_mask_extent], 180), 
                    np.unwrap(y[mlat_mask_extent], 90), c=midnight_colors, 
                     lw=midnight_width, linestyle=midnight_style, zorder=90,
                     transform=ccrs.PlateCarree())
    # Set Extent
    if projection == 'north' or projection == 'south':
        import matplotlib.path as mpath
        ax.set_extent([-180, 181, latlim[0], latlim[1]], crs=ccrs.PlateCarree())
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)
        ax.set_boundary(circle, transform=ax.transAxes)
    elif lonlim[0] != -180 and lonlim[1] != 180:
        ax.set_extent([lonlim[0], lonlim[1], latlim[0], latlim[1]], crs=ccrs.PlateCarree())#ccrs.PlateCarree())
    
    if 'fig' in locals():
        return fig, ax
    else:
        return ax
