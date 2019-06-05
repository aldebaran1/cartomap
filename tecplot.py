#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import h5py
from argparse import ArgumentParser
from datetime import datetime

def plotter(root: str = None,
            slide: int = None):

    f = h5py.File(root, 'r')
    im = f['GPSTEC']['im'][0:][0:][slide]
    im = np.transpose(im)
    lat = f['GPSTEC']['lat']
    lon = f['GPSTEC']['lon']
    t = f['GPSTEC']['time']

    # scale cmap
    cmax = np.max(list(filter(lambda x: ~np.isnan(x), np.reshape(im, 64800))))
    cmin = np.min(list(filter(lambda x: ~np.isnan(x), np.reshape(im, 64800))))
    minoff = 0  # offset
    maxoff = 0

    fig = plt.figure('VTEC ({})'.format(datetime.fromtimestamp(t[slide])), figsize=(12, 6))

    ax1 = plt.subplot(1, 2, 1, projection=ccrs.SouthPolarStereo())  # projection
    ax1.title.set_text('South Polar Stereographic')
    ax1.add_feature(cfeature.OCEAN, zorder=1)
    ax1.add_feature(cfeature.LAKES, zorder=1)
    ax1.add_feature(cfeature.RIVERS, zorder=1)
    ax1.add_feature(cfeature.LAND, zorder=1)
    ax1.add_feature(cfeature.BORDERS, zorder=3)
    ax1.add_feature(cfeature.COASTLINE, zorder=3)
    ax1.gridlines()
    im1 = ax1.pcolormesh(lon, lat, im, transform=ccrs.PlateCarree(), vmin=cmin + minoff, vmax=cmax - maxoff,
                          cmap='viridis', zorder=2)
    cb1 =fig.colorbar(im1, shrink=0.5)
    cb1.set_label('Total Electron Content [$10^{16} e^{-}/m^{2}$]')

    ax2 = plt.subplot(1, 2, 2, projection=ccrs.NorthPolarStereo())  # projection
    ax2.title.set_text('North Polar Stereographic')
    ax2.add_feature(cfeature.OCEAN, zorder=1)
    ax2.add_feature(cfeature.LAKES, zorder=1)
    ax2.add_feature(cfeature.RIVERS, zorder=1)
    ax2.add_feature(cfeature.LAND, zorder=1)
    ax2.add_feature(cfeature.BORDERS, zorder=3)
    ax2.add_feature(cfeature.COASTLINE, zorder=3)
    ax2.gridlines()
    im2 = ax2.pcolormesh(lon, lat, im, transform=ccrs.PlateCarree(), vmin=cmin + minoff, vmax=cmax - maxoff,
                          cmap='viridis', zorder=2)
    cb2 = fig.colorbar(im2, shrink=0.5)
    cb2.set_label('Total Electron Content [$10^{16}  e^{-}/m^{2}$]')

    # fig.subplots_adjust(bottom=0.05, top=0.95, left=0.04, right=0.95, wspace=0.02)
    fig.tight_layout()
    plt.show()


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('file', type=str, help='local file address')
    p.add_argument('-s', '--slide', type=int, help='slide number [0,239]', default=0)

    P = p.parse_args()

    plotter(root=P.file, slide=P.slide)
