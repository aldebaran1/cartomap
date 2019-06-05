#!/usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib.path as mpath
import numpy as np
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import h5py
from argparse import ArgumentParser
from datetime import datetime
import os
from glob import glob

months = {1: 'jan', 2: 'feb', 3: 'mar', 4: 'apr', 5: 'may', 6: 'jun',
          7: 'jul', 8: 'aug', 9: 'sep', 10: 'oct', 11: 'nov', 12: 'dec'}

projections = {'plate': [ccrs.PlateCarree(), 'Plate Carree'],
               'near': [ccrs.NearsidePerspective(), 'Nearside Perspective'],
               'polar': [ccrs.NorthPolarStereo(), 'Polar Stereo'],
               'mercator': [ccrs.Mercator(), 'Mercator'],
               'geostat': [ccrs.Geostationary(), 'Geostationary']}

cmaps = plt.colormaps()

def save(root: str = None,
         slide: str = None,
         proj: str = None,
         lim: float = None,
         cmap: str = None):

    f = h5py.File(root, 'r')
    lat = f['GPSTEC']['lat']
    lon = f['GPSTEC']['lon']
    t = f['GPSTEC']['time']

    if cmap not in cmaps:
        cmap = 'gist_ncar'

    if slide is not None:
        slide = int(slide)
        im = f['GPSTEC']['im'][0:][0:][slide]
        im = np.transpose(im)
        time = datetime.fromtimestamp(t[slide])

        # scale cmap
        cmax = np.max(list(filter(lambda x: ~np.isnan(x), np.reshape(im, 64800))))
        cmin = np.min(list(filter(lambda x: ~np.isnan(x), np.reshape(im, 64800))))
        minoff = 0  # offset
        maxoff = 0

        if proj == 'polar':
            fig = plt.figure()

            theta = np.linspace(0, 2 * np.pi, 100)
            center, radius = [0.5, 0.5], 0.5
            verts = np.vstack([np.sin(theta), np.cos(theta)]).T
            circle = mpath.Path(verts * radius + center)

            ax1 = plt.subplot(1, 2, 1, projection=ccrs.SouthPolarStereo())
            ax1.set_extent([-180, 180, -90, 0], ccrs.PlateCarree())
            ax1.title.set_text('South Polar Stereographic')
            ax1.add_feature(cfeature.OCEAN, zorder=1)
            ax1.add_feature(cfeature.LAKES, zorder=1)
            ax1.add_feature(cfeature.RIVERS, zorder=1)
            ax1.add_feature(cfeature.LAND, zorder=1)
            ax1.add_feature(cfeature.BORDERS, zorder=3)
            ax1.add_feature(cfeature.COASTLINE, zorder=3)
            ax1.gridlines()
            ax1.set_boundary(circle, transform=ax1.transAxes)

            im1 = ax1.pcolormesh(lon, lat, im, transform=ccrs.PlateCarree(), vmin=cmin + minoff, vmax=cmax - maxoff,
                                 cmap=cmap, zorder=2)

            ax2 = plt.subplot(1, 2, 2, projection=ccrs.NorthPolarStereo())
            ax2.set_extent([-180, 180, 90, 0], ccrs.PlateCarree())
            ax2.title.set_text('North Polar Stereographic')
            ax2.add_feature(cfeature.OCEAN, zorder=1)
            ax2.add_feature(cfeature.LAKES, zorder=1)
            ax2.add_feature(cfeature.RIVERS, zorder=1)
            ax2.add_feature(cfeature.LAND, zorder=1)
            ax2.add_feature(cfeature.BORDERS, zorder=3)
            ax2.add_feature(cfeature.COASTLINE, zorder=3)
            ax2.gridlines()
            ax2.set_boundary(circle, transform=ax2.transAxes)
            im2 = ax2.pcolormesh(lon, lat, im, transform=ccrs.PlateCarree(), vmin=cmin + minoff, vmax=cmax - maxoff,
                                 cmap=cmap, zorder=2)
            fig.subplots_adjust(right=0.8)
            cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
            fig.colorbar(im2, cax=cbar_ax, label='Total Electron Concentration [TECu]')

        elif proj is not 'polar':
            fig = plt.figure('TEC ({})'.format(datetime.fromtimestamp(t[slide])))
            ax = plt.subplot(1, 1, 1, projection=projections[proj][0])
            ax.title.set_text(projections[proj][1])
            ax.add_feature(cfeature.OCEAN, zorder=1)
            ax.add_feature(cfeature.LAKES, zorder=1)
            ax.add_feature(cfeature.RIVERS, zorder=1)
            ax.add_feature(cfeature.LAND, zorder=1)
            ax.add_feature(cfeature.BORDERS, zorder=3)
            ax.add_feature(cfeature.COASTLINE, zorder=3)
            ax.gridlines()
            imcm = ax.pcolormesh(lon, lat, im, transform=ccrs.PlateCarree(), vmin=cmin + minoff, vmax=cmax - maxoff,
                                 cmap=cmap, zorder=2)
            cb = fig.colorbar(imcm, shrink=0.5)
            cb.set_label('Total Electron Content [TECu]')
        print('Saving slide {}'.format(slide))
        os.mkdir(os.path.split(root)[0] + '\\{}{}'.format(months[time.month], time.day))
        folder = os.path.join(os.path.split(root)[0], '{}{}'.format(months[time.month], time.day))
        print(folder)

        # plt.savefig(os.path.join(folder, '{}.png'.format(slide)))
        figsav = plt.gcf()
        figsav.suptitle('{}'.format(time))
        figsav.set_size_inches((10, 5), forward=False)
        figsav.savefig(os.path.join(folder, '{}.png'.format(str(slide).zfill(3))), dpi=200)
        plt.close(fig)
        plt.close(figsav)
    else:
        t0 = datetime.fromtimestamp(t[0])
        os.mkdir(os.path.split(root)[0] + '\\{}{}'.format(months[t0.month], t0.day))
        folder = os.path.join(os.path.split(root)[0], '{}{}'.format(months[t0.month], t0.day))

        if proj == 'polar':

            theta = np.linspace(0, 2 * np.pi, 100)
            center, radius = [0.5, 0.5], 0.5
            verts = np.vstack([np.sin(theta), np.cos(theta)]).T
            circle = mpath.Path(verts * radius + center)

            for slide in list(range(len(t))):
                time = datetime.fromtimestamp(t[slide])
                if lim is not 0:
                    cmax = lim
                    cmin = 0.0
                else:
                    cmax = np.max(list(filter(lambda x: ~np.isnan(x), np.reshape(im, 64800))))
                    cmin = np.min(list(filter(lambda x: ~np.isnan(x), np.reshape(im, 64800))))

                fig = plt.figure()
                im = f['GPSTEC']['im'][0:][0:][slide]
                im = np.transpose(im)

                ax1 = plt.subplot(1, 2, 1, projection=ccrs.SouthPolarStereo())
                ax1.set_extent([-180, 180, -90, 0], ccrs.PlateCarree())
                ax1.title.set_text('South Polar Stereographic')
                ax1.add_feature(cfeature.OCEAN, zorder=1)
                ax1.add_feature(cfeature.LAKES, zorder=1)
                ax1.add_feature(cfeature.RIVERS, zorder=1)
                ax1.add_feature(cfeature.LAND, zorder=1)
                ax1.add_feature(cfeature.BORDERS, zorder=3)
                ax1.add_feature(cfeature.COASTLINE, zorder=3)
                ax1.gridlines()
                ax1.set_boundary(circle, transform=ax1.transAxes)

                ax1.pcolormesh(lon, lat, im, transform=ccrs.PlateCarree(), vmin=cmin, vmax=cmax,
                               cmap=cmap, zorder=2)

                ax2 = plt.subplot(1, 2, 2, projection=ccrs.NorthPolarStereo())
                ax2.set_extent([-180, 180, 90, 0], ccrs.PlateCarree())
                ax2.title.set_text('North Polar Stereographic')
                ax2.add_feature(cfeature.OCEAN, zorder=1)
                ax2.add_feature(cfeature.LAKES, zorder=1)
                ax2.add_feature(cfeature.RIVERS, zorder=1)
                ax2.add_feature(cfeature.LAND, zorder=1)
                ax2.add_feature(cfeature.BORDERS, zorder=3)
                ax2.add_feature(cfeature.COASTLINE, zorder=3)
                ax2.gridlines()
                ax2.set_boundary(circle, transform=ax2.transAxes)
                imcm = ax2.pcolormesh(lon, lat, im, transform=ccrs.PlateCarree(), vmin=cmin, vmax=cmax,
                                     cmap=cmap, zorder=2)

                fig.subplots_adjust(right=0.8)
                cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
                fig.colorbar(imcm, cax=cbar_ax, label='Total Electron Concentration [TECu]')

                print('Saving slide {}/{}...'.format(slide+1, len(t)))
                figsav = plt.gcf()
                figsav.suptitle('{}'.format(time))
                figsav.set_size_inches((10, 5), forward=False)
                figsav.savefig(os.path.join(folder, '{}.png'.format(str(slide).zfill(3))), dpi=200)
                plt.close(fig)
                plt.close(figsav)

            print(folder)
        else:
            for slide in list(range(len(t))):
                fig = plt.figure()
                im = f['GPSTEC']['im'][0:][0:][slide]
                im = np.transpose(im)
                if lim is not 0:
                    cmax = lim
                    cmin = 0.0
                else:
                    cmax = np.max(list(filter(lambda x: ~np.isnan(x), np.reshape(im, 64800))))
                    cmin = np.min(list(filter(lambda x: ~np.isnan(x), np.reshape(im, 64800))))

                ax = plt.subplot(1, 1, 1, projection=projections[proj][0])
                ax.title.set_text(projections[proj][1])
                ax.add_feature(cfeature.OCEAN, zorder=1)
                ax.add_feature(cfeature.LAKES, zorder=1)
                ax.add_feature(cfeature.RIVERS, zorder=1)
                ax.add_feature(cfeature.LAND, zorder=1)
                ax.add_feature(cfeature.BORDERS, zorder=3)
                ax.add_feature(cfeature.COASTLINE, zorder=3)
                ax.gridlines()
                im = ax.pcolormesh(lon, lat, im, transform=ccrs.PlateCarree(), vmin=cmin, vmax=cmax,
                                   cmap=cmap, zorder=2)
                cb = fig.colorbar(im, shrink=0.5)
                cb.set_label('Total Electron Content [TECu]')
                # fig.tight_layout()
                print('Saving slide {}/{}...'.format(slide+1, len(t)))
                # plt.savefig(os.path.join(folder, '{}.png'.format(slide)))

                figsav = plt.gcf()
                figsav.set_size_inches((10, 5), forward=False)
                figsav.savefig(os.path.join(folder, '{}.png'.format(slide)), dpi=200)
                plt.close(fig)
                plt.close(figsav)

            print(folder)


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('root', type=str, help='local address')
    p.add_argument('-s', '--slide', type=str, help='slide number [0,239]')
    # p.add_argument('-o', '--odir', type=str, help='directory to save images')
    p.add_argument('-p', '--proj', type=str, help='map projection - plate or polar', default='polar')
    p.add_argument('-l', '--lim', type=float, help='absolute limit of colorbar - 0 for no absolute', default=70)
    p.add_argument('-c', '--cmap', type=str, help='colormap', default=None)

    P = p.parse_args()

    root = P.root

    if os.path.splitext(root)[1] in ['.h5', '.hdf5']:
        save(root=P.root, slide=P.slide, proj=P.proj, lim=P.lim, cmap=P.cmap)

    else:
        flist = sorted(glob(os.path.split(root)[0] + '\\conv*.h5'))
        if len(flist) > 0:
            for file in flist:
                save(file, slide=P.slide, proj=P.proj, lim=P.lim, cmap=P.cmap)

