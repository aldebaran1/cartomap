#!/usr/bin/env python
import cartomap as cm
import matplotlib.pyplot as plt
import cartopy.crs as ccrs


cm.plotCartoMap(projection='north', terrain=True)

ny_lon, ny_lat = -74.00, 40.71
delhi_lon, delhi_lat = 77.23, 28.61

plt.plot([ny_lon, delhi_lon], [ny_lat, delhi_lat],
         color='blue', linewidth=2, marker='x',
         transform=ccrs.Geodetic()
         )

plt.show()



######################### Supported Arguments ##############################
# latlim, # Latitude limits
# lonlim, # Longitude limits
# parallels, # Specify parallels to draw
# meridians, # Specify meridians to draw
# figsize, # Define figure size -> plt.figure(figsize=figsize)
# projection, # Projection type, Look below
# title, # Figure title
# resolution, # As per CartoPy, three options are possible 110m, 50m or 10m
# states, # Draw states
# grid_linewidth, # Grid == meridians&parallels
# grid_color,
# grid_linestyle,
# terrain, # Orographic colormap, defaults as it comes with CartoPy
# background_color,
# border_color='k'. # Border=states and countries

# projections
# Sterographic as 'stereo',
# Mercator as 'merc',
# PlateCarree as 'plate',
# LambertConformal as 'lambert'
# Mollweide as 'mollweide'
# NorthPolarStereo as 'north'
# SouthPolarStereo as 'south'
