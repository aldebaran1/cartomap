# CartoMap

CartoPy based map utils

Cartomap is a simplified interface with an easy API for CartoPy library.

## Install

1. Install cartopy: due to the PROJ.4 prereqs, most users install by:
   ```sh
   conda install cartopy
   ```
   optionally, latest development version:
   ```sh
   git clone https://github.com/SciTools/Cartopy
   cd cartopy
   pip install -e .
   ```
2. Get and install cartomap
   ```sh
   git clone https://github.com/aldebaran1/cartomap.git
   cd cartomap
   pip install -e .
   ```

## API
Load and retun a map:
```python
import cartomap as cm

fig = cm.plotCartoMap(arguments)
```

Supported arguments list:
* latlim, # Latitude limits
* lonlim, # Longitude limits
* parallels, # Specify parallels to draw
* meridians, # Specify meridians to draw
* figsize, # Define figure size -> plt.figure(figsize=figsize)
* projection, # Projection type, Look below
* title, # Figure title
* resolution, # As pet CartoPy, three options are possible 110m, 50m or 10m
* states, # Draw states
* grid_linewidth, # Grid == meridians&parallels
* grid_color,
* grid_linestyle, 
* terrain, # Orographic colormap, defults as it comes with CartoPy
* background_color,
* border_color='k'. # Border=states and countries

Cartomap supports the following projections:
* Sterographic as 'stereo',
* Mercator as 'merc',
* PlateCarree as 'plate',
* LambertConformal as 'lambert'.

Gridlines are automatically computed, this interface includes additional routines to include Gridline ticks and labels for 
Sterographic and Lambert projection, which are not included in CartoPy as per today, ie, Cartopy v0.16.0
