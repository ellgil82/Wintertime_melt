# coding= ASCII
""" Script for plotting model domain and plotting point source data (such as automatic weather stations) on a map.

Dependencies:
- iris 1.11.0
- matplotlib 1.5.1
- numpy 1.10.4

Author: Ella Gilbert, 2018.

"""

# Import modules
import iris
import numpy as np
import matplotlib.pyplot as plt
import fnmatch
import os
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable

## Define functions
def rotate_data(var, lat_dim, lon_dim):
    ## Rotate projection
    #create numpy arrays of coordinates
    rotated_lat = var.coord('grid_latitude').points
    rotated_lon = var.coord('grid_longitude').points
    ## set up parameters for rotated projection
    pole_lon = var.coord('grid_longitude').coord_system.grid_north_pole_longitude
    pole_lat = var.coord('grid_latitude').coord_system.grid_north_pole_latitude
    #rotate projection
    real_lon, real_lat = iris.analysis.cartography.unrotate_pole(rotated_lon,rotated_lat, pole_lon, pole_lat)
    print ('\nunrotating pole...')
    lat = var.coord('grid_latitude')
    lon = var.coord('grid_longitude')
    lat = iris.coords.DimCoord(real_lat, standard_name='latitude',long_name="grid_latitude",var_name="lat",units=lat.units)
    lon= iris.coords.DimCoord(real_lon, standard_name='longitude',long_name="grid_longitude",var_name="lon",units=lon.units)
    var.remove_coord('grid_latitude')
    var.add_dim_coord(lat, data_dim=lat_dim)
    var.remove_coord('grid_longitude')
    var.add_dim_coord(lon, data_dim=lon_dim)
    return real_lon, real_lat

def load_vars(file):
	new_lsm = iris.load_cube('/data/clivarm/wip/ellgil82/May_2016/Re-runs/CS2/20160522T1200Z_Peninsula_km1p5_Smith_tnuc_pa000.pp', 'land_binary_mask')
	new_orog = iris.load_cube('/data/clivarm/wip/ellgil82/May_2016/Re-runs/CS2/20160522T1200Z_Peninsula_km1p5_Smith_tnuc_pa000.pp', 'surface_altitude')
	me_too = [new_lsm, new_orog]
	for j in me_too:
		real_lon, real_lat = rotate_data(j, 0, 1)
	return new_orog, new_lsm, real_lon, real_lat

## Set up plotting options
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Helvetica', 'Liberation sans', 'Tahoma', 'DejaVu Sans','Verdana']

f = '/data/clivarm/wip/ellgil82/May_2016/Re-runs/CS2/20160525T1200Z_Peninsula_km1p5_Smith_tnuc_pa012.pp'

orog, lsm, real_lon, real_lat = load_vars(f)

def draw_map():
    fig, ax = plt.subplots(figsize=(12,8), edgecolor = 'dimgrey')
    #ax = fig.add_subplot(111)
    m = Basemap(projection='stere', llcrnrlat= min(real_lat+0.2) , urcrnrlat= max(real_lat), llcrnrlon= min(real_lon), urcrnrlon= max(real_lon), lat_ts = -67, lon_0 = min(real_lon)+((max(real_lon)-min(real_lon))/2), lat_0 = min(real_lat)+((max(real_lat)-min(real_lat+0.2))/2))
    #ax.patch.set_facecolor('#566da0')
    ax.set_xlim(min(real_lon), max(real_lon))
    ax.set_ylim(min(real_lat +0.2), max(real_lat))
    xlon, ylat = np.meshgrid(real_lon,real_lat)
    lsm_masked = np.ma.masked_where(lsm.data==0, lsm.data)
    CI = m.scatter(-63.37105, -66.48272, s=100, marker='o', facecolor='#fb9a99', color = '#222222', latlon=True, zorder=10)
    #AWS14 = m.scatter(-61.03, -67.01, marker='o', s=100, facecolor='#fb9a99', color = '#222222',latlon=True, zorder=10)
    #AWS15 = m.scatter(-62.09, -67.34, marker='o', s=100, facecolor='#fb9a99', color = '#222222',  latlon=True, zorder=10)
    #ax.annotate(xy = (0.38, 0.48), s=' Cabinet \nInlet AWS', fontsize='26', color='#222222',  xycoords = 'axes fraction' )
    #ax.annotate(xy = (0.45, 0.53), s='AWS 14', fontsize='26', color='#222222',  xycoords = 'axes fraction' )
    #ax.annotate(xy=(0.37, 0.335), s='AWS 15', fontsize = '26', color = '#222222',  xycoords = 'axes fraction')
    m.contourf(xlon, ylat, lsm_masked, colors='w', latlon=True, zorder  =1)
    m.contour(xlon, ylat, lsm.data, levels = [0], colors='#222222', lw = 2, latlon=True, zorder=2)
    m.contour(xlon, ylat, orog.data, levels = [10], colors = '#222222', linewidth = 1.5, latlon= True, zorder= 3)
    merid = m.drawmeridians(meridians=np.arange(np.around(min(real_lon)),np.around(max(real_lon)),6), labels = [False,False, True,False,True,True], fontsize =30 , color = '#222222')
    parallels=np.arange(np.around(min(real_lat)),np.around(max(real_lat)),2)
    m.drawmapscale(max(real_lon)-2,min(real_lat)+0.6,-70,-67.54,200,'fancy','km', fontsize = 24, fillcolor1 = 'w', fillcolor2 = '#222222', fontcolor = '#222222')
    par = m.drawparallels(parallels, labels=[True, False,True,False],fontsize =30, color = '#222222' )
    m.drawmapboundary(color='dimgrey', linewidth=2, fill_color='#566da0', ax = ax)
    # Change colours of labels
    def set_colour(x, colour):
        for m in x:
            for t in x[m][1]:
                t.set_color(colour)
    set_colour(par, 'dimgrey')
    set_colour(merid, 'dimgrey')
    plt.savefig('/users/ellgil82/figures/Maps/AWS_loc_melt.eps')
    plt.savefig('/users/ellgil82/figures/Maps/AWS_loc_melt.png')
    plt.show()

draw_map()

