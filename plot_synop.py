import iris
import iris.plot as iplt
import iris.quickplot as qplt
import numpy as np
import matplotlib.pyplot as plt
import fnmatch
import os
import matplotlib.dates as mdates
import sys
reload(sys)
sys.getdefaultencoding()
import numpy.ma as ma
import cartopy.crs as ccrs
import matplotlib

#make sure Python is looking in the right place for files
os.chdir('/data/clivarm/wip/ellgil82/May_2016/Re-runs/CS2/')

def load_files():
    surf1 = '/data/clivarm/wip/ellgil82/May_2016/Compare/CS2_new/km4p0/20160522T1200Z_Peninsula_km4p0_smoothed_pa000.pp'
    surf2 = '/data/clivarm/wip/ellgil82/May_2016/Re-runs/CS2/20160522T1200Z_Peninsula_km4p0_Smith_tnuc_pa012.pp'
    print ('importing cubes...')
    T_air = iris.load_cube(surf2, 'air_temperature')
    T_surf = iris.load_cube(surf2, 'surface_temperature')
    lsm = iris.load_cube(surf1, 'land_binary_mask')
    orog = iris.load_cube(surf1, 'surface_altitude')
    P = iris.load_cube(surf1, 'surface_air_pressure')
    P.convert_units('hPa')
    T_air.convert_units('celsius')
    T_surf.convert_units('celsius')
    ## Iris v1.11 version
    u_wind = iris.load_cube(surf2, 'x_wind')
    v_wind = iris.load_cube(surf2, 'y_wind')
    v_wind = v_wind[:,1:,:]
    Var = [T_air, T_surf, u_wind, v_wind, P]
    ## Rotate projection
    #create numpy arrays of coordinates
    rotated_lat = P.coord('grid_latitude').points
    rotated_lon = P.coord('grid_longitude').points
    ## set up parameters for rotated projection
    pole_lon = 298.5
    pole_lat = 22.99
    #rotate projection
    real_lon, real_lat = iris.analysis.cartography.unrotate_pole(rotated_lon,rotated_lat, pole_lon, pole_lat)
    print ('\nunrotating pole...')
    lat = P.coord('grid_latitude')
    lon = P.coord('grid_longitude')
    lat = iris.coords.DimCoord(real_lat, standard_name='latitude',long_name="grid_latitude",var_name="lat",units=lat.units)
    lon= iris.coords.DimCoord(real_lon, standard_name='longitude',long_name="grid_longitude",var_name="lon",units=lon.units)
    for var in Var:
        var.remove_coord('grid_latitude')
        var.add_dim_coord(lat, data_dim=1)
        var.remove_coord('grid_longitude')
        var.add_dim_coord(lon, data_dim=2)
    lsm.remove_coord('grid_latitude')
    lsm.add_dim_coord(lat, data_dim=0)
    lsm.remove_coord('grid_longitude')
    lsm.add_dim_coord(lon, data_dim=1)
    P_ma = np.ma.masked_where(orog.data > 5, P.data)
    return u_wind[0,:,:], v_wind[0,:,:], T_air[0,:,:], T_surf[0,:,:], P_ma[0,:,:], lsm, real_lat, real_lon

u_wind, v_wind, T_air, T_surf, P, lsm, real_lat, real_lon = load_files()

def shiftedColorMap(cmap, min_val, max_val, name, var):
    epsilon = 0.001
    start, stop =  0.00, 1.0 #0.15, 0.85
    min_val, max_val =  min(0.0, min_val), max(0.0, max_val)
    midpoint = 1.0 - max_val / (max_val + abs(min_val))
    cdict = {'red': [], 'green': [], 'blue': [], 'alpha': []}
    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)
    # shifted index to match the data
    shift_index = np.hstack([np.linspace(0.0, midpoint, 128, endpoint=False), np.linspace(midpoint, 1.0, 129, endpoint=True)])
    for ri, si in zip(reg_index, shift_index):
        if abs(si - midpoint) < epsilon:
            r, g, b, a = cmap(0.5) # 0.5 = original midpoint.
        else:
            r, g, b, a = cmap(ri)
        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))
    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)
    return newcmap

## Caption:

def plot_synop():
    fig = plt.figure(figsize=(10, 12))
    ax = fig.add_axes([0.18, 0.25, 0.75, 0.63], frameon=False)#, projection=ccrs.PlateCarree())#
    ax.tick_params(which='both', axis='both', labelsize=34, labelcolor='dimgrey', pad=10, size=0, tick1On=False, tick2On=False)
    PlotLonMin = np.min(real_lon)
    PlotLonMax = np.max(real_lon)
    PlotLatMin = np.min(real_lat)
    PlotLatMax = np.max(real_lat)
    XTicks = np.linspace(PlotLonMin, PlotLonMax, 3)
    XTickLabels = [None] * len(XTicks)
    for i, XTick in enumerate(XTicks):
        if XTick < 0:
            XTickLabels[i] = '{:.0f}{:s}'.format(np.abs(XTick), '$^{\circ}$W')
        else:
            XTickLabels[i] = '{:.0f}{:s}'.format(np.abs(XTick), '$^{\circ}$E')
    plt.xticks(XTicks, XTickLabels)
    ax.set_xlim(PlotLonMin, PlotLonMax)
    ax.tick_params(which='both', pad=10, labelsize = 34, color = 'dimgrey')
    YTicks = np.linspace(PlotLatMin, PlotLatMax, 4)
    YTickLabels = [None] * len(YTicks)
    for i, YTick in enumerate(YTicks):
        if YTick < 0:
            YTickLabels[i] = '{:.0f}{:s}'.format(np.abs(YTick), '$^{\circ}$S')
        else:
            YTickLabels[i] = '{:.0f}{:s}'.format(np.abs(YTick), '$^{\circ}$N')
    plt.yticks(YTicks, YTickLabels)
    ax.set_ylim(PlotLatMin, PlotLatMax)
    plt.contour(lsm.coord('longitude').points, lsm.coord('latitude').points,lsm.data, colors='0.3', linewidths=2.5, zorder=3)#, transform=RH.coord('grid_longitude').coord_system.as_cartopy_projection())
    P_lev = plt.contour(lsm.coord('longitude').points, lsm.coord('latitude').points, P_ma, colors='0.3', linewidths=3, levels = range(960,1000,4), zorder=4)
    ax.clabel(P_lev, v = [960,968,976,982,990], inline=True, inline_spacing  = 3, fontsize = 28, fmt = '%1.0f')
    bwr_zero = shiftedColorMap(cmap=matplotlib.cm.bwr, min_val=-30, max_val=0, name='bwr_zero', var=T_air.data)
    c = ax.pcolormesh(real_lon, real_lat, T_air.data, cmap=bwr_zero, vmin=-30, vmax=0, zorder=1)  #
    CBarXTicks = [-30, -20, -10, 0]  # CLevs[np.arange(0,len(CLevs),int(np.ceil(len(CLevs)/5.)))]
    CBAxes = fig.add_axes([0.25, 0.15, 0.6, 0.03])
    CBar = plt.colorbar(c, cax=CBAxes, orientation='horizontal', ticks=CBarXTicks)  #
    CBar.set_label('1.5 m air temperature ($^{\circ}$C)', fontsize=34, labelpad=10, color = 'dimgrey')
    CBar.solids.set_edgecolor("face")
    CBar.outline.set_edgecolor('dimgrey')
    CBar.ax.tick_params(which='both', axis='both', labelsize=34, labelcolor='dimgrey', pad=10, size=0, tick1On=False, tick2On=False)
    CBar.outline.set_linewidth(2)
    x, y = np.meshgrid(real_lon, real_lat)
    q = ax.quiver(x[::25, ::25], y[::25, ::25], u_wind.data[::25, ::25], v_wind.data[::25, ::25], color = '0.5', pivot='middle', scale=100, zorder = 5)
    plt.quiverkey(q, 0.25, 0.9, 10, r'$10$ $m$ $s^{-1}$', labelpos='N', color = '0.3', labelcolor = 'dimgrey',
                  fontproperties={'size': '32',  'weight': 'bold'},
                  coordinates='figure', )
    plt.draw()
    plt.savefig('/users/ellgil82/figures/Wintertime melt/Re-runs/synop_cond_CS2.png', transparent = True)
    plt.savefig('/users/ellgil82/figures/Wintertime melt/Re-runs/synop_cond_CS2.eps', transparent=True)
    plt.show()

plot_synop()
