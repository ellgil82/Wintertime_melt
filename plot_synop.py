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
from tools import rotate_data

def load_files(case):
    if case == 'CS1':
        os.chdir('/data/clivarm/wip/ellgil82/May_2016/Re-runs/CS1/')  # path to data
        filepath = '/data/clivarm/wip/ellgil82/May_2016/Re-runs/CS1/'
        case_date = '20160507T1200Z'
    elif case == 'CS2':
        os.chdir('/data/clivarm/wip/ellgil82/May_2016/Re-runs/CS2/')
        filepath = '/data/clivarm/wip/ellgil82/May_2016/Re-runs/CS2/'
        case_date = '20160522T1200Z'
    os.chdir(filepath)
    surf = filepath+case_date+'_Peninsula_km4p0_ctrl_pa000.pp'
    print ('importing cubes...')
    T_air = iris.load_cube(surf, 'air_temperature')
    T_surf = iris.load_cube(surf, 'surface_temperature')
    lsm = iris.load_cube(surf, 'land_binary_mask')
    orog = iris.load_cube(surf, 'surface_altitude')
    P = iris.load_cube(surf, 'surface_air_pressure')
    P.convert_units('hPa')
    T_air.convert_units('celsius')
    T_surf.convert_units('celsius')
    ## Iris v1.11 version
    u_wind = iris.load_cube(surf, 'x_wind')
    v_wind = iris.load_cube(surf, 'y_wind')
    v_wind = v_wind[:,1:,:]
    ## Rotate projection
    for var in [T_air, T_surf, u_wind, v_wind, P]:
        real_lon, real_lat = rotate_data(var,1,2)
    for var in [lsm, orog]:
        real_lon, real_lat = rotate_data(var, 0, 1)
    P_ma = np.ma.masked_where(lsm.data == 1, P[0,:,:].data) #orog.data > 0 &
    return u_wind[0,:,:], v_wind[0,:,:], T_air[0,:,:], T_surf[0,:,:], P_ma, lsm, real_lat, real_lon, orog

u_wind, v_wind, T_air, T_surf, P, lsm, real_lat, real_lon, orog = load_files('CS1')

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
    plt.contour(lsm.coord('longitude').points, lsm.coord('latitude').points,lsm.data, colors='#757676', linewidths=2.5, zorder=2)#, transform=RH.coord('grid_longitude').coord_system.as_cartopy_projection())
    plt.contour(orog.coord('longitude').points, orog.coord('latitude').points, orog.data, colors='#757676', levels=[15], linewidths=2.5, zorder=3)
    P_lev = plt.contour(lsm.coord('longitude').points, lsm.coord('latitude').points, P, colors='#222222', linewidths=3, levels = range(960,1020,4), zorder=4)
    ax.clabel(P_lev, v = [960,968,976,982,990], inline=True, inline_spacing  = 3, fontsize = 28, fmt = '%1.0f')
    bwr_zero = shiftedColorMap(cmap=matplotlib.cm.bwr, min_val=-30, max_val=10, name='bwr_zero', var=T_air.data)
    c = ax.pcolormesh(real_lon, real_lat, T_air.data, cmap=bwr_zero, vmin=-30, vmax=10, zorder=1)  #
    CBarXTicks = [-30,  -10, 10]  # CLevs[np.arange(0,len(CLevs),int(np.ceil(len(CLevs)/5.)))]
    CBAxes = fig.add_axes([0.25, 0.15, 0.6, 0.03])
    CBar = plt.colorbar(c, cax=CBAxes, orientation='horizontal', ticks=CBarXTicks)  #
    CBar.set_label('1.5 m air temperature ($^{\circ}$C)', fontsize=34, labelpad=10, color = 'dimgrey')
    CBar.solids.set_edgecolor("face")
    CBar.outline.set_edgecolor('dimgrey')
    CBar.ax.tick_params(which='both', axis='both', labelsize=34, labelcolor='dimgrey', pad=10, size=0, tick1On=False, tick2On=False)
    CBar.outline.set_linewidth(2)
    x, y = np.meshgrid(real_lon, real_lat)
    q = ax.quiver(x[::25, ::25], y[::25, ::25], u_wind.data[::25, ::25], v_wind.data[::25, ::25], color = '#414345', pivot='middle', scale=100, zorder = 5)
    plt.quiverkey(q, 0.25, 0.9, 10, r'$10$ $m$ $s^{-1}$', labelpos='N', color = '#414345', labelcolor = '#414345',
                  fontproperties={'size': '32',  'weight': 'bold'},
                  coordinates='figure', )
    plt.draw()
    plt.savefig('/users/ellgil82/figures/Wintertime melt/Re-runs/synop_cond_CS1.png', transparent = True)
    plt.savefig('/users/ellgil82/figures/Wintertime melt/Re-runs/synop_cond_CS1.eps', transparent=True)
    plt.show()

plot_synop()

def plot_both():
    fig, axs = plt.subplots(1,2, figsize=(20, 12), frameon=False)
    lab_dict = {0: 'a', 1: 'b'}
    #ax = fig.add_axes([0.18, 0.25, 0.75, 0.63], frameon=False) # , projection=ccrs.PlateCarree())#
    case_list = ['CS1', 'CS2']
    for a in [0,1]:
        axs[a].spines['right'].set_visible(False)
        axs[a].spines['left'].set_visible(False)
        axs[a].spines['top'].set_visible(False)
        axs[a].spines['bottom'].set_visible(False)
        u_wind, v_wind, T_air, T_surf, P, lsm, real_lat, real_lon, orog = load_files(case_list[a])
        axs[a].contour(lsm.coord('longitude').points, lsm.coord('latitude').points, lsm.data, colors='#757676', linewidths=2.5,
                zorder=2)  # , transform=RH.coord('grid_longitude').coord_system.as_cartopy_projection())
        axs[a].contour(orog.coord('longitude').points, orog.coord('latitude').points, orog.data, colors='#757676', levels=[15],
                linewidths=2.5, zorder=3)
        P_lev = axs[a].contour(lsm.coord('longitude').points, lsm.coord('latitude').points, P, colors='#222222', linewidths=3,
                        levels=range(960, 1020, 4), zorder=4)
        axs[a].clabel(P_lev, v=[960, 968, 976, 982, 990], inline=True, inline_spacing=3, fontsize=28, fmt='%1.0f')
        bwr_zero = shiftedColorMap(cmap=matplotlib.cm.bwr, min_val=-30, max_val=10, name='bwr_zero', var=T_air.data)
        c = axs[a].pcolormesh(real_lon, real_lat, T_air.data, cmap=bwr_zero, vmin=-30, vmax=10, zorder=1)  #
        x, y = np.meshgrid(real_lon, real_lat)
        q = axs[a].quiver(x[::25, ::25], y[::25, ::25], u_wind.data[::25, ::25], v_wind.data[::25, ::25], color='#414345',
                      pivot='middle', scale=100, zorder=5)
        axs[a].tick_params(which='both', axis='both', labelsize=34, labelcolor='dimgrey', pad=10, size=0, tick1On=False, tick2On=False)
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
        plt.sca(axs[a])
        plt.xticks(XTicks, XTickLabels)
        axs[a].set_xlim(PlotLonMin, PlotLonMax)
        axs[a].tick_params(which='both', pad=10, labelsize=34, color='dimgrey')
        YTicks = np.linspace(PlotLatMin, PlotLatMax, 4)
        YTickLabels = [None] * len(YTicks)
        for i, YTick in enumerate(YTicks):
            if YTick < 0:
                YTickLabels[i] = '{:.0f}{:s}'.format(np.abs(YTick), '$^{\circ}$S')
            else:
                YTickLabels[i] = '{:.0f}{:s}'.format(np.abs(YTick), '$^{\circ}$N')
        plt.sca(axs[a])
        plt.yticks(YTicks, YTickLabels)
        axs[a].set_ylim(PlotLatMin, PlotLatMax)
        lab = axs[a].text(-80, -61.5, zorder=100,  s=lab_dict[a], fontsize=32, fontweight='bold', color='dimgrey')
        axs[a].set_title(case_list[a], fontsize=34, color='dimgrey')
    CBarXTicks = [-30, -10, 10]  # CLevs[np.arange(0,len(CLevs),int(np.ceil(len(CLevs)/5.)))]
    CBAxes = fig.add_axes([0.35, 0.15, 0.3, 0.03])
    CBar = plt.colorbar(c, cax=CBAxes, orientation='horizontal', ticks=CBarXTicks)  #
    CBar.set_label('1.5 m air temperature ($^{\circ}$C)', fontsize=34, labelpad=10, color='dimgrey')
    CBar.solids.set_edgecolor("face")
    CBar.outline.set_edgecolor('dimgrey')
    CBar.ax.tick_params(which='both', axis='both', labelsize=34, labelcolor='dimgrey', pad=10, size=0, tick1On=False,
                        tick2On=False)
    CBar.outline.set_linewidth(2)
    plt.sca(axs[1])
    plt.tick_params(axis = 'y', which=  'both', labelleft = 'off', labelright = 'on')
    #yaxis.set_label_coords(1.27, 0.5)
    plt.quiverkey(q, 0.51, 0.9, 10, r'$10$ $m$ $s^{-1}$', labelpos='N', color='#414345', labelcolor='#414345',
                  fontproperties={'size': '32', 'weight': 'bold'},
                  coordinates='figure', )
    plt.draw()
    plt.subplots_adjust(bottom = 0.25, top = 0.85)
    plt.savefig('/users/ellgil82/figures/Wintertime melt/Re-runs/synop_cond_both_cases.png', transparent=True)
    plt.savefig('/users/ellgil82/figures/Wintertime melt/Re-runs/synop_cond_both_cases.eps', transparent=True)
    plt.show()

plot_both()