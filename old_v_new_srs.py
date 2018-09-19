# coding= ASCII
""" Script for calculating differences between one UM run and another, using different ancillary files (land-sea mask
and orography) and plotting snapshot maps of near-surface temperature and wind vectors to visualise this.

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
import matplotlib.dates as mdates
from matplotlib.lines import Line2D
import sys
reload(sys)
sys.getdefaultencoding()
from matplotlib import rcParams
import matplotlib
import numpy.ma as ma

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

## Set up plotting options
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Helvetica', 'Liberation sans', 'Tahoma', 'DejaVu Sans','Verdana']

## Set-up cases
case = 'CS2' # string of case study in the format 'CS' + number, e.g. 'CS1'

# Make sure Python is looking in the right place for files
if case == 'CS1':
    os.chdir('/data/clivarm/wip/ellgil82/May_2016/Re-runs/CS1/') # path to data
    filepath = '/data/clivarm/wip/ellgil82/May_2016/Re-runs/CS1/'
    f_old = '/data/clivarm/wip/ellgil82/May_2016/Compare/CS1/km1p5/'
    f_new = '/data/clivarm/wip/ellgil82/May_2016/Re-runs/CS1/'
    AWS_idx = (25404,25812) # Indices of corresponding times in AWS data (hacky workaround)
    i = np.arange(16)
elif case == 'CS2':
    os.chdir('/data/clivarm/wip/ellgil82/May_2016/Re-runs/CS2/')
    filepath = '/data/clivarm/wip/ellgil82/May_2016/Re-runs/CS2/'
    f_old = '/data/clivarm/wip/ellgil82/May_2016/Compare/CS2/km1p5/'
    f_new = '/data/clivarm/wip/ellgil82/May_2016/Re-runs/CS2/'
    AWS_idx = (26124,26508)
    i = np.arange(16)

def find_gridbox(x, y, real_lat, real_lon): # Finds the indices of the inputted lat/lon coordinates
    global lon_index, lat_index
    lat_index = np.argmin((real_lat - x) ** 2)  # take whole array and subtract lat you want from
    lon_index = np.argmin((real_lon - y) ** 2)  # each point, then find the smallest difference
    return lon_index, lat_index

def construct_srs(var_name):
    i = np.arange(16)
    k = var_name.data
    series = []
    for j in i:
        a = k[:12, j]
        a = np.array(a)
        series = np.append(series, a)
    return series

def construct_Timesrs(var_name):
    i = np.arange(16)
    m = var_name
    x = []
    for j in i:
        b = m[j, :12]
        x = np.append(x, b)
    return x

def load_surf():
    old = []
    new = []
    for file in os.listdir(f_old):
        if fnmatch.fnmatch(file, '*km1p5_ctrl_pa012.pp'):
            old.append(file)
    for file in os.listdir(f_new):
        if fnmatch.fnmatch(file, '*km1p5_Smith_tnuc_pa012.pp'):
            new.append(file)
    print ('\n importing cubes...')
    os.chdir(f_old)
    Ts_old = iris.load_cube(old, 'surface_temperature')
    Ta_old = iris.load_cube(old, 'air_temperature')
    u_old = iris.load_cube(old, 'x_wind')
    v_old = iris.load_cube(old, 'y_wind')
    v_old = v_old[:,1:,:]
    RH_old = iris.load_cube(old, 'relative_humidity')
    Ts_old.convert_units('celsius')
    Ta_old.convert_units('celsius')
    Time = RH_old.coord('time')
    Time_old = Time.units.num2date(Time.points)
    #Time_old = construct_Timesrs(Time_old)
    os.chdir(f_new)
    Ta_new = iris.load_cube(new, 'air_temperature')
    Ts_new = iris.load_cube(new, 'surface_temperature')
    Ts_new.convert_units('celsius')
    Ta_new.convert_units('celsius')
    u_new = iris.load_cube(new, 'x_wind')
    v_new = iris.load_cube(new, 'y_wind')
    v_new = v_new[:,1:,:]
    RH_new = iris.load_cube(new, 'relative_humidity')
    old_var = [Ts_old,  Ta_old,  u_old, v_old,  RH_old]
    new_var = [Ts_new, Ta_new, u_new, v_new, RH_new]
    ## Rotate projection
    print ('\n rotating pole...')
    #create numpy arrays of coordinates
    rotated_lat = Ts_new.coord('grid_latitude').points
    rotated_lon = Ts_new.coord('grid_longitude').points
    ## set up parameters for rotated projection
    pole_lon = 298.5
    pole_lat = 22.99
    #rotate projection
    real_lon, real_lat = iris.analysis.cartography.unrotate_pole(rotated_lon,rotated_lat, pole_lon, pole_lat)
    print ('\nunrotating pole...')
    lat = Ts_new.coord('grid_latitude')
    lon = Ts_new.coord('grid_longitude')
    lat = iris.coords.DimCoord(real_lat, standard_name='latitude',long_name="grid_latitude",var_name="lat",units=lat.units)
    lon = iris.coords.DimCoord(real_lon, standard_name='longitude',long_name="grid_longitude",var_name="lon",units=lon.units)
    for var in old_var:
        var.remove_coord('grid_latitude')
        var.add_dim_coord(lat, data_dim=1)
        var.remove_coord('grid_longitude')
        var.add_dim_coord(lon, data_dim=2)
    for var in new_var:
        var.remove_coord('grid_latitude')
        var.add_dim_coord(lat, data_dim=1)
        var.remove_coord('grid_longitude')
        var.add_dim_coord(lon, data_dim=2)
    ## Find the nearest grid box to the latitude of interest
    print ('\n finding AWS...')
    lon_index, lat_index = find_gridbox(-66.48272, -63.37105, real_lat, real_lon)
    print('\n converting time units...')
    #convert units within iris
    Time = Ta_new.coord('time')
    Time_new = Time.units.num2date(Time.points)
    print('\n extracting time series from cubes...')
    #just one grid box
    ## Construct time series for variables with forecast_period and reference_time
    Ts_new = Ts_new[:,lat_index,lon_index].data
    Ta_new = Ta_new[:,lat_index,lon_index].data
    RH_new = RH_new[:,lat_index,lon_index].data
    U = u_new[:,lat_index,lon_index]
    V = v_new[:,lat_index,lon_index]
    v_CI = (V.data)
    u_CI = (U.data)
    sp_new = np.sqrt((u_CI**2)+(v_CI**2))
    os.chdir(f_old)
    Ts_old = Ts_old[:, lat_index, lon_index]
    Ta_old = Ta_old[:, lat_index, lon_index]
    RH_old = RH_old[:, lat_index,lon_index]
    U = u_old[:, lat_index,lon_index]
    V = v_old[:, lat_index,lon_index]
    #Ts_old = construct_srs(Ts_old)
    #Ta_old = construct_srs(Ta_old)
    #RH_old = construct_srs(RH_old)
    #U = construct_srs(U)
    #V = construct_srs(V)
    print ('\n calculating wind speed...')
    ##convert u and v wind to wind speed
    #convert to numpy array
    v_CI = (V)
    u_CI = (U)
    sp_old = np.sqrt((u_CI**2)+(v_CI**2))
    #sp_old = 'nope'
    return  Ts_new, Ts_old, Ta_new, Ta_old, RH_new, RH_old,  Time_new, Time_old, sp_new, sp_old

#Ts_new, Ts_old, Ta_new, Ta_old, RH_new, RH_old,  Time_new, Time_old, sp_new, sp_old = load_surf()

def load_AWS():
    print '\nimporting AWS observations...'
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from numpy import loadtxt
    #Load AWS data and create variables
    fAWS = '/data/clivarm/wip/ellgil82/May_2016/CS2/data/AWS_data.txt'
    AWS_srs = loadtxt (fAWS, skiprows=1)
    print '\nsubsetting for Case Study 2...'
    # Create subset for Case Study 2
    CS = AWS_srs[AWS_idx[0]:AWS_idx[1],:]
    AWS_time = loadtxt(fAWS, usecols=(0, 1, 2), skiprows=1)
    CS_time = AWS_time[AWS_idx[0]:AWS_idx[1],:]
    #make the times match the model data
    RH_CS = CS[:,5]
    Tair_CS = CS[:,4]
    Ts_CS = CS[:,-1]
    wind_CS = CS[:,8]
    Time_CS = CS[:,3]
    SW_CS = CS[:,11]
    LW_CS = CS[:,14]
    LH_CS = CS[:,-4]
    SH_CS = CS[:,-5]
    G_CS = CS[:,-3]
    melt_CS = CS[:,-2]
    print '\nconverting times...'
    # Convert times so that they can be plotted
    Time_list = np.empty(1,)
    from datetime import datetime, timedelta
    for each_day in Time_CS:
        epoch = datetime(datetime.now().year-3, 12, 31)
        result = epoch + timedelta(days=each_day)
        Time_list = np.append(Time_list, result)
    Time_list = np.array(Time_list)
    Time_list = np.delete(Time_list,[0])
    return Time_list, melt_CS, SH_CS, LH_CS, LW_CS, SW_CS, Time_CS, RH_CS, Ts_CS, Tair_CS, wind_CS,

#Time_list,melt_CS, SH_CS, LH_CS, LW_CS, SW_CS, Time_CS, RH_CS,Ts_CS, Tair_CS, wind_CS, = load_AWS()

# Old: '#6a3d9a', New: '#33a02c'

def surf_plot():
    fig, ax = plt.subplots(4,1,sharex= True, figsize=(18, 16))
    ax = ax.flatten()
    Time_list, melt_CS, SH_CS, LH_CS, LW_CS, SW_CS, Time_CS, RH_CS, Ts_CS, Tair_CS, wind_CS, = load_AWS()
    Ts_new, Ts_old, Ta_new, Ta_old, RH_new, RH_old, Time_new, Time_old, sp_new, sp_old= load_surf()
    days = mdates.DayLocator(interval=1)
    dayfmt = mdates.DateFormatter('%d %b')
    for axs in ax:
        axs.spines['top'].set_visible(False)
        plt.setp(axs.spines.values(), linewidth=2, color='dimgrey', )
        axs.set_xlim(Time_list[1], Time_list[-1])
        axs.tick_params(axis='both', which='both', labelsize=24, tick1On=False, tick2On=False, labelcolor='dimgrey', pad=10)
        [l.set_visible(False) for (w, l) in enumerate(axs.yaxis.get_ticklabels()) if w % 2 != 0]
        [l.set_visible(False) for (w, l) in enumerate(axs.xaxis.get_ticklabels()) if w % 2 != 0]
    # Surface temperature
    ax[0].plot(Time_list, Ts_CS, color='k', linewidth=2.5, label="Cabinet Inlet AWS")
    ax2 = ax[0].twiny()
    ax3 = ax2.twinx()
    ax[0].spines['right'].set_visible(False)
    ax2.plot(Time_new, Ts_new, linewidth=2.5, color='#33a02c',  label='Updated set-up')
    ax3.plot(Time_old, Ts_old, linewidth=2.5, color='#6a3d9a', label="Default set-up")
    ax2.axis('off')
    ax3.axis('off')
    ax2.set_xlim(Time_new[1], Time_new[-1])
    ax2.tick_params(axis='both', tick1On=False, tick2On=False)
    ax3.tick_params(axis='both', tick1On=False, tick2On=False)
    ax2.set_ylim(-30, 15)
    ax3.set_ylim(-30, 15)
    ax[0].set_ylabel('Surface \nTemperature \n($^{\circ}$C)', rotation=0, fontsize=24, color='dimgrey',
                     labelpad=80)
    ax[0].tick_params(axis='both', which='both', labelsize=24, tick1On=False, tick2On=False)
    plt.setp(ax[0].get_yticklabels()[-1], visible=False)
    plt.setp(ax[0].get_yticklabels()[-3], visible=False)
    # Air temperature
    ax[1].plot(Time_list, Tair_CS, color='k', linewidth=2.5, label="Cabinet Inlet AWS")
    ax2 = ax[1].twiny()
    ax3 = ax2.twinx()
    ax[1].spines['right'].set_visible(False)
    ax2.plot(Time_new, Ta_new, linewidth=2.5, color = '#33a02c', lw = 2.5, label = 'Updated set-up')
    ax3.plot(Time_old, Ta_old, linewidth = 2.5, color = '#6a3d9a', label = "Default set-up")
    ax2.axis('off')
    ax3.axis('off')
    ax2.set_xlim(Time_new[1], Time_new[-1])
    ax2.tick_params(axis='both', tick1On = False, tick2On = False)
    ax3.tick_params(axis='both', tick1On = False, tick2On = False)
    ax2.set_ylim(-30,15)
    ax3.set_ylim(-30,15)
    ax[1].set_ylabel('Near-surface \nAir Temperature \n($^{\circ}$C)', rotation=0, fontsize=24, color='dimgrey',labelpad=80)
    ax[1].tick_params(axis='both', which='both', labelsize=24, tick1On=False, tick2On=False)
    plt.setp(ax[1].get_yticklabels()[-1], visible=False)
    plt.setp(ax[1].get_yticklabels()[-3], visible=False)
    # Relative humidity
    ax[2].plot(Time_list, RH_CS, color='k', linewidth=2.5, label="Cabinet Inlet AWS")
    ax2 = ax[2].twiny()
    ax3 = ax2.twinx()
    ax[2].spines['right'].set_visible(False)
    ax2.plot(Time_new, RH_new, linewidth=2.5, color='#33a02c', lw=2.5, label='Updated set-up')
    ax3.plot(Time_old, RH_old, linewidth=2.5, color='#6a3d9a', label="Default set-up")
    ax2.yaxis.set_label_coords(-0.35, 0.5)
    ax2.set_xlim(Time_new[1], Time_new[-1])
    ax[2].set_ylim([0, 100])
    ax2.set_ylim(0,100)
    ax3.set_ylim(0,100)
    ax2.axis('off')
    ax3.axis('off')
    ax2.tick_params(axis='both', tick1On = False, tick2On = False)
    ax3.tick_params(axis='both', tick1On = False, tick2On = False)
    ax[2].xaxis.set_major_formatter(dayfmt)
    ax[2].set_ylabel('Relative \nHumidity \n(%)', rotation=0, fontsize=24, color = 'dimgrey', labelpad = 80)
    ax[2].tick_params(axis='both', which='both', labelsize=24, tick1On = False, tick2On = False)
    # Wind speed
    ax[3].plot(Time_list, wind_CS, color='k', linewidth=2.5, label="Cabinet Inlet AWS")
    ax2 = ax[3].twiny()
    ax3 = ax2.twinx()
    ax[3].spines['right'].set_visible(False)
    ax2.plot(Time_new, sp_new, linewidth=2.5, color='#33a02c', lw=2.5, label='Updated set-up')
    ax3.plot(Time_old, sp_old, linewidth=2.5, color='#6a3d9a', label="Default set-up")
    ax2.yaxis.set_label_coords(-0.35, 0.5)
    ax2.set_xlim(Time_new[1], Time_new[-1])
    ax[3].set_ylim([0, 40])
    ax2.set_ylim(0,40)
    ax3.set_ylim(0,40)
    ax2.axis('off')
    ax3.axis('off')
    ax2.tick_params(axis='both', tick1On = False, tick2On = False)
    ax3.tick_params(axis='both', tick1On = False, tick2On = False)
    ax[3].set_ylabel('Wind \nSpeed \n(m s$^{-1}$)', rotation=0, fontsize=24, color = 'dimgrey', labelpad = 80)
    ax[3].tick_params(axis='both', which='both', labelsize=24, tick1On = False, tick2On = False)
    ax[3].xaxis.set_major_formatter(dayfmt)
    plt.setp(ax[3].get_yticklabels()[-2], visible=False)
   # Legend
    lns = [Line2D([0], [0], color='k', linewidth=2.5),
           Line2D([0], [0], color='#33a02c', linewidth=2.5),
           Line2D([0], [0], color='#6a3d9a', linewidth=2.5)]
    labs = ['Observations from \nCabinet Inlet', 'Updated set-up', 'Default set-up']
    lgd = ax[3].legend(lns, labs, bbox_to_anchor=(0.88, 1.22), loc=2, fontsize=20)
    for ln in lgd.get_texts():
        plt.setp(ln, color='dimgrey')
    lgd.get_frame().set_linewidth(0.0)
    plt.subplots_adjust(wspace = 0.13, hspace = 0.13, top = 0.98, right = 0.9, left = 0.18, bottom = 0.08)
    plt.savefig('/users/ellgil82/figures/Wintertime melt/Re-runs/'+case+'_old_v_new_srs.png', transparent = True)
    plt.savefig('/users/ellgil82/figures/Wintertime melt/Re-runs/'+case+'_old_v_new_srs.eps', transparent = True)
    plt.show()

surf_plot()