# coding= ASCII
""" Script for calculating and plotting time series of meteorological variables and surface energy fluxes from model
data and automatic weather station observations.

This can be adapted to be include observations from any location within the model domain, and data can be manipulated
to include moving time-window averaging, mean statistics etc.

Dependencies:
- iris 1.11.0
- matplotlib 1.5.1
- numpy 1.10.4
- rotate_data.py script in same folder

Author: Ella Gilbert, 2017. Updated April 2018.

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
import numpy.ma as ma
import scipy
import pandas as pd
from tools import rotate_data, find_gridbox, compose_date, compose_time
import pandas as pd
from numpy import loadtxt

## Set-up cases
case = 'CS1' # string of case study in the format 'CS' + number, e.g. 'CS1'
#res = 'km4p0' # string of resolution to match filename, e.g. 'km4p0'

# Make sure Python is looking in the right place for files
if case == 'CS1':
    os.chdir('/data/clivarm/wip/ellgil82/May_2016/Re-runs/CS1/') # path to data
    filepath = '/data/clivarm/wip/ellgil82/May_2016/Re-runs/CS1/'
    case_start = '2016-05-08' # Must be in datetime format as a string, e.g. '2016-05-08'
    case_end = '2016-05-13' 
    #AWS_idx = (25404,25812) # Indices of corresponding times in AWS data (hacky workaround)
    res_list = [ 'km1p5', 'km4p0'] # List of model resolutions you want to process
elif case == 'CS2':
    os.chdir('/data/clivarm/wip/ellgil82/May_2016/Re-runs/CS2/')
    filepath = '/data/clivarm/wip/ellgil82/May_2016/Re-runs/CS2/'
    case_start = '2016-05-22' # Must be in datetime format as a string, e.g. '2016-05-08'
    case_end = '2016-05-30' 
    #AWS_idx = (26124,26508)
    res_list = ['km1p5', 'km4p0'] 

def load_AWS():
    print '\nimporting AWS observations...'
    # Load AWS data and create variables
    AWS_srs = np.genfromtxt ('/data/clivarm/wip/ellgil82/AWS/Cabinet_Inlet_AWS18.txt', names = True)
    AWS_srs = pd.DataFrame(AWS_srs) # Convert to pandas DataFrame this way because it loads in incorrectly using pd.from_csv
    print '\nsubsetting for Case Study 2...'
    AWS_srs.index = compose_date(AWS_srs['Year'], days=AWS_srs['day'])
    AWS_srs['Time'] = 24*(AWS_srs['Time'] - AWS_srs['day'])
    # Create subset for Case Study 
    case = AWS_srs.loc[case_start:case_end]
    print '\nconverting times...'
    # Convert times so that they can be plotted
    Time_list = np.array(Time_list)
    Time_list = np.delete(Time_list,[0])
    E = AWS_srs['LWnet'] + AWS_srs['SWnet'] + AWS_srs['Hlat'] + AWS_srs['Hsen'] - AWS_srs['Gs']
    var_dict = {
    'Time_list': Time_list, 
    'E': E,
    'AWS_obs': AWS_srs
    }
    return var_dict

AWS_var = load_AWS()

## Function to load in SEB data from UM. Make sure the file names point to the correct file stream where your variables are stored. 
## This can be adapted to use other formats, e.g. NetCDF, GRIB etc. (see Iris docs for further information: https://scitools.org.uk/iris/docs/latest/#)

def load_SEB(res):
    surf = []
    SEB = []
    for file in os.listdir(filepath):
        if fnmatch.fnmatch(file, '*%(res)s_*_pa012.pp' % locals()):
            surf.append(file)
        elif fnmatch.fnmatch(file, '*%(res)s_*_pb012.pp' % locals()):
            SEB.append(file)
    print ('\n importing cubes at %(res)s resolution...' % locals())
    os.chdir(filepath)
    print ('\n Downwelling shortwave...')
    SW_d = iris.load_cube(SEB, 'surface_downwelling_shortwave_flux_in_air')
    print('\n Downwelling longwave...')
    LW_d = iris.load_cube(SEB, 'surface_downwelling_longwave_flux')
    print('\n Net shortwave...')
    SW_n = iris.load_cube(SEB, 'surface_net_downward_shortwave_flux')
    print('\n Net longwave...')
    LW_n = iris.load_cube(SEB, 'surface_net_downward_longwave_flux')
    print('\n Latent heat...')
    LH = iris.load_cube(SEB, 'surface_upward_latent_heat_flux')
    print('\n Sensible heat...')
    SH = iris.load_cube(SEB, 'surface_upward_sensible_heat_flux')
    print('\n Surface temperature...')
    T_surf = iris.load_cube(surf, 'surface_temperature')
    T_surf.convert_units('celsius')
    Var = [SH, LH, LW_d, SW_d, LW_n, SW_n, T_surf]
    ## Rotate projection
    print ('\n rotating pole...')
    for var in Var:
        real_lon, real_lat = rotate_data(var, 1, 2)
    ## Find the nearest grid box to the latitude of interest
    print ('\n finding AWS...')
    lon_index, lat_index = find_gridbox(-66.48272, -63.37105, real_lat, real_lon)
    print('\n converting time units...')
    #convert units within iris
    Time = SH.coord('time')
    Time_srs = Time.units.num2date(Time.points)
    # Create Larsen mask to return values only a) on the ice shelf, b) orography is < 15 m
    orog = iris.load_cube(filepath + res + '_orog.pp', 'surface_altitude')
    lsm = iris.load_cube(filepath + res + '_lsm.pp', 'land_binary_mask')
    a = np.ones((orog.shape))
    orog = orog[:270, 95:240].data
    lsm = lsm[:270, 95:240].data
    b = np.zeros((270, 95))
    c = np.zeros((270, 160))
    d = np.zeros((130, 400))
    orog = np.hstack((b, orog))
    orog = np.hstack((orog, c))
    orog = np.vstack((orog, d))
    lsm = np.hstack((b, lsm))
    lsm = np.hstack((lsm, c))
    lsm = np.vstack((lsm, d))
    mask2d = np.ma.masked_where(orog.data > 15, a)
    Larsen_mask = np.ma.masked_where(lsm.data == 0, mask2d)
    Larsen_mask = np.broadcast_to(Larsen_mask == 1, T_surf.shape, subok=True)
    SH = np.ma.masked_array(SH.data, Larsen_mask.mask)
    LH = np.ma.masked_array(LH.data, Larsen_mask.mask)
    SW_d = np.ma.masked_array(SW_d.data, Larsen_mask.mask)
    LW_d = np.ma.masked_array(LW_d.data, Larsen_mask.mask)
    SW_n = np.ma.masked_array(SW_n.data, Larsen_mask.mask)
    LW_n = np.ma.masked_array(LW_n.data, Larsen_mask.mask)
    # Flip turbulent fluxes to match convention (positive = down)
    LH = 0 - LH
    SH = 0 - SH
    # Calculate 5th and 95th percentiles to give estimate of variability in time series
    print('\n calculating percentiles... (this may take some time)')
    percentiles = []
    for each_var in [SH, LH, SW_d, SW_n, LW_d, LW_n]:
        p95 = np.percentile(each_var, 95, axis=(1, 2))
        p5 = np.percentile(each_var, 5, axis=(1, 2))
        percentiles.append(p5)
        percentiles.append(p95)
    print('\n extracting time series from cubes...')
    # Just one grid box
    SW_d = SW_d[:,lat_index,lon_index].data
    LW_d = LW_d[:,lat_index,lon_index].data
    SW_n = SW_n[:,lat_index,lon_index].data
    LW_n = LW_n[:,lat_index,lon_index].data
    LH = LH[:,lat_index,lon_index].data
    SH = SH[:,lat_index,lon_index].data
    T_surf = T_surf[:,lat_index, lon_index].data
    print('\n constructing %(res)s series...' % locals())
    print ('\n making melt variable...')
    # Calculate total SEB (without Gs, which is unavailable in the UM)
    E = SW_n + LW_n + LH + SH
    # Create melt variable
    # Create masked array when Ts<0
    melt_masked = np.ma.masked_where(T_surf<-0.025, E)
    melt = melt_masked.clip(min=0)
    melt_forced = np.ma.masked_where(AWS_var['Ts']<-0.025, E)
    var_dict = {
    'Time_srs': Time_srs, 
    'Ts': T_surf, 
    'SW_n': SW_n, 
    'SW_d': SW_d, 
    'LW_n': LW_n, 
    'LW_d': LW_d, 
    'SH': SH, 
    'LH': LH, 
    'melt': melt_masked, 
    'melt_forced': melt_forced, 
    'percentiles': percentiles, 
    'E': E}
    return var_dict

## Function, as above, to load in surface data from the UM.

def load_surf(res):
    surf = []
    for file in os.listdir(filepath):
        if fnmatch.fnmatch(file, '*%(res)s_*_pa012.pp' % locals()):
            surf.append(file)
    print ('\n importing cubes...')
    os.chdir(filepath)
    T_surf = iris.load_cube(surf, 'surface_temperature')
    T_air = iris.load_cube(surf, 'air_temperature')
    orog = iris.load_cube(filepath+res+'_orog.pp', 'surface_altitude') # This assumes your orography and land-sea mask are stored separately,
    lsm = iris.load_cube(filepath + res + '_lsm.pp', 'land_binary_mask') # but can be adapted to read in from one of your file streams.
    T_surf.convert_units('celsius')
    T_air.convert_units('celsius')
    RH = iris.load_cube(surf, 'relative_humidity')
    ## Iris v1.11 version
    u_wind = iris.load_cube(surf, 'x_wind')
    v_wind = iris.load_cube(surf, 'y_wind')
    v_wind = v_wind[:,1:,:]
    Var = [T_surf, T_air, RH, u_wind, v_wind]
    ## Rotate projection
    print ('\n rotating pole...')
    for var in Var:
        real_lon, real_lat = rotate_data(var, 1,2)
    ## Find the nearest grid box to the latitude of interest
    print ('\n finding AWS...')
    lon_index, lat_index = find_gridbox(-66.48272, -63.37105, real_lat, real_lon)
    print('\n converting time units...')
    #convert units within iris
    Time = T_surf.coord('time')
    Time_srs = Time.units.num2date(Time.points)
    print ('\n calculating wind speed...')
    ##convert u and v wind to wind speed
    #convert to numpy array
    v_CI = (v_wind.data)
    u_CI = (u_wind.data)
    sp_srs = np.sqrt((u_CI**2)+(v_CI**2))
    # Create Larsen mask !! Is this for 4 km or 1.5 ??
    a = np.ones((orog.shape))
    orog = orog[:270,95:240].data
    lsm = lsm[:270,95:240].data
    b =  np.zeros((270,95))
    c = np.zeros((270,160))
    d = np.zeros((130,400))
    orog = np.hstack((b, orog))
    orog = np.hstack((orog, c))
    orog = np.vstack((orog,d))
    lsm = np.hstack((b, lsm))
    lsm = np.hstack((lsm, c))
    lsm = np.vstack((lsm,d))
    mask2d = np.ma.masked_where(orog.data > 15, a)
    Larsen_mask = np.ma.masked_where(lsm.data == 0, mask2d)
    Larsen_mask = np.broadcast_to(Larsen_mask == 1, T_surf.shape, subok =True)
    T_surf = np.ma.masked_array(T_surf.data, Larsen_mask.mask)
    RH = np.ma.masked_array(RH.data, Larsen_mask.mask)
    T_air = np.ma.masked_array(T_air.data, Larsen_mask.mask)
    sp_srs = np.ma.masked_array(sp_srs, Larsen_mask.mask)
    # Calculate 5th and 95th percentiles to give estimate of variability in time series
    print ('\n calculating percentiles... (this may take some time)')
    percentiles = []
    for each_var in [T_surf, T_air, RH, sp_srs]:
        p95 = np.percentile(each_var, 95, axis = (1,2))
        p5 = np.percentile(each_var, 5, axis = (1,2))
        percentiles.append(p5)
        percentiles.append(p95)
    print('\n extracting time series from cubes...')
    # just one grid box
    T_surf = T_surf[:,lat_index,lon_index].data
    T_air = T_air[:,lat_index,lon_index].data
    RH = RH[:,lat_index,lon_index].data
    RH[RH>100] = 100
    sp_srs = sp_srs[:,lat_index, lon_index]
    print('\n constructing %(res)s series...' % locals())
    var_dict = {
    'sp_srs': sp_srs, 
    'Ts': T_surf, 
    'T_air': T_air, 
    'RH': RH, 
    'Time_srs': Time_srs, 
    'percentiles': percentiles}
    return var_dict 

SEB_1p5 = load_SEB('km1p5')
surf_1p5 = load_surf('km1p5')
# SEB_4p0 = load_SEB('km4p0')
# surf_4p0 = load_surf('km4p0')

total_SEB_obs = SH_CS + LH_CS + SWn_CS + LWn_CS

## =========================================== SENSITIVITY TESTING ================================================== ##

## Does the temperature bias improve if we choose a more representative grid point?

def sens_test(res):
    surf = []
    for file in os.listdir(filepath):
        if fnmatch.fnmatch(file, '*%(res)s_*_pa012.pp' % locals()):
            surf.append(file)
    print ('\n importing cubes...')
    os.chdir(filepath)
    T_surf = iris.load_cube(surf, 'surface_temperature')
    T_air = iris.load_cube(surf, 'air_temperature')
    orog = iris.load_cube(filepath+res+'_orog.pp', 'surface_altitude')
    lsm = iris.load_cube(filepath + res + '_lsm.pp', 'land_binary_mask')
    T_surf.convert_units('celsius')
    T_air.convert_units('celsius')
    RH = iris.load_cube(surf, 'relative_humidity')
    ## Iris v1.11 version
    u_wind = iris.load_cube(surf, 'x_wind')
    v_wind = iris.load_cube(surf, 'y_wind')
    v_wind = v_wind[:,1:,:]
    Var = [T_surf, T_air, RH, u_wind, v_wind]
    ## Rotate projection
    print ('\n rotating pole...')
    #create numpy arrays of coordinates
    rotated_lat = RH.coord('grid_latitude').points
    rotated_lon = RH.coord('grid_longitude').points
    ## set up parameters for rotated projection
    pole_lon = 298.5
    pole_lat = 22.99
    #rotate projection
    real_lon, real_lat = iris.analysis.cartography.unrotate_pole(rotated_lon,rotated_lat, pole_lon, pole_lat)
    print ('\nunrotating pole...')
    lat = RH.coord('grid_latitude')
    lon = RH.coord('grid_longitude')
    lat = iris.coords.DimCoord(real_lat, standard_name='latitude',long_name="grid_latitude",var_name="lat",units=lat.units)
    lon = iris.coords.DimCoord(real_lon, standard_name='longitude',long_name="grid_longitude",var_name="lon",units=lon.units)
    for var in Var:
        var.remove_coord('grid_latitude')
        var.add_dim_coord(lat, data_dim=1)
        var.remove_coord('grid_longitude')
        var.add_dim_coord(lon, data_dim=2)
    ## Find the nearest grid box to the latitude of interest
    print ('\n finding AWS...')
    lon_index, lat_index = find_gridbox(-66.48272, -63.37105, real_lat, real_lon)
    print('\n converting time units...')
    #convert units within iris
    Time = T_surf.coord('time')
    Time_srs = Time.units.num2date(Time.points)
    #convert to numpy array
    v_CI = (v_wind.data)
    u_CI = (u_wind.data)
    sp_srs = np.sqrt((u_CI**2)+(v_CI**2))
    Ts_subset = T_surf[:, lat_index, 140 ]
    Ta_subset = T_air[:, lat_index, 140 ]
    wind_subset = sp_srs[:, lat_index, 140 ]
    return Ts_subset, Ta_subset, wind_subset

#Ts_subset, Ta_subset, wind_subset = sens_test('km1p5')

#Ts_mean = np.mean(Ts_subset.data, axis = (1,2))
#Ta_mean = np.mean(Ta_subset.data, axis = (1,2))

#plt.plot(wind_CS, color = 'k')
#plt.plot(wind_subset, color ='r')
#plt.plot(sp_srs, color ='b')
#plt.show()

## ============================================ COMPUTE STATISTICS ================================================== ##

## Set up plotting options
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Helvetica', 'Liberation sans', 'Tahoma', 'DejaVu Sans',
                               'Verdana']

#Create empty pandas dataframe to put data in
def correl_plot():
    Time_list, melt_CS, SH_CS, LH_CS, LWd_CS, SWd_CS, LWn_CS, SWn_CS, Time_CS, RH_CS, Ts_CS, Tair_CS, wind_CS = load_AWS()
    SH, LH, T_surf, Time_srs, melt, SW_n, LW_n, LW_d, SW_d, percentiles_SEB, melt_forced, E = load_SEB('km1p5')
    sp_srs, T_surf, T_air, RH, Time_srs, percentiles_surf = load_surf('km1p5')
    R_net = SW_n + LW_n
    Rnet_CS = SWn_CS + LWn_CS
    fig, ax = plt.subplots(4,2, figsize = (16,28))
    ax2 = ax[:,1]
    ax = ax.flatten()
    ax2.flatten()
    lab_dict = {0: 'a', 1: 'b', 2: 'c', 3: 'd', 4: 'e', 5: 'f', 6: 'g', 7: 'h'}
    for axs in ax:
        axs.spines['top'].set_visible(False)
        axs.spines['right'].set_visible(False)
        plt.setp(axs.spines.values(), linewidth=2, color='dimgrey', )
        axs.set(adjustable='box-forced', aspect='equal')
        axs.tick_params(axis='both', which='both', labelsize=24, tick1On=False, tick2On=False, labelcolor='dimgrey', pad=10)
        [l.set_visible(False) for (w, l) in enumerate(axs.yaxis.get_ticklabels()) if w % 2 != 0]
        [l.set_visible(False) for (w, l) in enumerate(axs.xaxis.get_ticklabels()) if w % 2 != 0]
        axs.yaxis.set_label_coords(-0.4, 0.5)
    for axs in ax2:
        axs.spines['left'].set_visible(False)
        axs.spines['right'].set_visible(True)
        axs.yaxis.set_label_position('right')
        axs.yaxis.set_ticks_position('right')
        axs.tick_params(axis='both', which='both', labelsize=24, tick1On=False, tick2On=False, labelcolor='dimgrey', pad=10)
        [l.set_visible(False) for (w, l) in enumerate(axs.yaxis.get_ticklabels()) if w % 2 != 0]
        axs.yaxis.set_label_coords(1.45, 0.57)
    plot = 0
    surf_met_mod = [T_surf, T_air, RH, sp_srs, SW_d, LW_d, R_net, melt ]
    surf_met_obs = [Ts_CS, Tair_CS, RH_CS, wind_CS, SWd_CS, LWd_CS, Rnet_CS, melt_CS]
    titles = ['$T_S$', '$T_{air}$', '\nRelative \nHumidity', '\nWind speed', '$SW_\downarrow$',  '$LW_\downarrow$', '$R_{net}$', 'melt']
    from itertools import chain
    for i in range(len(surf_met_mod)):
        slope, intercept, r2, p, sterr = scipy.stats.linregress(surf_met_obs[i], surf_met_mod[i])
        if p <= 0.01:
            ax[plot].text(0.75, 0.9, horizontalalignment='right', verticalalignment='top',
                          s='r$^{2}$ = %s' % np.round(r2, decimals=2), fontweight = 'bold', transform=ax[plot].transAxes, size=24,
                          color='dimgrey')
        else:
            ax[plot].text(0.75, 0.9, horizontalalignment='right', verticalalignment='top',
                          s='r$^{2}$ = %s' % np.round(r2, decimals=2), transform=ax[plot].transAxes, size=24,
                          color='dimgrey')
        ax[plot].scatter(surf_met_obs[i], surf_met_mod[i], color = '#f68080', s = 50)
        ax[plot].set_xlim(min(chain(surf_met_obs[i],surf_met_mod[i])), max(chain(surf_met_obs[i], surf_met_mod[i])))
        ax[plot].set_ylim(min(chain(surf_met_obs[i], surf_met_mod[i])), max(chain(surf_met_obs[i], surf_met_mod[i])))
        ax[plot].plot(ax[plot].get_xlim(), ax[plot].get_ylim(), ls="--", c = 'k', alpha = 0.8)
         #'r$^{2}$ = %s' % r2,
        #ax[plot].text(0.5, 1.1, s = titles[i], transform = ax[plot].transAxes)
        ax[plot].set_xlabel('Observed %s' % titles[i], size = 24, color = 'dimgrey', rotation = 0, labelpad = 10)
        ax[plot].set_ylabel('Modelled %s' % titles[i], size = 24, color = 'dimgrey', rotation =0, labelpad= 80)
        lab = ax[plot].text(0.1, 0.85, transform = ax[plot].transAxes, s=lab_dict[plot], fontsize=32, fontweight='bold', color='dimgrey')
        plot = plot +1
    plt.subplots_adjust(top = 0.98, hspace = 0.1, bottom = 0.05, wspace = 0.15, left = 0.2, right = 0.8)
    plt.setp(ax[2].get_xticklabels()[-1], visible=False)
    plt.setp(ax[5].get_xticklabels()[-2], visible=False)
    plt.setp(ax[6].get_xticklabels()[-2], visible=False)
    plt.setp(ax[1].get_xticklabels()[-3], visible=False)
    plt.setp(ax[2].get_xticklabels()[-3], visible=False)
    plt.setp(ax[2].get_yticklabels()[-1], visible=False)
    plt.setp(ax[5].get_yticklabels()[-2], visible=False)
    plt.setp(ax[6].get_yticklabels()[-2], visible=False)
    plt.setp(ax[1].get_yticklabels()[-3], visible=False)
    plt.setp(ax[2].get_yticklabels()[-3], visible=False)
    plt.savefig('/users/ellgil82/figures/Wintertime melt/Re-runs/correlations.png', transparent=True )
    plt.savefig('/users/ellgil82/figures/Wintertime melt/Re-runs/correlations.eps', transparent=True)
    plt.savefig('/users/ellgil82/figures/Wintertime melt/Re-runs/correlations.pdf', transparent=True)
    plt.show()

#correl_plot()
#bias_srs()

def calc_bias(res):
    Time_list, melt_CS, SH_CS, LH_CS, LWd_CS, SWd_CS, LWn_CS, SWn_CS, Time_CS, RH_CS, Ts_CS, Tair_CS, wind_CS, = load_AWS()
    sp_srs, T_surf, T_air, RH, Time_srs, percentiles_surf = load_surf(res)
    SH, LH, T_surf, Time_srs, melt, SW_n, LW_n, LW_d, SW_d, percentiles_SEB, melt_forced, E = load_SEB(res)
    R_net = SW_n + LW_n
    Rnet_CS = SWn_CS + LWn_CS
    total_SEB_obs_CS = SH_CS + LH_CS + SWn_CS + LWn_CS
    # Calculate bias of time series
    # Forecast error
    surf_obs = [RH_CS, wind_CS, Tair_CS, Ts_CS, melt_CS, SH_CS, LH_CS, LWn_CS, SWn_CS, SWd_CS, LWd_CS, total_SEB_obs_CS, melt_CS]
    surf_mod = [RH, sp_srs, T_air, T_surf, melt, SH, LH, LW_n, SW_n, SW_d, LW_d, E, melt_forced]
    mean_obs = []
    mean_mod = []
    bias = []
    errors = []
    for i in np.arange(len(surf_obs)):
        b = surf_mod[i] - surf_obs[i]
        errors.append(b)
        mean_obs.append(np.mean(surf_obs[i]))
        mean_mod.append(np.mean(surf_mod[i]))
        bias.append(mean_mod[i] - mean_obs[i])
    mean_error = np.mean(errors, axis=1)
    mean_sq_error = np.mean(np.power(errors, 2), axis=1)
    rmse = np.sqrt(mean_sq_error)
    idx = ['RH', 'wind', 'Tair', 'Ts', 'melt','SH', 'LH', 'LWn', 'SWn', 'SWd', 'LWd', 'total', 'melt forced']
    df = pd.DataFrame(index = idx)
    df['obs mean'] = pd.Series(mean_obs, index = idx)
    df['mod mean'] = pd.Series(mean_mod, index = idx)
    df['bias'] =pd.Series(bias, index=idx)
    df['rmse'] = pd.Series(rmse, index = idx)
    df['% RMSE'] = ( df['rmse']/df['obs mean'] ) * 100
    for i in range(len(surf_mod)):
        slope, intercept, r2, p, sterr = scipy.stats.linregress(surf_obs[i], surf_mod[i])
        print(idx[i])
        print('\nr2 = %s\n' % r2)
    print('RMSE/bias = \n\n\n')
    print(df)

#calc_bias('km1p5')

## ================================================ PLOTTING ======================================================== ##

## Set up plotting options
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Helvetica', 'Liberation sans', 'Tahoma', 'DejaVu Sans',
                               'Verdana']
#SH_srs, LH_srs,  T_surf, Time_srs, melt, SW_n, LW_n, LW_d, SW_d, percentiles_SEB = load_SEB('km1p5')
#sp_srs, T_surf, T_air, RH, Time_srs, percentiles_surf = load_surf('km1p5')

def SEB_plot():
    fig, ax = plt.subplots(2,2,sharex= True, sharey = True, figsize=(18, 12))
    ax = ax.flatten()
    Time_list, melt_CS, SH_CS, LH_CS,  LWd_CS, SWd_CS, LW_CS, SW_CS, Time_CS, RH_CS, Ts_CS, Tair_CS, wind_CS, = load_AWS()
    col_dict = {'km0p5': '#33a02c', 'km1p5': '#f68080', 'km4p0': '#1f78b4'}
    days = mdates.DayLocator(interval=1)
    dayfmt = mdates.DateFormatter('%d %b')
    for axs in ax:
        axs.spines['top'].set_visible(False)
        plt.setp(axs.spines.values(), linewidth=2, color='dimgrey', )
        axs.set_xlim(Time_list[1], Time_list[-1])
        axs.tick_params(axis='both', which='both', labelsize=24, tick1On=False, tick2On=False, labelcolor='dimgrey', pad=10)
        [l.set_visible(False) for (w, l) in enumerate(axs.yaxis.get_ticklabels()) if w % 2 != 0]
        [l.set_visible(False) for (w, l) in enumerate(axs.xaxis.get_ticklabels()) if w % 2 != 0]
        # Plot each res in turn
    res_list = ['km4p0', 'km1p5' ]
    for r in res_list:
        SH_srs, LH_srs, T_surf, Time_srs, melt, SW_n, LW_n, LW_d, SW_d, percentiles_SEB, melt_forced, E = load_SEB(r)
       # Shortwave
        ax2 = ax[0].twiny()
        obs = ax[0].plot(Time_list, SWd_CS, color='k', linewidth=2.5, label="Cabinet Inlet AWS", zorder = 3)
        ax[0].spines['right'].set_visible(False)
        ax2.plot(Time_srs, SW_d, linewidth=2.5, color=col_dict[r], label='*%(r)s UM output for Cabinet Inlet' % locals(), zorder = 4)
        ax2.axis('off')
        ax2.set_xlim(Time_srs[1], Time_srs[-1])
        ax2.tick_params(axis='both', tick1On = False, tick2On = False)
        ax2.set_ylim([-200,400])
        ax[0].set_ylim([-200, 400])
        ax[0].set_ylabel('Downwelling \nShortwave \nRadiation \n(W m$^{-2}$)', rotation=0, fontsize=24, color = 'dimgrey')
        ax[0].yaxis.set_label_coords(-0.4, 0.5)
        ax[0].tick_params(axis='both', which='both', labelsize=24, tick1On = False, tick2On = False)
        ax[0].text(x=Time_srs[25], y=330, fontweight='bold', s='a', fontsize=30, color='dimgrey', zorder=5)
        # Longwave
        ax2 = ax[1].twiny()
        obs = ax[1].plot(Time_list, LWd_CS, color='k', linewidth=2.5, label="Cabinet Inlet AWS", zorder =3)
        ax[1].spines['left'].set_visible(False)
        ax2.plot(Time_srs, LW_d, linewidth=2.5, color=col_dict[r], label='*%(r)s UM output for Cabinet Inlet' % locals(), zorder = 4)
        ax[1].set_ylabel('Downwelling \nLongwave \nRadiation \n(W m$^{-2}$)',  rotation=0, fontsize = 24,  color = 'dimgrey')
        ax[1].yaxis.set_label_coords(1.3, 0.5)
        ax[1].spines['right'].set_visible(False)
        ax2.set_xlim(Time_srs[1], Time_srs[-1])
        ax[1].tick_params(axis='x', tick1On=False, tick2On=False)
        ax2.axis('off')
        ax2.tick_params(axis='both', tick1On=False, tick2On=False)
        ax2.tick_params(axis='x', tick1On=False, tick2On=False)
        ax[1].set_ylim(-200, 400)
        ax2.set_ylim(-200, 400)
        ax[1].text(x=Time_srs[25], y=330, fontweight='bold', s='b', fontsize=30, color='dimgrey', zorder=5)
        # Sensible Heat
        ax2 = ax[2].twiny()
        ax[2].spines['right'].set_visible(False)
        obs = ax[2].plot(Time_list, SH_CS, color='k', linewidth=2.5, label="Cabinet Inlet AWS", zorder= 3)
        ax2.plot(Time_srs, SH_srs, linewidth=2.5, color=col_dict[r], label='*%(r)s UM output for Cabinet Inlet' % locals(), zorder = 4)
        ax[2].set_ylabel('Sensible \nHeat Flux \n(W m$^{-2}$)',  rotation=0, fontsize = 24,  color = 'dimgrey')
        ax[2].yaxis.set_label_coords(-0.4, 0.5)
        ax[2].tick_params(axis='both', which='both', labelsize=24, tick1On = False, tick2On = False )
        ax2.axis('off')
        ax2.set_xlim(Time_srs[1], Time_srs[-1])
        ax2.tick_params(axis='both', tick1On=False, tick2On=False)
        ax[2].set_ylim(-200, 400)
        ax2.set_ylim(-200, 400)
        ax[2].text(x=Time_srs[25], y=330, fontweight='bold', s='c', fontsize=30, color='dimgrey', zorder=5)
        # Latent Heat
        ax2 = ax[3].twiny()
        obs = ax[3].plot(Time_list, LH_CS, color='k', linewidth=2.5, label="Cabinet Inlet AWS", zorder = 3)
        ax2.plot(Time_srs, LH_srs, linewidth=2.5, color=col_dict[r], label='*%(r)s UM output for Cabinet Inlet' % locals(), zorder = 4)
        ax[3].set_ylabel('Latent \nHeat Flux \n(W m$^{-2}$)',  rotation=0, fontsize = 24,  color = 'dimgrey')
        ax[3].spines['left'].set_visible(False)
        ax[3].spines['right'].set_visible(False)
        ax[3].yaxis.set_label_coords(1.3, 0.5)
        ax[3].tick_params(axis='x', which='both', labelsize=24 )
        ax2.axis('off')
        ax2.set_xlim(Time_srs[1], Time_srs[-1])
        ax[3].tick_params(axis='both', tick1On=False, tick2On=False)
        ax2.tick_params(axis='both', tick1On=False, tick2On=False)
        ax[3].set_ylim(-200, 400)
        ax2.set_ylim(-200, 400)
        ax[3].text(x=Time_srs[25], y=330, fontweight='bold', s='d', fontsize=30, color='dimgrey', zorder=5)
        ax[3].fill_between(Time_srs, percentiles_SEB[2], percentiles_SEB[3], facecolor=col_dict[r], alpha=0.4, zorder=2)
        ax[2].fill_between(Time_srs, percentiles_SEB[0], percentiles_SEB[1], facecolor=col_dict[r], alpha=0.4, zorder=2)
        ax[1].fill_between(Time_srs, percentiles_SEB[8], percentiles_SEB[9], facecolor=col_dict[r], alpha=0.4, zorder=2)
        ax[0].fill_between(Time_srs, percentiles_SEB[4], percentiles_SEB[5], facecolor=col_dict[r], alpha=0.4, zorder=2)
    ax[2].xaxis.set_major_formatter(dayfmt)
    ax[3].xaxis.set_major_formatter(dayfmt)
    #Legend
    lns = [Line2D([0],[0], color='k', linewidth = 2.5)]
    labs = ['Observations from Cabinet Inlet']
    for r in res_list:
        lns.append(Line2D([0],[0], color=col_dict[r], linewidth = 2.5))
        labs.append(r[2]+'.'+r[4]+' km UM output for Cabinet Inlet')
    lgd = ax[1].legend(lns, labs, bbox_to_anchor=(0.55, 1.1), loc=2, fontsize=20)
    frame = lgd.get_frame()
    frame.set_facecolor('white')
    for ln in lgd.get_texts():
        plt.setp(ln, color='dimgrey')
    lgd.get_frame().set_linewidth(0.0)
    plt.subplots_adjust(wspace = 0.05, hspace = 0.05, top = 0.95, right = 0.8, left = 0.2, bottom = 0.08)
    plt.savefig('/users/ellgil82/figures/Wintertime melt/Re-runs/SEB_'+case+'_with_percentiles.png', transparent = True)
    plt.savefig('/users/ellgil82/figures/Wintertime melt/Re-runs/SEB_'+case+'_with_percentiles.eps', transparent = True)
    plt.savefig('/users/ellgil82/figures/Wintertime melt/Re-runs/SEB_' + case + '_with_percentiles.pdf',transparent=True)
    plt.show()

def surf_plot():
    fig, ax = plt.subplots(2,2,sharex= True, figsize=(22, 12))
    ax = ax.flatten()
    Time_list, melt_CS, SH_CS, LH_CS,  LWd_CS, SWd_CS,  LW_CS, SW_CS, Time_CS, RH_CS, Ts_CS, Tair_CS, wind_CS, = load_AWS()
    col_dict = {'km0p5': '#33a02c', 'km1p5': '#f68080', 'km4p0': '#1f78b4'}
    days = mdates.DayLocator(interval=1)
    dayfmt = mdates.DateFormatter('%d %b')
    for axs in ax:
        axs.spines['top'].set_visible(False)
        plt.setp(axs.spines.values(), linewidth=2, color='dimgrey', )
        axs.set_xlim(Time_list[1], Time_list[-1])
        axs.tick_params(axis='both', which='both', labelsize=24, tick1On=False, tick2On=False, labelcolor='dimgrey', pad=10)
        [l.set_visible(False) for (w, l) in enumerate(axs.yaxis.get_ticklabels()) if w % 2 != 0]
        [l.set_visible(False) for (w, l) in enumerate(axs.xaxis.get_ticklabels()) if w % 2 != 0]
    # Plot each res in turn
    res_list = [ 'km4p0', 'km1p5']
    for r in res_list:
        sp_srs, T_surf, T_air, RH, Time_srs, percentiles_surf = load_surf(r)
       # RH first
        obs = ax[0].plot(Time_list, RH_CS, color='k', linewidth=2.5, label="Cabinet Inlet AWS")
        ax2 = ax[0].twiny()
        ax[0].spines['right'].set_visible(False)
        RH[RH > 100] = 100
        ax2.plot(Time_srs, RH, linewidth=2.5, color=col_dict[r], label='*%(r)s UM output for Cabinet Inlet' % locals())
        ax2.axis('off')
        ax2.set_xlim(Time_srs[1], Time_srs[-1])
        ax[0].set_xlim(Time_srs[1], Time_srs[-1])
        ax2.tick_params(axis='both', tick1On = False, tick2On = False)
        ax2.set_ylim([0,100])
        ax2.yaxis.set_label_coords(-0.35, 0.5)
        ax[0].set_ylim([0, 100])
        ax[0].set_ylabel('Relative \nHumidity \n(%)', rotation=0, fontsize=24, color = 'dimgrey', labelpad = 80)
        ax[0].tick_params(axis='both', which='both', labelsize=24, tick1On = False, tick2On = False)
        ax[0].text(x=Time_srs[15], y= 90, fontweight='bold', s = 'a', fontsize = 30, color = 'dimgrey', zorder = 5)
        # Wind speed
        ax2 = ax[1].twiny()
        obs = ax[1].plot(Time_list, wind_CS, color='k', linewidth=2.5, label="Cabinet Inlet AWS")
        ax[1].spines['left'].set_visible(False)
        ax2.plot(Time_srs, sp_srs, linewidth=2.5, color=col_dict[r], label='*%(r)s UM output for Cabinet Inlet' % locals())
        ax[1].set_ylabel('Wind Speed \n(m s$^{-1}$)',  rotation=0, fontsize = 24,  color = 'dimgrey', labelpad = 80)
        ax[1].yaxis.set_label_coords(1.25, 0.5)
        ax[1].yaxis.tick_right()
        [l.set_visible(False) for (w, l) in enumerate(ax[1].yaxis.get_ticklabels()) if w % 2 != 0]
        ax2.set_xlim(Time_srs[1], Time_srs[-1])
        ax[1].set_xlim(Time_srs[1], Time_srs[-1])
        ax[1].tick_params(axis='both', tick1On=False, tick2On=False)
        ax[1].tick_params(axis='x', tick1On=False, tick2On=False)
        ax2.axis('off')
        ax2.tick_params(axis='both', tick1On=False, tick2On=False)
        ax2.tick_params(axis='x', tick1On=False, tick2On=False)
        ax[1].set_ylim(0, 40)
        ax2.set_ylim(0, 40)
        ax[1].text(x=Time_srs[15], y=35, fontweight='bold', s='b', fontsize=30, color='dimgrey', zorder=5)
        #Temperature
        ax2 = ax[2].twiny()
        ax[2].spines['right'].set_visible(False)
        obs = ax[2].plot(Time_list, Tair_CS, color='k', linewidth=2.5, label="Cabinet Inlet AWS")
        ax2.plot(Time_srs, T_air,linewidth=2.5, color=col_dict[r], label='*%(r)s UM output for Cabinet Inlet' % locals())
        ax[2].set_ylabel('Near-surface \nAir Temperature \n($^{\circ}$C)',  rotation=0, fontsize = 24,  color = 'dimgrey', labelpad = 80)
        ax[2].tick_params(axis='both', which='both', labelsize=24, tick1On = False, tick2On = False )
        ax2.axis('off')
        ax[2].yaxis.set_label_coords(-0.1, 0.5)
        ax2.set_xlim(Time_srs[1], Time_srs[-1])
        ax[2].set_xlim(Time_srs[1], Time_srs[-1])
        ax2.tick_params(axis='both', tick1On=False, tick2On=False)
        ax[2].set_ylim(-30, 15)
        ax2.set_ylim(-30, 15)
        ax[2].text(x=Time_srs[15], y=10, fontweight='bold', s='c', fontsize=30, color='dimgrey', zorder=5)
        plt.setp(ax[2].get_yticklabels()[-1], visible=False)
        plt.setp(ax[2].get_yticklabels()[-3], visible=False)
        #Temperature
        ax2 = ax[3].twiny()
        obs = ax[3].plot(Time_list, Ts_CS, color='k', linewidth=2.5, label="Cabinet Inlet AWS")
        ax2.plot(Time_srs, T_surf,linewidth=2.5, color=col_dict[r], label='*%(r)s UM output for Cabinet Inlet' % locals())
        ax[3].set_ylabel('Surface \nTemperature\n($^{\circ}$C)',  rotation=0, fontsize = 24,  color = 'dimgrey',labelpad = 80)
        ax[3].spines['left'].set_visible(False)
        ax[3].yaxis.set_label_coords(1.25, 0.5)
        ax[3].yaxis.tick_right()
        ax[3].tick_params(axis='both', which='both', labelsize=24 )
        ax2.axis('off')
        ax2.set_xlim(Time_srs[1], Time_srs[-1])
        ax[3].set_xlim(Time_srs[1], Time_srs[-1])
        ax[3].tick_params(axis='both', tick1On=False, tick2On=False)
        ax2.tick_params(axis='both', tick1On=False, tick2On=False)
        [l.set_visible(False) for (w, l) in enumerate(ax[3].yaxis.get_ticklabels()) if w % 2 != 0]
        ax[3].set_ylim(-30, 15)
        ax2.set_ylim(-30, 15)
        ax[3].text(x=Time_srs[15], y=10, fontweight='bold', s='d', fontsize=30, color='dimgrey', zorder=5)
        plt.setp(ax[3].get_yticklabels()[-3], visible=False)
        plt.setp(ax[3].get_yticklabels()[-1], visible=False)
        plt.setp(ax[3].get_yticklabels()[-5], visible=False)
        ax[3].fill_between(Time_srs, percentiles_surf[0], percentiles_surf[1], facecolor=col_dict[r], alpha=0.4)
        ax[2].fill_between(Time_srs, percentiles_surf[2], percentiles_surf[3], facecolor=col_dict[r], alpha=0.4)
        ax[0].fill_between(Time_srs, percentiles_surf[4], percentiles_surf[5], facecolor=col_dict[r], alpha=0.4)
        ax[1].fill_between(Time_srs, percentiles_surf[6], percentiles_surf[7], facecolor=col_dict[r], alpha=0.4)
    ax[2].yaxis.set_label_coords(-0.3, 0.5)
    ax[2].xaxis.set_major_formatter(dayfmt)
    ax[3].xaxis.set_major_formatter(dayfmt)
    # Legend
    lns = [Line2D([0],[0], color='k', linewidth = 2.5)]
    labs = ['Observations from Cabinet Inlet']
    for r in res_list:
        lns.append(Line2D([0],[0], color=col_dict[r], linewidth = 2.5))
        labs.append(r[2]+'.'+r[4]+' km UM output for Cabinet Inlet')
    lgd = ax[1].legend(lns, labs, bbox_to_anchor=(0.55, 1.1), loc=2, fontsize=20)
    frame = lgd.get_frame()
    frame.set_facecolor('white')
    for ln in lgd.get_texts():
        plt.setp(ln, color='dimgrey')
    lgd.get_frame().set_linewidth(0.0)
    plt.subplots_adjust(wspace = 0.05, hspace = 0.05, top = 0.95, right = 0.85, left = 0.16, bottom = 0.08)
    plt.savefig('/users/ellgil82/figures/Wintertime melt/Re-runs/surface_met_'+case+'_with_range.png', transparent = True)
    plt.savefig('/users/ellgil82/figures/Wintertime melt/Re-runs/surface_met_'+case+'_with_range.eps', transparent = True)
    plt.savefig('/users/ellgil82/figures/Wintertime melt/Re-runs/surface_met_'+case+'_with_range.pdf', transparent = True)
    plt.show()

print('\nplotting surface vars....')
#surf_plot()

print('\nplotting SEBs....')
#SEB_plot()

def total_SEB_model():
    SH_srs, LH_srs, T_surf, Time_srs, melt, SW_n, LW_n, LW_d, SW_d = load_SEB('km1p5')
    fig, ax = plt.subplots(1,1, figsize = (18,8), frameon= False)
    days = mdates.DayLocator(interval=1)
    dayfmt = mdates.DateFormatter('%b %d')
    plt.setp(ax.spines.values(), linewidth=2, color = 'dimgrey')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.plot(Time_srs, SW_n, color = '#6fb0d2', lw = 2.5, label = 'Net shortwave flux')
    ax.plot(Time_srs, LW_n, color = '#86ad63', lw = 2.5, label = 'Net longwave flux')
    ax.plot(Time_srs, SH_srs, color = '#1f78b4', lw = 2.5, label = 'Sensible heat flux')
    ax.plot(Time_srs, LH_srs, color = '#33a02c', lw = 2.5, label = 'Latent heat flux')
    ax.plot(Time_srs, melt, color = '#f68080', lw = 2.5, label = 'Melt flux')
    lgd = plt.legend(fontsize = 18, frameon = False)
    for ln in lgd.get_texts():
        plt.setp(ln, color = 'dimgrey')
    ax.xaxis.set_major_formatter(dayfmt)
    ax.text(x=Time_srs[6], y=350, s='b', fontsize=32, fontweight='bold', color='dimgrey')
    ax.tick_params(axis='both', which='both', labelsize=24, tick1On = False, tick2On=False, labelcolor = 'dimgrey', pad = 10)
    ax.axhline(y=0, xmin = 0, xmax = 1, linestyle = '--', linewidth = 1)
    ax.set_ylabel('Energy flux \n (W m$^{-2}$)',  rotation = 0, fontsize = 28, labelpad = 70, color='dimgrey')
    plt.subplots_adjust(left = 0.2, bottom = 0.1, right = 0.95)
    [l.set_visible(False) for (w,l) in enumerate(ax.yaxis.get_ticklabels()) if w % 2 != 0]
    [l.set_visible(False) for (w,l) in enumerate(ax.xaxis.get_ticklabels()) if w % 2 != 0]
    plt.savefig('/users/ellgil82/figures/Wintertime melt/Re-runs/CS1_model_SEB_km1p5.eps', transparent=True )
    plt.savefig('/users/ellgil82/figures/Wintertime melt/Re-runs/CS1_model_SEB_km1p5.png', transparent=True )
    plt.show()

def total_SEB_obs():
    Time_list, melt_CS, SH_CS, LH_CS, LWd_CS, SWd_CS, LWn_CS, SWn_CS, Time_CS, RH_CS, Ts_CS, Tair_CS, wind_CS = load_AWS()
    fig, ax = plt.subplots(1,1, figsize = (18,8), frameon= False)
    days = mdates.DayLocator(interval=1)
    dayfmt = mdates.DateFormatter('%b %d')
    Time_list[0] = Time_list[0].replace(second=1)
    plt.setp(ax.spines.values(), linewidth=2, color = 'dimgrey')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.text(x=Time_srs[6], y=350, s='a', fontsize=32, fontweight='bold', color='dimgrey')
    ax.set_ylim(-200,400)
    ax.plot(Time_list[:360], SWn_CS[:360], color = '#6fb0d2', lw = 2.5, label = 'Net shortwave \nradiation flux')
    ax.plot(Time_list[:360], LWn_CS[:360], color = '#86ad63', lw = 2.5, label = 'Net longwave \nradiation flux')
    ax.plot(Time_list[:360], SH_CS[:360], color = '#1f78b4', lw = 2.5, label = 'Sensible heat flux')
    ax.plot(Time_list[:360], LH_CS[:360], color = '#33a02c', lw = 2.5, label = 'Latent heat flux')
    ax.plot(Time_list[:360], melt_CS[:360], color = '#f68080', lw = 2.5, label = 'Melt flux')
    lgd = plt.legend(fontsize = 18, frameon = False)
    for ln in lgd.get_texts():
        plt.setp(ln, color = 'dimgrey')
    ax.xaxis.set_major_formatter(dayfmt)
    ax.tick_params(axis='both', which='both', labelsize=24, tick1On = False, tick2On=False, labelcolor = 'dimgrey', pad = 10)
    ax.axhline(y=0, xmin = 0, xmax = 1, linestyle = '--', linewidth = 1)
    ax.set_ylabel('Energy flux \n (W m$^{-2}$)',  rotation = 0, fontsize = 28, labelpad = 70, color='dimgrey')
    plt.subplots_adjust(left = 0.2, bottom = 0.1, right = 0.95)
    [l.set_visible(False) for (w,l) in enumerate(ax.yaxis.get_ticklabels()) if w % 2 != 0]
    [l.set_visible(False) for (w,l) in enumerate(ax.xaxis.get_ticklabels()) if w % 2 != 0]
    plt.savefig('/users/ellgil82/figures/Wintertime melt/Re-runs/CS1_obs_SEB.eps', transparent=True )
    plt.savefig('/users/ellgil82/figures/Wintertime melt/Re-runs/CS1_obs_SEB.png', transparent=True )
    plt.show()

def total_SEB():
    fig, axs = plt.subplots(2, 1, figsize=(18, 18), sharex = True, frameon=False)
    axs = axs.flatten()
    days = mdates.DayLocator(interval=1)
    dayfmt = mdates.DateFormatter('%b %d')
    Time_list, melt_CS, SH_CS, LH_CS, LWd_CS, SWd_CS, LWn_CS, SWn_CS, Time_CS, RH_CS, Ts_CS, Tair_CS, wind_CS = load_AWS()
    SH, LH, T_surf, Time_srs, melt, SW_n, LW_n, LW_d, SW_d, percentiles_SEB, melt_forced, E = load_SEB('km1p5')
    Time_list[0] = Time_list[0].replace(second=1)
    for ax in axs:
        plt.setp(ax.spines.values(), linewidth=2, color='dimgrey')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_xlim(Time_list[0], Time_list[-1])
        ax.xaxis.set_major_formatter(dayfmt)
        ax.tick_params(axis='both', which='both', labelsize=24, tick1On=False, tick2On=False, labelcolor='dimgrey', pad=10)
        ax.axhline(y=0, xmin=0, xmax=1, linestyle='--', linewidth=1)
        ax.set_ylabel('Energy flux \n (W m$^{-2}$)', rotation=0, fontsize=28, labelpad=70, color='dimgrey')
        [l.set_visible(False) for (w, l) in enumerate(ax.yaxis.get_ticklabels()) if w % 2 != 0]
        [l.set_visible(False) for (w, l) in enumerate(ax.xaxis.get_ticklabels()) if w % 2 != 0]
        ax.set_ylim(-200, 400)
    axs[0].plot(Time_list, SWn_CS, color='#6fb0d2', lw=2.5, label='Net shortwave flux')
    axs[0].plot(Time_list, LWn_CS, color='#86ad63', lw=2.5, label='Net longwave flux')
    axs[0].plot(Time_list, SH_CS, color='#1f78b4', lw=2.5, label='Sensible heat flux')
    axs[0].plot(Time_list, LH_CS, color='#33a02c', lw=2.5, label='Latent heat flux')
    axs[0].plot(Time_list, melt_CS, color='#f68080', lw=2.5, label='Melt flux')
    axs[1].plot(Time_srs, SW_n, color='#6fb0d2', lw=2.5, label='Net shortwave flux')
    axs[1].plot(Time_srs, LW_n, color='#86ad63', lw=2.5, label='Net longwave flux')
    axs[1].plot(Time_srs, SH, color='#1f78b4', lw=2.5, label='Sensible heat flux')
    axs[1].plot(Time_srs, LH, color='#33a02c', lw=2.5, label='Latent heat flux')
    axs[1].plot(Time_srs, melt, color='#f68080', lw=2.5, label='Melt flux')
    lgd = plt.legend(fontsize=24, bbox_to_anchor = (0.9, 1.2))
    for ln in lgd.get_texts():
        plt.setp(ln, color='dimgrey')
    lgd.get_frame().set_linewidth(0.0)
    axs[0].text(x=Time_list[10], y=320, s='a', fontsize=32, fontweight='bold', color='dimgrey')
    axs[1].text(x=Time_srs[10], y=320, s='b', fontsize=32, fontweight='bold', color='dimgrey')
    plt.subplots_adjust(left=0.2, bottom=0.1, right=0.95, hspace = 0.1, wspace = 0.1)
    plt.savefig('/users/ellgil82/figures/Wintertime melt/Re-runs/'+case+'total_SEB.eps', transparent = True)
    plt.savefig('/users/ellgil82/figures/Wintertime melt/Re-runs/'+case+'total_SEB.png', transparent=True)
    plt.savefig('/users/ellgil82/figures/Wintertime melt/Re-runs/' + case + 'total_SEB.pdf', transparent=True)
    plt.show()

total_SEB()


def bias_colour():
    fig, ax = plt.subplots(2, 2, figsize=(14, 12))
    ax = ax.flatten()
    Time_list, melt_CS, SH_CS, LH_CS, LW_CS, SW_CS, Time_CS, RH_CS, Ts_CS, Tair_CS, wind_CS, = load_AWS()
    col_dict = {'km0p5': '#33a02c', 'km1p5': '#f68080', 'km4p0': '#1f78b4'}
    CbAx = fig.add_axes([0.9, 0.3, 0.03, 0.4])  # x0, y0, width, height
    for axs in ax:
        axs.spines['top'].set_visible(False)
        plt.setp(axs.spines.values(), linewidth=2, color='dimgrey', )
        axs.tick_params(axis='both', which='both', labelsize=24, tick1On=False, tick2On=False, labelcolor='dimgrey', pad=10)
        [l.set_visible(False) for (w, l) in enumerate(axs.yaxis.get_ticklabels()) if w % 2 != 0]
        [l.set_visible(False) for (w, l) in enumerate(axs.xaxis.get_ticklabels()) if w % 2 != 0]
    # Plot each res in turn
    sp_srs, T_surf, T_air, RH_srs, Time_srs = load_surf('km1p5')
    # RH
    fit = np.polyfit(RH_srs, RH_CS, deg=1)
    ax[0].spines['right'].set_visible(False)
    ax[0].scatter(RH_srs, RH_CS, c = mdates.date2num(Time_list), cmap = 'GnBu', vmin = mdates.date2num(Time_list[0])-2, vmax = mdates.date2num(Time_list[-1])+2, marker = '+', lw = 2, s = 200)
    ax[0].plot(RH_srs, fit[0] * RH_srs + fit[1], color = 'k', lw = 1)
    ax[0].set_xlim(50, 110)
    ax[0].set_ylim(30, 105)
    ax[0].set_ylabel('Relative \nHumidity \n(%)', rotation=0, fontsize=24, color='dimgrey')
    ax[0].tick_params(axis='both', which='both', labelsize=24, tick1On=False, tick2On=False)
    # Wind speed
    fit = np.polyfit(sp_srs, wind_CS, deg=1)
    ax[1].scatter(sp_srs, wind_CS, c = mdates.date2num(Time_list), cmap = 'GnBu', vmin = mdates.date2num(Time_list[0])-2, vmax = mdates.date2num(Time_list[-1])+2, marker = '+', lw = 2, s = 200)
    ax[1].plot(sp_srs, fit[0] * sp_srs + fit[1], color = 'k', lw = 1)
    ax[1].set_xlim(0, 30)
    ax[1].set_ylim(0, 22)
    ax[1].set_ylabel('Wind \nSpeed \n(m s${-1}$)', rotation=0, fontsize=24, color='dimgrey')
    ax[1].tick_params(axis='both', which='both', labelsize=24, tick1On=False, tick2On=False)
    # Air temp
    fit = np.polyfit(T_air, Tair_CS, deg=1)
    ax[2].scatter(T_air, Tair_CS, c = mdates.date2num(Time_list), cmap = 'GnBu', vmin = mdates.date2num(Time_list[0])-2, vmax = mdates.date2num(Time_list[-1])+2, marker = '+', lw = 2, s = 200)
    ax[2].plot(T_air, fit[0] * T_air + fit[1], color = 'k', lw = 1)
    ax[2].spines['right'].set_visible(False)
    ax[2].set_xlim(-30, 8)
    ax[2].set_ylim(-30, 8)
    ax[2].set_ylabel('Near-surface \nAir Temperature \n($^{\circ}$C)', rotation=0, fontsize=24, color='dimgrey')
    ax[2].tick_params(axis='both', which='both', labelsize=24, tick1On=False, tick2On=False)
    # Surf temp
    fit = np.polyfit(T_surf, Ts_CS, deg=1)
    scatter = ax[3].scatter(T_surf, Ts_CS, c = mdates.date2num(Time_list), cmap = 'GnBu',vmin = mdates.date2num(Time_list[0])-2, vmax = mdates.date2num(Time_list[-1])+2, marker = '+', lw = 2, s = 200)
    ax[3].plot(T_surf, fit[0] * T_surf + fit[1], color = 'k', lw = 1)
    ax[3].set_xlim(-30, 1)
    ax[3].set_ylim(-30, 1)
    ax[3].set_ylabel('Surface \nTemperature \n($^{\circ}$C)', rotation=0, fontsize=24, color='dimgrey')
    ax[3].tick_params(axis='both', which='both', labelsize=24, tick1On=False, tick2On=False)
    cb = plt.colorbar(scatter, format = mdates.DateFormatter('%d %b'), cax = CbAx)
    plt.subplots_adjust(left = 0.1, bottom = 0.08, top = 0.92, right = 0.88)
    plt.show()

def bias_srs():
    fig, ax = plt.subplots(2, 2, figsize=(18, 12))
    ax = ax.flatten()
    Time_list, melt_CS, SH_CS, LH_CS, LW_CS, SW_CS, Time_CS, RH_CS, Ts_CS, Tair_CS, wind_CS, = load_AWS()
    col_dict = {'km0p5': '#33a02c', 'km1p5': '#f68080', 'km4p0': '#1f78b4'}
    for axs in ax:
        axs.spines['top'].set_visible(False)
        plt.setp(axs.spines.values(), linewidth=2, color='dimgrey', )
        axs.tick_params(axis='both', which='both', labelsize=24, tick1On=False, tick2On=False, labelcolor='dimgrey',pad = 10)
        [l.set_visible(False) for (w, l) in enumerate(axs.yaxis.get_ticklabels()) if w % 2 != 0]
        [l.set_visible(False) for (w, l) in enumerate(axs.xaxis.get_ticklabels()) if w % 2 != 0]
    # Plot each res in turn
    res_list = ['km1p5']
    for r in res_list:
        sp_srs, T_surf, T_air, RH_srs, Time_srs = load_surf(r)
        # RH
        one_to_one = np.polyfit([0,1],[0,1], deg=1)
        fit = np.polyfit(RH_srs, RH_CS, deg=1)
        ax[0].spines['right'].set_visible(False)
        ax[0].scatter(RH_srs, RH_CS, c = col_dict[r], marker = '+', lw = 2, s = 200)
        ax[0].plot(RH_srs, fit[0] * RH_srs + fit[1], color = 'k', lw = 1)
        ax[0].plot(RH_srs, one_to_one[0] * RH_srs + one_to_one[1], color='dimgrey', lw=1, linestyle='--')
        ax[0].set_xlim(40, 100)
        ax[0].set_ylim(40, 100)
        ax[0].set_ylabel('Relative \nHumidity \n(%)', rotation=0, fontsize=24, color='dimgrey', labelpad = 80)
        ax[0].tick_params(axis='both', which='both', labelsize=24, tick1On=False, tick2On=False)
        # Wind speed
        fit = np.polyfit(sp_srs, wind_CS, deg=1)
        ax[1].spines['left'].set_visible(False)
        ax[1].yaxis.set_label_coords(1.25, 0.5)
        ax[1].scatter(sp_srs, wind_CS, c = col_dict[r], marker = '+', lw = 2, s = 200)
        ax[1].plot(sp_srs, fit[0] * sp_srs + fit[1], color = 'k', lw = 1)
        ax[1].plot(sp_srs, one_to_one[0] * sp_srs + one_to_one[1], color='dimgrey', lw=1, linestyle='--')
        ax[1].set_xlim(0, 30)
        ax[1].set_ylim(0, 30)
        ax[1].set_ylabel('Wind \nSpeed \n(m s${-1}$)', rotation=0, fontsize=24, color='dimgrey', labelpad = 80)
        ax[1].tick_params(axis='both', which='both', labelsize=24, tick1On=False, tick2On=False)
        # Air temp
        fit = np.polyfit(T_air, Tair_CS, deg=1)
        ax[2].scatter(T_air, Tair_CS, c = col_dict[r], marker = '+', lw = 2, s = 200)
        ax[2].plot(T_air, fit[0] * T_air + fit[1], color = 'k', lw = 1)
        ax[2].plot(T_air, one_to_one[0] * T_air + one_to_one[1], color='dimgrey', lw=1, linestyle='--')
        ax[2].spines['right'].set_visible(False)
        ax[1].yaxis.tick_right()
        ax[2].set_xlim(-30, 10)
        ax[2].set_ylim(-30, 10)
        ax[2].set_ylabel('Near-surface \nAir Temperature \n($^{\circ}$C)', rotation=0, fontsize=24, color='dimgrey', labelpad = 80)
        ax[2].tick_params(axis='both', which='both', labelsize=24, tick1On=False, tick2On=False)
        # Surf temp
        fit = np.polyfit(T_surf, Ts_CS, deg=1)
        ax[3].spines['left'].set_visible(False)
        ax[3].yaxis.set_label_coords(1.25, 0.5)
        ax[3].scatter(T_surf, Ts_CS, c = col_dict[r], marker = '+', lw = 2, s = 200)
        ax[3].plot(T_surf, fit[0] * T_surf + fit[1], color = 'k', lw = 1)
        ax[3].plot(T_surf, one_to_one[0] * T_surf + one_to_one[1], color='dimgrey', lw=1, linestyle='--')
        ax[3].set_xlim(-30, 1)
        ax[3].set_ylim(-30, 1)
        ax[3].set_ylabel('Surface \nTemperature \n($^{\circ}$C)', rotation=0, fontsize=24, color='dimgrey', labelpad = 120)
        ax[3].tick_params(axis='both', which='both', labelsize=24, tick1On=False, tick2On=False)
        ax[3].yaxis.tick_right()
    plt.subplots_adjust(left = 0.2, bottom = 0.08, top = 0.92, right = 0.8, wspace = 0.08,)
    plt.show()


'''
Time_list, melt_CS, SH_CS, LH_CS, LWd_CS, SWd_CS, LWn_CS, SWn_CS, Time_CS, RH_CS, Ts_CS, Tair_CS, wind_CS,  = load_AWS()
sp_srs, T_surf, T_air, RH, Time_srs, = load_surf('km4p0')
SH_srs, LH_srs, T_surf, Time_srs, melt, SW_n, LW_n, LW_d, SW_d = load_SEB('km4p0')

#SH_srs, LH_srs, LW_d, SW_d,  T_surf, Time_srs, melt = load_SEB(r)

#bias_srs()

total_SEB_obs_CS = SH_CS + LH_CS + SWn_CS + LWn_CS
total_SEB_model_CS = SH_srs + LH_srs + LW_n + SW_n

# Calculate bias of time series
# Forecast error
surf_obs = [RH_CS, wind_CS, Tair_CS, Ts_CS, melt_CS, SH_CS, LH_CS, LWn_CS, SWn_CS, SWd_CS, LWd_CS, total_SEB_obs_CS]
surf_mod = [RH, sp_srs, T_air, T_surf, melt, SH_srs, LH_srs, LW_n, SW_n, SW_d, LW_d, total_SEB_model_CS]
mean_obs = []
mean_mod = []
bias = []
errors = []
for i in np.arange(len(surf_obs)):
    b = surf_mod[i] - surf_obs[i]
    errors.append(b)
    mean_obs.append(np.mean(surf_obs[i]))
    mean_mod.append(np.mean(surf_mod[i]))
    bias.append(mean_mod[i] - mean_obs[i])

mean_error = np.mean(errors, axis=1)
mean_sq_error = np.mean(np.power(errors, 2), axis=1)
rmse = np.sqrt(mean_sq_error)

idx = ['RH', 'wind', 'Tair', 'Ts', 'melt','SH', 'LH', 'LWn', 'SWn', 'SWd', 'LWd', 'total']
import pandas as pd
df = pd.DataFrame(index = idx)
df['obs mean'] = pd.Series(mean_obs, index = idx)
df['mod mean'] = pd.Series(mean_mod, index = idx)
df['bias'] =pd.Series(bias, index=idx)
df['rmse'] = pd.Series(rmse, index = idx)
df['% RMSE'] = ( df['rmse']/df['obs mean'] ) * 100

print df '''