""" Script for calculating differences between default MetUM ancillary files for the land-sea mask
and orography, and for plotting snapshot maps of near-surface temperature and wind vectors to visualise the effect of this.

N.B. Default UM orography and coastlines are based on the GLOBE dataset (NGDC, 1999) at 1 km resolution and are based on
data collected in 1993, so predates the collapse of the Larsen A and B ice shelves.

The updated land-sea mask is based on the SCAR Antarctic Digital Database coastline, version 7.0 (released January 2016
and available at https://www.add.scar.org/). The orography ancillary file is based on the Ohio State University RAMP
200 m resolution Antarctic Digital Elevation Model (DEM) (Hongxing, 1999).

References:
Hongxing, L. (1999). Development of an Antarctic Digital Elevation model. Technical report. Byrd Polar Research Center,
 	Ohio State University, page 157.
NGDC (National Geophysical Data Center) (1999). Global Land One-kilometer Base Elevation (GLOBE) v.1.

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
sys.path.append('/users/ellgil82/scripts/Tools/')
from rotate_data import rotate_data
from divg_temp_colourmap import shiftedColorMap

def load_vars(file):
	T_air = iris.load_cube(file, 'air_temperature')
	T_surf = iris.load_cube(file, 'surface_temperature')
	T_air.convert_units('celsius')
	T_surf.convert_units('celsius')
	u = iris.load_cube(file, 'x_wind')
	v = iris.load_cube(file, 'y_wind')
	v = v[:,:400]
	old_lsm = iris.load_cube('/data/clivarm/wip/ellgil82/May_2016/Compare/CS1/km1p5/20160525T1200Z_Peninsula_km1p5_ctrl_pa000.pp', 'land_binary_mask')
	old_orog = iris.load_cube('/data/clivarm/wip/ellgil82/May_2016/Compare/CS1/km1p5/20160525T1200Z_Peninsula_km1p5_ctrl_pa000.pp', 'surface_altitude')
	new_lsm = iris.load_cube('/data/clivarm/wip/ellgil82/May_2016/Re-runs/CS2/20160522T1200Z_Peninsula_km1p5_Smith_tnuc_pa000.pp', 'land_binary_mask')
	new_orog = iris.load_cube('/data/clivarm/wip/ellgil82/May_2016/Re-runs/CS2/20160522T1200Z_Peninsula_km1p5_Smith_tnuc_pa000.pp', 'surface_altitude')
	rotate_me = [T_air, T_surf, u, v]
	for i in rotate_me:
		real_lon, real_lat = rotate_data(i, 1, 2)
	me_too = [new_lsm, new_orog, old_lsm, old_orog]
	for j in me_too:
		real_lon, real_lat = rotate_data(j, 0, 1)
	return T_air[0,:,:], T_surf[0,:,:], u[0,:,:], v[0,:,:], new_orog, new_lsm, old_orog, old_lsm,  real_lon, real_lat

## Set up plotting options
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Segoe UI', 'Helvetica', 'Liberation sans', 'Tahoma', 'DejaVu Sans','Verdana']

## Caption: Near-surface (1.5 m) air temperatures and 10 m wind speeds during non-foehn (left column; panels a and c)
## and foehn conditions (right column; panels b and d) during CS2 using the default MetUM orography and land-sea mask
## (top row; panels a and b) and updated orography and land-sea mask (bottom row; panels c and d). Non-foehn panels show
## hourly mean conditions for the period 22 May 23:00 until 23 May 00:00 UTC and foehn panels show hourly mean
## conditions for the period 26 May 11:00 until 12:00 UTC. Filled contours show 1.5 m air temperatures, while vectors
## show 10 m wind speeds. A scale vector of 10 m s$^{-1}$ is given at top centre for reference.

def plot_difs():
	''' Function to plot surface temperature and wind vectors during foehn/non-foehn conditions using both default and
	updated ancillary files'''
	# Define files from which to import data
	new_non_foehn = '/data/clivarm/wip/ellgil82/May_2016/Re-runs/CS2/20160522T1200Z_Peninsula_km1p5_ctrl_pa000.pp'
	new_foehn = '/data/clivarm/wip/ellgil82/May_2016/Re-runs/CS2/20160525T1200Z_Peninsula_km1p5_ctrl_pa000.pp'
	old_non_foehn = '/data/clivarm/wip/ellgil82/May_2016/Compare/CS2/km1p5/20160522T1200Z_Peninsula_km1p5_ctrl_pa012.pp'
	old_foehn = '/data/clivarm/wip/ellgil82/May_2016/Compare/CS2/km1p5/20160525T1200Z_Peninsula_km1p5_ctrl_pa012.pp'
	files = [old_non_foehn, old_foehn, new_non_foehn, new_foehn]
	# Set up figure
	fig, ax = plt.subplots(2,2, figsize = (20,20))
	ax = ax.flatten()
	cbaxes = fig.add_axes([0.35, 0.15, 0.5, 0.02])
	plot = 0
	lab_dict = {0: 'a', 1: 'b', 2: 'c', 3:'d'}
	for axs in ax:
		axs.spines['right'].set_visible(False)
		axs.spines['left'].set_visible(False)
		axs.spines['top'].set_visible(False)
		axs.spines['bottom'].set_visible(False)
		axs.tick_params(axis='both', which='both', length = 0, labelbottom='off', labelleft='off')
	for file in files:
		# Load_everything
		T_air, T_surf, u, v, new_orog, new_lsm, old_orog, old_lsm, real_lon, real_lat = load_vars(file)
		# Re-centre colour map around zero
		bwr_zero = shiftedColorMap(cmap = matplotlib.cm.bwr, min_val = -35, max_val =15, name = 'bwr_zero', var = T_air.data,start = 0, stop = 1)
		# Plot air temperatures
		x, y = np.meshgrid(T_air.coord('longitude').points,T_air.coord('latitude').points)
		c = ax[plot].pcolormesh(T_air.coord('longitude').points, T_air.coord('latitude').points, T_air.data, cmap = bwr_zero, vmin = -35, vmax = 15, zorder = 1)
		# Plot AWS location
		ax[plot].plot(-63.37105, -66.48272, markersize=20, marker='o', color='lime', zorder = 2)
		# Plot wind vectors
		q1 = ax[plot].quiver(x[::60, ::60], y[::60, ::60], u.data[::60, ::60], v.data[::60, ::60], pivot='middle', scale=100, zorder = 4)
		# Adjust axes limits
		ax[plot].set_xlim(min(T_air.coord('longitude').points),max(T_air.coord('longitude').points))
		ax[plot].set_ylim(min(T_air.coord('latitude').points),max(T_air.coord('latitude').points))
		ax[plot].tick_params(axis='both', which='both', length = 0, labelbottom='off', labelleft='off')
		# Add sub-plot label
		lab = ax[plot].text(0.1, 0.85, transform=ax[plot].transAxes, s=lab_dict[plot], fontsize=32, fontweight='bold', color='dimgrey', zorder = 5)
		plot = plot+1
	# Set up colourbar
	cb = plt.colorbar(c, cax = cbaxes, orientation = 'horizontal')
	[l.set_visible(False) for (w, l) in enumerate(cb.ax.xaxis.get_ticklabels()) if w % 2 != 0]
	cb.ax.set_xlabel('1.5 m air temperature ($^{\circ}$C)', fontsize=30, labelpad = 20, color = 'dimgrey')
	cb.outline.set_edgecolor('dimgrey')
	cb.outline.set_linewidth(2)
	cb.ax.tick_params(labelsize= 30, labelcolor='dimgrey', pad = 10)
	# Add key for vectors
	plt.quiverkey(q1,0.56, 0.95, 10, r'$10$ $m$ $s^{-1}$', labelcolor = 'dimgrey', color = 'dimgrey',  labelpos='N', fontproperties={'size':'30', 'weight':'bold'},coordinates='figure')
	# Add LSM and orography contours into correct axes
	for axs in [ax[0], ax[1]]:
		axs.contour(T_air.coord('longitude').points, T_air.coord('latitude').points, old_lsm.data, colors='dimgrey', linewidth=3, zorder = 3)
		axs.contour(T_air.coord('longitude').points, T_air.coord('latitude').points, old_lsm.data, levels = [100], colors='dimgrey', linewidth=2, zorder = 3)
	for axs in [ax[2], ax[3]]:
		axs.contour(T_air.coord('longitude').points, T_air.coord('latitude').points, new_lsm.data, colors='dimgrey', linewidth=3, zorder = 3)
		axs.contour(T_air.coord('longitude').points, T_air.coord('latitude').points, new_lsm.data, levels=[100], colors='dimgrey', linewidth=2, zorder=3)
	# Add titles
	ax[0].set_title('Non-foehn conditions', fontsize = 30, y = 1.05, color = 'dimgrey')
	ax[1].set_title('Foehn conditions', fontsize = 30,  y = 1.05, rotation = 0, color = 'dimgrey')
	ax[0].set_ylabel('Default set-up', fontsize = 30, x = -0.5, y = 0.5, rotation = 0, color = 'dimgrey', labelpad = 120)
	ax[2].set_ylabel('Updated set-up', fontsize = 30, x = -0.4, y = 0.5, rotation = 0, color = 'dimgrey', labelpad=120)
	plt.subplots_adjust(bottom=0.25, top=0.9, left=0.2, right=0.95, wspace=0.08, hspace = 0.08)
	plt.savefig('/users/ellgil82/figures/Wintertime melt/Re-runs/new_v_old.png', transparent = True)
	plt.savefig('/users/ellgil82/figures/Wintertime melt/Re-runs/new_v_old.eps', transparent = True)
	plt.show()

#plot_difs()

## Caption: Difference between default and updated surface elevation (orography) and land-sea mask (LSM) ancillary files
## used in MetUM runs. The default orography and LSM of the 1.5 km resolution Antarctic Peninsula model domain are based
## on the 1 km resolution GLOBE dataset (NGDC, 1999). The updated land-sea mask is derived from the Scientific Committee
## for Antarctic Research (SCAR) Antarctic Digital Database, and updated orography is based on the Ohio State University
## RAMP 200 m Digital Elevation Model (Hongxing, 1999), produced by Tony Phillips at the British Antarctic Survey. Blue
## colours indicate where a) ice is removed from the updated LSM compared to the default, and b) where surface elevation
## is lower than in the default ancillaries, while red colours indicate the opposite. The location of the Cabinet Inlet
## AWS is shown in pink.

def orog_dif_plot():
	# Load necessary files
	old_orog = iris.load_cube('/data/clivarm/wip/ellgil82/new_ancils/km1p5/orog/orog_original.nc', 'OROGRAPHY (/STRAT LOWER BC)')[0,0,:,:]
	new_orog = iris.load_cube('/data/clivarm/wip/ellgil82/new_ancils/km1p5/orog/new_orog_smoothed.nc', 'Height')[0, 0,:, :]
	old_lsm = iris.load_cube('/data/clivarm/wip/ellgil82/new_ancils/km1p5/lsm/lsm_original.nc', 'LAND MASK (No halo) (LAND=TRUE)')[0, 0, :, :]
	new_lsm = iris.load_cube('/data/clivarm/wip/ellgil82/new_ancils/km1p5/lsm/new_mask.nc', 'LAND MASK (No halo) (LAND=TRUE)')[0,0,:,:]
	# Set up figure
	fig, ax = plt.subplots(1, 1, figsize=(10,11.5))
	ax.spines['right'].set_visible(False)
	ax.spines['left'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.tick_params(axis='both', which='both', length=0, labelbottom='off', labelleft='off')
	# Plot Cabinet Inlet AWS
	cube = iris.load_cube('/data/clivarm/wip/ellgil82/May_2016/Re-runs/CS2/20160522T1200Z_Peninsula_km1p5_ctrl_pa000.pp','surface_altitude')
	real_lon, real_lat = rotate_data(cube, 0, 1)
	ax.plot(-63.37105, -66.48272, markersize=15, marker='o', color='#f68080', zorder=10)
	# Calculate differences
	orog_dif = new_orog.data - old_orog.data
	lsm_dif = new_lsm.data - old_lsm.data
	# Mask data where no difference is seen
	orog_dif = ma.masked_where((lsm_dif == 0) & (new_lsm.data == 0), orog_dif)
	lsm_dif = ma.masked_where(lsm_dif == 0, lsm_dif)
	# Truncate colormap to minimise visual impact of one or two extreme values
	squished_bwr = shiftedColorMap(cmap=matplotlib.cm.bwr, min_val=-800, max_val=800, name='squished_bwr', var=orog_dif, start = .15, stop = .85)
	# Plot differences between old and new orography and LSM
	c = ax.pcolormesh(real_lon, real_lat, orog_dif, cmap='squished_bwr', vmin=-800, vmax=800, zorder = 1)
	lsm = ax.contourf(real_lon, real_lat, lsm_dif, cmap = 'bwr', vmax = 1, vmin = -1, zorder = 2)
	# Add new LSM and 25 m orography contour
	ax.contour(real_lon, real_lat, new_lsm.data, lw=3, colors='dimgrey', zorder = 3)
	ax.contour(real_lon, real_lat, new_orog.data, lw = 2, levels = [100], colors = 'dimgrey', zorder = 4)
	# Set up colour bar
	cbaxes = fig.add_axes([0.22, 0.12, 0.56, 0.03])
	cbticks = np.linspace(-800, 800, 4)
	cbticklabs = [-800, 0, 800]
	cb = plt.colorbar(c, cax=cbaxes, orientation='horizontal', ticks=cbticks)
	cb.set_ticks(cbticks, cbticklabs)
	cb.ax.set_xlabel('Surface elevation difference (m)', fontsize=30, labelpad=20, color='dimgrey')
	cb.ax.text(-0.3, 2.2, 'Area removed \nfrom new LSM', fontsize = 30, color = 'dimgrey')
	cb.ax.text(0.78, 2.2, 'Area added \nto new LSM', fontsize = 30, color = 'dimgrey')
	cb.outline.set_edgecolor('dimgrey')
	cb.outline.set_linewidth(2)
	cb.solids.set_edgecolor('face')
	cb.ax.tick_params(labelsize=30, tick1On=False, tick2On=False, labelcolor='dimgrey', pad=10)
	[l.set_visible(False) for (w, l) in enumerate(cb.ax.xaxis.get_ticklabels()) if w % 3 != 0]
	plt.subplots_adjust(bottom=0.27, left=0.11, right=0.89, top=0.95, hspace=0.05)
	plt.savefig('/users/ellgil82/figures/new_ancils/orog_difs_km1p5.png', transparent = True)
	plt.savefig('/users/ellgil82/figures/new_ancils/orog_difs_km1p5.eps', transparent = True)
	plt.savefig('/users/ellgil82/figures/new_ancils/orog_difs_km1p5.pdf', transparent=True)
	plt.show()

orog_dif_plot()