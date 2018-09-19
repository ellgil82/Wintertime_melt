''' Tools for processing Unified Model data

Author: Ella Gilbert, 2018.
'''

## Function to rotate iris cubes plotted on a rotated pole grid to regular lats/lons

def rotate_data(var, lat_dim, lon_dim):
    ## Rotate projection
    # create numpy arrays of coordinates
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

## Function to find indices of the inputted lat/lon coordinates to place real-world coordinates on model grid

def find_gridbox(x, y, real_lat, real_lon): 
    global lon_index, lat_index
    lat_index = np.argmin((real_lat - x) ** 2)  # take whole array and subtract lat you want from
    lon_index = np.argmin((real_lon - y) ** 2)  # each point, then find the smallest difference
    return lon_index, lat_index