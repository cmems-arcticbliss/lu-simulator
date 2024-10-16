
# Definition of NetCDF format for grid, mask and input fields

xdim='x'
ydim='y'
lonvar='nav_lon'
latvar='nav_lat'

# Special value for masked variables in mask file
mask_spval=0.

# Special value for masked variables in output file
output_spval=0.

# Definition of internal file format
member='member'
location_index='location_index'
dx='dx'
dy='dy'
theta='theta'

# Type of staggered grid (upright or downleft)
stag_type='upright'


# name of grid cell size in x direction in grid/mesh file
gridxname='e1t'
gridyname='e2t'
cfactor=1e-3   # conversion factor from initial unit as given in input grid file and desired unit in damping file  (example: 1e-3 will convert from meters (grid cell size as given in grid file) to kilometers as we want to have the multiplicative factor in the damping file
