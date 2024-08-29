from netCDF4 import Dataset
import netcdf_format as ncf
import numpy as np

# Shape of grid array
grid_shape = None

# Identifier for NetCDF file with perturbations
ncpert_file = None
# Identifier for dx and dy variables in perturbation file
nc_dx = None
nc_dy = None

# Identifier for NetCDF file with geophysical variables
ncvar_file = None

# Identifier for NetCDF file with rotations (for diagnostics)
ncrotation_file = None
# Identifier for dx and dy variables in perturbation file
nc_theta = None

# Read grid from grid file
def read_grid(grid_file):
  global grid_shape

  # Open NetCDF file
  ncgrd_file = Dataset(grid_file, 'r')
  ncgrd_lon = ncgrd_file.variables[ncf.lonvar]
  ncgrd_lat = ncgrd_file.variables[ncf.latvar]

  # Save grid shape for future use
  grid_shape = ncgrd_lon.shape

  # Read grid array and
  # flatten grid locations into 1D vector
  coord_array = ncgrd_lon[:]
  lon1d = coord_array.astype(np.float64).flatten()
  coord_array = ncgrd_lat[:]
  lat1d = coord_array.astype(np.float64).flatten()

  # Close NetCDF file
  ncgrd_file.close()

  return lon1d,lat1d

# Save grid in file
def save_grid(grid_file,lon1d,lat1d):

  # Open NetCDF file
  ncgrd_file = Dataset(grid_file, 'w', format='NETCDF4')
  ncgrd_file.createDimension(ncf.location_index, lon1d.shape[0])
  ncgrd_lon = ncgrd_file.createVariable(ncf.lonvar, lon1d.dtype, (ncf.location_index))
  ncgrd_lat = ncgrd_file.createVariable(ncf.latvar, lat1d.dtype, (ncf.location_index))

  # Flatten grid locations into 1D vector
  ncgrd_lon[:] = lon1d
  ncgrd_lat[:] = lat1d

  # Close NetCDF file
  ncgrd_file.close()

# Read perturbation corresponding to one single ensemble member
def read_perturbation(input_file,imember,multiple_files=False):

  # Open NetCDF file
  nc_file = Dataset(input_file, 'r')

  nc_dx = nc_file.variables[ncf.dx]
  nc_dy = nc_file.variables[ncf.dy]

  # Load variables for requested ensemble member
  if multiple_files:
    dx = nc_dx[imember,:]
    dy = nc_dy[imember,:]
  else:
    dx = nc_dx[imember,:,:]
    dy = nc_dy[imember,:,:]

  # Close NetCDF file
  nc_file.close()

  return dx, dy

# Create perturbation file
def create_perturbation_file(output_file,ensemble_size,multiple_files=False,vector_size=1):
  global ncpert_file, nc_dx, nc_dy, grid_shape

  # Open NetCDF file
  ncpert_file = Dataset(output_file,'w', format='NETCDF4')

  if multiple_files:
    # Create dimensions and variables
    ncpert_file.createDimension(ncf.member, ensemble_size)
    ncpert_file.createDimension(ncf.location_index, vector_size)
    nc_dx = ncpert_file.createVariable(ncf.dx, np.double, (ncf.member,ncf.location_index))
    nc_dy = ncpert_file.createVariable(ncf.dy, np.double, (ncf.member,ncf.location_index))
  else:
    # Create dimensions and variables
    ncpert_file.createDimension(ncf.member, ensemble_size)
    ncpert_file.createDimension(ncf.xdim, grid_shape[1])
    ncpert_file.createDimension(ncf.ydim, grid_shape[0])
    nc_dx = ncpert_file.createVariable(ncf.dx, np.double, (ncf.member,ncf.ydim,ncf.xdim))
    nc_dy = ncpert_file.createVariable(ncf.dy, np.double, (ncf.member,ncf.ydim,ncf.xdim))

# Write perturbation corresponding to one single ensemble member
def write_perturbation(imember,dx,dy,multiple_files=False):
  global ncpert_file, nc_dx, nc_dy, grid_shape

  if multiple_files:
    # Save perturbation arrays
    nc_dx[imember,:] = dx
    nc_dy[imember,:] = dy
  else:
    # Save perturbation arrays
    nc_dx[imember,:,:] = dx.reshape(grid_shape)
    nc_dy[imember,:,:] = dy.reshape(grid_shape)

# Close perturbation file
def close_perturbation_file():
  global ncpert_file

  ncpert_file.close()

# Create rotation file
def create_rotation_file(output_file,ensemble_size):
  global ncrotation_file, nc_theta, grid_shape

  # Open NetCDF file
  ncrotation_file = Dataset(output_file,'w', format='NETCDF4')

  # Create dimensions and variables
  ncrotation_file.createDimension(ncf.member, ensemble_size)
  ncrotation_file.createDimension(ncf.xdim, grid_shape[1])
  ncrotation_file.createDimension(ncf.ydim, grid_shape[0])
  nc_theta = ncrotation_file.createVariable(ncf.theta, np.double, (ncf.member,ncf.ydim,ncf.xdim))

# Write rotation angle corresponding to one single ensemble member
def write_rotation(imember,theta):
  global ncrotation_file, nc_theta, grid_shape

  nc_theta[imember,:,:] = theta.reshape(grid_shape)

# Close rotation file
def close_rotation_file():
  global ncrotation_file

  ncrotation_file.close()

# Read geophysical variable from input file
def read_variable(input_file,varname):

  # Open NetCDF file
  nc_file = Dataset(input_file, 'r')
  nc_var = nc_file.variables[varname]

  # Load variables from file
  variable = nc_var[:]

  # Close NetCDF file
  nc_file.close()

  return variable

# Open file with geophysical variables (for overwriting in it)
def open_variable_file(output_file):
  global ncvar_file

  # Open NetCDF file
  ncvar_file = Dataset(output_file,'r+', format='NETCDF4')

# Write geophysical variable in output file
def write_variable(variable,varname):
  global ncvar_file

  # Get variable identifier from variable name
  nc_var = ncvar_file.variables[varname]
  # Save variable in file
  nc_var[:] = variable

# Close file with geophysical variables
def close_variable_file():
  global ncvar_file

  ncvar_file.close()

