import numpy as np
from netCDF4 import Dataset

# Create a 100x100 grid for x and y
grid_size = 100
x = np.linspace(0, 1, grid_size)
y = np.linspace(0, 1, grid_size)
x_grid, y_grid = np.meshgrid(x, y)

# Define the scalar function f(x, y) = y - x
f = y_grid - x_grid
#f = np.floor( 10 *  ( y_grid - x_grid ) ) / 10
#f = y_grid

# Define the constant vector field u = (1, 1)
u = np.ones((grid_size, grid_size, 2))  # Last dimension is 2 for the vector components

# Create a new NetCDF file
nc_filename = 'test.nc'
with Dataset(nc_filename, 'w', format='NETCDF4') as ncfile:
    
    # Create the dimensions
    ncfile.createDimension('x', grid_size)
    ncfile.createDimension('y', grid_size)

    # Create the coordinate variables
    x_var = ncfile.createVariable('nav_lon', 'f4', ('x','y'))
    y_var = ncfile.createVariable('nav_lat', 'f4', ('x','y'))

    # Create the scalar function variable f(x, y)
    f_var = ncfile.createVariable('f', 'f4', ('x', 'y'))

    # Create the vector field variable u(x, y)
    u_var = ncfile.createVariable('u', 'f4', ('x', 'y'))
    v_var = ncfile.createVariable('v', 'f4', ('x', 'y'))

    # Assign data to the coordinate variables
    x_var[:,:] = x_grid
    y_var[:,:] = y_grid

    # Assign data to the scalar function variable
    f_var[:, :] = f

    # Assign data to the vector field variable
    u_var[:, :] = u[:, :, 0]
    v_var[:, :] = u[:, :, 1]

    # Optional: Add units or descriptions
    x_var.units = 'nondimensional'
    y_var.units = 'nondimensional'
    f_var.units = 'nondimensional'
    u_var.units = 'nondimensional'
    v_var.units = 'nondimensional'
    f_var.description = 'Scalar function f(x, y) = y - x'
    u_var.description = 'Constant vector field u(x, y) = 1'
    v_var.description = 'Constant vector field v(x, y) = 1'

print(f"NetCDF file '{nc_filename}' created successfully.")

