import numpy as np
from netCDF4 import Dataset

# Create a 100x100 grid for x and y
grid_size = 100
x = np.linspace(0, 1, grid_size)
y = np.linspace(0, 1, grid_size)
x_grid, y_grid = np.meshgrid(x, y)

# Define a scalar function f(x, y) = y - x
f = y_grid - x_grid

# Define a constant vector field u(x, y) = (1, 1)
u = np.ones((grid_size, grid_size, 2))  # Last dimension is 2 for the vector components

# Define a constant tensor field t(x, y) = ( [a, c] , [c, b] )
# (symmetric, positive definite)
a = 1.2 ; b = 1 ; c = 1
t = np.zeros((grid_size, grid_size, 2, 2))
for i in range(grid_size):
    for j in range(grid_size):
        t[i,j,:,:] = np.array([[a, c],
                               [c, b]])

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

    # Create the vector field variable u(x, y)
    txx_var = ncfile.createVariable('txx', 'f4', ('x', 'y'))
    txy_var = ncfile.createVariable('txy', 'f4', ('x', 'y'))
    tyy_var = ncfile.createVariable('tyy', 'f4', ('x', 'y'))

    # Assign data to the coordinate variables
    x_var[:,:] = x_grid
    y_var[:,:] = y_grid

    # Assign data to the scalar function variable
    f_var[:,:] = f

    # Assign data to the vector field variable
    u_var[:,:] = u[:, :, 0]
    v_var[:,:] = u[:, :, 1]

    # Assign data to the vector field variable
    txx_var[:,:] = t[:, :, 0, 0]
    txy_var[:,:] = t[:, :, 0, 1]
    tyy_var[:,:] = t[:, :, 1, 1]

    # Add  descriptions
    f_var.description = 'Scalar function f(x, y)'
    u_var.description = 'Constant vector field u(x, y)'
    v_var.description = 'Constant vector field v(x, y)'
    txx_var.description = 'Constant tensor field txx(x, y)'
    txy_var.description = 'Constant tensor field txy(x, y)'
    tyy_var.description = 'Constant tensor field txy(y, y)'

print(f"NetCDF file '{nc_filename}' created successfully.")

