import numpy as np
import argparse
from netCDF4 import Dataset
import matplotlib.pyplot as plt

# Parse command-line arguments
parser = argparse.ArgumentParser(prog='plot_test_case',description='Plot output from test case')
parser.add_argument('-i','--input_file', type=str, required=True, help='name of input file')
args = parser.parse_args()

# Open NetCDF file
input_file = args.input_file
nc_file = Dataset(input_file, 'r')
nc_x = nc_file.variables['nav_lon']
nc_y = nc_file.variables['nav_lat']
nc_f = nc_file.variables['f']
nc_u = nc_file.variables['u']
nc_v = nc_file.variables['v']

# Load variables from file
x = nc_x[:]
y = nc_y[:]
f = nc_f[:]
u = nc_u[:]
v = nc_v[:]

# Close NetCDF file
nc_file.close()

# Plotting
plt.figure(figsize=(8, 8))

# Display the scalar field with imshow
plt.imshow(f, extent=[0, 1, 0, 1], origin='lower', cmap='cividis', alpha=0.7)
#plt.imshow(f, extent=[0, 1, 0, 1], origin='lower', alpha=0.7)
plt.contour(x, y, f, colors='white', levels=20, linewidths=1.5, linestyles='solid')

# Overlay the vector field with quiver
plt.quiver(x[2::5,2::5], y[2::5,2::5], u[2::5,2::5], v[2::5,2::5], color='black',pivot='middle',scale=20)

# Add labels and title
plt.xlabel('x')
plt.ylabel('y')
plt.title('Scalar Field with Vector Field')

# Show the plot
plt.colorbar(label='Scalar value',shrink=0.7) # Add color bar to show scalar values
#plt.show()

plt.savefig(input_file+'.png', dpi=100, bbox_inches='tight',pad_inches=0.01)
