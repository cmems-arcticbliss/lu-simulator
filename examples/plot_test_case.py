import numpy as np
import argparse
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse


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
nc_txx = nc_file.variables['txx']
nc_txy = nc_file.variables['txy']
nc_tyy = nc_file.variables['tyy']

# Load variables from file
x = nc_x[:]
y = nc_y[:]
f = nc_f[:]
u = nc_u[:]
v = nc_v[:]
txx = nc_txx[:]
txy = nc_txy[:]
tyy = nc_tyy[:]

# Close NetCDF file
nc_file.close()

# Get dimensions
idim = x.shape[0]
jdim = x.shape[1]

# Plotting
fig, ax = plt.subplots(figsize=(8, 8))

# Display the scalar field with colors
plt.imshow(f, extent=[0, 1, 0, 1], origin='lower', cmap='cividis', alpha=0.7,vmin=-1,vmax=1)
plt.colorbar(label='Scalar value',cmap='cividis',shrink=0.7) # Add color bar to show scalar values

# Display the scalar field with white contours
plt.contour(x, y, f, colors='white', levels=20, linewidths=1.5, linestyles='solid')

# Display the vector field with black arrows
plt.quiver(x[2::5,2::5], y[2::5,2::5], u[2::5,2::5], v[2::5,2::5], color='black',pivot='middle',scale=20)

# Display each tensor as an ellipse
for i in range(9,idim-1,10):
  for j in range(9,jdim-1,10):
    tensor = np.array([[txx[i,j], txy[i,j]],
                       [txy[i,j], tyy[i,j]]])

    # Compute eigenvalues and eigenvectors
    eigenvalues, eigenvectors = np.linalg.eigh(tensor)

    # The angle of the ellipse (orientation)
    angle = np.degrees(np.arctan2(*eigenvectors[:, 0][::-1]))

    # Width and height of the ellipse are proportional to the square root of the eigenvalues
    plot_scale = 15
    width, height = np.sqrt(eigenvalues) / plot_scale

    # Position
    ellipse = Ellipse(xy=(x[i, j], y[i, j]), width=width, height=height,
                      angle=angle, edgecolor='blue', facecolor='none')
    ax.add_patch(ellipse)

# Add labels and title
plt.title('Scalar Field with Vector Field and Tensor field')

# Save figure in file
plt.savefig(input_file+'.png', dpi=100, bbox_inches='tight',pad_inches=0.01)
