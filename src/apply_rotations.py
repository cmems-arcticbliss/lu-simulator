"""
Apply rotation to vectors and tensors to follow the perturbation of the coordinates

Available functions:
 - rotate_vector: rotate vector to follow the perturbation of the coordinates
 - rotate_tensor: rotate symmetric tensor to follow the perturbation of the coordinates

Module parameters:

"""
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import util_grid

def rotate_vector(u,v,dx,dy):
  """
  Rotate vectors to follow the perturbation of the coordinates

  Args:
  u,v : components of the vector
  dx : perturbation along x coordinate
  dy : perturbation along y coordinate

  Returns:
  u,v are modified in place
  """

  # Dimension of input variable
  ndim = u.ndim

  # Calculate the angle of rotation at U points
  theta_u = compute_rotation_angle(dx,dy,grid_type='U')

  # Calculate the angle of rotation at V points
  theta_v = compute_rotation_angle(dx,dy,grid_type='V')

  # Rotate vector for every 2d slice
  if ndim == 2 :
    rotate_vector_2d(u,v,theta_u,theta_v)
  elif ndim == 3 :
    for k in range(u.shape[0]):
      u2d = u[k,:,:]
      v2d = v[k,:,:]
      rotate_vector_2d(u2d,v2d,theta_u,theta_v)
      u[k,:,:] = u2d
      v[k,:,:] = v2d
  elif ndim == 4 :
    for l in range(u.shape[0]):
      for k in range(u.shape[1]):
        u2d = u[l,k,:,:]
        v2d = v[l,k,:,:]
        rotate_vector_2d(u2d,v2d,theta_u,theta_v)
        u[l,k,:,:] = u2d
        v[l,k,:,:] = v2d
  else:
    raise ValueError("Bad variable dimension")

def rotate_vector_2d(u,v,theta_u,theta_v):
  """
  Rotate vectors to follow the perturbation of the coordinates

  Args:
  u,v : components of the vector
  theta_u, theta_v : angle of rotation at U and V points

  Returns:
  u,v are modified in place
  """

  # Define working arrays
  u_v = np.zeros_like(u)
  v_u = np.zeros_like(v)

  # Compute u at V points
  u_v[1:,1:] = ( u[0:-1, 0:-1] + u[1:, 0:-1] + u[0:-1, 1:] + u[1:, 1:] ) / 4
  u_v[0,:] = u[0,:] ; u_v[:,0] = u[:,0]

  # Compute v at U points
  v_u[1:,1:] = ( v[0:-1, 0:-1] + v[1:, 0:-1] + v[0:-1, 1:] + v[1:, 1:] ) / 4
  v_u[0,:] = v[0,:] ; v_u[:,0] = v[:,0]

  # Apply rotation angle to input vector
  u =   u   * np.cos(theta_u) +  v_u * np.sin(theta_u)
  v = - u_v * np.sin(theta_v) +  v   * np.cos(theta_v)

def rotate_tensor(txx,txy,tyy,dx,dy,grid_type='T'):
  """
  Rotate symmetric tensor to follow the perturbation of the coordinates

  Args:
  txx,txy,tyy : components of the tensor
  dx : perturbation along x coordinate
  dy : perturbation along y coordinate
  grid_type : type of staggered grid (T, U, V or F)

  Returns:
  rtxx,rtxy,rtyy : transformed components of the tensor
  """

  # Dimension of input variable
  ndim = txx.ndim

  # Calculate the angle of rotation at appropriate grid points
  theta = compute_rotation_angle(dx,dy,grid_type)

  # Rotate tensor for every 2d slice
  if ndim == 2 :
    rtxx, rtxy, rtyy = rotate_tensor_2d(txx,txy,tyy,theta)
  elif ndim == 3 :
    for k in range(txx.shape[0]):
      txx2d = txx[k,:,:]
      txy2d = txy[k,:,:]
      tyy2d = tyy[k,:,:]
      rtxx, rtxy, rtyy = rotate_tensor_2d(txx2d,txy2d,tyy2d,theta)
      rtxx[k,:,:] = txx2d
      rtxy[k,:,:] = txy2d
      rtyy[k,:,:] = tyy2d
  elif ndim == 4 :
    for l in range(txx.shape[0]):
      for k in range(txx.shape[1]):
        txx2d = txx[l,k,:,:]
        txy2d = txy[l,k,:,:]
        tyy2d = tyy[l,k,:,:]
        rtxx, rtxy, rtyy = rotate_tensor_2d(txx2d,txy2d,tyy2d,theta)
        rtxx2d = txx[l,k,:,:]
        rtxy2d = txy[l,k,:,:]
        rtyy2d = tyy[l,k,:,:]
  else:
    raise ValueError("Bad variable dimension")

  return rtxx, rtxy, rtyy

def rotate_tensor_2d(txx,txy,tyy,theta):
  """
  Rotate symmetric tensor to follow the perturbation of the coordinates

  Args:
  txx,txy,tyy : components of the tensor
  theta : angle of rotation at T points

  Returns:
  rtxx,rtxy,rtyy : transformed components of the tensor
  """

  # Apply rotation angle to input tensor
  cos_theta = np.cos(theta) ; sin_theta = np.sin(theta)
  rtxx = txx * cos_theta**2 + 2 * txy * cos_theta * sin_theta + tyy * sin_theta**2
  rtxy = ( tyy - txx ) * sin_theta * cos_theta + txy * ( cos_theta**2 - sin_theta**2 )
  rtyy = txx * sin_theta**2 - 2 * txy * cos_theta * sin_theta + tyy * cos_theta**2

  return rtxx, rtxy, rtyy

def compute_rotation_angle(dx,dy,grid_type='T'):
  """
  Compute rotation angle associated to perturbation

  Args:
  dx : perturbation along x coordinate
  dy : perturbation along y coordinate
  grid_type : type of staggered grid (T, U, V or F)

  Returns:
  theta : angle of rotation
  """

  # Initialize rotation angle
  theta = np.zeros_like(dx)

  # Set shifts to use for each type of grid
  if grid_type == 'T':
    ish = 2 ; jsh = 2
  elif grid_type == 'U':
    ish = 1 ; jsh = 2
  elif grid_type == 'V':
    ish = 2 ; jsh = 1
  elif grid_type == 'F':
    ish = 1 ; jsh = 1
  else:
    raise ValueError("Bad type of staggered grid (T, U , V or F)")

  # We remove translation by computing difference of dx and dy
  # (dx and dy are always assumed to correspond to grid T)
  ilo = ish - 1 ; jlo = jsh - 1
  dx_x = ( dx[ish:,jlo:-1] - dx[:-ish,jlo:-1] ) / ish
  dx_y = ( dx[ilo:-1,jsh:] - dx[ilo:-1,:-jsh] ) / jsh
  dy_x = ( dy[ish:,jlo:-1] - dy[:-ish,jlo:-1] ) / ish
  dy_y = ( dy[ilo:-1,jsh:] - dy[ilo:-1,:-jsh] ) / jsh

  # The transformation of the coordinates is then T = I + dx'/dx
  # where T is the combination of a rotation and deformation: T = R D.
  # The rotation angle is computed to obtain a symmetric matrix D from T:
  # tan(theta) = [ (Tyx - Txy)/2 ] / [ 1 + (Txx + Tyy)/2 ]
  trace = ( dx_x + dy_y ) / 2
  asymm = ( dy_x - dx_y ) / 2

  theta[ilo:-1,jlo:-1] = np.arctan2( asymm, 1 + trace )

  # Extrapolate to boundaries
  theta[:ilo,:] = theta[ilo,:][np.newaxis, :]
  theta[-1,:]   = theta[-2,:]
  theta[:,:jlo] = theta[:,jlo][:, np.newaxis]
  theta[:,-1]   = theta[:,-2]

  return theta

if __name__ == "__main__":
  import argparse
  import netcdf_io as ncio
  import netcdf_format as ncf

  # Parse command-line arguments
  parser = argparse.ArgumentParser(prog='apply_rotations',description='Compute rotation angle from perturbation')
  parser.add_argument('-m',      '--sample_size',      type=int,   required=True,  help='sample size')
  parser.add_argument('-i',      '--input_file',       type=str,   required=True,  help='name of input file with perturbations')
  parser.add_argument('-grid',   '--grid_file',        type=str,   required=True,  help='name of grid file')
  parser.add_argument('-std',    '--perturbation_std', type=float, required=True,  help='standard deviation of perturbations')
  parser.add_argument('-o',      '--output_file',      type=str,   required=True,  help='name of output file with perturbed fields')
  args = parser.parse_args()

  # Read grid from grid file
  lon1d, lat1d = ncio.read_grid(args.grid_file)

  # Create output file
  ncio.create_rotation_file(args.output_file,args.sample_size)

  # Loop on sample/ensemble size
  for imember in range(args.sample_size):
    print('  Computing rotation angle for member: ',imember)

    # Get perturbation of grid locations
    dx, dy = ncio.read_perturbation(args.input_file,imember)

    # Scale perturbation with standard deviation
    dx = dx * args.perturbation_std
    dy = dy * args.perturbation_std

    # Compute rotation angle
    theta = compute_rotation_angle(dx,dy)

    # Save perturbation in output file
    print('  Rotation angle stored in file: ',args.output_file)
    ncio.write_rotation(imember,theta)

  # Close output file
  ncio.close_rotation_file()

