"""
Grid interpolation utilities

Available functions:
 - interpolate_from_grid_T : interpolated from grid T to scattered grid

Module parameters:
 - stag_type: type of staggered grid ('upright' or 'downleft')

"""
import numpy as np

# Type of staggered grid (upright or downleft)
stag_type='upright'

def interpolate_from_grid_T(phi,grid_type):
  """
  Interpolate from grid T to scattered grid

  Args:
  phi : input field on grid T
  grid_type : type of staggered grid (T, U, V or F)

  Returns:
  phi_interp : angle of rotation
  """
  global stag_type

  # Initialize output array
  phi_interp = np.zeros_like(phi)

  # Case of up-right staggered grid (like NEMO)
  if stag_type == 'upright' :
    if grid_type == 'U':
      phi_interp[:-1,:] = (phi[:-1, :] + phi[1:, :]) / 2
      phi_interp[-1,:] = phi[-1,:]
    elif grid_type == 'V':
      phi_interp[:,:-1] = (phi[:, :-1] + phi[:, 1:]) / 2
      phi_interp[:,-1] = phi[:,-1]
    elif grid_type == 'F':
      phi_interp[:-1,:-1] = ( phi[:-1, :-1] + phi[1:, :-1] + phi[:-1, 1:] + phi[1:, 1:] ) / 4
      phi_interp[-1,:] = phi[-1,:]
      phi_interp[:,-1] = phi[:,-1]
    else:
      raise ValueError("Bad type of staggered grid (T, U , V or F)")
  # Case of down-left staggered grid (like CROCO)
  elif stag_type == 'downleft' :
    if grid_type == 'U':
      phi_interp[1:,:] = (phi[0:, :] + phi[1:, :]) / 2
      phi_interp[0,:] = phi[0,:]
    elif grid_type == 'V':
      phi_interp[:,1:] = (phi[:, 0:] + phi[:, 1:]) / 2
      phi_interp[:,0] = phi[:,0]
    elif grid_type == 'F':
      phi_interp[1:,1:] = ( phi[0:, 0:] + phi[1:, 0:] + phi[0:, 1:] + phi[1:, 1:] ) / 4
      phi_interp[0,:] = phi[0,:]
      phi_interp[:,0] = phi[:,0]
    else:
      raise ValueError("Bad type of staggered grid (upright or downleft)")

  return phi_interp

