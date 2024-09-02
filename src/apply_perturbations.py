"""
Apply random location perturbation to geophysical variables

Available functions:
 - apply_perturbation: apply a given perturbation to a given variable

Module parameters:

"""
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import apply_rotations as app_rot
import util_grid

no_extrapolation=False  # No extrapolation beyond grid limits

def apply_perturbation(reference,dx,dy,grid_type='T'):
  """
  Apply a given perturbation to a given variable

  Args:
  reference : reference variable to be perturbed
  dx : perturbation along x coordinate
  dy : perturbation along y coordinate
  grid_type : type of staggered grid (T, U, V or F)

  Returns:
  perturbed: perturbed variable
  """
  global no_extrapolation

  # Dimension of input variable
  ndim = reference.ndim

  # Define horizontal grid coordinates
  x = np.arange(0, reference.shape[ndim-1], dtype=float)
  y = np.arange(0, reference.shape[ndim-2], dtype=float)

  # Define reference output grid (to be perturbed)
  xx, yy = np.meshgrid(x, y)

  # Interpolate perturbation on staggered grid (if not T)
  # and shift location of grid points
  if grid_type != 'T':
    dxi = util_grid.interpolate_from_grid_T(dx,grid_type)
    dyi = util_grid.interpolate_from_grid_T(dy,grid_type)
    points = (yy+dyi, xx+dxi)
  else:
    points = (yy+dy, xx+dx)

  # Initialize unmasked field
  perturbed = np.empty_like(reference)

  # Interpolate at new locations separately for every 2d slice of the array
  if ndim == 2 :
    slice2d = reference[:,:]
    if no_extrapolation:
      interpolator = RegularGridInterpolator((y,x), slice2d, method='linear', bounds_error=False)
    else:
      interpolator = RegularGridInterpolator((y,x), slice2d, method='linear', bounds_error=False, fill_value=None)
    perturbed = interpolator(points)
  elif ndim == 3 :
    for k in range(reference.shape[0]):
      slice2d = reference[k,:,:]
      if no_extrapolation:
        interpolator = RegularGridInterpolator((y,x), slice2d, method='linear', bounds_error=False, fill_value=None)
      else:
        interpolator = RegularGridInterpolator((y,x), slice2d, method='linear', bounds_error=False, fill_value=None)
      perturbed[k,:,:] = interpolator(points)
  elif ndim == 4 :
    for l in range(reference.shape[0]):
      for k in range(reference.shape[1]):
        slice2d = reference[l,k,:,:]
        if no_extrapolation:
          interpolator = RegularGridInterpolator((y,x), slice2d, method='linear', bounds_error=False, fill_value=None)
        else:
          interpolator = RegularGridInterpolator((y,x), slice2d, method='linear', bounds_error=False, fill_value=None)
        perturbed[l,k,:,:] = interpolator(points)
  else:
    raise ValueError("Bad variable dimension")

  return perturbed

if __name__ == "__main__":
  import os
  import shutil
  import argparse
  import netcdf_io as ncio
  import netcdf_format as ncf
  import vars_def as vdef

  # Set type of staggered grid
  util_grid.stag_type=ncf.stag_type

  # Parse command-line arguments
  parser = argparse.ArgumentParser(prog='apply_perturbations',description='Apply perturbations to geophysical fields')
  parser.add_argument('-m',      '--sample_size',      type=int,   required=True,  help='sample size')
  parser.add_argument('-i',      '--input_file',       type=str,   required=True,  help='name of input file with perturbations')
  parser.add_argument('-std',    '--perturbation_std', type=float, required=True,  help='standard deviation of perturbations')
  parser.add_argument('-grid',   '--grid_file',        type=str,   required=True,  help='name of grid file')
  parser.add_argument('-ref',    '--reference_file',   type=str,   required=True,  help='name of input file with reference fields')
  parser.add_argument('-v',      '--varlist_file',     type=str,   required=True,  help='file with list of variables to perturb')
  parser.add_argument('-o',      '--output_file',      type=str,   required=True,  help='name of output file with perturbed fields')
  parser.add_argument('-mask',   '--mask_file',        type=str,   required=False, help='name of mask file')
  parser.add_argument('-factor', '--damping_factor',   type=str,   required=False, help='file with damping factor to apply to perturbations')
  parser.add_argument('-N',      '--nproc',            type=int,   required=False, help='number of processors to use with MPI')
  args = parser.parse_args()

  # MPI parallelization settings
  use_mpi = args.nproc is not None
  if use_mpi:
    from mpi4py import MPI
    mpicomm = MPI.COMM_WORLD
    mpirank = mpicomm.Get_rank()
    mpisize = mpicomm.Get_size()
    if mpisize != args.nproc:
      raise ValueError("Inconsistent number of processors")
  else:
    mpisize=1
    mpirank=0

  # Get list of variables from file
  vdef.read_vars(args.varlist_file)

  # Loop on sample/ensemble size
  for imember in range(mpirank, args.sample_size, mpisize):
    print('Applying perturbation for member: ',imember)

    # Get perturbation of grid locations
    dx, dy = ncio.read_perturbation(args.input_file,imember)

    # Scale perturbation with standard deviation
    dx = dx * args.perturbation_std
    dy = dy * args.perturbation_std

    # Define output file name with ensemble member index as a prefix
    output_file=f'{imember+1:0>3}' + args.output_file
    # Check existence of output file
    # If not, copy reference file to output file
    if not os.path.exists(output_file):
      shutil.copy(args.reference_file, output_file)
    # Open output file for this member
    ncio.open_variable_file(output_file)

    # Initialize vector and tensor component index
    component = 0

    # Loop on variables
    for ivar, variable in enumerate(vdef.var_list):
      varname = variable['name']
      print('Applying perturbation for variable: ',varname)

      # Read this variable from reference file
      reference_field = ncio.read_variable(args.reference_file,varname)

      # Read and apply damping factor for this variable if any
      if args.damping_factor is not None:
        damping_factor = ncio.read_variable(args.damping_factor,varname)
        dx = dx * damping_factor
        dy = dy * damping_factor

      # Apply perturbation to this variable
      perturbed_field = apply_perturbation(reference_field,dx,dy,grid_type=variable['grid_type'])
      # Warning: variables may need proper remasking before writing

      if variable['tensor_type'] == 'scalar' :
        # Read and apply mask for this variable if any
        if args.mask_file is not None:
          varmask = ncio.read_variable(args.mask_file,varname)
          perturbed_field[np.where(varmask==ncf.mask_spval)] = ncf.output_spval
        # Save scalar variable to file
        ncio.write_variable(perturbed_field,varname)
      elif variable['tensor_type'] == 'vector' :
        component = component + 1
        if component == 1 :
          u_name = variable['name'] ; u = perturbed_field ; u_type = variable['grid_type']
        else:
          v_name = variable['name'] ; v = perturbed_field ; v_type = variable['grid_type']
          component = 0
          # Rotate vector once all components have been perturbed
          # (vector components are assumed on grid U and V)
          if u_type == v_type:
            # vector components are on the same grid
            app_rot.rotate_vector(u,v,dx,dy,grid_type=u_type)
          else:
            # vector components are assumed on grid U and V
            if u_type != 'U':
              raise ValueError("Bad grid type of vector components")
            if v_type != 'V':
              raise ValueError("Bad grid type of vector components")
            app_rot.rotate_vector(u,v,dx,dy,grid_type='UV')
          # Read and apply mask for this variable if any
          if args.mask_file is not None:
            varmask = ncio.read_variable(args.mask_file,u_name)
            u[np.where(varmask==ncf.mask_spval)] = ncf.output_spval
            varmask = ncio.read_variable(args.mask_file,v_name)
            v[np.where(varmask==ncf.mask_spval)] = ncf.output_spval
          # Save vector components to file
          ncio.write_variable(u,u_name)
          ncio.write_variable(v,v_name)
      elif variable['tensor_type'] == 'tensor' :
        component = component + 1
        if component == 1 :
          txx_name = variable['name'] ; txx = perturbed_field
        elif component == 2 :
          txy_name = variable['name'] ; txy = perturbed_field
        else:
          tyy_name = variable['name'] ; tyy = perturbed_field
          component = 0
          # Rotate tensor once all components have been perturbed
          # (all tensor components are assumed on the same grid)
          app_rot.rotate_tensor(txx,txy,tyy,dx,dy,grid_type=variable['grid_type'])
          # Read and apply mask for this variable if any
          if args.mask_file is not None:
            varmask = ncio.read_variable(args.mask_file,txx_name)
            txx[np.where(varmask==ncf.mask_spval)] = ncf.output_spval
            varmask = ncio.read_variable(args.mask_file,txy_name)
            txy[np.where(varmask==ncf.mask_spval)] = ncf.output_spval
            varmask = ncio.read_variable(args.mask_file,tyy_name)
            tyy[np.where(varmask==ncf.mask_spval)] = ncf.output_spval
          # Save tensor components to file
          ncio.write_variable(txx,txx_name)
          ncio.write_variable(txy,txy_name)
          ncio.write_variable(tyy,tyy_name)
      else:
        raise ValueError("Bad variable tensor type")

    # Close output file for this member
    ncio.close_variable_file()
