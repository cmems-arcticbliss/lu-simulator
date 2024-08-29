import numpy as np
import pyensdam as edam

# Compute damping factor from mask array
def compute_damping_factor(mask,mask_spval,length_scale,maxrange=1000,window=5):

  # Define unmasking parameters
  edam.interpolation.unmask_spval = 0.
  edam.interpolation.unmask_max = maxrange
  edam.interpolation.unmask_window = window

  # Define length scale of damping factor buffer region
  edam.interpolation.unmask_damping = length_scale

  # Set factor to 1 on masked values and 0 elsewhere
  damping_factor = np.zeros_like(mask)                # do not damp
  damping_factor[np.where(mask==mask_spval)] = 1  # damp on mask

  # Extrapolate inside the domain to generate the damping factor
  ndim = mask.ndim
  if ndim == 2 :
    slice2d = damping_factor[:,:]
    edam.interpolation.unmask2D(sliced2d)
    damping_factor[:,:] = slice2d
  elif ndim == 3 :
    for k in range(mask.shape[0]):
      slice2d = damping_factor[k,:,:]
      edam.interpolation.unmask2D(slice2d)
      damping_factor[k,:,:] = slice2d
  elif ndim == 4 :
    for l in range(mask.shape[0]):
      for k in range(mask.shape[1]):
        slice2d = damping_factor[l,k,:,:]
        edam.interpolation.unmask2D(slice2d)
        damping_factor[l,k,:,:] = slice2d
  else:
    raise ValueError("Bad variable dimension")

  # Reverse -> from 0 to 1
  damping_factor = 1 - damping_factor

  return damping_factor

if __name__ == "__main__":
  import os
  import shutil
  import argparse
  import netcdf_io as ncio
  import netcdf_format as ncf
  import vars_def as vdef

  # Parse command-line arguments
  parser = argparse.ArgumentParser(prog='generate_damping_factor',description='Generate damping factor to apply to perturbations')
  parser.add_argument('-mask', '--mask_file',     type=str,   required=True,  help='name of mask file')
  parser.add_argument('-v',    '--varlist_file',  type=str,   required=True,  help='file with list of variables to perturb')
  parser.add_argument('-o',    '--output_file',   type=str,   required=True,  help='name of output file with damping factor')
  parser.add_argument('-l',    '--length_scale',  type=float, required=True,  help='width of the buffer zone along the mask')
  args = parser.parse_args()

  # Get list of variables from file
  vdef.read_vars(args.varlist_file)

  # Check existence of output file
  # If not, copy mask file to output file
  if not os.path.exists(args.output_file):
    shutil.copy(args.mask_file, args.output_file)

  # Open output file
  ncio.open_variable_file(args.output_file)

  # Loop on variables
  for ivar, variable in enumerate(vdef.var_list):
    varname = variable['name']

    # Get mask from mask file
    mask = ncio.read_variable(args.mask_file,varname)

    # Compute damping factor
    damping_factor = compute_damping_factor(mask,ncf.mask_spval,args.length_scale)

    # Write damping_factor in file
    ncio.write_variable(damping_factor,varname)

  # Close output file
  ncio.close_variable_file()

