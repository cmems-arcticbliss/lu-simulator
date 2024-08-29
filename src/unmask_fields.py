import numpy as np
import pyensdam as edam

# Unmask field according to mask
def unmask_field(reference,maxrange=1000,window=5,spval=-9999.):

  # Define unmasking parameters
  edam.interpolation.unmask_spval = spval
  edam.interpolation.unmask_max = maxrange
  edam.interpolation.unmask_window = window

  # Initialize unmasked field
  unmasked_field = np.empty_like(reference)

  # Apply unmasking algorithm separately on every 2d slice of the array
  ndim = reference.ndim
  if ndim == 2 :
    slice2d = reference[:,:]
    edam.interpolation.unmask2D(slice2d)
    unmasked_field[:,:] = slice2d
  elif ndim == 3 :
    for k in range(reference.shape[0]):
      slice2d = reference[k,:,:]
      edam.interpolation.unmask2D(slice2d)
      unmasked_field[k,:,:] = slice2d
  elif ndim == 4 :
    for l in range(reference.shape[0]):
      for k in range(reference.shape[1]):
        slice2d = reference[l,k,:,:]
        edam.interpolation.unmask2D(slice2d)
        unmasked_field[l,k,:,:] = slice2d
  else:
    raise ValueError("Bad variable dimension")

  return unmasked_field

# Apply the module to NetCDF files
if __name__ == "__main__":
  import os
  import shutil
  import argparse
  import netcdf_io as ncio
  import netcdf_format as ncf
  import vars_def as vdef

  # Parse command-line arguments
  parser = argparse.ArgumentParser(prog='unmask_field',description='Unmask reference field')
  parser.add_argument('-ref',  '--reference_file', type=str,   required=True,  help='name of input file with reference fields')
  parser.add_argument('-v',    '--varlist_file',   type=str,   required=True,  help='file with list of variables to perturb')
  parser.add_argument('-mask', '--mask_file',      type=str,   required=True,  help='name of mask file')
  parser.add_argument('-o',    '--output_file',    type=str,   required=True,  help='name of output file with unmasked reference field')
  args = parser.parse_args()

  # Get list of variables from file
  vdef.read_vars(args.varlist_file)

  # Check existence of output file
  # If not, copy reference file to output file
  if not os.path.exists(args.output_file):
    shutil.copy(args.reference_file, args.output_file)

  # Open output file
  ncio.open_variable_file(args.output_file)

  # Loop on variables
  for ivar, variable in enumerate(vdef.var_list):
    varname = variable['name']

    # Get array from reference file
    reference = ncio.read_variable(args.reference_file,varname)

    # Get mask from mask file
    varmask = ncio.read_variable(args.mask_file,varname)

    # Set reference field to special value on land
    reference[np.where(varmask==ncf.mask_spval)] = -9999.

    # Unmask reference field
    unmasked_array = unmask_field(reference)

    # Write unmasked reference field
    ncio.write_variable(unmasked_array,varname)

  # Close output file
  ncio.close_variable_file()
