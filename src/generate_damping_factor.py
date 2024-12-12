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

def resetmask(factor,mask,mask_spval):
    """ Reset factor to 1 on masked values and 0 elsewhere"""
    factor[np.where(mask==mask_spval)] = 0  # damp on mask
    return factor

# Function to create the 1-2-1 smoothing kernel
def get_smoothing_kernel():
    """Returns a normalized 1-2-1 smoothing kernel."""
    kernel = np.array([[1, 2, 1],
                       [2, 4, 2],
                       [1, 2, 1]]) / 16.0
    return kernel

# Function to apply the smoothing filter multiple times on a 2D slice
def apply_smoothing(data, kernel, num_iterations):
    """Apply smoothing filter num_iterations times on 2D data."""
    smoothed_data = data
    for _ in range(num_iterations):
        smoothed_data = scipy.ndimage.convolve(smoothed_data, kernel)
    return smoothed_data

# Function to apply smoothing on 2D slices for multi-dimensional arrays
def smooth_multidimensional(data, kernel, num_iterations):
    """
    Applies smoothing on the last two dimensions (y,x) of a 2D or 3D (or higher) array.
    """
    if data.ndim == 2:
        # If the array is 2D, apply smoothing directly
        return apply_smoothing(data, kernel, num_iterations)
    else:
        # If the array is 3D or higher, apply smoothing over the last two dimensions (y,x)
        smoothed_data = np.empty_like(data)
        
        # Iterate over all other dimensions except the last two (y,x)
        for idx in np.ndindex(data.shape[:-2]):
            # Select the y,x slice
            slice_2d = data[idx]
            smoothed_data[idx] = apply_smoothing(slice_2d, kernel, num_iterations)
        
        return smoothed_data

if __name__ == "__main__":
    import os
    import shutil
    import argparse
    import netcdf_io as ncio
    import netcdf_format as ncf
    import vars_def as vdef
    import numpy as np
    import scipy.ndimage
    
    # Parse command-line arguments
    parser = argparse.ArgumentParser(prog='generate_damping_factor', description='Generate damping factor to apply to perturbations')
    parser.add_argument('-mask', '--mask_file', type=str, required=True, help='name of mask file')
    parser.add_argument('-v', '--varlist_file', type=str, required=True, help='file with list of variables to perturb')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='name of output file with damping factor')
    parser.add_argument('-l', '--length_scale', type=float, required=True, help='width of the buffer zone along the mask')
    parser.add_argument('-s', '--smoothing_iterations', type=int, default=0, help='Number of smoothing iterations (optional)')
    parser.add_argument('-kmscaling', '--kmscaling', type=str,required=False,  help='Mesh file to get the size of the grid cells if scaling in km is required (optional)')
    args = parser.parse_args()

    # Get the smoothing kernel (can be reused across variables)
    kernel = get_smoothing_kernel()

    # Get list of variables from file
    vdef.read_vars(args.varlist_file)

    # Check existence of output file, if not, copy mask file to output file
    if not os.path.exists(args.output_file):
        shutil.copy(args.mask_file, args.output_file)

    # Open output file
    ncio.open_variable_file(args.output_file)

    # Read gridcell size from mesh file if required
    if args.kmscaling is not None:
        gridx = ncio.read_variable(args.kmscaling, ncf.gridxname)
        gridy = ncio.read_variable(args.kmscaling, ncf.gridyname)
        gridsize = 0.5*gridx + 0.5*gridy   # average grid size in meters
        gridsize = np.squeeze(gridsize*ncf.cfactor)   # unit conversion factor if needed
        #print(gridsize.shape)

    # Loop on variables
    for ivar, variable in enumerate(vdef.var_list):
        varname = variable['name']

        # Get mask from mask file
        mask = ncio.read_variable(args.mask_file, varname)

        # Compute damping factor
        damping_factor = compute_damping_factor(mask, ncf.mask_spval, args.length_scale)

        # Apply smoothing on damping factor if iterations are specified
        if args.smoothing_iterations > 0:
            print(f"Applying {args.smoothing_iterations} smoothing iterations to variable {varname}")
            damping_factor = smooth_multidimensional(damping_factor, kernel, args.smoothing_iterations)
            # Reset factor to 1 on masked values
            damping_factor = resetmask(damping_factor,mask,ncf.mask_spval)

        if args.kmscaling is not None:
            damping_factor = damping_factor / gridsize   # should handle different shapes of damping_factor automatically (?)

        # Write damping_factor in file
        ncio.write_variable(damping_factor, varname)

    # Close output file
    ncio.close_variable_file()

