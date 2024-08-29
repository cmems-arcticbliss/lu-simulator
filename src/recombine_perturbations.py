import argparse
import numpy as np
import netcdf_io as ncio

if __name__ == "__main__":
  # Parse command-line arguments
  parser = argparse.ArgumentParser(prog='recombine_perturbations',description='Recombine perturbations after parallel computation')
  parser.add_argument('-N',    '--nproc',             type=int, required=True, help='number of processors to use with MPI')
  parser.add_argument('-i',    '--perturbation_file', type=str, required=True, help='name of perturbation file')
  parser.add_argument('-grid', '--grid_file',         type=str, required=True, help='name of grid file')
  parser.add_argument('-m',    '--sample_size',       type=int, required=True, help='sample size')
  args = parser.parse_args()

  # MPI settings
  mpisize = args.nproc

  # Read grid from grid file
  lon1d, lat1d = read_grid(args.grid_file)

  # Create output file
  ncio.create_perturbation_file(args.perturbation_file,args.sample_size)

  # Loop on ensemble members
  for imember in args.sample_size:
    # Initialize global member
    dxglo = np.zeros(ncio.grid_shape,dtype=np.double)
    dyglo = np.zeros(ncio.grid_shape,dtype=np.double)
    x_ind, y_ind = np.where(~np.isnan(dxglo))

    # Reconstruct ensemble member by recombining contribution from processors
    for mpirank in range(mpisize):
      input_file = args.perturbation_file + f'{mpirank:0>4}.nc'
      dx, dy = ncio.read_perturbation(input_file,imember)

      dxglo[x_ind[mpirank::mpisize],y_ind[mpirank::mpisize]] = dx
      dyglo[x_ind[mpirank::mpisize],y_ind[mpirank::mpisize]] = dy

    # Save member in output file
    print('Save recombined member: ',imember)
    ncio.write_perturbation(imember,dxglo,dyglo)

  # Close output file
  ncio.close_perturbation_file()
