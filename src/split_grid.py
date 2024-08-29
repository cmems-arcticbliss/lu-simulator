import argparse
import numpy as np
import netcdf_io as ncio

if __name__ == "__main__":
  # Parse command-line arguments
  parser = argparse.ArgumentParser(prog='split_grid',description='Split grid for parallelization')
  parser.add_argument('-N',    '--nproc',     type=int,   required=True, help='number of processors to use with MPI')
  parser.add_argument('-grid', '--grid_file', type=str,   required=True, help='name of grid file')
  args = parser.parse_args()

  # read grid from grid file
  lon1dglo, lat1dglo = ncio.read_grid(args.grid_file)

  # MPI settings
  mpisize = args.nproc

  # Loop on processors (this code is not parallelized)
  for mpirank in range(mpisize):
    # Extract piece of grid to be used by processor mpirank
    lon1d = np.ascontiguousarray(lon1dglo[mpirank::mpisize])
    lat1d = np.ascontiguousarray(lat1dglo[mpirank::mpisize])
    # Save it to a specific file
    out_grid_file = args.grid_file + f'{mpirank:0>4}'
    ncio.save_grid(out_grid_file,lon1d,lat1d)
