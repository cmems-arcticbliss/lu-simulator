"""
Sample random location perturbation with given spectrum

Available functions:
 - sample_perturbation_spherical: in sphericat coordinates

Module parameters:
 - lmin : minimum degree of the spherical harmonics used in the computations
 - lmax : maximum degree of the spherical harmonics used in the computations
 - lcut : degree of the spherical harmonics defining the length scale of the perturbation spectrum 

"""
import numpy as np
import pyensdam as edam

# Set default values of module attributes
lmin=0       # minimum degree of the spherical harmonics used in the computations
lmax=60      # maximum degree of the spherical harmonics used in the computations
lcut=12      # degree of the spherical harmonics defining the length scale of the perturbation spectrum 

nwnbr=100    # Number of wave numbers in the definition of the spectrum (in cartesian or grid coordinates)
nharm=100    # Number of harmonic functions to superpose (in cartesian or grid coordinates)

# Callback function defining the power spectrum in spherical coordinates
# of the perturbations (in spherical coordinates)
# as a function of degree l and order m of the spherical harmonics
def power_spectrum(l,m):
  global lmax, lcut

  # Power spectrum
  power = 1. / ( 1. + l*l/(lcut*lcut) )

  # Normalize the spectrum
  norm = 0.
  for ll in range(0,lmax+1):
    norm = norm + 1. / ( 1. + ll*ll/(lcut*lcut) )
  power = power / norm

  # Scale to account for the multiplicity of each degree
  power = power / ( 1. + 2. * l )

  return power

# Sample random perturbations at given grid locations in spherical coordinates
def sample_perturbation_spherical(lon1d,lat1d,power_spectrum):
  """
  Sample random perturbation with required spectrum at given grid locations

  Args:
  lon1d : [double array]: longitudes of the output random field (1D)
  lat1d : [double array]: latitudes of the output random field (1D)

  Returns:
  perturbation: perturbation at given grid locations
  """
  global lmin, lmax

  perturbation = edam.random.field2s_sample(lon1d,lat1d,power_spectrum,lmin,lmax)

  return perturbation

# Sample random perturbations at given grid locations in cartesian coordinates
def sample_perturbation_cartesian(lon1d,lat1d,length_scale):
  """
  Sample random perturbation with required spectrum at given grid locations

  Args:
  lon1d : [double array]: longitudes of the output random field (1D)
  lat1d : [double array]: latitudes of the output random field (1D)

  Returns:
  perturbation: perturbation at given grid locations
  """
  global nharm

  # Set Gaussian spectrum with given length scale
  dk = 4 / ( nwnbr * length_scale )
  spct_freq = np.arange(dk, dk + nwnbr * dk, dk, dtype=np.double)
  spct_power = np.exp( - spct_freq * spct_freq )

  # Initialize definition of the spectrum
  edam.random.field2d_init(spct_freq,spct_power)

  # Compute random perturbation
  perturbation = edam.random.field2d_sample(lon1d.reshape(-1,1),lat1d.reshape(-1,1),nharm)

  return perturbation[0]

# Sample random perturbations at given grid locations in grid coordinates
def sample_perturbation_grid(grid_shape,length_scale):
  """
  Sample random perturbation with required spectrum at given grid locations

  Args:
  grid_shape :   length scale of perturbation (in grid points)
  length_scale : length scale of perturbation (in grid points)

  Returns:
  perturbation: perturbation at given grid locations
  """
  global nharm

  # Set Gaussian spectrum with given length scale
  dk = 4 / ( nwnbr * length_scale )
  spct_freq = np.arange(dk, dk + nwnbr * dk, dk, dtype=np.double)
  spct_power = np.exp( - spct_freq * spct_freq )

  # Initialize definition of the spectrum
  edam.random.field2d_init(spct_freq,spct_power)

  # Set grid coordinates
  indices = np.indices(grid_shape)
  lon2d = indices[0].astype(np.float64)  # Row indices as double (float64)
  lat2d = indices[1].astype(np.float64)  # Column indices as double (float64)

  # Compute random perturbation
  perturbation = edam.random.field2d_sample(lon2d,lat2d,nharm)

  return perturbation.reshape(grid_shape)

# Apply the module to NetCDF files
if __name__ == "__main__":
  import argparse
  import netcdf_io as ncio

  # Parse command-line arguments
  parser = argparse.ArgumentParser(prog='sample_perturbations',description='Sample perturbations of the coordinates')
  parser.add_argument('-m',    '--sample_size',  type=int,   required=True,  help='sample size')
  parser.add_argument('-l',    '--length_scale', type=float, required=True,  help='correlation length scale')
  parser.add_argument('-grid', '--grid_file',    type=str,   required=True,  help='name of grid file')
  parser.add_argument('-o',    '--output_file',  type=str,   required=True,  help='name of output file with perturbations')
  parser.add_argument('-c',    '--cartesian',    action='store_true', required=False, help='use cartesian coordinates')
  parser.add_argument('-s',    '--spherical',    action='store_true', required=False, help='use spherical coordinates')
  parser.add_argument('-seed', '--seed_index',   type=int,   required=False, help='index of seed to use in the random number generator')
  parser.add_argument('-N',    '--nproc',        type=int,   required=False, help='number of processors to use with MPI')
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
    grid_file=args.grid_file + f'{mpirank:0>4}'
    output_file=args.output_file + f'{mpirank:0>4}'
  else:
    mpisize=1
    mpirank=0
    grid_file=args.grid_file
    output_file=args.output_file

  # Seed random number generator
  # (reproducible for a given seed index, default=0)
  if args.seed_index is not None:
    edam.random.seed(args.seed_index)
  else:
    edam.random.seed(0)

  # Read grid from grid file
  lon1d, lat1d = ncio.read_grid(args.grid_file)

  # Create output file
  ncio.create_perturbation_file(output_file,args.sample_size,multiple_files=use_mpi,vector_size=lon1d.size)

  # Loop on sample/ensemble size
  for imember in range(args.sample_size):
    if mpirank == 0:
      print('  Generating sample member:',imember)

    # Generate new random perturbation
    if args.spherical:
      lcut = 1 + int( 360. / ( 2 * np.pi * args.length_scale ) ) ; lmax = 5 * lcut
      dx = sample_perturbation_spherical(lon1d,lat1d,power_spectrum)
      dy = sample_perturbation_spherical(lon1d,lat1d,power_spectrum)
    elif args.cartesian:
      dx = sample_perturbation_cartesian(lon1d,lat1d,args.length_scale)
      dy = sample_perturbation_cartesian(lon1d,lat1d,args.length_scale)
    else:
      dx = sample_perturbation_grid(ncio.grid_shape,args.length_scale)
      dy = sample_perturbation_grid(ncio.grid_shape,args.length_scale)

    # Save perturbation in output file
    if mpirank == 0:
      print('  Perturbation stored in file: ',output_file)
    ncio.write_perturbation(imember,dx,dy,multiple_files=use_mpi)

  # Close output file
  ncio.close_perturbation_file()
