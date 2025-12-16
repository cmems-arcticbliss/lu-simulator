# lu-simulator : location uncertainty simulator

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17896830.svg)](https://doi.org/10.5281/zenodo.17896830)

A tool to simulate horizontal location uncertainty in 2D or 3D geophysical fields.

This is done by a random displacement perturbations applied to the geophysical field, then re-interpolated on the original grid.

The code is developed in the framework of the [Arctic-BLISS project](https://github.com/cmems-arcticbliss)
to study the sensitivity of sea-ice prediction systems to small perturbations
of the initial condition (predictability experiment). 

Further possible applications include:

- simulate initial uncertainty in geophysical ensemble forecasting systems,

- enhance the ensemble spread in ensemble data assimilation systems,

- augment the size of the ensemble (or of the ML catalogue)
  with additional possible states of the system.

As compared to the classic method, which produces Gaussian perturbations of the model variables
in a specific modal subspace (EOFs, singular vectors, bred vectors,...),
this method has the originality of producing non-Gaussian perturbations
of the geoohysical fields through a Gaussian perturbation of the coordinates.
This can be especially useful in applications where the system
is very sensitive to small perturbations in the location of structures.

See the example directory for an illustration.

### Credits:

Developed by Jean-Michel Brankart, IR CNRS @IGE in Grenoble, France in collaboration with S. Leroux, Datlas, Grenoble.

### Requirements

This package makes use of:

- the [pyensdam library](https://github.com/brankart/ensdam/),

- the scipy interpolation tool: scipy.interpolate.

### Basic usage

This tool is primarily intended to ba applied on NetCDF files (as inputs and outputs).
To apply it on data in memory, just import and apply the functions embedded in the python code.

In this basic usage section (see below for generalizations), it is assumed that:

- the x and y coordinates are the grid indices,

- all grid points are valid (there is no land mask),

- the NetCDF files follow the standard of the [NEMO ocean model](https://www.nemo-ocean.eu/),

- everything is performed on one single processor.

##### Sample perturbations of the coordinates

```
sample_perturbations.py
  -m sample_size
  -l correlation_length_scale
  -grid grid_file.nc
  -o perturbation_file.nc
```

In this simple case, the correlation length scale is given in grid points,
and only the size of the grid is used from the grid file.

To seed the random number generator, use the optional argument:

```
  -seed seed_index
```

where `seed_index` is the index of a (pseudo-)randomly generated seed.
For a given index (default=0), the same seed will be used, which makes
the generation of the perturbations fully reproducible.

##### Apply perturbations to the geophysical fields

```
apply_perturbations.py
  -m sample_size
  -i perturbation_file.nc
  -std perturbation_std
  -grid grid_file.nc
  -ref reference_fields.nc
  -v listvar.txt
  -o perturbed_fields.nc
```

The reference file contains the geophysical fields that must be perturbed.
The standard deviation of the perturbation (assumed here homogeneous)
is given by `perturbation_std` (in grid points in this simplified case).
The list of variables to perturb is given by a text file with the following format:

```
var1 scalar T
var2 scalar F
...
var4 vector U
var5 vector V
...
var7 tensor T
var8 tensor T
var9 tensor T
...
varn scalar F
```

In this file, `varn` are the names of the variables in the reference NetCDF file.
The second column gives the type of each variable (scalar, vector or tensor),
and the third column gives the grid type (T, U, V or F).
Vectors must come by groups of 2 successive variables,
and tensors by groups of 3 successive variables.
For vectors, we must provide the components x and y, and
for tensors, we must provide the components xx, yy, and xy
(in that order, assuming symmetry of the tensor).

### About the perturbation

The perturbed of a variable at grid point $(x_i,y_i)$
is obtained by interpolating the reference field at location $(x_i+\delta x_i,y_i+\delta y_i)$,
where $(\delta x_i,\delta y_i)$ is a random perturbation of the location.
All variables at a given grid point are perturbed using the same $(\delta x_i,\delta y_i)$,
so that they all come from the same location in the reference field.
For staggered grid, the perturbation $(\delta x_i,\delta y_i)$ is interpolated
from grid T to the required grid.

At a given location, the perturbation $(\delta x_i,\delta y_i)$
is an isotropic normal random vector with a user-defined standard deviation (`perturbation_std`).
It is horizontally correlated according to a user-defined correlation length scale (`correlation_length_scale`).
By default, the horizontal correlation structure is Gaussian,
but other options are possible (see below).

Since $(\delta x,\delta y)$ is not constant, the perturbation of the scalar fields
can locally be viewed as a combination of a translation, a rotation and a deformation.
Because of the deformation, the local structure of the reference field can be modified,
as for instance the magnitude of the gradient, or the orthogonality bewteen isolines.
The typical magnitude of the local rotation and deformation is the ratio
between the standard deviation and the correlation length scale
($\rho=$ `perturbation_std` / `correlation_length_scale`).
Choose it much smaller than 1 if you want to keep the transformation regular.
This code is primarily intended for use with scalar fields and with small value of $\rho$.

Regarding vectors and tensors, we start by applying the perturbation to each component as explained above.
But this is usually not enough to maintain a sufficient local coherence of the state of the system.
As explained above, everything cannot be preserved, and the appropriate method
to adjust vectors and tensors depends on the application.
We have included only a few in the code:

1. Transform vectors and tensors with the inverse Jacobian of the change of coordinates.
In this case, vectors and tensors are transformed to follow the change of coordinates,
and their invariants are modified, which may not always be what is expected
to simulate uncertainty in their location.

2. Transform vectors and tensors with the inverse Jacobian and renormalize
to restore their invariants as they were in the reference state.
This is the method used by default, unless otherwise specified.

3. View the inverse Jacobian as a combination of a rotation and a deformation
($J^{-1}= R D$) and apply only the rotation to vectors and tensors.

The first two solutions requires that the transformation is regular, but not necessarily the third one.
As already said, we see that we cannot conserve everything, and that
location uncertainty inevitably involves uncertainty
in the relative orientation of the structures.

### Use world coordinates

To compute the perturbations with world coordinates rather than grid indices,
use the options:

`  -c`, for cartesian coordinates, or  
`  -s`, for spherical coordinates.

With these options, the location of the grid points are read from the grid file.
The unit of length scales to provide (`correlation_length_scale` and `perturbation_std`)
is then the same unit as the coordinates in the grid file.

In all cases, the spectrum of the perturbations can be modified by the option:

```
  -spct spectrum_type
```

where `spectrum_type` can be either 'gaussian', 'k-2' or 'k-4'.
It is easy to add a new option in the code.

### Account for a land mask

- On the one hand, before interpolating in the reference field, it is necessary
to extrapolate it over the mask to avoid interpolating using missing values.
This can be done with the tool:

```
unmask_fields.py
  -ref reference_fields.nc
  -v listvar.txt
  -mask mask.nc
  -o unmasked_reference_fields.nc
```

The unmasked reference field must then be used instead of the reference field
in the application of the perturbations.

- On the other hand, to avoid problems with extrapolated values, one may want
to reduce the amplitude of the perturbation to zero along the land mask.
The reduction factor can be generated from the mask file with the tool:

```
generate_damping_factor.py
  -mask mask.nc
  -v listvar.txt
  -o damping_factor.nc
  -l damping_length_scale
```

where `damping_length_scale` is the typical size of the buffer zone along the mask.
The resulting damping factor approximately follows the formula:
$\mbox{factor} = 1 - \exp( -x^2 / l^2 )$, where the distances x and l are expressed
in number of grid points from the coast line.

Then the damping factor must be used as an additional (optional) argument
in the application of the perturbations
(with `apply_perturbations.py`):

```
  -factor damping_factor.nc
```

More generally, this option can also be used to modulate spatially
the amplitude of the perturbations.

In addition, to re-apply the mask on the perturbed field, we can use
another optional argument in the application of the perturbations
(with `apply_perturbations.py`):

```
  -mask mask.nc
```

### Modify the format of the NetCDF files

To change the default convention for the NetCDF files (which follows the NEMO standard),
just modify the settings in the file `netcdf_format.py`.
It is possible to modify the name of dimensions, dimension variables,
grid variables and mask variables.

### MPI parallelization

For large size problems, the cost of the sampling of the perturbations
with `sample_perturbations.py` can be substantial
(especially in spherical coordinates and for a short correlation length scale),
and it can be useful to run it on more than one processor. This requires:

- splitting the grid file in pieces (one per processor):

```
split_grid.py
  -N nproc
  -grid grid.nc
```

- running `sample_perturbations.py` in parallel, with the additional option:

```
  -N nproc
```

- recombining the output perturbations, with:

```
recombine_perturbations.py
  -N nproc
  -i perturbation_file.nc
  -grid grid.nc
  -m sample_size
```

The gain in clock time should be roughly proportional to the number of processors used.

On the other hand, `apply_perturbations.py` can also be parallelized
(with the same optional parameter) by running different members on different processors.
In this case, be sure that the number of processors divide the number of members
so that each processor computes the same number of members.
