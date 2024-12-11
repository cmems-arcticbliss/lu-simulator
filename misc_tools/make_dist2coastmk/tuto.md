_Last update:_ 2024-12-11

__Goal:__ my notes to explain how to install the lu-simulator on Jean Zay@IDRIS (and also the required ensdam package) and how to use the script `/src/generate_damping_factor.py`
 to define a mask at xxx distance (in grid points) from the coast. Note that this is not the primary goal of this script (initial goal is to compute a damping coefficient from 0 over land to 1 over ocean-sea-ice with a decreasing coefficient from 0 to 1 until xx distance from thecoast (in grid points).

### 1. Installer ensdam:
```
git clone https://github.com/brankart/ensdam.git

module purge
module load gcc/9.1.0
module load cmake
module load intel-all
module load gcc/9.1.0
module load hdf5/1.10.5-mpi
module load netcdf/4.7.2-mpi
module load netcdf-fortran/4.5.2-mpi
module unload python
module load climate_science

chmod u+x compile.bash
./compile.bash

# test with examples:
cd example/python
python example_random.py

```

### 2. Installer lu-simulator
```
git clone https://github.com/cmems-arcticbliss/lu-simulator.git

## NB: si nouveau terminal, toujours reloader module python : load climate_science
```

### 3. Utiliser le script ci dessous:

* Go to subdirectory: `lu-simulator/misc_tools/make_dist2coastmk`
```
cd lu-simulator/misc_tools/make_dist2coastmk

vi mkmask.sh
```
  
```bash
#!/bin/bash
# input parameters
lgrid=15 # in grid points. Distance to the coast at which damping decreases 
sm=3 # number of iterations of the 1-2-1 smoother
vargridinfo="vargridinfonew.asc"
mask="mask.nc"

# output file name
prefix_ofile="damping_factor"
suffix='pourLaurine'

python ../src/generate_damping_factor.py -l ${lgrid} -v ${vargridinfo} -mask ${mask} -s ${sm}  -o ${prefix_ofile}_lgrid${lgrid}_sm${sm}_${suffix}.nc
        
```

* `lgrid` is the number of grid point to set according to what you want.

* to run: `./mkmask.sh`
