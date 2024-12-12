_Last update:_ 2024-12-12

_Author:_ Stephanie Leroux, Datlas.

__Goal:__ my notes to explain how to install the lu-simulator on Jean Zay@IDRIS (and also the required ensdam package).

---

### 1. Install  and compile ensdam package (on Jean-Zay@IDRIS HPC): 

* Do it once for all:
```
# clone repository
git clone https://github.com/brankart/ensdam.git

# install modules to get at minimum fortran, cmake and python with cython
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
# note on Jean Zay, once fortran, cmake, etc are installed via their modules, no need to specify their path in the above script, the machine "knows" where to find them. Might not be the case on other machines.

# test with examples:
cd example/python
python example_random.py
```
* Note: Each time you will use the lu simulator (which needs some ensdam libraries) you'll need to lead the same python module as the one you have used to compile ensdam once.

### 2. Install the lu-simulator
```
git clone https://github.com/cmems-arcticbliss/lu-simulator.git
```
