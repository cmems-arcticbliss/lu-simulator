#!/bin/bash

# Lu SIMULATOR
ludir="/lustre/fswork/projects/rech/cli/regi915/DEVGIT/ArcticBLISS-all/lu-simulator"

# input parameters
lgrid=15 # in grid points. Distance to the coast at which damping decreases 
sm=3 # number of iterations of the 1-2-1 smoother
vargridinfo="varlist.asc"
mask="reducedmask.nc"

# output file name
prefix_ofile="damping_factor"
suffix='pourLaurine'

python ${ludir}/src/generate_damping_factor.py -l ${lgrid} -v ${vargridinfo} -mask ${mask} -s ${sm}  -o ${prefix_ofile}_lgrid${lgrid}_sm${sm}_${suffix}.nc

