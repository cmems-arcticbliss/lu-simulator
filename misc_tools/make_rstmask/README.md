The lu generator can be applied  to perturb a NEMO-SI3 restart file for sea ice. This notebook explains how to create a file that will contain the correct mask for each variable in the restart file. The lu generator takes this restart-mask file as an input as well as the actual restart file.

## 1. Rebuild the restart file (if needed)
You might need to use the "rebuild" tools from NEMO to rebuild your restart file over the entire domain.

Example how to do so (rebuilds restart decomposed on 38 subdomains):
````
rebuild_nemo -d 2 -x 492 -y 566 NANUK4_ICE_ABL-ABLEVP903_00007080_restart_ice  38
````


## 2. Read the restart file and sort out variables according to some pre-defined rules

Usage:
```
# load python with xarray (on jean jay: module load climate_science)

cd path/to/restart/file/directory/and/mesh/mask/file

# sort out variables from the restart file according to what we can guess (this will have to be modified when we know for sure!)
./sortvarsrst.py NANUK4_ICE_ABL-ABLBBM903_00007080_restart_ice.nc vargridinfo.asc
```

In this example, the predefined rules to sort out the variables are based on LB's indications. The variables are sorted based on the following criteria:
* Variables with dimensions strictly smaller than 3 are added to `varlist_skip`.
* Variables 'Uv_sub', 'Vv_sub’, 'Uu_sub', 'Vu_sub' are added to `varlist_skip` .
* Variables ending with 't' are added to `varlist_tmask`.
* Variables ending with 'f' are added to `varlist_fmask`.
* Variable names 'uVice', 'v_ice' are added to `varlist_vmask`.
* Variable names 'vUice', 'u_ice' are added to `varlist_umask`.
* Variables starting with 'sx' or 'sy' and not finishing by 't' or 'f' are added to `varlist_tmask`
* All other variables are added to the `varlist_tmask`.

An output text file is then written containing for each variable its type (scalar,vector,tensor) and its grid type (U,V,T,F) depending on the above criteria.
This text file is shared in the current directory.

## 3. Create a masked copy of the restart file
Create a copy of the restart file  where each variable will be replaced by its mask. The corresponding mask   is  read from a separate NetCDF file, and the corresponding mask type for each variable is determined based on information  read from a text file that can be prepared manually or using the `sortvarsrst.py` script.
If some variables from the restart file are not listed in the text file and if dimension larger than 2, the values are replaced by ones where the variable is not zero and by zeroes elsewhere (i.e. poor man's mask).
Usage:
```
# load python with xarray (on jean jay: module load climate_science)

cd path/to/restart/file/directory/and/mesh/mask/file

# load python with xarray
module load climate_science

# copy restart file and repalce variables by their mask with correct dimensions
./mkrstmask.py mesh_mask_NANUK4_L31_4.2.nc NANUK4_ICE_ABL-ABLBBM903_00007080_restart_ice.nc NANUK4_ICE_ABL-ABLBBM903_00007080_restart_ice_mask.nc vargridinfo.asc 
```

Note that the created files has added a  `_Fillvalue` attribute to each variable compared to the original restart files. It  should not interfere with the usage the lu generator makes of the mask file (?).

In the end, if you aim to use the text file as input for the lu simultor, you must check that each vector components appear as vx then vy on the next line. You must also check tensor components order: xx, yy, xy. If you have not applied the official mask for given variables (e.g. 'Uv_sub', 'Vv_sub’, 'Uu_sub', 'Vu_sub') but still want to perturbed them with the lu generator, you need to add them in the list.
