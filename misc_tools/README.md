The lu generator can be applied  to perturb a NEMO-SI3 restart file for sea ice. This notebook explains how to create a file that will contain the correct mask for each variable in the restart file. The lu generator takes this restart-mask file as an input as well as the actual restart file.

## 1. Rebuild the restart file (if needed)
You might need to use the "rebuild" tools from NEMO to rebuild your restart file over the entire domain.

Example how to do so (rebuilds restart decomposed on 38 subdomains):
````
rebuild_nemo -d 2 -x 492 -y 566 NANUK4_ICE_ABL-ABLEVP903_00007080_restart_ice  38
````


## 2. Read the restart file and sort out variables according to some pre-defined rules
The goal is to create an output text file that will contain four lists of the variables names, depending on their corresponding mask type. The first line has the variable names to which we need to apply the tmask, second line for fmask, third line for umask, fourth line for vmask. 

Usage:
```
# load python with xarray (on jean jay: module load climate_science)

cd path/to/restart/file/directory/and/mesh/mask/file

# sort out variables from the restart file according to what we can guess (this will have to be modified when we know for sure!)
./sortvarsrst.py NANUK4_ICE_ABL-ABLBBM_restart_ice_mask.nc varlistall.asc
```

Note that at this stage, the predefined rules to sort out the variables are not all correct because for some variables, we don't know yet to which mask type they correspond. THIS CODE WILL NEED TO BE CORRECTED before final use. 
The variables are sorted based on the following criteria:
* Variables with dimensions strictly smaller than 3 are added to `varlist_skip`.
* Variables ending with 't' are added to `varlist_tmask`.
* Variables ending with 'f' are added to `varlist_fmask`.
* Variables starting with 'u' or 'U' are added to `varlist_umask`.
* Variables starting with 'v' or 'V' are added to `varlist_vmask`.
* Variables starting with 'sx' or 'sy' are added to `varlist_moments`.
* All other variables are added to `varlist_remain`.

Only the lists corresponding to t,f,u,v masks are written in the text file. But all the lists are displayed on the screen for your information.

An example of what the output text file looks like is:
```
dmgt,sgm11t,sgm12t,sgm22t
dmgf,sgm11f,sgm12f,sgm22f
Uu_sub,Uv_sub,uVice,u_ice
Vu_sub,Vv_sub,vUice,v_i,v_ice,v_il,v_ip,v_s                                            
```
This file doesn't need to contain all the variables. If not in the list, the variables will be treated as before (i.e. non-zero values replaced by ones, and everything else set to zero).

## 3. Create a masked copy of the restart file
Create a copy of the restart file  where each variable will be replaced by its mask according to the text file as an argument that lists the variables for each mask type.

Usage:
```
# load python with xarray (on jean jay: module load climate_science)

cd path/to/restart/file/directory/and/mesh/mask/file

# load python with xarray
module load climate_science

# copy restart file and repalce variables by their mask with correct dimensions
./mkrstmask.py mesh_mask_NANUK4_L31_4.2.nc NANUK4_ICE_ABL-ABLBBM903_00007080_restart_ice.nc NANUK4_ICE_ABL-ABLBBM903_00007080_restart_ice_mask.nc varlistall.asc 
```

Note that the created files has added a  `_Fillvalue` attribute to each variable compared to the original restart files. It  should not interfere with the usage the lu generator makes of the mask file (?).
