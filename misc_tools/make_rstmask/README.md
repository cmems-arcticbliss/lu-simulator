The lu generator can be applied  to perturb a NEMO-SI3 restart file for sea ice. This notebook explains how to create a file that will contain the correct mask for each variable in the restart file. The lu generator takes this restart-mask file as an input as well as the actual restart file.

## 1. Rebuild the restart file (if needed)
You might need to use the "rebuild" tools from NEMO to rebuild your restart file over the entire domain.

Example how to do so (rebuilds restart decomposed on 38 subdomains):
````
rebuild_nemo -d 2 -x 492 -y 566 NANUK4_ICE_ABL-ABLEVP903_00007080_restart_ice  38
````


## 2. Read the restart file and sort out variables according to some pre-defined rules
The goal is to create two output text files:
* one that will contain four lists of the variables names, depending on their corresponding mask type. The first line has the variable names to which we need to apply the tmask, second line for fmask, third line for umask, fourth line for vmask. 

Usage:
```
# load python with xarray (on jean jay: module load climate_science)

cd path/to/restart/file/directory/and/mesh/mask/file

# sort out variables from the restart file according to what we can guess (this will have to be modified when we know for sure!)
./sortvarsrst.py NANUK4_ICE_ABL-ABLBBM903_00007080_restart_ice.nc varlistall_LB.asc
```

In this example, the predefined rules to sort out the variables are based on LB's indications. The variables are sorted based on the following criteria:
* Variables with dimensions strictly smaller than 3 are added to `varlist_skip`.
* Variables ending with 't' are added to `varlist_tmask`.
* Variables ending with 'f' are added to `varlist_fmask`.
* Variable names 'uVice', 'Uv_sub', 'Vv_sub', and 'v_ice' are added to `varlist_vmask`.
* Variable names 'vUice', 'Uu_sub', 'Vu_sub', and 'u_ice' are added to `varlist_umask`.
* Variables starting with 'sx' or 'sy' and not finishing by 't' or 'f' are added to `varlist_tmask`
* All other variables are added to the `varlist_tmask`.

The lists corresponding to t,f,u,v masks are written in the text file. An example of how the output text file looks like is:
```
a_i,a_ip,dmgt,e_i_l01,e_i_l02,e_s_l01,e_s_l02,oa_i,sgm11t,sgm12t,sgm22t,snwice_mass,snwice_mass_b,sv_i,sx1mdt,sxa,sxage,sxap,sxc0_l01,sxc0_l02,sxdd1t,sxdd2t,sxdd3t,sxe_l01,sxe_l02,sxice,sxsal,sxsn,sxvl,sxvp,sxx1mdt,sxxa,sxxage,sxxap,sxxc0_l01,sxxc0_l02,sxxdd1t,sxxdd2t,sxxdd3t,sxxe_l01,sxxe_l02,sxxice,sxxsal,sxxsn,sxxvl,sxxvp,sxy1mdt,sxya,sxyage,sxyap,sxyc0_l01,sxyc0_l02,sxydd1t,sxydd2t,sxydd3t,sxye_l01,sxye_l02,sxyice,sxysal,sxysn,sxyvl,sxyvp,sy1mdt,sya,syage,syap,syc0_l01,syc0_l02,sydd1t,sydd2t,sydd3t,sye_l01,sye_l02,syice,sysal,sysn,syvl,syvp,syy1mdt,syya,syyage,syyap,syyc0_l01,syyc0_l02,syydd1t,syydd2t,syydd3t,syye_l01,syye_l02,syyice,syysal,syysn,syyvl,syyvp,t_su,v_i,v_il,v_ip,v_s
dmgf,sgm11f,sgm12f,sgm22f,sx1mdf,sxdd1f,sxdd2f,sxdd3f,sxx1mdf,sxxdd1f,sxxdd2f,sxxdd3f,sxy1mdf,sxydd1f,sxydd2f,sxydd3f,sy1mdf,sydd1f,sydd2f,sydd3f,syy1mdf,syydd1f,syydd2f,syydd3f
Uu_sub,Vu_sub,u_ice,vUice
Uv_sub,Vv_sub,uVice,v_ice                                           
```
This text file is shared in the current directory.

## 3. Create a masked copy of the restart file
Create a copy of the restart file  where each variable will be replaced by its mask according to the text file as an argument that lists the variables for each mask type.

Usage:
```
# load python with xarray (on jean jay: module load climate_science)

cd path/to/restart/file/directory/and/mesh/mask/file

# load python with xarray
module load climate_science

# copy restart file and repalce variables by their mask with correct dimensions
./mkrstmask.py mesh_mask_NANUK4_L31_4.2.nc NANUK4_ICE_ABL-ABLBBM903_00007080_restart_ice.nc NANUK4_ICE_ABL-ABLBBM903_00007080_restart_ice_mask.nc varlistall_LB.asc 
```

Note that the created files has added a  `_Fillvalue` attribute to each variable compared to the original restart files. It  should not interfere with the usage the lu generator makes of the mask file (?).
