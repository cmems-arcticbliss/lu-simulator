#!/usr/bin/env python

import xarray as xr
import argparse

def create_rst_mask(mask_file, data_file, output_file, lists_file):
    """
    Copy and replaces values in a NEMO restart file (NetCDF) with the  masks corresponding to each variable. The masks are read from another NetCDF file.

    The function applies masks to the variables in the data file based on predefined lists read from a text file that can be prepared by the `sortvarsrst.py` script.
    For those variables not in the list (because we don't know its mask type yet, and if dimension larger than 2, the values are replaced by ones where the variable is not zero and by zeroes elsewhere.

    Parameters
    ----------
    mask_file : str
        Path to the NetCDF file containing the mask variables (tmask, fmask, umask, vmask).
    data_file : str
        Path to the NEMO restart file (NetCDF) to be copied and modified with mask values.
    output_file : str
        Path to the output NetCDF file where the modified restart data will be saved.
    lists_file : str
        Path to the text file containing variable lists for t, f, u, v masks.

    Notes
    -----
    This function expects the mask file to contain the variables `tmask`, `fmask`, `umask`, and `vmask` 
    and the lists file to contain the variables categorized into different masks in this order: tmask,fmask,umask,vmask.
    """
    # Load the mask dataset
    mask_ds = xr.open_dataset(mask_file)

    # Load the data dataset
    data_ds = xr.open_dataset(data_file)

    # Create a new dataset to store modifications
    #data_ds_mod = data_ds.copy(deep=True)

    # Read the variable lists from the text file
    var_lists = {}
    with open(lists_file, 'r') as f:
        var_lists['tmask'] = f.readline().strip().split(',')
        var_lists['fmask'] = f.readline().strip().split(',')
        var_lists['umask'] = f.readline().strip().split(',')
        var_lists['vmask'] = f.readline().strip().split(',')

    # Loop through all variables in the data dataset
    for var_name in data_ds.data_vars:
        # Access the variable data from the dataset
        var_data = data_ds[var_name]

        # Determine which mask to use based on which list the variable belongs to
        if var_name in var_lists['tmask']:
            mask = mask_ds['tmask'][0, 0, :, :].squeeze()
            mask_label = "t"
        elif var_name in var_lists['fmask']:
            mask = mask_ds['fmask'][0, 0, :, :].squeeze()
            mask_label = "f"
        elif var_name in var_lists['umask']:
            mask = mask_ds['umask'][0, 0, :, :].squeeze()
            mask_label = "u"
        elif var_name in var_lists['vmask']:
            mask = mask_ds['vmask'][0, 0, :, :]
            mask_label = "v"
        else:
            if len(var_data.dims) > 2:
                print(f"!!!!!!!!! Tricky variable... {var_name}: we don't know yet its real mask type. As a temporary trick its values will  be replaced by ones where it's initially non-zero and zeroes elsewhere).")
                # Update the dataset with modified data
                var_data = xr.where(var_data != 0, 1., 0.)
                data_ds[var_name] = var_data
                continue
            else:
                print(f"*** Skipping this variable... {var_name} as it does not need to be modified in the file")
                continue
                

        # Remove coordinate metadata if needed (otherwise it appears in final output file)
        mask = mask.drop_vars(['nav_lev', 'time_counter'], errors='ignore')

        print(f"==== Dealing with... {var_name} and {mask_label}mask.")
        
        # Broadcast mask to the shape of the current variable
        mask_broadcasted = xr.broadcast(mask, var_data)[0]

        # Reorder mask_broadcasted dimensions to match the order of the dimensions of var_data
        mask_broadcasted = mask_broadcasted.transpose(*var_data.dims)

        # Replace values in the variable where mask == 1
        var_data = xr.where(mask_broadcasted == 1, 1., 0.)
        
        # Update the dataset with modified data
        data_ds[var_name] = var_data
        # Remove specific attributes after assignment
        
        if '_FillValue' in data_ds[var_name].attrs:
            del data_ds[var_name].attrs['_FillValue']
        
    # Save the modified dataset to a new file
    data_ds.to_netcdf(output_file)
    print(f"Modified dataset saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Replace values in a NetCDF dataset based on various masks from another NetCDF file.")
    parser.add_argument('mask_file', type=str, help='Path to the NetCDF file containing the mask variables.')
    parser.add_argument('data_file', type=str, help='Path to the restart file (NetCDF) to copy and replace by masks corresponding to all variables.')
    parser.add_argument('output_file', type=str, help='Path to the output NetCDF file where the mask data will be saved.')
    parser.add_argument('lists_file', type=str, help='Path to the text file containing variable lists.')

    args = parser.parse_args()
    create_rst_mask(args.mask_file, args.data_file, args.output_file, args.lists_file)
