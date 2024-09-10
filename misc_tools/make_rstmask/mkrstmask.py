#!/usr/bin/env python

import xarray as xr
import argparse

def read_and_categorize(file_path):
    """
    Reads a text file containing variable names, their type (scalar,vector,tensor) and grid type (T,F,U,V)  and categorizes them into lists based on the grid type. 

    Args:
        file_path (str): Path to the input text file.

    Returns:
        dict: A dictionary with grid types as keys and lists of variable names as values.
    """
    # Initialize the dictionary with empty lists for each grid type
    var_lists = {
        'tmask': [],
        'fmask': [],
        'umask': [],
        'vmask': []
    }

    # Read the file and categorize each variable
    with open(file_path, 'r') as file:
        for line in file:
            # Strip any leading/trailing whitespace and split by space
            line = line.strip()
            if not line:
                continue
            
            parts = line.split()
            if len(parts) != 3:
                print(f"Skipping malformed line: {line}")
                continue
            
            var_name, var_type, grid = parts
            
            # Append the variable to the appropriate list based on the grid type
            if grid == 'T':
                var_lists['tmask'].append(var_name)
            elif grid == 'F':
                var_lists['fmask'].append(var_name)
            elif grid == 'U':
                var_lists['umask'].append(var_name)
            elif grid == 'V':
                var_lists['vmask'].append(var_name)
            else:
                print(f"Unknown grid type '{grid}' for variable '{var_name}'")

    return var_lists


def create_rst_mask(mask_file, data_file, output_file, lists_file):
    """
    Copy and replaces values in a NEMO restart file (NetCDF) with the  mask corresponding to each variable. The corresponding mask   is  read from a separate NetCDF file, and the corresponding mask type for each variable is determined based on information  read from a text file that can be prepared manually or using the `sortvarsrst.py` script.
    If some variables from the restart file are not listed in the text file and if dimension larger than 2, the values are replaced by ones where the variable is not zero and by zeroes elsewhere.

    Parameters
    ----------
    mask_file : str
        Path to the NetCDF file containing the mask variables (tmask, fmask, umask, vmask).
    data_file : str
        Path to the NEMO restart file (NetCDF) to be copied and modified with mask values.
    output_file : str
        Path to the text file containing the grid type info. Format: one line per variable: varname, variable type (scalar,vector,tensor), variable grid type (T,F,U,V) 

    Notes
    -----
    This function expects the mask file to contain the variables `tmask`, `fmask`, `umask`, and `vmask`.
    """
    # Load the mask dataset
    mask_ds = xr.open_dataset(mask_file)

    # Load the data dataset
    data_ds = xr.open_dataset(data_file)

    # Read the variable lists from the text file
    var_lists = read_and_categorize(lists_file)

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
                print(f"!!!!!!!!! Tricky variable... {var_name}: This variable was not in the list.")
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
