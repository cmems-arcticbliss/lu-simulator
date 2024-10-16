#!/usr/bin/env python

import xarray as xr
from pprint import pprint
import argparse

def determine_type(var_name):
    """
    Determines whether the variable is a scalar, vector, or tensor in a NEMO restart file (BBM rheology).

    Args:
        var_name (str): The name of the variable.

    Returns:
        str: 'vector', 'tensor', or 'scalar'.
    """
    # Check for tensor
    if var_name.startswith('sxx') or var_name.startswith('sxy') or var_name.startswith('syy'):
        return 'tensor'
    # Check for tensor based on ending and specific characters
    if (var_name.endswith('t') or var_name.endswith('f')) and len(var_name) >= 3:
        if var_name[-3] in {'1', '2'} and var_name[-2] in {'1', '2'}:
            return 'tensor'

    # Check for vector
    # List of exceptions for variables starting with 'v'
    excluded_vector_vars = {'v_i', 'v_ip', 'v_il', 'v_s'}
    
    if var_name.lower().startswith('u'):
        return 'vector'
    elif var_name.lower().startswith('v') and var_name not in excluded_vector_vars:
        return 'vector'
    elif (var_name.startswith('sx') or var_name.startswith('sy')) and not (var_name[2:3] == 'x' or var_name[2:3] == 'y'):
        return 'vector'    

    # Otherwise, it's a scalar
    return 'scalar'


def sort_variables_rst(file_path, output_file_path, verbose=False):
    """
    Sorts variables from a NEMO restart file in NetCDF format into different lists based on their names and writes this information in an output text file.

    The variables are sorted based on the following criteria:
    - Variables with only one dimension or dimensions (y, x) are added to `varlist_skip`.
    - Variables ending with 't' are added to `varlist_tmask`.
    - Variables ending with 'f' are added to `varlist_fmask`.
    - Variables starting with 'sx' or 'sy' and not finishing by 't' or 'f' are added to `varlist_tmask` and to `varlist_moments`.
    - Variable names 'uVice', 'v_ice' are added to `varlist_vmask`.
    - Variable names 'vUice', 'u_ice' are added to `varlist_umask`.
    - All other variables are added to the `varlist_tmask` and `varlist_remain`.
    
    Based on these lists, a text file is written,  containaing variable name, type (scalar, vector, tensor), and grid type (T, U, V, F) for each non-skiped variable.

    Args:
        file_path (str): Path to the input NetCDF file.
        output_file_path (str): Path to the  output text file where is written the variable name, type, and grid type.
        
    """
    # Load the dataset
    ds = xr.open_dataset(file_path)

    # Initialize lists
    varlist_tmask = []
    varlist_fmask = []
    varlist_umask = []
    varlist_vmask = []
    varlist_skip = []
    varlist_remain = []
    varlist_moments = []

    # List of specific variable names
    specific_vmask = ['uVice',  'v_ice']
    specific_umask = ['vUice',  'u_ice']

    # Iterate over each variable in the dataset
    for var_name in ds.data_vars:
        var_data = ds[var_name]

        # Check if the variable has only one dimension or dimensions (y, x)
        if len(var_data.dims) < 2 or (len(var_data.dims) == 2 and 'y' in var_data.dims and 'x' in var_data.dims):
            varlist_skip.append(var_name)


        # Check if the variable ends with 't'
        elif var_name.endswith('t'):
            varlist_tmask.append(var_name)

        # Check if the variable ends with 'f'
        elif var_name.endswith('f'):
            varlist_fmask.append(var_name)

        # Check if the variable starts with 'sx' or 'sy' and does not end with 't' or 'f'
        elif (var_name.startswith('sx') or var_name.startswith('sy')) and not (var_name.endswith('t') or var_name.endswith('f')):
            varlist_tmask.append(var_name)
            varlist_moments.append(var_name)

        # Check if the variable is one of the specific vmask variables
        elif var_name in specific_vmask:
            varlist_vmask.append(var_name)

        # Check if the variable is one of the specific umask variables
        elif var_name in specific_umask:
            varlist_umask.append(var_name)

        # Otherwise, put it in the remain list and the tlist
        else:
            varlist_remain.append(var_name)
            varlist_tmask.append(var_name)

    # Sort the  u,v lists only
    varlist_umask.sort()
    varlist_vmask.sort()
   

    # Write to the  output file 
    with open(output_file_path, 'w') as f2:
        for var_name in varlist_tmask:
            if var_name not in varlist_skip:
                f2.write(f"{var_name} {determine_type(var_name)} T\n")
        for var_name in varlist_fmask:
            if var_name not in varlist_skip:
                f2.write(f"{var_name} {determine_type(var_name)} F\n")
        for var_name in varlist_umask:
            if var_name not in varlist_skip:
                f2.write(f"{var_name} {determine_type(var_name)} U\n")
        for var_name in varlist_vmask:
            if var_name not in varlist_skip:
                f2.write(f"{var_name} {determine_type(var_name)} V\n")

    if verbose:
        # Print the lists to the terminal
        print("\nVariables with only 1 dimension or (y, x) (varlist_skip):")
        pprint(varlist_skip)
    
        print("\nVariables ending with 'f' (varlist_fmask):")
        pprint(varlist_fmask)
    
        print("\nVariables in the umask list (varlist_umask):")
        pprint(varlist_umask)
    
        print("\nVariables in the vmask list (varlist_vmask):")
        pprint(varlist_vmask)
    
        print("\nVariables ending with 't' (varlist_tmask):")
        pprint(varlist_tmask)

    # Check if the total number of variables is consistent
    total_vars = len(ds.data_vars)
    sum_vars = (len(varlist_tmask) + len(varlist_fmask) + len(varlist_umask) +
                len(varlist_vmask) + len(varlist_skip))

    print(f"Total number of variables: {total_vars}")
    print(f"Sum of variables in lists: {sum_vars}")

    if total_vars == sum_vars:
        print("The total number of variables matches the sum of variables in each list.")
    else:
        print("Warning: The total number of variables does not match the sum of variables in each list.")

    print(f"Variable grid and type information has been written to this file {output_file_path}.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sort variables in a NetCDF file into different lists and write to text files.")
    parser.add_argument('input_file', type=str, help='Path to the input NetCDF file.')
    parser.add_argument('output_file', type=str, help='Path to the  output text file.')
    parser.add_argument('--verbose', action='store_true', help='Print all lists.')
    args = parser.parse_args()

    sort_variables_rst(args.input_file, args.output_file,  verbose=args.verbose)
