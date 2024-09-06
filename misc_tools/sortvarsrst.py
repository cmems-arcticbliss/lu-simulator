#!/usr/bin/env python

import xarray as xr
from pprint import pprint
import argparse

def sort_variables_rst(file_path, output_file_path, include_all_lists=False):
    """
    Sorts variables from a NEMO restart file in NetCDF format into different lists based on their names and writes these lists to a text file.

    The variables are sorted based on the following criteria:
    - Variables with only one dimension or dimensions (y, x) are added to `varlist_skip`.
    - Variables ending with 't' are added to `varlist_tmask`.
    - Variables ending with 'f' are added to `varlist_fmask`.
    - Variables starting with 'u' or 'U' are added to `varlist_umask`.
    - Variables starting with 'v' or 'V' are added to `varlist_vmask`.
    - Variables starting with 'sx' or 'sy' are added to `varlist_moments`.
    - All other variables are added to `varlist_remain`.

    A test is applied to check if all variables have been sorted out in a list.

    Args:
        file_path (str): Path to the input NetCDF file.
        output_file_path (str): Path to the output text file where the sorted variable lists will be written.
        include_all_lists (bool): If True, includes the 'moments', 'skip', and 'remain' lists in the output file. Default is False.
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

    # Iterate over each variable in the dataset
    for var_name in ds.data_vars:
        var_data = ds[var_name]

        # Check if the variable has only one dimension or dimensions (y, x)
        if len(var_data.dims) < 2 or (len(var_data.dims) == 2 and 'y' in var_data.dims and 'x' in var_data.dims):
            varlist_skip.append(var_name)

        # Check if the variable starts with 'sx' or 'sy'
        elif var_name.startswith('sx') or var_name.startswith('sy'):
            varlist_moments.append(var_name)
        
        # Check if the variable ends with 't'
        elif var_name.endswith('t'):
            varlist_tmask.append(var_name)
        
        # Check if the variable ends with 'f'
        elif var_name.endswith('f'):
            varlist_fmask.append(var_name)
        
        # Check if the variable starts with 'u' or 'U'
        elif var_name.startswith(('u', 'U')):
            varlist_umask.append(var_name)
        
        # Check if the variable starts with 'v' or 'V'
        elif var_name.startswith(('v', 'V')):
            varlist_vmask.append(var_name)
        
        # Otherwise, put it in the remain list
        else:
            varlist_remain.append(var_name)

    # Sort the lists for better readability
    varlist_tmask.sort()
    varlist_fmask.sort()
    varlist_umask.sort()
    varlist_vmask.sort()
    varlist_skip.sort()
    varlist_remain.sort()
    varlist_moments.sort()

    # Prepare the output file path
    with open(output_file_path, 'w') as f:
        f.write(','.join(varlist_tmask) + '\n')
        f.write(','.join(varlist_fmask) + '\n')
        f.write(','.join(varlist_umask) + '\n')
        f.write(','.join(varlist_vmask) + '\n')
        if include_all_lists:
            f.write(','.join(varlist_moments) + '\n')
            f.write(','.join(varlist_skip) + '\n')
            f.write(','.join(varlist_remain) + '\n')


    # Print the lists to the terminal
    print("\nVariables with only 1 dimension or (y, x) (varlist_skip):")
    pprint(varlist_skip)

    print("\nVariables ending with 't' (varlist_tmask):")
    pprint(varlist_tmask)

    print("\nVariables ending with 'f' (varlist_fmask):")
    pprint(varlist_fmask)

    print("\nVariables starting with 'u' (varlist_umask):")
    pprint(varlist_umask)

    print("\nVariables starting with 'v' (varlist_vmask):")
    pprint(varlist_vmask)

    print("\nVariables starting with 'sx' or 'sy' (varlist_moments):")
    pprint(varlist_moments)

    print("\nVariables not sorted into any list (varlist_remain):")
    pprint(varlist_remain)

    # Check if the total number of variables is consistent
    total_vars = len(ds.data_vars)
    sum_vars = (len(varlist_tmask) + len(varlist_fmask) + len(varlist_umask) +
                len(varlist_vmask) + len(varlist_skip) + len(varlist_remain) +
                len(varlist_moments))

    print(f"Total number of variables: {total_vars}")
    print(f"Sum of variables in lists: {sum_vars}")

    if total_vars == sum_vars:
        print("The total number of variables matches the sum of variables in each list.")
    else:
        print("Warning: The total number of variables does not match the sum of variables in each list.")

    print(f"Lists have been written to the text file: {output_file_path} in the current directory.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sort variables in a NetCDF file into different lists and write to a text file.")
    parser.add_argument('input_file', type=str, help='Path to the input NetCDF file.')
    parser.add_argument('output_file', type=str, help='Path to the output text file.')
    parser.add_argument('-alllists', action='store_true', help='Include the "moments", "skip", and "remain" lists in the output file.')
    args = parser.parse_args()

    sort_variables_rst(args.input_file, args.output_file, include_all_lists=args.alllists)
