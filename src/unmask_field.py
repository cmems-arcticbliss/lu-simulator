import argparse
import numpy as np

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(prog='unmask_field',description='Unmask reference field')
    parser.add_argument('-ref',  '--reference_file', type=str,   required=True,  help='name of input file with reference fields')
    parser.add_argument('-v',    '--variable_list',  type=str,   required=True,  help='file with list of varoables to perturb')
    parser.add_argument('-mask', '--mask_file',      type=str,   required=True,  help='name of mask file')
    parser.add_argument('-o',    '--output_file',    type=str,   required=True,  help='name of output file with unmasked reference field')
    args = parser.parse_args()


