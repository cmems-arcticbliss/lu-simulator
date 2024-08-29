
# Initialize a list to hold the variables
var_list = []

# Define the attribute names
attribute_names = ['name', 'tensor_type', 'grid_type']

def read_vars(varlist_file):

  # Read the text file
  with open(varlist_file, 'r') as file:
    for line in file:
      # Split the line by spaces
      attributes = line.strip().split()
      # Create a dictionary for each variable
      variable = dict(zip(attribute_names, attributes))
      var_list.append(variable)


if __name__ == "__main__":
  import argparse

  # Parse command-line arguments
  parser = argparse.ArgumentParser(prog='vars_def',description='Read variable file')
  parser.add_argument('-v', '--varlist_file', type=str, required=True, help='file with list of variables')
  args = parser.parse_args()

  # Read variable file
  read_vars(args.varlist_file)

  # Print the variables
  for variable in var_list:
    print(variable) 

  # Print name of all variables
  for ivar, variable in enumerate(var_list):
    print(variable['name']) 

  # Print name of first variable
  print(var_list[0]['name']) 
