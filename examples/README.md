# lu-simulator : examples


Prepare the data for the example with:

```
python prepare_test_case.py
```


Sample perturbations of the coordinates with:

```
python sample_perturbations.py -m 2 -l 20. -grid test.nc -o test_perturbation.nc
```


Apply perturbations to the reference field (scalar and vector) with:

```
python apply_perturbations.py -m 2 -std 3. -grid test.nc -i test_perturbation.nc -ref test.nc -v test_varlist.txt -o test_perturbed.nc
```


Plot reference and perturbed field with

```
python plot_test_case.py -i input_file
```

This should provide figures comparable to those provided here as an illustration.
The scalar field is in color, with white isolines;
the vector field is shown with black arrows;
the tensor field is shown with blue ellipses.

Reference field:

![plot](test.nc.png)

Perturbed field, member 1:

![plot](001test_perturbed.nc.png)

Perturbed field, member 2:

![plot](002test_perturbed.nc.png)
