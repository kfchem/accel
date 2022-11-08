<p align="center">
  <img src="./images/logo.svg" alt="ACCeL" width="200px">
</p>

#
ACCeL is a python package that enables batch processing of conformational isomers, making it easy to create, run, and analyze files for computational chemistry programs and automate calculation processes.
## Installation
```
pip install accel
```
## Examples
```
from accel import Box
Box("*.log").read_atoms().read_energy().energy_limit().rmsd_limit().write_input("template_file.inp")
```
## How to use "Box"
![Box](./images/box_figure.png)

#### `add(contents)`
> Adds files to the Box. Contents can be specified flexibly as directory name, filename, wildcard, iterative object, etc.
#### `read_atoms(filetype=None, **options)`
> Read coordinate information and other information from files.
#### `check_end(filetype=None, **options)`
> Verify that output files are terminated successfully.
#### `read_energy(filetype=None, **options)`
> Read energy values from output files.
#### `read_correction(filetype=None, **options)`
> Read thermodynamic correction values from output files.
#### `check_freq(filetype=None, **options)`
> Read frequencis from files and check the number of imaginary frequencies.
#### `run(filetype=None, **options)`
> Execute files and wait for them to complete.
#### `submit(filetype=None, **options)`
> Execute all files immediately.
#### `calc_free_energy(filetype=None, **options)`
> Calculate the free energy using loaded values.
#### `labeling(separator="_", index_list=[0, 1])`
> Set the label name (isomer name other than the conformational isomer) automatically. By default, a conformation with the name `KEF20958_a_256` will have the label `KEF20958_a`. Therefore, it is recommended that the file name handled in ACCeL be BaseName_IsomerName_Number.
#### `set_data(key, value=None)`
> Set data in batches.
#### `set_state(flag=True)`
> Set state at once.
#### `set_label(label="")`
> Set label names in batches.
#### `zero_fill(digit=3, separator="_", position=3)`
> Fill in the numbers in the name with zeros. By default, `KEF20958_a_2` is set to `KEF20958_a_002`.
#### `count(comment="")`
> Count up the number of active systems with their respective labels.
#### `export_data(filepath)`
> Output information on the conformers to a CSV file.
#### `energy_limit(threshold=3.0, max_limit=None, in_label=True)`
> Enable only stable conformation (default 3.0 kcal/mol) per label.
#### `calc_rel_energy(in_label=True)`
> Convert energy values to relative values.
#### `calc_distribution(in_label=True, temperature=298.15)`
> The Boltzmann distribution is calculated according to energy values. By default, it is calculated for each label.
#### `calc_energy(keys=[], unit=Units.kcal_mol)`
> The data specified for `keys` are added together and the energy value is calculated.
#### `copy_files(directory, change_path=False, suffix=None)`
> Duplicate the file in `directory`.
#### `search(directory=None, existing_check=True, suffix=None)`
> Search for files with the same name in `directory`.
#### `read_xyz()`
> Read coordinate information, etc. from XYZ files (.xyz).
#### `write_xyz(directory=None, change_path=True, centering=True)`
> Export XYZ files.
#### `read_mol()`
> Read coordinate information, etc. from MOL files (.mol, .sdf, .sd).
#### `write_mol(directory=None, change_path=True, centering=True)`
> Export MOL files.
#### `write_input(template, directory=None, change_path=True)`
> Create an input file. The following keywords in the template will be replaced.
> #NAME# -> name
> #LABEL# -> label name
> #AXYZ# -> Element name + XYZ list
> #ATOMS# -> number of atoms
> #CHG# -> charge
> #MULT# -> multiplicity
> #ENRGY# -> energy value
> #PATH# -> file path
> #DATA[]# -> data value
#### `calc_bonds(cov_scaling=1.1, vdw_scaling=1.0)`
> Embed bonding information. Necessary to calculate the symmetric information.
#### `calc_symm(calc_all=True)`
> Embed symmetry information; required for RMSD calculations.
#### `rmsd_limit(threshold=0.01, all_combinations_of_confs=False, redundant_check=3, all_perturbation_of_rotamers=False)`
> Disable one conformation among the conformations whose RMSD is less than a certain value. In other words, it removes the identical conformations. The default value is 0.01 angstrom.
#### `map_numbers(reference_box=None)`
> Correct misnumbering of the same substituent (e.g., two H's in methylene) between conformations.
#### `calc_length(number_a, number_b, key="")`
> Calculate the distance between atoms.
#### `calc_dihedral(number_a, number_b, number_c, number_d, key="")`
> Calculate the dihedral angle between atoms.
#### `calc_angle(number_a, number_b, number_c, key="")`
> Calculate the angle between atoms.
#### `modify_length(number_a, number_b, target, fix_a=False, fix_b=False, numbers_along_with_a=[], numbers_along_with_b=[])`
> Correct the distance between atoms.
#### `convert_to_mirror(centering=True)`
> Convert to a mirror image.
#### `only_minimum(in_label=True)`
> Enable only the most stable conformation. The default is to calculate with the respective label.
#### `get() -> Systems`
> Returns the Systems of the active coordination group.
#### `get_average(keys=[], keys_for_atoms=[]) -> Systems`
> Return Systems for the locus group with energy-weighted average data for each label.
#### `duplicate() -> Box`
> Returns a duplicate Box.

