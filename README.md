<p align="center">
  <img src="./images/logo.png" alt="ACCeL Logo" width="200px">
</p>

# ACCeL

**ACCeL** (Automated Conformer Calculations and Exploration Library) is a Python package designed to streamline the processing of conformational ensembles in computational chemistry. It provides an intuitive and flexible framework for managing large sets of molecular structures, automating workflows such as input generation, energy analysis, structure filtering, and file export. ACCeL supports various file formats and is well suited for processing quantum chemical outputs in a reproducible and scalable way.

## 🔧 Installation

Install ACCeL using pip:

```bash
pip install accel
```

## 🚀 Quick Start

```python
from accel import Box

Box("*.log")\
    .read_atoms()\
    .read_energy()\
    .energy_limit()\
    .rmsd_limit()\
    .write_input("template_file.inp")
```

This script performs the following:

- Reads atomic coordinates and energy values from `.log` files.
- Filters conformers based on energy and RMSD thresholds.
- Generates input files using a specified template.

## 📦 The `Box` Class

The `Box` class serves as the core component of ACCeL, managing collections of molecular structures, typically conformers and their variants. Key functionalities include:

- **File Management**: Add files or directories, read atomic coordinates, energies, and thermodynamic corrections.
- **Validation**: Check for successful completion of calculations and analyze vibrational frequencies.
- **Data Processing**: Calculate free energies, relative energies, and Boltzmann distributions.
- **Structure Handling**: Filter conformers by energy or RMSD thresholds, generate mirror images, and modify molecular geometries.
- **Input/Output**: Create input files from templates, export data to CSV, and read/write XYZ or MOL files.

## 🖼️ Visual Overview

<p align="center">
  <img src="./images/box_figure.png" alt="Box Class Overview" width="600px">
</p>

## 📊 Method Overview

| Method                                                                                                                  | Description                                                                                                                                                                                                                                                          |
| ----------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `add(contents)`                                                                                                         | Add files or directories into the Box                                                                                                                                                                                                                                |
| `read_atoms(filetype=None, **options)`                                                                                  | Read atomic coordinates and relevant information                                                                                                                                                                                                                     |
| `check_end(filetype=None, **options)`                                                                                   | Check whether output files finished successfully                                                                                                                                                                                                                     |
| `read_energy(filetype=None, **options)`                                                                                 | Extract energy values from output files                                                                                                                                                                                                                              |
| `read_correction(filetype=None, **options)`                                                                             | Read thermodynamic correction values                                                                                                                                                                                                                                 |
| `check_freq(filetype=None, **options)`                                                                                  | Analyze vibrational frequencies and imaginary modes                                                                                                                                                                                                                  |
| `run(filetype=None, **options)`                                                                                         | Run calculations and wait for completion                                                                                                                                                                                                                             |
| `submit(filetype=None, **options)`                                                                                      | Submit all calculations immediately without waiting                                                                                                                                                                                                                  |
| `calc_free_energy(filetype=None, **options)`                                                                            | Calculate Gibbs free energy by combining electronic energy and thermodynamic correction values                                                                                                                                                                       |
| `labeling(separator="_", index_list=[0, 1])`                                                                            | Automatically set labels based on filename patterns. For example, a file named `KEF20958_a_256.log` will be labeled as `KEF20958_a` if `separator="_"` and `index_list=[0, 1]`. This helps group conformers under common labels for energy comparison and filtering. |
| `set_data(key, value=None)`                                                                                             | Batch-set arbitrary metadata                                                                                                                                                                                                                                         |
| `set_state(flag=True)`                                                                                                  | Batch-enable or disable entries                                                                                                                                                                                                                                      |
| `set_label(label="")`                                                                                                   | Assign label names in bulk                                                                                                                                                                                                                                           |
| `zero_fill(digit=3, separator="_", position=3)`                                                                         | Zero-pad numeric suffixes in names. For example, `KEF20958_a_2` becomes `KEF20958_a_002` when `digit=3`, `separator="_"`, and `position=3`. This ensures proper alphanumeric sorting of filenames.                                                                   |
| `count(comment="")`                                                                                                     | Count the number of active entries by label                                                                                                                                                                                                                          |
| `export_data(filepath)`                                                                                                 | Output data to CSV file                                                                                                                                                                                                                                              |
| `energy_limit(threshold=3.0, max_limit=None, in_label=True)`                                                            | Filter conformers by relative energy threshold                                                                                                                                                                                                                       |
| `calc_rel_energy(in_label=True)`                                                                                        | Convert energy values to relative values within each label group                                                                                                                                                                                                     |
| `calc_distribution(in_label=True, temperature=298.15)`                                                                  | Calculate Boltzmann-weighted populations                                                                                                                                                                                                                             |
| `calc_energy(keys=[], unit=Units.kcal_mol)`                                                                             | Calculate total energy using specified keys                                                                                                                                                                                                                          |
| `copy_files(directory, change_path=False, suffix=None)`                                                                 | Copy files to specified directory                                                                                                                                                                                                                                    |
| `search(directory=None, existing_check=True, suffix=None)`                                                              | Search for files with matching names                                                                                                                                                                                                                                 |
| `read_xyz()`                                                                                                            | Read XYZ coordinate files                                                                                                                                                                                                                                            |
| `write_xyz(directory=None, change_path=True, centering=True)`                                                           | Export XYZ files                                                                                                                                                                                                                                                     |
| `read_mol()`                                                                                                            | Read MOL-format structure files                                                                                                                                                                                                                                      |
| `write_mol(directory=None, change_path=True, centering=True)`                                                           | Export MOL files                                                                                                                                                                                                                                                     |
| `write_input(template, directory=None, change_path=True)`                                                               | Generate input files using template variables                                                                                                                                                                                                                        |
| `calc_bonds(cov_scaling=1.1, vdw_scaling=1.0)`                                                                          | Calculate bonding information for symmetry calculations                                                                                                                                                                                                              |
| `calc_symm(calc_all=True)`                                                                                              | Calculate symmetry data required for RMSD                                                                                                                                                                                                                            |
| `rmsd_limit(threshold=0.01, all_combinations_of_confs=False, redundant_check=3, all_perturbation_of_rotamers=False)`    | Filter out redundant structures based on RMSD                                                                                                                                                                                                                        |
| `map_numbers(reference_box=None)`                                                                                       | Align atom numbering between conformers to ensure consistent indexing, particularly for symmetric atoms (e.g., hydrogens in methylene groups). Useful for accurate comparison of geometries and RMSD calculations.                                                   |
| `calc_length(number_a, number_b, key="")`                                                                               | Calculate interatomic distances                                                                                                                                                                                                                                      |
| `calc_dihedral(number_a, number_b, number_c, number_d, key="")`                                                         | Calculate dihedral angles                                                                                                                                                                                                                                            |
| `calc_angle(number_a, number_b, number_c, key="")`                                                                      | Calculate bond angles                                                                                                                                                                                                                                                |
| `modify_length(number_a, number_b, target, fix_a=False, fix_b=False, numbers_along_with_a=[], numbers_along_with_b=[])` | Modify bond lengths between atoms                                                                                                                                                                                                                                    |
| `convert_to_mirror(centering=True)`                                                                                     | Generate mirror images of structures                                                                                                                                                                                                                                 |
| `only_minimum(in_label=True)`                                                                                           | Keep only the lowest-energy conformer per label                                                                                                                                                                                                                      |
| `get()`                                                                                                                 | Return a `Systems` object representing active conformers. This object behaves similarly to a list but includes additional methods for structure handling and analysis.                                                                                               |
| `get_average(keys=[], keys_for_atoms=[])`                                                                               | Return a `Systems` object containing Boltzmann-averaged data per label                                                                                                                                                                                               |
| `duplicate()`                                                                                                           | Return a copy of the Box object                                                                                                                                                                                                                                      |

## 📝 Template Keywords

When generating input files, the following placeholders in your template will be replaced:

- `#NAME#`: File name without extension
- `#LABEL#`: Label name
- `#AXYZ#`: Atomic symbols with XYZ coordinates
- `#ATOMS#`: Number of atoms
- `#CHG#`: Molecular charge
- `#MULT#`: Spin multiplicity
- `#ENRGY#`: Energy value
- `#PATH#`: File path
- `#DATA[]#`: Custom data values (e.g., `#DATA[rmsd]#` to insert the RMSD value if previously set)

## 📄 License

This project is licensed under the MIT License.

## 🗂️ Project Information

- **Repository**: [github.com/kfchem/accel](https://github.com/kfchem/accel)
- **Author**: Keisuke Fukaya
- **License**: MIT
- **Python Compatibility**: Python 3.9+
- **Keywords**: conformer, quantum chemistry, automation, batch processing

## 📬 Contact

For bug reports or feature requests, please open an issue on the [GitHub repository](https://github.com/kfchem/accel/issues).

For other inquiries (e.g., academic or collaboration-related), feel free to contact the author via email: [kfukaya@pu-toyama.ac.jp](mailto:kfukaya@pu-toyama.ac.jp)

Thank you for your interest in ACCeL!
