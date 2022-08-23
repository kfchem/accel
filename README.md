<p align="center">
  <img src="./images/logo.svg" alt="ACCeL" width="150px">
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
# How to use ACCeL
## figure of Box
![Box](./images/box_figure.png)

## methods of Box
#### `add(contents)`
> Adds files to the Box. Contents can be specified flexibly as directory name, filename, wildcard, iterative object, etc.
#### `read_atoms(filetype=None, **options)`
> Read coordinate information and other information from files.