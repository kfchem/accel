ACCeL: Automated Conformer Calculations and Exploration Library
================================================================

ACCeL is a Python package designed to streamline the processing of conformational ensembles in computational chemistry.

It provides a simple and flexible framework for managing large sets of molecular structures, especially conformers and their variants. ACCeL supports automated workflows such as:

- Reading and organizing output files from quantum chemistry calculations
- Extracting atomic coordinates, energies, and thermodynamic corrections
- Filtering structures based on energy thresholds or RMSD values
- Generating input files using user-defined templates
- Labeling and organizing conformers for downstream processing
- Exporting data to CSV or coordinate formats

The core object, ``Box``, enables method-chained operations to efficiently process molecular datasets. ACCeL is particularly useful for researchers working with large conformer sets and quantum chemical outputs.

Features
--------

- Automated file handling and structure filtering
- Energy-based and structural deduplication
- Labeling and grouping by file name
- Boltzmann population and free energy calculations
- Input/output support for XYZ and MOL formats
- Extensible plugin system

Project Information
-------------------

- Repository: https://github.com/kfchem/accel
- License: MIT
- Python version: 3.9+
