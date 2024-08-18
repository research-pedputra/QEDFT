# Quantum Espresso Input File Generators

This Python script provides classes for generating input files for Quantum Espresso calculations based on crystal structure data from CIF files. It includes the following classes:

## `SCFNSCFGenerator`

Generates input files for both Self-Consistent Field (SCF) and Non-Self-Consistent Field (NSCF) calculations. 

- **Constructor**:
  - `cif_file_path`: Path to the CIF file containing the crystal structure.
  - `output_dir`: Directory to save the generated input files.
  - `calculation_type`: Type of calculation (`'scf'` or `'nscf'`).

- **Methods**:
  - `get_input()`: Generates the Quantum Espresso input file based on the specified calculation type.

## `PDOSGenerator`

Generates input files for Partial Density of States (PDOS) calculations.

- **Constructor**:
  - `cif_file_path`: Path to the CIF file containing the crystal structure.
  - `output_dir`: Directory to save the generated input files.

- **Methods**:
  - `get_input()`: Generates the PDOS input file.

## `BandGenerator`

Generates input files for band structure calculations and extracts the number of bands from SCF output.

- **Constructor**:
  - `cif_file_path`: Path to the CIF file containing the crystal structure.
  - `output_dir`: Directory to save the generated input files.

- **Methods**:
  - `get_input(kpoints=None)`: Generates the band structure input file. Optionally specify k-points if not using automatic k-point grid.

## `BandsGenerator`

Generates input files specifically for band structure calculations.

- **Constructor**:
  - `cif_file_path`: Path to the CIF file containing the crystal structure.
  - `output_dir`: Directory to save the generated input files.

- **Methods**:
  - `get_input()`: Generates the band structure input file.

### Usage

To use these classes, initialize them with the appropriate CIF file path and output directory, and call the `get_input()` method to generate the input files for Quantum Espresso.

Ensure that the `ase` library is installed and that you have set up your Quantum Espresso pseudopotential files and output directories correctly.

#### Import the class

```python
from qeinputmaker import SCFNSCFGenerator
```

#### Initialize for file input making
```python
file = classname('path/to/cif_file.cif', 'path/to/output_dir', 'typeofinputfile')
```
#### Generate the file input
```python
file.get_input()
```
