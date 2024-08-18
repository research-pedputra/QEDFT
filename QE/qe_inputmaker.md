## SCF Input File Generator

This script generates a Self-Consistent Field (SCF) input file for Quantum ESPRESSO based on a CIF file. It reads the atomic structure from the CIF file and generates the necessary input parameters for SCF calculations.

## Explanation

- **Imports**: The script imports `read` from `ase.io` to read CIF files.
  
- **`SCFGenerator` Class**: Manages the creation of the SCF input file.

  - **`__init__` Method**: Initializes with a CIF file path and output directory. Reads the atomic structure and sets the output directory.
  
  - **`_get_unique_atom_symbols` Method**: Returns a sorted string of unique atomic symbols.
  
  - **`_write_atomic_species` Method**: Writes the `ATOMIC_SPECIES` section, including element types, masses, and pseudopotential file names.
  
  - **`_write_atomic_positions` Method**: Writes atomic positions in the `ATOMIC_POSITIONS` section.
  
  - **`_write_cell_parameters` Method**: Writes lattice vectors in the `CELL_PARAMETERS` section.
  
  - **`get_input` Method**: Creates and writes the SCF input file, including sections for control parameters, system details, and electron settings. It also calls methods to write atomic species, positions, and cell parameters, and sets up k-points for the automatic grid.

```python
from ase.io import read

class SCFGenerator:
    def __init__(self, cif_file_path, output_dir):
        self.atoms = read(cif_file_path)
        self.output_dir = output_dir

    def _get_unique_atom_symbols(self):
        return ''.join(sorted(set(self.atoms.get_chemical_symbols())))

    def _write_atomic_species(self, f):
        unique_atom_types = set(self.atoms.get_chemical_symbols())
        f.write("ATOMIC_SPECIES\n")
        for symbol in unique_atom_types:
            f.write(f" {symbol} {self.atoms.get_masses()[self.atoms.get_chemical_symbols().index(symbol)]} {symbol}.UPF\n")

    def _write_atomic_positions(self, f):
        f.write("ATOMIC_POSITIONS angstrom\n")
        for symbol, position in zip(self.atoms.get_chemical_symbols(), self.atoms.get_positions()):
            f.write(f"  {symbol}  {position[0]:.8f}  {position[1]:.8f}  {position[2]:.8f}\n")

    def _write_cell_parameters(self, f):
        f.write("CELL_PARAMETERS angstrom\n")
        for vector in self.atoms.cell:
            f.write(f"  {vector[0]:.8f}  {vector[1]:.8f}  {vector[2]:.8f}\n")

    def get_input(self):
        output_file_path = f"{self.output_dir}/scf_{self.atoms.symbols}.in"
        with open(output_file_path, 'w') as f:
            f.write("&CONTROL\n")
            f.write("  calculation = 'scf'\n")
            f.write("  etot_conv_thr =   5.0000000000d-05\n")
            f.write("  forc_conv_thr =   1.0000000000d-04\n")
            f.write(f"  outdir = './out_{self.atoms.symbols}/'\n")
            f.write(f"  prefix = '{self.atoms.symbols}'\n")
            f.write("  pseudo_dir = './pseudo/'\n")
            f.write("  tprnfor = .true.\n")
            f.write("  tstress = .true.\n")
            f.write("  verbosity = 'high'\n")
            f.write("/\n")
            f.write("&SYSTEM\n")
            f.write("  degauss = 1.4699723600d-02\n")
            f.write("  ecutrho = 3.2000000000d+02\n")
            f.write("  ecutwfc = 4.0000000000d_01\n")
            f.write("  ibrav = 0\n")
            f.write(f"  nat = {len(self.atoms)}\n")
            f.write(f"  ntyp = {len(set(self.atoms.get_chemical_symbols()))}\n")
            f.write("  occupations = 'smearing'\n")
            f.write("  smearing = 'gaussian'\n")
            f.write("/\n")
            f.write("&ELECTRONS\n")
            f.write("  conv_thr = 1.0e-8\n")
            f.write("  electron_maxstep = 80\n")
            f.write("  mixing_beta = 0.4\n")
            f.write("/\n")

            self._write_atomic_species(f)
            self._write_atomic_positions(f)
            
            f.write("K_POINTS automatic\n")
            f.write("3 3 3 0 0 0\n")
            self._write_cell_parameters(f)
```


