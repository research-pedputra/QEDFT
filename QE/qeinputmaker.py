from ase.io import read
import re

class SCFNSCFGenerator:
    def __init__(self, cif_file_path, output_dir, calculation_type='scf'):
        self.atoms = read(cif_file_path)
        self.output_dir = output_dir
        self.calculation_type = calculation_type

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
        if self.calculation_type not in ['scf', 'nscf']:
            raise ValueError("Invalid calculation type. Choose 'scf' or 'nscf'.")
        
        output_file_path = f"{self.output_dir}/{self.calculation_type}_{self.atoms.symbols}.in"
        with open(output_file_path, 'w') as f:
            f.write("&CONTROL\n")
            f.write(f"  calculation = '{self.calculation_type}'\n")
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
            
            if self.calculation_type == 'scf':
                f.write("K_POINTS automatic\n")
                f.write("3 3 3 0 0 0\n")
            elif self.calculation_type == 'nscf':
                f.write("K_POINTS automatic\n")
                f.write("12 12 12 0 0 0\n")
            
            self._write_cell_parameters(f)

class PDOSGenerator:
    def __init__(self, cif_file_path, output_dir):
        self.atoms = read(cif_file_path)
        self.output_dir = output_dir

    def get_input(self):
        output_file_path = f"{self.output_dir}/PDOS_{self.atoms.symbols}.in"
        with open(output_file_path, 'w') as f:
            f.write("&PROJWFC\n")
            f.write(f"  outdir = './out_{self.atoms.symbols}/'\n")
            f.write(f"  prefix = '{self.atoms.symbols}'\n")
            f.write("  degauss = 0.01\n")
            f.write("/\n")

class BandGenerator:
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

    def _get_nbnd_from_scf_output(self):
        # Read SCF output file
        scf_output_file = f"./scf_{self.atoms.symbols}.out"
        with open(scf_output_file, 'r') as f:
            for line in f:
                match = re.match(r'\s*number\s+of\s+electrons\s+=\s+(\d+\.?\d*)', line)
                if match:
                    nbnd = int(float(match.group(1)))
                    return nbnd
        raise ValueError("Could not find 'number of electrons' in SCF output file.")

    def get_input(self, kpoints=None):
        output_file_path = f"{self.output_dir}/band_{self.atoms.symbols}.in"
        with open(output_file_path, 'w') as f:
            f.write("&CONTROL\n")
            f.write("  calculation = 'bands'\n")
            f.write("  etot_conv_thr =   5.0000000000d-05\n")
            f.write("  forc_conv_thr =   1.0000000000d-04\n")
            f.write(f"  outdir = './out_scf_{self.atoms.symbols}/'\n")
            f.write(f"  prefix = '{self.atoms.symbols}'\n")
            f.write("  pseudo_dir = './pseudo/'\n")
            f.write("  tprnfor = .true.\n")
            f.write("  tstress = .true.\n")
            f.write("  verbosity = 'high'\n")
            f.write("/\n")
            f.write("&SYSTEM\n")
            f.write("  degauss =   1.4699723600d-02\n")
            f.write("  ecutrho =   2.5000000000d+02\n")
            f.write("  ecutwfc =   6.0000000000d+01\n")
            f.write("  ibrav = 0\n")
            f.write(f"  nat = {len(self.atoms)}\n")
            f.write(f"  nbnd = {self._get_nbnd_from_scf_output()}\n")
            f.write("  nosym = .false.\n")
            f.write(f"  ntyp = {len(set(self.atoms.get_chemical_symbols()))}\n")
            f.write("  occupations = 'smearing'\n")
            f.write("  smearing = 'gaussian'\n")
            f.write("  nspin=2,\n")
            f.write("  starting_magnetization(3)= 0.6\n")
            f.write("  starting_magnetization(4)= 0.6\n")
            f.write("/\n")
            f.write("&ELECTRONS\n")
            f.write("  conv_thr =   1.0000000000d-09\n")
            f.write("  electron_maxstep = 80\n")
            f.write("  mixing_beta =   4.0000000000d-01\n")
            f.write("/\n")
            
            self._write_atomic_species(f)
            self._write_atomic_positions(f)
            
            if kpoints is None:
                f.write("K_POINTS automatic\n")
                f.write("12 12 12 0 0 0\n")
            else:
                f.write("K_POINTS {crystal_b}\n")
                f.write(f"{len(kpoints)}\n")
                for kpoint in kpoints:
                    f.write(f"{kpoint}\n")

            self._write_cell_parameters(f)

class BandsGenerator:
    def __init__(self, cif_file_path, output_dir):
        self.atoms = read(cif_file_path)
        self.output_dir = output_dir

    def get_input(self):
        output_file_path = f"{self.output_dir}/{self.atoms.symbols}.bands.in"
        with open(output_file_path, 'w') as f:
            f.write("&bands\n")
            f.write(f"  outdir = './out_{self.atoms.symbols}/'\n")
            f.write(f"  prefix = '{self.atoms.symbols}'\n")
            f.write(f"  filband = '{self.atoms.symbols}.bands'\n")
            f.write("/\n")


