# QEDFT Toolkit

A practical Python toolkit for Quantum ESPRESSO and XRD workflows.

## Structure

- `myqetools/` (inside `QE/`) → high-level API (recommended entry point)
- `QE/` → existing QE utilities + wrapper package
- `XRD/` → existing XRD processing utilities (still imported by the wrapper)

## What is new

`QE/myqetools` now includes a **Quantum ESPRESSO binary wrapper** that provides:

1. Input generation wrappers (`scf`, `nscf`, `bands`, `projwfc`).
2. A unified `run()` interface for QE executables.
3. A registry of supported QE binaries (`supported_binaries()`).

## Install dependencies

```bash
pip install numpy matplotlib scipy ase PyQt5
```

## Quick usage

### 1) Use one API for QE input generation + execution

```python
from QE.myqetools import QuantumEspressoAPI

qe = QuantumEspressoAPI(qe_bin_dir="/path/to/qe/bin", mpiexec="mpirun", mpi_nprocs=8)

# generate input
scf = qe.scf_or_nscf("NiO.cif", "./inputs", calculation_type="scf")
scf.get_input()

# inspect available binaries
print(qe.supported_binaries())

# run pw.x
qe.run("pw", input_file="./inputs/scf_NiO.in", output_file="./outputs/scf.out")

# run bands.x
qe.run("bands", input_file="./inputs/NiO.bands.in", output_file="./outputs/bands.out")
```

### 2) Dry-run command construction

```python
from QE.myqetools import QERunner

runner = QERunner(qe_bin_dir="/path/to/qe/bin")
result = runner.run("projwfc", input_file="projwfc.in", dry_run=True)
print(" ".join(result.command))
```

### 3) Keep using analysis wrappers

```python
import numpy as np
from QE.myqetools import BandAnalysis, PDOSAnalysis

bands = BandAnalysis().load(np.loadtxt("sample.bands.gnu"), 5.2, [("G", 0.0), ("X", 1.0)])
bands.plot()

pd = PDOSAnalysis()
data = pd.read({"O_p": ("blue", "O_p.dat")}, fermi_energy=3.5)
pd.plot(data, xlim=(-5, 5), ylim=(-20, 20))
```

## Supported QE binaries

The binary wrapper supports shorthand aliases for many QE executables, including:

- `pw`, `cp`, `pp`, `dos`, `bands`, `projwfc`, `ph`, `q2r`, `matdyn`, `dynmat`, `neb`
- `pw2wannier90`, `pw2bgw`, `pw2gw`, `pw2casino`, `pw2critic`
- `turbo_lanczos`, `turbo_davidson`, `turbo_eels`, `turbo_spectrum`
- `epsilon`, `xspectra`, `epw`, `hp`, `ld1`, `atomic`, `cppp`, `upfconv`

You can always call `qe.run("some_binary.x", ...)` for explicit executables.
