# QEDFT Toolkit (QE API)

This repository now exposes a **broad Quantum ESPRESSO API** under `QE/myqetools`, designed to cover nearly all QE capabilities through a unified execution layer.

## Location

- API package: `QE/myqetools/`
- Core QE utilities: `QE/*.py`
- XRD helpers used by API: `XRD/`

## Main import

```python
from QE.myqetools import QuantumEspressoAPI, QERunner
```

## What “covering QE capability” means here

The API provides:

1. **Input-generation wrappers** for SCF/NSCF/bands/PDOS input files from CIF.
2. **Full binary runner** that can run any executable name (`"something.x"`) and also ships with extensive shorthand aliases.
3. **Workflow execution** (`run_many`) for chained QE jobs.
4. **Capability grouping** (`supported_groups`) and runtime discovery (`discover_available_binaries`).

## Quick examples

### A) Build and run common PW workflow

```python
from QE.myqetools import QuantumEspressoAPI

qe = QuantumEspressoAPI(qe_bin_dir="/opt/qe/bin", mpiexec="mpirun", mpi_nprocs=16)

qe.scf_or_nscf("sample.cif", "./inputs", calculation_type="scf").get_input()
qe.run("pw", input_file="./inputs/scf_sample.in", output_file="./outputs/scf.out")
```

### B) Discover what is available on your machine

```python
from QE.myqetools import QERunner

runner = QERunner(qe_bin_dir="/opt/qe/bin")
print(runner.supported_groups())
print(runner.discover_available_binaries())
```

### C) Run a full chain (SCF -> NSCF -> bands.x -> projwfc.x)

```python
from QE.myqetools import QuantumEspressoAPI

qe = QuantumEspressoAPI(qe_bin_dir="/opt/qe/bin", mpiexec="mpirun", mpi_nprocs=8)

steps = [
    {"binary": "pw", "input_file": "scf.in", "output_file": "scf.out"},
    {"binary": "pw", "input_file": "nscf.in", "output_file": "nscf.out"},
    {"binary": "bands", "input_file": "bands_pp.in", "output_file": "bands_pp.out"},
    {"binary": "projwfc", "input_file": "projwfc.in", "output_file": "projwfc.out"},
]
qe.run_many(steps)
```

## Coverage areas (aliases)

- **PW / postproc**: `pw`, `cp`, `pp`, `dos`, `bands`, `projwfc`, `average`, `plan_avg`, `plotband`, `plotproj`, `plotrho`, `sumpdos`, `open_grid`
- **Phonons**: `ph`, `q2r`, `matdyn`, `dynmat`, `lambda`
- **NEB**: `neb`, `path_interpolation`
- **Interfaces/converters**: `pw2wannier90`, `wannier_ham`, `pw2bgw`, `pw2gw`, `pw2casino`, `pw2critic`, `bands_unfold`, `kcw`, `kcw_pp`
- **TDDFPT**: `turbo_lanczos`, `turbo_davidson`, `turbo_eels`, `turbo_spectrum`
- **Spectroscopy**: `gww_fit`, `gww_pp`, `epsilon`, `xspectra`
- **Transport/response**: `epw`, `hp`
- **Pseudo/atomic tools**: `ld1`, `atomic`, `cppp`, `upfconv`, `virtual_v2`

If your build has an executable not listed, run it directly:

```python
qe.run("my_custom_qe_binary.x", input_file="input.in", output_file="out.log")
```
