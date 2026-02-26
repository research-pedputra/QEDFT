# myqetools API Guide

This guide explains the modernized API in this repository, with emphasis on the new Quantum ESPRESSO binary wrapper.

## Main entry points

```python
from QE.myqetools import (
    QuantumEspressoAPI,
    QERunner,
    QEInputFactory,
    BandAnalysis,
    PDOSAnalysis,
    ChargeDensityAnalysis,
    StructureTools,
    XRDAnalysis,
)
```

## 1) Quantum ESPRESSO wrapper (inputs + binaries)

### `QuantumEspressoAPI`

A unified façade for QE workflows.

- `scf_or_nscf(cif, outdir, calculation_type="scf")`
- `pdos_input(cif, outdir)`
- `band_input(cif, outdir)`
- `bands_post_input(cif, outdir)`
- `supported_binaries()`
- `run(binary, input_file=..., output_file=..., npool=..., extra_args=...)`

Example:

```python
from QE.myqetools import QuantumEspressoAPI

qe = QuantumEspressoAPI(qe_bin_dir="/opt/qe/bin", mpiexec="mpirun", mpi_nprocs=16)

# make SCF input
qe.scf_or_nscf("sample.cif", "./inputs", "scf").get_input()

# execute pw.x
qe.run("pw", input_file="./inputs/scf_sample.in", output_file="./outputs/scf.out")

# execute projected DOS
qe.run("projwfc", input_file="./inputs/projwfc.in", output_file="./outputs/projwfc.out")
```

## 2) Low-level binary execution only

### `QERunner`

If you only want command execution without the extra façade:

```python
from QE.myqetools import QERunner

runner = QERunner(qe_bin_dir="/opt/qe/bin")
print(runner.supported_binaries())

# dry run
preview = runner.run("bands", input_file="bands.in", dry_run=True)
print(" ".join(preview.command))
```

## 3) Existing analysis wrappers

### `BandAnalysis`

- `load(data, fermi_energy, kpoints)`
- `plot()`
- `plot_band(band_number, overlay=True)`
- `plot_shifted(start_band_number, shift_value)`
- `combine_with_pdos(...)`

### `PDOSAnalysis`

- `read(files, fermi_energy)`
- `plot(aopdos_data, xlim=None, ylim=None)`

### `ChargeDensityAnalysis`

- initialize with `cube_files`, `dimensions`, `scale_factor`
- call `run(output_file)`

### `StructureTools`

- `qe_to_xyz(input_file, output_file)`
- `html_view(file_path)`

### `XRDAnalysis`

- `subtract_background(x, y, tol=1)`
- `find_peaks(...)`
- `plot_folder(folder_path)`

## QE binaries coverage

Supported shorthand aliases include major executables:

`pw`, `cp`, `pp`, `dos`, `bands`, `projwfc`, `ph`, `q2r`, `matdyn`, `dynmat`, `neb`,
`pw2wannier90`, `pw2bgw`, `pw2gw`, `pw2casino`, `pw2critic`, `turbo_lanczos`,
`turbo_davidson`, `turbo_eels`, `turbo_spectrum`, `epsilon`, `xspectra`, `epw`, `hp`,
`ld1`, `atomic`, `cppp`, `upfconv`.

You can also pass an explicit binary name such as `"pw.x"` directly.
