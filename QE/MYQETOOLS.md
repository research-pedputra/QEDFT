# myqetools API Guide (QE coverage update)

This guide documents the updated QE API that now covers broader QE capability domains and workflow orchestration.

## Entry points

```python
from QE.myqetools import (
    QuantumEspressoAPI,
    QERunner,
    QERunResult,
    QE_BINARY_ALIASES,
    QE_BINARY_GROUPS,
    QEInputFactory,
    BandAnalysis,
    PDOSAnalysis,
    ChargeDensityAnalysis,
    StructureTools,
    XRDAnalysis,
)
```

## QuantumEspressoAPI

High-level façade around input generation + execution.

### Added execution helpers

- `supported_binaries()`
- `supported_groups()`
- `discover_available_binaries()`
- `run(...)`
- `run_many(steps)`
- `run_group(group, input_map, output_dir=".")`

### Example

```python
from QE.myqetools import QuantumEspressoAPI

qe = QuantumEspressoAPI(qe_bin_dir="/opt/qe/bin", mpiexec="mpirun", mpi_nprocs=8)

qe.scf_or_nscf("sys.cif", "./in", "scf").get_input()
qe.run("pw", input_file="./in/scf_sys.in", output_file="./out/scf.out", npool=2)

qe.run_many([
    {"binary": "pw", "input_file": "./in/scf_sys.in", "output_file": "./out/scf.out"},
    {"binary": "pw", "input_file": "./in/nscf_sys.in", "output_file": "./out/nscf.out"},
    {"binary": "bands", "input_file": "./in/sys.bands.in", "output_file": "./out/bands.out"},
])
```

## QERunner

Low-level runner for direct control.

```python
from QE.myqetools import QERunner

runner = QERunner(qe_bin_dir="/opt/qe/bin")
print(runner.supported_groups())
print(runner.discover_available_binaries())
preview = runner.run("ph", input_file="ph.in", dry_run=True)
print(" ".join(preview.command))
```

Group execution shortcut:

```python
runner.run_group(
    group="interfaces",
    input_map={
        "pw2wannier90": "pw2wan.in",
        "wannier_ham": "wham.in",
    },
    output_dir="./if_out",
)
```

## Capability domains

- PW/post-processing
- Phonon
- NEB
- Interfaces/converters
- TDDFPT
- Spectroscopy
- Transport/response
- Pseudo/atomic

## Analysis wrappers retained

The following wrappers remain available and unchanged in intent:

- `BandAnalysis`
- `PDOSAnalysis`
- `ChargeDensityAnalysis`
- `StructureTools`
- `XRDAnalysis`


## Band unfolding (Python-first)

Use the Python unfolding path (bands_unfold-like) by default:

```python
import numpy as np
from QE.myqetools import QuantumEspressoAPI

qe = QuantumEspressoAPI()

# synthetic or parsed supercell data
k_sc = np.random.rand(200, 3)
E_sc = np.random.normal(loc=0.0, scale=2.0, size=(200, 8))
S = np.diag([2, 2, 1])

unfolded = qe.run_python_unfold(
    kpoints_supercell=k_sc,
    energies_supercell=E_sc,
    supercell_matrix=S,
    sigma_k=0.02,
    sigma_e=0.08,
)

fig, ax, im = qe.unfolding().plot_intensity_notebook(unfolded, fermi_energy=0.0)
qe.unfolding().save_unfolded_python(unfolded, "./out/unfolded_python.dat")
```

Optional legacy path (if desired):

```python
qe.run_bands_unfold(input_file="./in/unfold.in", output_file="./out/unfold.log")
```
