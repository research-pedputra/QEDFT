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
