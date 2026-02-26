# myqetools API Guide (QE coverage)

This guide documents the QE API designed to cover broad Quantum ESPRESSO functionality.

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

## 1) QuantumEspressoAPI (high-level façade)

### Responsibilities

- Generate common QE inputs from CIF (`scf`, `nscf`, `bands`, `projwfc` helpers).
- Run QE binaries via one normalized `run()` method.
- Run chained workflows via `run_many()`.
- Expose capability maps via `supported_binaries()` and `supported_groups()`.
- Detect available executables at runtime with `discover_available_binaries()`.

### Typical use

```python
from QE.myqetools import QuantumEspressoAPI

qe = QuantumEspressoAPI(qe_bin_dir="/opt/qe/bin", mpiexec="mpirun", mpi_nprocs=32)

# Input generation
qe.scf_or_nscf("sys.cif", "./in", "scf").get_input()

# Binary execution
qe.run("pw", input_file="./in/scf_sys.in", output_file="./out/scf.out", npool=4)

# Multi-step workflow
qe.run_many([
    {"binary": "pw", "input_file": "./in/scf_sys.in", "output_file": "./out/scf.out"},
    {"binary": "pw", "input_file": "./in/nscf_sys.in", "output_file": "./out/nscf.out"},
    {"binary": "bands", "input_file": "./in/sys.bands.in", "output_file": "./out/bands.out"},
])
```

## 2) QERunner (low-level execution engine)

Use `QERunner` when you only need execution, not the full façade.

```python
from QE.myqetools import QERunner

runner = QERunner(qe_bin_dir="/opt/qe/bin")
print(runner.supported_groups())
print(runner.discover_available_binaries())

preview = runner.run("ph", input_file="ph.in", dry_run=True)
print(" ".join(preview.command))
```

## 3) QE capability domains supported by aliases

- PW/post-processing
- Phonons
- NEB
- Interfaces/converters
- TDDFPT
- Spectroscopy
- Transport/response
- Pseudopotential/atomic tools

All aliases are exposed in `QE_BINARY_ALIASES`. Domain grouping is exposed in `QE_BINARY_GROUPS`.

## 4) Additional analysis wrappers

These wrap existing repository utilities:

- `BandAnalysis`
- `PDOSAnalysis`
- `ChargeDensityAnalysis`
- `StructureTools`
- `XRDAnalysis`

Use them directly for plotting/analysis while using `QuantumEspressoAPI` for execution orchestration.
