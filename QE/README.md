# QEDFT Toolkit (QE API)

This repository exposes a broad Quantum ESPRESSO API under `QE/myqetools`, intended to cover most practical QE capabilities while preserving the original utility modules.

## Location

- API package: `QE/myqetools/`
- Core QE utilities: `QE/*.py`
- XRD helpers consumed by API: `XRD/`

## Main import

```python
from QE.myqetools import QuantumEspressoAPI, QERunner
```

## Coverage in this API

The wrapper includes:

1. Input-generation wrappers (`scf`, `nscf`, `bands`, `pdos`).
2. Binary execution for many QE tools through aliases.
3. Grouped capability metadata (`supported_groups`).
4. Runtime binary discovery (`discover_available_binaries`).
5. Multi-step execution (`run_many`) and group execution (`run_group`).

## Quick usage

### 1) Generate inputs + run a binary

```python
from QE.myqetools import QuantumEspressoAPI

qe = QuantumEspressoAPI(qe_bin_dir="/opt/qe/bin", mpiexec="mpirun", mpi_nprocs=16)
qe.scf_or_nscf("sample.cif", "./inputs", calculation_type="scf").get_input()
qe.run("pw", input_file="./inputs/scf_sample.in", output_file="./outputs/scf.out")
```

### 2) Inspect groups and discover available binaries

```python
from QE.myqetools import QERunner

runner = QERunner(qe_bin_dir="/opt/qe/bin")
print(runner.supported_groups())
print(runner.discover_available_binaries())
```

### 3) Run multiple steps

```python
from QE.myqetools import QuantumEspressoAPI

qe = QuantumEspressoAPI(qe_bin_dir="/opt/qe/bin")
steps = [
    {"binary": "pw", "input_file": "scf.in", "output_file": "scf.out"},
    {"binary": "pw", "input_file": "nscf.in", "output_file": "nscf.out"},
    {"binary": "bands", "input_file": "bands_pp.in", "output_file": "bands_pp.out"},
]
qe.run_many(steps)
```

### 4) Run a capability group with alias->input mapping

```python
from QE.myqetools import QERunner

runner = QERunner(qe_bin_dir="/opt/qe/bin")
runner.run_group(
    group="phonon",
    input_map={
        "ph": "ph.in",
        "q2r": "q2r.in",
        "matdyn": "matdyn.in",
    },
    output_dir="./phonon_out",
)
```

## QE capability groups exposed

- `pw`
- `phonon`
- `neb`
- `interfaces`
- `tddfpt`
- `spectroscopy`
- `transport`
- `pseudo`

The aliases can be read directly from `QE_BINARY_ALIASES`, and non-listed binaries can still be run by explicit name:

```python
qe.run("my_custom_qe_binary.x", input_file="input.in", output_file="out.log")
```
